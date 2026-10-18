#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

class Quad4Hourglass
{
public:
  using Nodes = std::array<std::array<double, 2>, 4>;
  using Mode = std::array<double, 4>;
  using Matrix = std::array<std::array<double, 4>, 4>;

  static Mode mode(const Nodes & nodes)
  {
    constexpr Mode seed = {{1.0, -1.0, 1.0, -1.0}};
    double normal[3][4] = {};
    for (unsigned int i = 0; i < 4; ++i)
    {
      const double basis[3] = {1.0, nodes[i][0], nodes[i][1]};
      for (unsigned int row = 0; row < 3; ++row)
      {
        normal[row][3] += basis[row] * seed[i];
        for (unsigned int col = 0; col < 3; ++col)
          normal[row][col] += basis[row] * basis[col];
      }
    }

    for (unsigned int pivot = 0; pivot < 3; ++pivot)
    {
      unsigned int best = pivot;
      for (unsigned int row = pivot + 1; row < 3; ++row)
        if (std::abs(normal[row][pivot]) > std::abs(normal[best][pivot]))
          best = row;
      if (std::abs(normal[best][pivot]) < 1e-14)
        throw std::runtime_error("degenerate QUAD4 affine basis");
      for (unsigned int col = pivot; col < 4; ++col)
        std::swap(normal[pivot][col], normal[best][col]);
      const double scale = normal[pivot][pivot];
      for (unsigned int col = pivot; col < 4; ++col)
        normal[pivot][col] /= scale;
      for (unsigned int row = 0; row < 3; ++row)
        if (row != pivot)
        {
          const double factor = normal[row][pivot];
          for (unsigned int col = pivot; col < 4; ++col)
            normal[row][col] -= factor * normal[pivot][col];
        }
    }

    Mode gamma;
    double norm = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
      gamma[i] = seed[i] - normal[0][3] - normal[1][3] * nodes[i][0] -
                 normal[2][3] * nodes[i][1];
      norm += gamma[i] * gamma[i];
    }
    if (norm < 1e-24)
      throw std::runtime_error("degenerate QUAD4 hourglass mode");
    norm = std::sqrt(norm);
    for (auto & value : gamma)
      value /= norm;
    return gamma;
  }

  static double area(const Nodes & nodes)
  {
    double twice_area = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
      const unsigned int next = (i + 1) % 4;
      twice_area += nodes[i][0] * nodes[next][1] - nodes[next][0] * nodes[i][1];
    }
    const double value = 0.5 * std::abs(twice_area);
    if (value < 1e-14)
      throw std::runtime_error("degenerate QUAD4 area");
    return value;
  }

  static double characteristicLengthSquared(const Nodes & nodes)
  {
    double edge_square_sum = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
      const unsigned int next = (i + 1) % 4;
      const double dx = nodes[next][0] - nodes[i][0];
      const double dy = nodes[next][1] - nodes[i][1];
      edge_square_sum += dx * dx + dy * dy;
    }
    const double value = 0.25 * edge_square_sum;
    if (value < 1e-14)
      throw std::runtime_error("degenerate QUAD4 edge length");
    return value;
  }

  // Per-unit-area coefficient for the energy-checked amplitude h^T u / 4.
  static double gradientRawStiffnessPerShear(const Nodes & nodes)
  {
    constexpr std::array<std::array<double, 2>, 4> dn =
        {{{{-0.25, -0.25}}, {{0.25, -0.25}}, {{0.25, 0.25}}, {{-0.25, 0.25}}}};
    constexpr Mode seed = {{1.0, -1.0, 1.0, -1.0}};
    double jac[2][2] = {};
    std::array<double, 2> seed_moment = {};
    for (unsigned int i = 0; i < 4; ++i)
      for (unsigned int a = 0; a < 2; ++a)
      {
        const double x = nodes[i][a] - nodes[0][a];
        seed_moment[a] += x * seed[i];
        for (unsigned int b = 0; b < 2; ++b)
          jac[a][b] += x * dn[i][b];
      }
    const double det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    if (!std::isfinite(det) || det <= 1e-14)
      throw std::runtime_error("degenerate or inverted QUAD4 center Jacobian");
    double metric = 0.0, norm_squared = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
      const double bx = (dn[i][0] * jac[1][1] - dn[i][1] * jac[1][0]) / det;
      const double by = (-dn[i][0] * jac[0][1] + dn[i][1] * jac[0][0]) / det;
      const double h = seed[i] - bx * seed_moment[0] - by * seed_moment[1];
      metric += bx * bx + by * by;
      norm_squared += h * h;
    }
    return (0.005 / 16.0) * norm_squared * metric;
  }

  template <typename Scalar>
  static Scalar gradientRawResidualPerShear(const Nodes & nodes,
                                           double alpha,
                                           double dt,
                                           double zeta,
                                           double history_correction,
                                           const Scalar & current,
                                           double old,
                                           double older)
  {
    Scalar elastic = current, damping = 0.0;
    if (dt > 0)
    {
      elastic = (1.0 + alpha) * current - alpha * old;
      damping = current * ((1.0 + alpha) * zeta / dt) -
                old * ((1.0 + 2.0 * alpha) * zeta / dt) + older * (alpha * zeta / dt) +
                history_correction * (older - old);
    }
    // Direct Kref+Clegacy avoids subtracting nearly equal elastic terms (~400:1).
    return gradientRawStiffnessPerShear(nodes) * elastic +
           damping / characteristicLengthSquared(nodes);
  }

  static Matrix consistentMassCorrection(const Nodes & nodes)
  {
    Matrix correction = {};
    constexpr double gauss = 0.57735026918962576451;
    for (const double xi : {-gauss, gauss})
      for (const double eta : {-gauss, gauss})
      {
        const std::array<double, 4> shape = {{0.25 * (1.0 - xi) * (1.0 - eta),
                                              0.25 * (1.0 + xi) * (1.0 - eta),
                                              0.25 * (1.0 + xi) * (1.0 + eta),
                                              0.25 * (1.0 - xi) * (1.0 + eta)}};
        const std::array<std::array<double, 2>, 4> derivative =
            {{{{-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)}},
              {{0.25 * (1.0 - eta), -0.25 * (1.0 + xi)}},
              {{0.25 * (1.0 + eta), 0.25 * (1.0 + xi)}},
              {{-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)}}}};
        double dx_dxi = 0.0;
        double dx_deta = 0.0;
        double dy_dxi = 0.0;
        double dy_deta = 0.0;
        for (unsigned int i = 0; i < 4; ++i)
        {
          dx_dxi += nodes[i][0] * derivative[i][0];
          dx_deta += nodes[i][0] * derivative[i][1];
          dy_dxi += nodes[i][1] * derivative[i][0];
          dy_deta += nodes[i][1] * derivative[i][1];
        }
        const double determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;
        if (determinant <= 1e-14)
          throw std::runtime_error("degenerate or inverted QUAD4 Jacobian");
        for (unsigned int i = 0; i < 4; ++i)
          for (unsigned int j = 0; j < 4; ++j)
            correction[i][j] += shape[i] * shape[j] * determinant;
      }

    const double reduced_entry = area(nodes) / 16.0;
    for (auto & row : correction)
      for (auto & value : row)
        value -= reduced_entry;
    return correction;
  }
};
