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
};
