#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

// Reference-geometry center mode and Q4-Q1 compensation. No material state/update.
namespace DS04GSP1
{
using Nodes = std::array<std::array<double, 2>, 4>;
using Matrix = std::array<std::array<double, 4>, 4>;
inline auto derivatives(double xi, double eta)
{
  return Nodes{{{{-(1-eta)/4, -(1-xi)/4}}, {{(1-eta)/4, -(1+xi)/4}},
                {{(1+eta)/4, (1+xi)/4}}, {{-(1+eta)/4, (1-xi)/4}}}};
}
inline std::array<double, 4> shape(double xi, double eta)
{
  return {{(1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4,
           (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4}};
}
template <typename Scalar>
Scalar determinant(const std::array<std::array<Scalar, 2>, 4> & nodes, double xi, double eta)
{
  const auto d = derivatives(xi,eta);
  Scalar a=0.,b=0.,c=0.,e=0.;
  for (unsigned n=0;n<4;++n)
  {
    const auto x=nodes[n][0]-nodes[0][0], y=nodes[n][1]-nodes[0][1];
    a+=x*d[n][0]; b+=x*d[n][1]; c+=y*d[n][0]; e+=y*d[n][1];
  }
  return a*e-b*c;
}
struct Geometry
{
  Nodes gradient;
  std::array<double,4> gamma;
  Matrix mass_correction;
  double area, stiffness_per_shear;
};
inline Geometry geometry(const Nodes & x)
{
  Geometry r{};
  for (const auto & n:x)
    for (double z:n)
      if (!std::isfinite(z)) throw std::invalid_argument("SP1 nonfinite geometry");
  for (double xi:{-1.,0.,1.})
    for (double eta:{-1.,0.,1.})
      if (!(determinant(x,xi,eta)>0)) throw std::invalid_argument("SP1 inverted/degenerate QUAD4");
  const auto d=derivatives(0,0);
  double a=0,b=0,c=0,e=0;
  std::array<double,2> moment{};
  const std::array<double,4> seed{{1,-1,1,-1}};
  for(unsigned n=0;n<4;++n)
  {
    const double xx=x[n][0]-x[0][0], yy=x[n][1]-x[0][1];
    a+=xx*d[n][0]; b+=xx*d[n][1]; c+=yy*d[n][0]; e+=yy*d[n][1];
    moment[0]+=xx*seed[n]; moment[1]+=yy*seed[n];
  }
  const double det=a*e-b*c; r.area=4*det;
  double metric=0;
  for(unsigned n=0;n<4;++n)
  {
    r.gradient[n]={{(d[n][0]*e-d[n][1]*c)/det,(-d[n][0]*b+d[n][1]*a)/det}};
    r.gamma[n]=seed[n]-r.gradient[n][0]*moment[0]-r.gradient[n][1]*moment[1];
    metric+=r.gradient[n][0]*r.gradient[n][0]+r.gradient[n][1]*r.gradient[n][1];
  }
  r.stiffness_per_shear=.005*r.area*metric/16.;
  constexpr double g=.57735026918962576451;
  for(double xi:{-g,g})
    for(double eta:{-g,g})
    {
      const auto N=shape(xi,eta); const double volume=determinant(x,xi,eta);
      for(unsigned i=0;i<4;++i)
        for(unsigned j=0;j<4;++j) r.mass_correction[i][j]+=N[i]*N[j]*volume;
    }
  for(auto & row:r.mass_correction)
    for(auto & v:row) v-=r.area/16.;
  return r;
}
template <typename Scalar>
std::array<Scalar,4> gravityCorrection(const std::array<std::array<Scalar,2>,4> & x)
{
  std::array<Scalar,4> f{{0.,0.,0.,0.}};
  constexpr double g=.57735026918962576451;
  for(double xi:{-g,g})
    for(double eta:{-g,g})
    {
      const auto N=shape(xi,eta); const Scalar v=determinant(x,xi,eta);
      if (!(v>0.)) throw std::invalid_argument("SP1 invalid current gravity geometry");
      for(unsigned i=0;i<4;++i) f[i]+=N[i]*v;
    }
  const Scalar one=determinant(x,0,0);
  for(auto & v:f) v-=one;
  return f;
}
template <typename Scalar>
Scalar hourglass(const Geometry & g, const std::array<Scalar,4> & u, unsigned i, double mu)
{
  // Remove translation before the modal dot product; affine/rotation remains null.
  Scalar q=0.;
  for(unsigned j=0;j<4;++j) q+=g.gamma[j]*(u[j]-u[0]);
  return mu*g.stiffness_per_shear*g.gamma[i]*q;
}
template <typename Scalar>
Scalar trialAcceleration(const Scalar & u, double old, double v, double a, double h, double beta)
{
  if (!(h>0 && beta>0)) throw std::invalid_argument("SP1 invalid Newmark interval");
  return ((u-old)/(h*h)-v/h-(.5-beta)*a)/beta;
}
struct Weights { double current, old, older; };
inline Weights history(double alpha, double zeta, double h)
{
  if (!(h>0)) throw std::invalid_argument("SP1 invalid history interval");
  return {1+alpha+(1+alpha)*zeta/h,-alpha-(1+2*alpha)*zeta/h,alpha*zeta/h};
}
} // namespace DS04GSP1
