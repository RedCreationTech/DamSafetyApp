#include "AqueductShellSection.h"
#include "MooseTypes.h"
#include "gtest/gtest.h"

namespace
{
constexpr double E = 200e9, nu = .3, fy = 400e6, fu = 696e6, epmax = .1;
using V = std::array<double, 4>;
auto step(std::array<double, 3> de, V old = {}, double ep = 0.)
{
  return AqueductShell::planeStress(de, old, ep, E, nu, fy, fu, epmax);
}
}

TEST(AqueductShellSection, ElasticPlaneStress)
{
  auto s = step({1e-4, 0, 0});
  EXPECT_NEAR(s.value[0], E / (1 - nu*nu) * 1e-4, 1e-5);
  EXPECT_NEAR(s.value[1], nu * s.value[0], 1e-5);
  EXPECT_NEAR(s.value[2], 0, 1e-6);
  EXPECT_DOUBLE_EQ(s.plastic_strain, 0);
}

TEST(AqueductShellSection, UniaxialPlasticAndUnload)
{
  // Prescribed transverse strain corresponding to zero transverse stress in
  // uniaxial J2 loading, with the source's linear hardening segment.
  const double H = (fu-fy)/epmax;
  const double e = .01, ep = (E*e-fy)/(E+H), sigma = fy+H*ep;
  auto s = step({e, -nu*sigma/E-.5*ep, 0});
  EXPECT_NEAR(s.value[0], sigma, 1e-3);
  EXPECT_NEAR(s.value[1], 0, 1e-3);
  EXPECT_NEAR(s.value[2], 0, 1e-3);
  EXPECT_NEAR(s.plastic_strain, ep, 1e-12);
  auto u = step({-1e-4, nu*1e-4, 0}, s.value, s.plastic_strain);
  EXPECT_NEAR(u.value[0], sigma-E*1e-4, 1e-3);
  EXPECT_NEAR(u.plastic_strain, ep, 1e-12);
}

TEST(AqueductShellSection, ConstantExtrapolationBeyondLastTableRow)
{
  const double e = .2, ep = e - fu/E;
  auto s = step({e, -nu*fu/E-.5*ep, 0});
  EXPECT_NEAR(s.value[0], fu, 1e-3);
  EXPECT_NEAR(s.value[1], 0, 1e-3);
  EXPECT_NEAR(s.plastic_strain, ep, 1e-12);
}

TEST(AqueductShellSection, PlasticShearAndSymmetry)
{
  auto s = step({0,0,.01});
  auto r = step({0,0,-.01});
  const double equivalent = std::sqrt(3.) * s.value[3];
  EXPECT_NEAR(equivalent, fy+(fu-fy)/epmax*s.plastic_strain, 1e-3);
  EXPECT_NEAR(r.value[3], -s.value[3], 1e-6);
  EXPECT_NEAR(s.value[0], 0, 1e-3);
  EXPECT_NEAR(s.value[2], 0, 1e-3);
}

TEST(AqueductShellSection, SPOSIntegrationAndParallelAxis)
{
  auto z = AqueductShell::layerPositions(.1,.5);
  auto m = AqueductShell::moments(.1,.5);
  double i0=0,i1=0,i2=0;
  for (unsigned i=0;i<3;++i)
  {
    i0 += .1*AqueductShell::layerWeights[i];
    i1 += .1*AqueductShell::layerWeights[i]*z[i];
    i2 += .1*AqueductShell::layerWeights[i]*z[i]*z[i];
  }
  EXPECT_NEAR(i0,m.zero,1e-16);
  EXPECT_NEAR(i1,m.first,1e-16);
  EXPECT_NEAR(i2,m.second,1e-16);
  EXPECT_DOUBLE_EQ(z[2],0);
  EXPECT_DOUBLE_EQ(z[0],-.1);
}

TEST(AqueductShellSection, ADPlasticTangentMatchesCentralDifference)
{
  std::array<ADReal,3> de = {ADReal(.007),ADReal(-.001),ADReal(.002)};
  for(unsigned j=0;j<3;++j)
    Moose::derivInsert(de[j].derivatives(),j,1.);
  auto s=AqueductShell::planeStress(de,V{},0.,E,nu,fy,fu,epmax);
  for(unsigned j=0;j<3;++j)
  {
    std::array<double,3> plus={.007,-.001,.002},minus=plus;
    plus[j]+=1e-7;minus[j]-=1e-7;
    auto p=step(plus),m=step(minus);
    for(unsigned i=0;i<4;++i)
    {
      double fd=(p.value[i]-m.value[i])/2e-7;
      EXPECT_NEAR(s.value[i].derivatives()[j],fd,5e-7*E);
    }
  }
}

TEST(AqueductShellSection, InvalidMaterialIsRejected)
{
  EXPECT_THROW(AqueductShell::planeStress(std::array<double,3>{},V{},0.,E,.5,fy,fu,epmax),std::invalid_argument);
}

TEST(AqueductShellSection, SixRigidBodyModesHaveZeroStrain)
{
  const AqueductShell::Quad q={{{1,2,3},{3,2,3},{3,5,3},{1,5,3}}};
  const AqueductShell::Vec rotation={.2,-.3,.4};
  auto B=AqueductShell::bMatrix(q,{1.4,3.1,3});
  std::array<double,24> dofs{};
  for(unsigned i=0;i<4;++i)
  {
    auto u=AqueductShell::cross(rotation,q[i]);
    for(unsigned j=0;j<3;++j){dofs[6*i+j]=u[j]+j+1;dofs[6*i+3+j]=rotation[j];}
  }
  for(const auto & row:B)
  {
    double strain=0;
    for(unsigned j=0;j<24;++j) strain+=row[j]*dofs[j];
    EXPECT_NEAR(strain,0,1e-14);
  }
}

TEST(AqueductShellSection, MITCTyingPassesConstantCurvatureWithoutShearLocking)
{
  const AqueductShell::Quad q={{{0,0,0},{2,0,0},{2,3,0},{0,3,0}}};
  std::array<double,24> dofs{};
  const double curvature=.001;
  for(unsigned i=0;i<4;++i)
  {
    dofs[6*i+2]=-.5*curvature*q[i][0]*q[i][0];
    dofs[6*i+4]=curvature*q[i][0];
  }
  for(double x:{.2,1.7})for(double y:{.3,2.8})
  {
    auto B=AqueductShell::bMatrix(q,{x,y,0});
    for(unsigned i=0;i<9;++i)
    {
      double strain=0;
      for(unsigned j=0;j<24;++j) strain+=B[i][j]*dofs[j];
      EXPECT_NEAR(strain,i==3?curvature:0,1e-14);
    }
  }
}

TEST(AqueductShellSection, DrillingModeAndGeometryRestrictionsAreExplicit)
{
  AqueductShell::Quad q={{{0,0,0},{2,0,0},{2,3,0},{0,3,0}}};
  auto B=AqueductShell::bMatrix(q,{.6,1.1,0});
  double drill=0;for(unsigned i=0;i<4;++i) drill+=B[8][6*i+5];
  EXPECT_NEAR(drill,1,1e-14);
  q[2][2]=.1;
  EXPECT_THROW(AqueductShell::frame(q),std::invalid_argument);
}
