#include "IntegratedGroundMotion.h"
#include "Quad4Hourglass.h"
#include "gtest/gtest.h"
#include <limits>

TEST(IntegratedGroundMotion, ConstantAccelerationWithInitialState)
{
  const IntegratedGroundMotion motion({0,1,3},{2,2,2},4,5);
  for (double t : {0.0,0.37,1.0,1.8,3.0})
  {
    const auto s=motion.at(t);
    EXPECT_NEAR(s.displacement,4+5*t+t*t,1e-13);
    EXPECT_NEAR(s.velocity,5+2*t,1e-13);
    EXPECT_DOUBLE_EQ(s.acceleration,2);
  }
}
TEST(IntegratedGroundMotion, LinearAccelerationIndependentOfSamplingPartition)
{
  const IntegratedGroundMotion coarse({0,2},{-0.004,1.996});
  const IntegratedGroundMotion fine({0,.07,.13,.51,1,2},{-.004,.066,.126,.506,.996,1.996});
  for(double t : {0.0,.013,.07,.131,.333,1.0,1.719,2.0})
  {
    const auto a=coarse.at(t), b=fine.at(t);
    EXPECT_NEAR(a.displacement,-.002*t*t+t*t*t/6,1e-14);
    EXPECT_NEAR(a.velocity,-.004*t+t*t/2,1e-14);
    EXPECT_NEAR(a.displacement,b.displacement,1e-14);
    EXPECT_NEAR(a.velocity,b.velocity,1e-14);
    EXPECT_NEAR(a.acceleration,b.acceleration,1e-14);
  }
}
TEST(IntegratedGroundMotion, CrossKnotTriangleHasKnownIntegral)
{
  const IntegratedGroundMotion motion({0,1,2},{0,1,0});
  const auto s=motion.at(2);
  EXPECT_NEAR(s.velocity,1,1e-14);
  EXPECT_NEAR(s.displacement,1,1e-14);
  EXPECT_NEAR(motion.at(1).displacement,1.0/6,1e-14);
  EXPECT_NEAR(motion.at(1).velocity,.5,1e-14);
}
TEST(IntegratedGroundMotion, QueryOrderAndRejectedTrialCannotChangeMotion)
{
  const IntegratedGroundMotion motion({0,.07,.13,.2},{-.004,.2,-.1,.05});
  const auto before=motion.at(.08);
  for(double t : {.2,.1,.09,.13,0.0,.2,.01}) (void)motion.at(t);
  const auto after=motion.at(.08);
  EXPECT_DOUBLE_EQ(before.displacement,after.displacement);
  EXPECT_DOUBLE_EQ(before.velocity,after.velocity);
  EXPECT_DOUBLE_EQ(before.acceleration,after.acceleration);
}
TEST(IntegratedGroundMotion, RejectsMalformedInputAndOutsideDomain)
{
  EXPECT_THROW(IntegratedGroundMotion({0},{1}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,0},{1,2}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({1,0},{1,2}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,1},{1}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,1},{1,std::numeric_limits<double>::infinity()}),std::invalid_argument);
  const IntegratedGroundMotion motion({0,1},{1,2});
  EXPECT_THROW(motion.at(-.001),std::out_of_range);
  EXPECT_THROW(motion.at(1.001),std::out_of_range);
}
TEST(IntegratedGroundMotion, DistortedQuadTranslationNeedsMassCorrection)
{
  // Trapezoid J=(3-eta)/8. Analytic int Ni*J is [5,5,4,4]/12.
  const Quad4Hourglass::Nodes nodes={{{{0,0}},{{2,0}},{{1,1}},{{0,1}}}};
  const auto correction=Quad4Hourglass::consistentMassCorrection(nodes);
  const std::array<double,4> exact={{5.0/12,5.0/12,1.0/3,1.0/3}};
  for(unsigned i=0;i<4;++i)
  {
    double row=Quad4Hourglass::area(nodes)/4;
    for(double v:correction[i]) row+=v;
    EXPECT_NEAR(row,exact[i],1e-14);
    EXPECT_GT(std::abs(row-Quad4Hourglass::area(nodes)/4),.04);
  }
}
TEST(IntegratedGroundMotion, MassAndHHTDampingCoordinateTransformation)
{
  const Quad4Hourglass::Nodes nodes={{{{0,0}},{{2,0}},{{1,1}},{{0,1}}}};
  auto mass=Quad4Hourglass::consistentMassCorrection(nodes);
  const double alpha=-.05,eta=.7,Ag=.2,Vg=.13,Vg_old=.09;
  const std::array<double,4> aq={{0,.2,-.3,.1}},vq={{0,.8,.3,-.2}},vold={{0,.7,.1,-.1}};
  for(unsigned i=0;i<4;++i)
  {
    double absolute=0,relative=0,load=0;
    for(unsigned j=0;j<4;++j)
    {
      mass[i][j]+=Quad4Hourglass::area(nodes)/16;
      absolute+=mass[i][j]*(aq[j]+Ag+eta*((1+alpha)*(vq[j]+Vg)-alpha*(vold[j]+Vg_old)));
      relative+=mass[i][j]*(aq[j]+eta*((1+alpha)*vq[j]-alpha*vold[j]));
      load+=mass[i][j]*(Ag+eta*((1+alpha)*Vg-alpha*Vg_old));
    }
    EXPECT_NEAR(absolute,relative+load,1e-14);
  }
}
TEST(IntegratedGroundMotion, HourglassDoesNotResistRigidGroundTranslation)
{
  const Quad4Hourglass::Nodes nodes={{{{0,0}},{{2,0}},{{1,1}},{{0,1}}}};
  const auto mode=Quad4Hourglass::mode(nodes);
  double amplitude=0;
  for(double v:mode) amplitude+=v*1.817e-3;
  EXPECT_NEAR(amplitude,0,1e-16);
}
