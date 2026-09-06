#include "CDPBbarContraction.h"
#include "gtest/gtest.h"
#include <cmath>
struct TestC
{
  double x[3][3][3][3]={};
  double operator()(unsigned i,unsigned j,unsigned k,unsigned l) const { return x[i][j][k][l]; }
};
static TestC tensor(bool coupled)
{
  TestC c;
  const double p[3][3]={{2,.3,-.2},{.3,1,.6},{-.2,.6,3}};
  const double q[3][3]={{1,-.4,.1},{-.4,4,.2},{.1,.2,2}};
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
    for(unsigned k=0;k<3;++k) for(unsigned l=0;l<3;++l)
      c.x[i][j][k][l]=3*(i==j)*(k==l)+2*((i==k)*(j==l)+(i==l)*(j==k))
                       +(coupled ? p[i][j]*q[k][l]:0);
  return c;
}
static double residual(const TestC & c,unsigned a,unsigned b,double u)
{
  const double g[3]={.7,-.3,.8},h[3]={-.4,.9,.2};
  const double ag[3]={.1,.2,-.1},ah[3]={.3,-.2,.4};
  double stress[3][3]={};
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
    for(unsigned k=0;k<3;++k) for(unsigned l=0;l<3;++l)
    {
      const double eps=u*(.5*((k==b ? h[l]:0)+(l==b ? h[k]:0))
                            +(k==l ? (ah[b]-h[b])/3:0));
      stress[i][j]+=c(i,j,k,l)*eps;
    }
  double r=0,trace=0;
  for(unsigned j=0;j<3;++j) { r+=stress[a][j]*g[j];trace+=stress[j][j]; }
  return r+trace*(ag[a]-g[a])/3;
}
static double legacy(const TestC & c,unsigned a,unsigned b)
{
  const double g[3]={.7,-.3,.8},h[3]={-.4,.9,.2};
  const double ag[3]={.1,.2,-.1},ah[3]={.3,-.2,.4};
  double v=0;
  for(unsigned j=0;j<3;++j) for(unsigned l=0;l<3;++l) v+=g[j]*c(a,j,b,l)*h[l];
  for(unsigned k=0;k<3;++k)
  {
    for(unsigned i=0;i<3;++i) v+=c(i,i,k,k)*(ag[a]-g[a])*(ah[b]-h[b])/9;
    v+=c(a,a,k,k)*g[a]*(ah[b]-h[b])/3;
  }
  for(unsigned i=0;i<3;++i) for(unsigned l=0;l<3;++l)
    v+=c(i,i,b,l)*h[l]*(ag[a]-g[a])/3*(a==b && l!=b ? 2:1);
  return v;
}
TEST(CDPBbar, FullContractionMatchesResidualFiniteDifference)
{
  const auto c=tensor(true);
  for(unsigned a=0;a<3;++a) for(unsigned b=0;b<3;++b)
  {
    const double fd=(residual(c,a,b,.2+1e-6)-residual(c,a,b,.2-1e-6))/2e-6;
    EXPECT_NEAR(CDPBbar::contract(c,a,b,{.7,-.3,.8},{-.4,.9,.2},
                  std::array<double,3>{.1,.2,-.1}[a],std::array<double,3>{.3,-.2,.4}[b]),fd,1e-8);
  }
}
TEST(CDPBbar, IsotropicElasticityPreservesLegacy)
{
  const auto c=tensor(false);
  for(unsigned a=0;a<3;++a) for(unsigned b=0;b<3;++b)
    EXPECT_NEAR(CDPBbar::contract(c,a,b,{.7,-.3,.8},{-.4,.9,.2},
                  std::array<double,3>{.1,.2,-.1}[a],std::array<double,3>{.3,-.2,.4}[b]),legacy(c,a,b),1e-12);
}
TEST(CDPBbar, LegacyMissesCoupledResidualDerivative)
{
  const auto c=tensor(true);
  for(unsigned a=0;a<3;++a) for(unsigned b=0;b<3;++b)
    EXPECT_GT(std::abs(legacy(c,a,b)-(residual(c,a,b,1)-residual(c,a,b,0))),1e-4);
}
