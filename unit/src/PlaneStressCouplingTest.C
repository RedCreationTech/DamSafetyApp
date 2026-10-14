#include "PlaneStressCoupling.h"
#include "RankFourTensor.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(PlaneStressCoupling, FrozenTensorResidualActualADMatchesBothCrossBlocks)
{
  RankFourTensor C;
  C(0,0,2,2) = 3.; C(1,1,2,2) = 7.;
  C(2,2,0,0) = 13.; C(2,2,1,1) = 19.;
  C(0,1,2,2) = C(1,0,2,2) = 23.;
  C(2,2,0,1) = C(2,2,1,0) = -11.; // Deliberately not major symmetric.
  const RealGradient gt(.37, -.21, 0.), gp(.43, .29, 0.);
  const double phi=.25, test=.6;
  ADReal ux=0., uy=0., z=0.;
  Moose::derivInsert(ux.derivatives(),0,1.);
  Moose::derivInsert(uy.derivatives(),1,1.);
  Moose::derivInsert(z.derivatives(),2,1.);
  ADReal e[3][3] = {};
  e[0][0]=ux*gp(0); e[1][1]=uy*gp(1); e[2][2]=z*phi;
  e[0][1]=e[1][0]=.5*(ux*gp(1)+uy*gp(0));
  ADReal stress[3][3] = {};
  for (unsigned int i=0;i<3;++i)
    for (unsigned int j=0;j<3;++j)
      for (unsigned int k=0;k<3;++k)
        for (unsigned int l=0;l<3;++l)
          stress[i][j] += C(i,j,k,l)*e[k][l];
  for (double alpha : {0., -.05})
    for (unsigned int c=0;c<2;++c)
    {
      const ADReal ru=(1+alpha)*(stress[c][0]*gt(0)+stress[c][1]*gt(1));
      const ADReal rz=stress[2][2]*test;
      const double uz0=(1+alpha)*C(c,c,2,2)*gt(c)*phi;
      const double zu0=C(2,2,c,c)*test*gp(c);
      EXPECT_NEAR(ru.derivatives()[2], uz0+PlaneStressCoupling::displacementShearColumn(C,c,gt,phi,alpha),1e-14);
      EXPECT_NEAR(rz.derivatives()[c], zu0+PlaneStressCoupling::thicknessShearRow(C,c,test,gp),1e-14);
      EXPECT_GT(std::abs(ru.derivatives()[2]-uz0), .1);
      EXPECT_GT(std::abs(rz.derivatives()[c]-zu0), .1);
    }
}

TEST(PlaneStressCoupling, IsotropicElasticityNeedsNoShearCompletion)
{
  RankFourTensor C;
  for (unsigned int i=0;i<3;++i)
    for (unsigned int j=0;j<3;++j)
      for (unsigned int k=0;k<3;++k)
        for (unsigned int l=0;l<3;++l)
          C(i,j,k,l)=17.*(i==j)*(k==l)+13.*((i==k&&j==l)+(i==l&&j==k));
  for (unsigned int c=0;c<2;++c)
  {
    EXPECT_EQ(PlaneStressCoupling::displacementShearColumn(C,c,RealGradient(.3,.7,0.),.25,-.05),0.);
    EXPECT_EQ(PlaneStressCoupling::thicknessShearRow(C,c,.25,RealGradient(.3,.7,0.)),0.);
  }
}
