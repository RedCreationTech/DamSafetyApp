#include "AbaqusCDPLocalIntegrator.h"
#include "FrozenPrincipalFixture.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

namespace
{
using Tensor = AbaqusCDPFormula::SymmetricTensor;
using Matrix = std::array<std::array<long double,3>,3>;
struct Spectrum { std::array<double,3> values; Matrix vectors; std::array<Tensor,3> gradients; double backward,orthogonal; };

// TEST ONLY: a scaled symmetric Jacobi solver, independent of production trig/AD.
Spectrum stableSpectrum(const Tensor & stress)
{
  long double scale=0;
  for (auto x:stress) scale=std::max(scale,std::abs(static_cast<long double>(x)));
  if (scale==0) scale=1;
  Matrix a={{{stress[0]/scale,stress[3]/scale,stress[5]/scale},
             {stress[3]/scale,stress[1]/scale,stress[4]/scale},
             {stress[5]/scale,stress[4]/scale,stress[2]/scale}}};
  const auto original=a;
  Matrix q={{{1,0,0},{0,1,0},{0,0,1}}};
  for (unsigned iteration=0;iteration<80;++iteration)
  {
    unsigned p=0,r=1;
    for (unsigned i=0;i<3;++i)
      for (unsigned j=i+1;j<3;++j)
        if (std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
    if (std::abs(a[p][r])<1e-30L) break;
    const long double tau=(a[r][r]-a[p][p])/(2*a[p][r]);
    const long double t=std::copysign(1.L,tau)/(std::abs(tau)+std::hypot(1.L,tau));
    const long double c=1/std::sqrt(1+t*t),s=t*c;
    const long double off=a[p][r];
    a[p][p]-=t*off; a[r][r]+=t*off; a[p][r]=a[r][p]=0;
    for (unsigned k=0;k<3;++k)
    {
      if (k!=p && k!=r)
      {
        const long double ap=a[k][p],ar=a[k][r];
        a[k][p]=a[p][k]=c*ap-s*ar; a[k][r]=a[r][k]=s*ap+c*ar;
      }
      const long double qp=q[k][p],qr=q[k][r];
      q[k][p]=c*qp-s*qr; q[k][r]=s*qp+c*qr;
    }
  }
  std::array<unsigned,3> order={{0,1,2}};
  std::sort(order.begin(),order.end(),[&](unsigned i,unsigned j){return a[i][i]<a[j][j];});
  Spectrum result={}; long double backward=0,orthogonal=0;
  for(unsigned k=0;k<3;++k)
  {
    auto n=order[k]; result.values[k]=static_cast<double>(a[n][n]*scale);
    for(unsigned i=0;i<3;++i) result.vectors[i][k]=q[i][n];
    const long double x=q[0][n],y=q[1][n],z=q[2][n];
    result.gradients[k]={{static_cast<double>(x*x),static_cast<double>(y*y),static_cast<double>(z*z),
                         static_cast<double>(2*x*y),static_cast<double>(2*y*z),static_cast<double>(2*x*z)}};
    for(unsigned i=0;i<3;++i)
    {
      long double v=-a[n][n]*q[i][n];
      for(unsigned j=0;j<3;++j) v+=original[i][j]*q[j][n];
      backward=std::max(backward,std::abs(v));
    }
  }
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
  {
    long double v=-(i==j?1.L:0.L);
    for(unsigned k=0;k<3;++k) v+=result.vectors[k][i]*result.vectors[k][j];
    orthogonal=std::max(orthogonal,std::abs(v));
  }
  result.backward=static_cast<double>(backward); result.orthogonal=static_cast<double>(orthogonal);
  return result;
}
}

TEST(AbaqusCDPFrozenPrincipal, ReplaysActualSavedUnknownWithoutIntegration)
{
  const std::string dir="test/tests/cdp_material_table/data/";
  const CDPMaterialTable table(dir+"compression_hardening.csv",dir+"compression_damage.csv",
                               dir+"tension_stiffening.csv",dir+"tension_damage.csv",3.04e10);
  const AbaqusCDPLocalIntegrator::Parameters p={3.04e10,.2,36.31,.1,1.16,.667,40,1e-9,1e-7,1e-6};
  const AbaqusCDPLocalIntegrator integrator(table,p);
  unsigned residuals=0,matrices=0;
  for(const auto & f:FrozenPrincipalFixture::states)
  {
    SCOPED_TRACE(f.name);
    const AbaqusCDPLocalIntegrator::State old={f.plastic,f.kt,f.kc};
    const auto result=integrator.fixedStateDiagnostic(f.unknown,f.target,old);
    EXPECT_EQ(result.stress_scale,f.stress_scale); EXPECT_EQ(result.strain_scale,f.strain_scale);
    for(unsigned i=0;i<9;++i) {EXPECT_EQ(result.residual[i],f.residual[i]);++residuals;}
    if(f.mode)
    {
      const auto & j=f.mode==1?result.automatic_differentiation:result.finite_difference;
      for(unsigned i=0;i<9;++i) for(unsigned k=0;k<9;++k) EXPECT_EQ(j[i][k],f.jacobian[9*i+k]);
      ++matrices;
    }
    EXPECT_EQ(old.plastic_strain,f.plastic); EXPECT_EQ(old.tensile_equivalent_plastic_strain,f.kt);
    EXPECT_EQ(old.compressive_equivalent_plastic_strain,f.kc);
  }
  std::cout<<"REPLAY,"<<std::size(FrozenPrincipalFixture::states)<<","<<residuals<<","<<matrices<<"\n";
}

TEST(AbaqusCDPFrozenPrincipal, StableEigenpairsAndProjectorsMatchNinetyDigitReference)
{
  std::cout<<std::setprecision(17);
  double max_original_value=0,max_original_gradient=0,max_stable_gradient=0;
  for(const auto & f:FrozenPrincipalFixture::spectra)
  {
    SCOPED_TRACE(f.name);
    const auto stable=stableSpectrum(f.stress);
    const auto original=AbaqusCDPLocalIntegrator::principalStressDiagnostic(f.stress);
    const auto ordinary=AbaqusCDPFormula::stressInvariants(f.stress);
    double scale=0,ve=0,ge=0,ov=0,og=0;
    for(auto x:f.stress) scale=std::max(scale,std::abs(x));
    const double bound=128*std::numeric_limits<double>::epsilon()*std::max(scale,1e-300);
    EXPECT_LE(stable.backward,1e-12); EXPECT_LE(stable.orthogonal,1e-12);
    for(unsigned k=0;k<3;++k)
    {
      ve=std::max(ve,std::abs(stable.values[k]-f.values[k]));
      ov=std::max(ov,std::abs(original.values[k]-f.values[k]));
      EXPECT_LE(std::abs(stable.values[k]-f.values[k]),bound);
      // Plain and Dual value paths should agree; diagnose any independent discrepancy.
      EXPECT_DOUBLE_EQ(original.values[k],ordinary.principal_stress[k]);
      if(f.distinct) for(unsigned i=0;i<6;++i)
      {
        ge=std::max(ge,std::abs(stable.gradients[k][i]-f.gradients[6*k+i]));
        og=std::max(og,std::abs(original.gradient[k][i]-f.gradients[6*k+i]));
      }
    }
    if(f.distinct) EXPECT_LE(ge,1e-8);
    if(std::string(f.name).find("actual-failed")==0)
    {max_original_value=std::max(max_original_value,ov);max_original_gradient=std::max(max_original_gradient,og);}
    max_stable_gradient=std::max(max_stable_gradient,ge);
    std::cout<<"SPECTRUM,"<<f.name<<","<<ve<<","<<ge<<","<<ov<<","<<og<<","<<stable.backward<<","<<stable.orthogonal<<"\n";
  }
  EXPECT_GT(max_original_value,.01); EXPECT_GT(max_original_gradient,1e-3);
  std::cout<<"MAXIMUM,"<<max_original_value<<","<<max_original_gradient<<","<<max_stable_gradient<<"\n";
}

TEST(AbaqusCDPFrozenPrincipal, OrdinaryDirectionalDifferenceConfirmsShearConvention)
{
  const Tensor stress={{3,-4,8,2,-1,4}},direction={{.2,-.3,.7,.4,-.6,.1}};
  const double h=1e-5; const auto center=stableSpectrum(stress);
  Tensor plus=stress,minus=stress;
  for(unsigned i=0;i<6;++i){plus[i]+=h*direction[i];minus[i]-=h*direction[i];}
  const auto a=stableSpectrum(plus),b=stableSpectrum(minus);
  for(unsigned k=0;k<3;++k)
  {
    double projected=0;for(unsigned i=0;i<6;++i)projected+=center.gradients[k][i]*direction[i];
    EXPECT_NEAR((a.values[k]-b.values[k])/(2*h),projected,2e-8);
  }
}

TEST(AbaqusCDPFrozenPrincipal, ExactRepeatedRootHasDifferentOneSidedDerivatives)
{
  const Tensor stress={{-14999999,-14999999,2,-15000001,0,0}};
  Tensor plus=stress,minus=stress;const double h=.01;plus[2]+=h;minus[2]-=h;
  const auto center=stableSpectrum(stress),a=stableSpectrum(plus),b=stableSpectrum(minus);
  EXPECT_NEAR(center.values[0],-3e7,1e-8);EXPECT_NEAR(center.values[1],2,1e-8);EXPECT_NEAR(center.values[2],2,1e-8);
  EXPECT_NEAR((a.values[2]-center.values[2])/h,1,1e-6);
  EXPECT_NEAR((center.values[2]-b.values[2])/h,0,1e-6);
  EXPECT_NEAR((a.values[1]-center.values[1])/h,0,1e-6);
  EXPECT_NEAR((center.values[1]-b.values[1])/h,1,1e-6);
  // The repeated eigenspace projector is unique even when its individual vectors are not.
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
  {
    const auto v=center.vectors;
    const long double block=v[i][1]*v[j][1]+v[i][2]*v[j][2];
    const long double expected=i==2&&j==2?1:(i<2&&j<2?(i==j?.5L:-.5L):0.L);
    EXPECT_NEAR(static_cast<double>(block),static_cast<double>(expected),1e-12);
  }
}
