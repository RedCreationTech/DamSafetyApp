#include "DS04GFalseRotationTangent.h"
#include "DS04GSP1.h"
#include "AbaqusCDPSubstepIntegrator.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <fstream>

namespace
{
using Six = std::array<double, 6>;
using State = AbaqusCDPStateIntegrator::State;
constexpr double E = 2.97915e10, nu = .2;
RankTwoTensor tensor(const Six &v)
{
  RankTwoTensor t;
  const unsigned a[] = {0,1,2,0,1,0}, b[] = {0,1,2,1,2,2};
  for (unsigned k = 0; k < 6; ++k) t(a[k],b[k]) = t(b[k],a[k]) = v[k];
  return t;
}
Six six(const RankTwoTensor &t) { return {{t(0,0),t(1,1),t(2,2),t(0,1),t(1,2),t(0,2)}}; }
RankFourTensor tangent(const AbaqusCDPSubstepIntegrator::TangentMatrix &m)
{
  RankFourTensor D;
  const unsigned a[] = {0,1,2,0,1,0}, b[] = {0,1,2,1,2,2};
  for (unsigned c = 0; c < 6; ++c)
    for (unsigned r = 0; r < 6; ++r)
    {
      const Real v = m[c][r] * (a[c] == b[c] ? 1 : .5);
      D(a[r],b[r],a[c],b[c]) = D(b[r],a[r],a[c],b[c]) =
          D(a[r],b[r],b[c],a[c]) = D(b[r],a[r],b[c],a[c]) = v;
    }
  return D;
}
RankTwoTensor FhatFromRaw(const RankTwoTensor &t)
{
  const Real theta = .5 * std::atan2(2*t(0,1),t(0,0)-t(1,1));
  const Real m = (t(0,0)+t(1,1))/2, d = std::hypot((t(0,0)-t(1,1))/2,t(0,1));
  const Real a = 1/std::sqrt(2-std::sqrt(1+4*(m+d)));
  const Real b = 1/std::sqrt(2-std::sqrt(1+4*(m-d)));
  const Real c = std::cos(theta), s = std::sin(theta);
  RankTwoTensor f(RankTwoTensor::initIdentity);
  f(0,0)=a*c*c+b*s*s;f(1,1)=a*s*s+b*c*c;f(0,1)=f(1,0)=(a-b)*c*s;
  return f;
}
DS04GSP1::Nodes nodes(bool distorted)
{
  return distorted ? DS04GSP1::Nodes{{{{41.0744375,30.2808984}},{{40.5477539,31.1235957}},
                                    {{39.777457,30.5690957}},{{40.3079492,29.7433555}}}}
                   : DS04GSP1::Nodes{{{{0,0}},{{1,0}},{{1,1}},{{0,1}}}};
}
RankTwoTensor Ffrom(const DS04GSP1::Nodes &u, const DS04GSP1::Geometry &g)
{
  RankTwoTensor F(RankTwoTensor::initIdentity);
  for (unsigned n=0;n<4;++n)
    for (unsigned c=0;c<2;++c)
      for (unsigned j=0;j<2;++j) F(c,j)+=u[n][c]*g.gradient[n][j];
  return F;
}
template <class A> void read(std::istream &in,A &a) { for (auto &v:a)in>>v; }
}

TEST(DS04GFalseRotationTangent, IndependentTaylorDerivativeAtFiniteOldAndRigidPose)
{
  for (Real theta : {0., .3, -.4})
  {
    RankTwoTensor F(RankTwoTensor::initIdentity), old(RankTwoTensor::initIdentity);
    F(0,0)=std::cos(theta);F(0,1)=-std::sin(theta);
    F(1,0)=std::sin(theta);F(1,1)=std::cos(theta);
    old(0,0)=1.002;old(1,1)=.998;old(0,1)=.001;
    const auto hat=F*old.inverse();
    for (unsigned c=0;c<2;++c)
      for (unsigned j=0;j<2;++j)
      {
        RankTwoTensor L;L(c,j)=1;
        const auto d=L*F;
        const auto analytic=DS04GFalseRotation::incrementDerivative(hat,L*hat);
        for (Real h : {1e-5,1e-6,1e-7})
        {
          const auto fd=(DS04GFalseRotation::increment((F+h*d)*old.inverse())-
                         DS04GFalseRotation::increment((F-h*d)*old.inverse()))/(2*h);
          EXPECT_LT((fd-analytic).L2norm()/std::max(analytic.L2norm(),1e-12),1e-6);
        }
      }
  }
}

TEST(DS04GFalseRotationTangent, SixRealMaterialInputsCompleteForceChainTwoGeometries)
{
  const std::string path="unit/data/ds04g-sp1/";
  CDPMaterialTable tables(path+"compression_hardening.csv",path+"compression_damage.csv",
                          path+"tension_stiffening.csv",path+"tension_damage.csv",E);
  AbaqusCDPLocalIntegrator local(tables,{E,nu,30.,.1,1.16,.6667});
  AbaqusCDPStateIntegrator core(local,{0,1,.0005});
  AbaqusCDPSubstepIntegrator sub(core,{256,2.5e-5,1e-8});
  for (const auto *id : {"P01","P02","P03","P04","P05","P06"})
  {
    std::ifstream in(path+id+".txt");std::string name;Six s0,s1,expected_stress;
    std::array<double,16> state;double dt;
    in>>name;read(in,s0);read(in,s1);read(in,state);in>>dt;read(in,expected_stress);ASSERT_TRUE(bool(in));
    State old;
    std::copy_n(state.begin(),6,old.backbone.plastic_strain.begin());
    old.backbone.tensile_equivalent_plastic_strain=state[6];
    old.backbone.compressive_equivalent_plastic_strain=state[7];
    std::copy_n(state.begin()+8,6,old.viscous_plastic_strain.begin());
    old.viscous_tension_damage=state[14];old.viscous_compression_damage=state[15];
    const auto baseline=sub.integrateLinearized(s0,s1,dt,old);
    const auto D=tangent(baseline.algorithmic_tangent);
    Six raw;for (unsigned k=0;k<6;++k)raw[k]=s1[k]-s0[k];
    const auto hat=FhatFromRaw(tensor(raw));
    EXPECT_LT((DS04GFalseRotation::increment(hat)-tensor(raw)).L2norm(),1e-14);
    RankTwoTensor Fold(RankTwoTensor::initIdentity);Fold(0,0)=1.002;Fold(1,1)=.998;Fold(0,1)=.001;
    const auto target=hat*Fold;
    const auto weights=DS04GSP1::history(-.05,.003604,dt);
    // Real constitutive inputs embedded into synthetic finite old geometry.
    // Old/older stresses below are fixed test loads, not reconstructed dam history.
    const RankTwoTensor sold=tensor(expected_stress), solder=.8*sold;
    for (bool distorted : {false,true})
    {
      auto X=nodes(distorted);
      std::array<double,2> center{};
      for (const auto &n:X)for(unsigned c=0;c<2;++c)center[c]+=n[c]/4;
      for (auto &n:X)for(unsigned c=0;c<2;++c)n[c]-=center[c];
      const auto g0=DS04GSP1::geometry(X);DS04GSP1::Nodes u{};
      for (unsigned n=0;n<4;++n)
        for (unsigned c=0;c<2;++c)
          for (unsigned j=0;j<2;++j)u[n][c]+=(target(c,j)-(c==j))*X[n][j];
      auto evaluate=[&](const DS04GSP1::Nodes &v)
      {
        auto x=X;for(unsigned n=0;n<4;++n)for(unsigned c=0;c<2;++c)x[n][c]+=v[n][c];
        const auto g=DS04GSP1::geometry(x);const auto F=Ffrom(v,g0);
        Six strain=six(DS04GFalseRotation::increment(F*Fold.inverse()));
        for(unsigned k=0;k<6;++k)strain[k]+=s0[k];
        const auto r=sub.integrateLinearized(s0,strain,dt,old);
        EXPECT_EQ(r.result.accepted_substeps,baseline.result.accepted_substeps);
        EXPECT_EQ(r.result.final_result.backbone.active_branch,
                  baseline.result.final_result.backbone.active_branch);
        const auto sigma=weights.current*tensor(r.result.final_result.cauchy_stress)+
                           weights.old*sold+weights.older*solder;
        std::array<double,8> force{};
        for(unsigned n=0;n<4;++n)for(unsigned c=0;c<2;++c)
          for(unsigned b=0;b<2;++b)force[2*n+c]+=g.area*sigma(c,b)*g.gradient[n][b];
        return force;
      };
      auto x=X;for(unsigned n=0;n<4;++n)for(unsigned c=0;c<2;++c)x[n][c]+=u[n][c];
      const auto g=DS04GSP1::geometry(x);const auto F=Ffrom(u,g0);
      const auto effective=weights.current*tensor(baseline.result.final_result.cauchy_stress)+
                             weights.old*sold+weights.older*solder;
      unsigned passed=0;
      for (Real h : {1e-8,3e-9,1e-9})
      {
        double delta2=0,ref2=0;
        for(unsigned j=0;j<8;++j)
        {
          auto p=u,m=u;p[j/2][j%2]+=h;m[j/2][j%2]-=h;
          const auto fp=evaluate(p),fm=evaluate(m);
          for(unsigned i=0;i<8;++i)
          {
            const auto parts=DS04GFalseRotation::contract(F,Fold,D,effective,weights.current,
              i%2,j%2,RealGradient(g.gradient[i/2][0],g.gradient[i/2][1],0),
              RealGradient(g.gradient[j/2][0],g.gradient[j/2][1],0));
            const double a=g.area*(parts.material+parts.geometry),fd=(fp[i]-fm[i])/(2*h);
            delta2+=(a-fd)*(a-fd);ref2+=fd*fd;
          }
        }
        const double error=std::sqrt(delta2/std::max(ref2,1e-24));
        std::cout<<"DS04G_FULL_CHAIN "<<id<<" distorted="<<distorted<<" h="<<h<<" relative="<<error<<'\n';
        passed+=error<=1e-6;
      }
      EXPECT_GE(passed,2u)<<id<<" distorted="<<distorted;
    }
  }
}

TEST(DS04GFalseRotationTangent, ManufacturedN8N16N64MaterialAndMappingDerivatives)
{
  const std::string path="unit/data/ds04g-sp1/";
  CDPMaterialTable tables(path+"compression_hardening.csv",path+"compression_damage.csv",
                          path+"tension_stiffening.csv",path+"tension_damage.csv",E);
  AbaqusCDPLocalIntegrator local(tables,{E,nu,30.,.1,1.16,.6667});
  AbaqusCDPStateIntegrator core(local,{0,1,.0005});
  AbaqusCDPSubstepIntegrator sub(core,{256,2.5e-5,1e-8});
  for (unsigned mode=0;mode<3;++mode)
  {
    RankTwoTensor F(RankTwoTensor::initIdentity);
    if(mode==0){F(0,0)+=1e-4;F(0,1)=2e-4;F(1,1)-=1e-4;}
    if(mode==1){F(0,0)+=2e-4;F(1,1)-=4e-5;}
    if(mode==2){F(0,0)-=8e-4;F(1,1)+=1.6e-4;}
    const Six zero{};const State old;
    const auto base=sub.integrateLinearized(zero,six(DS04GFalseRotation::increment(F)),.002,old);
    const auto D=tangent(base.algorithmic_tangent);
    unsigned passed=0;
    for (Real h : {1e-7,3e-8,1e-8})
    {
      double delta2=0,ref2=0;
      for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)
      {
        RankTwoTensor L;L(i,j)=1;
        const auto plus=sub.integrateLinearized(zero,six(DS04GFalseRotation::increment(F+h*L*F)),.002,old);
        const auto minus=sub.integrateLinearized(zero,six(DS04GFalseRotation::increment(F-h*L*F)),.002,old);
        EXPECT_EQ(plus.result.accepted_substeps,base.result.accepted_substeps);
        EXPECT_EQ(minus.result.accepted_substeps,base.result.accepted_substeps);
        EXPECT_EQ(plus.result.final_result.backbone.active_branch,base.result.final_result.backbone.active_branch);
        EXPECT_EQ(minus.result.final_result.backbone.active_branch,base.result.final_result.backbone.active_branch);
        const auto analytic=D*DS04GFalseRotation::incrementDerivative(F,L*F);
        const auto fd=(tensor(plus.result.final_result.cauchy_stress)-
                       tensor(minus.result.final_result.cauchy_stress))/(2*h);
        delta2+=(analytic-fd).L2norm()*(analytic-fd).L2norm();ref2+=fd.L2norm()*fd.L2norm();
      }
      const double error=std::sqrt(delta2/std::max(ref2,1e-24));
      std::cout<<"DS04G_PARTITION_CHAIN mode="<<mode<<" N="<<base.result.accepted_substeps
               <<" h="<<h<<" relative="<<error<<'\n';
      passed+=error<=1e-6;
    }
    EXPECT_GE(passed,2u)<<"mode="<<mode;
  }
}
