#include "AbaqusCDPSubstepIntegrator.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <iostream>

namespace
{
using Core = AbaqusCDPStateIntegrator;
using State = Core::State;
using Six = Core::SymmetricTensor;
using Packed = std::array<double, Core::state_size>;
const std::string data = "unit/data/ds04g-sp1/";
constexpr double E = 2.97915e10;
CDPMaterialTable tables()
{
  return CDPMaterialTable(data+"compression_hardening.csv", data+"compression_damage.csv",
                         data+"tension_stiffening.csv", data+"tension_damage.csv", E);
}
Packed pack(const State & s)
{
  Packed a{};
  std::copy(s.backbone.plastic_strain.begin(),s.backbone.plastic_strain.end(),a.begin());
  a[6]=s.backbone.tensile_equivalent_plastic_strain;
  a[7]=s.backbone.compressive_equivalent_plastic_strain;
  std::copy(s.viscous_plastic_strain.begin(),s.viscous_plastic_strain.end(),a.begin()+8);
  a[14]=s.viscous_tension_damage;a[15]=s.viscous_compression_damage;
  a[16]=s.viscous_combined_damage;return a;
}
State unpack(const Packed & a)
{
  State s;
  std::copy_n(a.begin(),6,s.backbone.plastic_strain.begin());
  s.backbone.tensile_equivalent_plastic_strain=a[6];
  s.backbone.compressive_equivalent_plastic_strain=a[7];
  std::copy_n(a.begin()+8,6,s.viscous_plastic_strain.begin());
  s.viscous_tension_damage=a[14];s.viscous_compression_damage=a[15];
  s.viscous_combined_damage=a[16];return s;
}
template <class A> void read(std::istream & in,A &a){for(auto &v:a)in>>v;}
double oldCombination(const AbaqusCDPLocalIntegrator &local,const Six &strain,const State &old)
{
  Six elastic{};for(unsigned i=0;i<6;++i)elastic[i]=strain[i]-old.viscous_plastic_strain[i];
  return AbaqusCDPFormula::combineDamage(old.viscous_compression_damage,old.viscous_tension_damage,
      0,1,AbaqusCDPFormula::stressInvariants(local.elasticStress(elastic)).tension_weight).damage;
}
}

TEST(CombinedDamageHistory, SixRealInputsKeepBackboneBranchesAndCompleteTangent)
{
  const auto table=tables();AbaqusCDPLocalIntegrator local(table,{E,.2,30.,.1,1.16,.6667});
  Core off(local,{0,1,.0005});Core on(local,{0,1,.0005,1e-12,true});
  for(const auto *id:{"P01","P02","P03","P04","P05","P06"})
  {
    std::ifstream in(data+id+".txt");std::string name;Six s0,s1;std::array<double,16> old16;double dt;
    in>>name;read(in,s0);read(in,s1);read(in,old16);in>>dt;ASSERT_TRUE(bool(in));
    Packed old_values{};std::copy(old16.begin(),old16.end(),old_values.begin());
    auto old=unpack(old_values);old.viscous_combined_damage=oldCombination(local,s0,old);
    const auto snapshot=pack(old);
    for(unsigned n:{1u,2u})
    {
      double delta=0.;for(unsigned i=0;i<6;++i)delta=std::max(delta,std::abs(s1[i]-s0[i]));
      const double limit=n==1?0.:delta/1.5;
      AbaqusCDPSubstepIntegrator a(off,{256,limit,1e-8}),b(on,{256,limit,1e-8});
      const auto base=b.integrateLinearized(s0,s1,dt,old),legacy=a.integrateLinearized(s0,s1,dt,old);
      const auto p=pack(base.result.final_result.state),l=pack(legacy.result.final_result.state);
      for(unsigned i=0;i<16;++i)EXPECT_DOUBLE_EQ(p[i],l[i]);
      EXPECT_EQ(base.result.accepted_substeps,legacy.result.accepted_substeps);
      unsigned passed=0;
      for(double h:{1e-8,3e-9,1e-9})
      {
        double d2=0.,r2=0.;
        for(unsigned j=0;j<6;++j)
        {
          auto plus=s1,minus=s1;plus[j]+=h;minus[j]-=h;
          const auto pp=b.integrate(s0,plus,dt,old),mm=b.integrate(s0,minus,dt,old);
          EXPECT_EQ(pp.accepted_substeps,base.result.accepted_substeps);
          EXPECT_EQ(mm.accepted_substeps,base.result.accepted_substeps);
          for(unsigned i=0;i<6;++i)
          {
            const double fd=(pp.final_result.cauchy_stress[i]-mm.final_result.cauchy_stress[i])/(2*h);
            const double d=base.algorithmic_tangent[j][i]-fd;d2+=d*d;r2+=fd*fd;
          }
        }
        const double rel=std::sqrt(d2/std::max(r2,1e-24));
        std::cout<<"SCALAR_REAL_FD "<<id<<" N="<<base.result.accepted_substeps<<" h="<<h<<" rel="<<rel<<'\n';
        passed+=rel<1e-6;
      }
      EXPECT_GE(passed,2u)<<id<<" N="<<n;
      EXPECT_EQ(pack(old),snapshot);
    }
  }
}

TEST(CombinedDamageHistory, Full23By23TransitionIncludesOldScalarHistory)
{
  const auto table=tables();AbaqusCDPLocalIntegrator local(table,{E,.2,30.,.1,1.16,.6667});
  Core core(local,{0,1,.0005,1e-12,true});
  State old;old.backbone.tensile_equivalent_plastic_strain=.0002134;
  old.backbone.compressive_equivalent_plastic_strain=.0001345;
  old.viscous_tension_damage=.1;old.viscous_compression_damage=.001;
  old.viscous_combined_damage=.4;
  const Six strain={1e-5,-2e-5,-1e-5,3e-6,0.,0.};
  const auto base=core.integrateLinearized(strain,.002,old);
  ASSERT_FALSE(base.result.backbone.plastic);
  for(unsigned input=0;input<Core::transition_size;++input)
  {
    const double h=input>=20?1e-6:1e-9;
    auto sp=strain,sm=strain;auto op=pack(old),om=op;
    if(input<6){sp[input]+=h;sm[input]-=h;}
    else{op[input-6]+=h;om[input-6]-=h;}
    const auto rp=core.integrate(sp,.002,unpack(op)),rm=core.integrate(sm,.002,unpack(om));
    auto p=pack(rp.state),m=pack(rm.state);
    for(unsigned output=0;output<Core::transition_size;++output)
    {
      const double fd=output<6?(rp.cauchy_stress[output]-rm.cauchy_stress[output])/(2*h):
                                  (p[output-6]-m[output-6])/(2*h);
      const double analytic=base.derivative[input][output];
      const double scale=output<6?E:1.;
      EXPECT_NEAR(analytic,fd,1e-6*std::max({1.,std::abs(analytic),std::abs(fd),scale}))
          <<"input="<<input<<" output="<<output;
    }
  }
  EXPECT_NEAR(base.derivative[22][22],.2,1e-15);
}

TEST(CombinedDamageHistory, RecoveryDecreasesScalarAndKeepsBranchIrreversibility)
{
  const auto table=tables();AbaqusCDPLocalIntegrator local(table,{E,.2,30.,.1,1.16,.6667});
  Core on(local,{0,1,.0005,1e-12,true});State old;
  old.backbone.tensile_equivalent_plastic_strain=.0002134;
  old.viscous_tension_damage=local.materialResponse(CDPMaterialTable::Branch::TENSION,.0002134).damage.value;
  old.viscous_combined_damage=old.viscous_tension_damage;
  const auto snap=pack(old);const Six compression={-1e-6,-1e-6,-1e-6,0,0,0};
  const auto closed=on.integrate(compression,.002,old);
  ASSERT_FALSE(closed.backbone.plastic);
  EXPECT_NEAR(closed.damage.damage,.2*old.viscous_combined_damage,1e-14);
  EXPECT_DOUBLE_EQ(closed.state.viscous_tension_damage,old.viscous_tension_damage);
  EXPECT_EQ(pack(old),snap);
  const auto hold=on.integrate(compression,0.,old);
  EXPECT_DOUBLE_EQ(hold.state.viscous_combined_damage,old.viscous_combined_damage);
  Core mu0(local,{0,1,0.,1e-12,true});
  EXPECT_DOUBLE_EQ(mu0.integrate(compression,.002,old).damage.damage,0.);
  old.viscous_combined_damage=std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(on.integrate(compression,.002,old),std::runtime_error);
  EXPECT_THROW(on.integrateLinearized(compression,.002,old),std::runtime_error);
}

TEST(CombinedDamageHistory, SubstepChainFailureRollbackAndCleanReplay)
{
  const auto table=tables();AbaqusCDPLocalIntegrator local(table,{E,.2,30.,.1,1.16,.6667});
  Core core(local,{0,1,.0005,1e-12,true});AbaqusCDPSubstepIntegrator sub(core,{16,2.5e-5,1e-8});
  const Six end={.0002,-.00005,0,0,0,0};State old;const auto snapshot=pack(old);
  const auto normal=sub.integrate({},end,.002,old);
  ASSERT_EQ(normal.accepted_substeps,8u);
  State manual;Core::Result last{};
  for(unsigned i=1;i<=8;++i){Six strain{};for(unsigned j=0;j<6;++j)strain[j]=end[j]*i/8.;last=core.integrate(strain,.002/8,manual);manual=last.state;}
  EXPECT_EQ(pack(normal.final_result.state),pack(manual));
  EXPECT_EQ(normal.final_result.cauchy_stress,last.cauchy_stress);
  auto limited=AbaqusCDPLocalIntegrator::Parameters{E,.2,30.,.1,1.16,.6667};limited.maximum_iterations=0;
  AbaqusCDPLocalIntegrator failing(table,limited);Core failcore(failing,{0,1,.0005,1e-12,true});
  AbaqusCDPSubstepIntegrator fail(failcore,{4,0.,1e-8});
  EXPECT_THROW(fail.integrate({},end,.002,old),std::runtime_error);
  EXPECT_THROW(fail.integrateLinearized({},end,.002,old),std::runtime_error);
  EXPECT_EQ(pack(old),snapshot);
  const auto replay=sub.integrate({},end,.002,old);
  EXPECT_EQ(pack(replay.final_result.state),pack(normal.final_result.state));
}
