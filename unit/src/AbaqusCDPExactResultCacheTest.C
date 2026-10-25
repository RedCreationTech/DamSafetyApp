#include "AbaqusCDPExactResultCache.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <vector>

namespace
{
using Cache=AbaqusCDPExactResultCache;
std::uint64_t raw(double x){std::uint64_t b;std::memcpy(&b,&x,8);return b;}
std::vector<std::uint64_t> fingerprint(const Cache::Result & r)
{
  std::vector<std::uint64_t> values;
  const auto add=[&](double x){values.push_back(raw(x));};
  const auto tensor=[&](const auto & x){for(double v:x)add(v);};
  const auto & f=r.result.final_result;const auto & b=f.backbone;
  tensor(b.effective_stress);tensor(b.state.plastic_strain);
  add(b.state.tensile_equivalent_plastic_strain);add(b.state.compressive_equivalent_plastic_strain);
  for(double x:{double(b.active_branch),double(b.plastic),double(b.iterations),double(b.jacobian_fallbacks),
      double(b.automatic_jacobian_evaluations),double(b.finite_difference_jacobian_evaluations),
      double(b.local_factorizations),double(b.local_backsolves),b.residual_norm,b.trial_yield,b.final_yield,
      b.plastic_multiplier,b.backbone_tension_damage,b.backbone_compression_damage})add(x);
  tensor(f.viscous_effective_stress);tensor(f.cauchy_stress);
  tensor(f.state.backbone.plastic_strain);tensor(f.state.viscous_plastic_strain);
  for(double x:{f.state.backbone.tensile_equivalent_plastic_strain,
      f.state.backbone.compressive_equivalent_plastic_strain,f.state.viscous_tension_damage,
      f.state.viscous_compression_damage,f.damage.tension_weight,f.damage.tensile_recovery_factor,
      f.damage.compressive_recovery_factor,f.damage.stiffness_factor,f.damage.damage,
      double(f.rate_independent),f.dt_over_relaxation_time,f.plastic_strain_lag_norm,
      f.tension_damage_lag,f.compression_damage_lag})add(x);
  const auto & s=r.result;
  for(double x:{double(s.accepted_substeps),double(s.cutback_count),double(s.attempted_partitions),
      double(s.total_local_iterations),double(s.total_jacobian_fallbacks),
      double(s.total_automatic_jacobian_evaluations),double(s.total_finite_difference_jacobian_evaluations),
      double(s.total_local_factorizations),double(s.total_local_backsolves),double(s.proactively_partitioned)})add(x);
  for(const auto & column:r.algorithmic_tangent)tensor(column);
  return values;
}
Cache::Key key(){return Cache::key(8,2,{},{},.0025,{});}
}

TEST(AbaqusCDPExactResultCache, All29InputComponentsAndPointIdentityAreRequired)
{
  Cache c(4);Cache::Result result{};auto k=key();c.store(k,result);
  ASSERT_NE(c.find(k,true),nullptr);
  for(std::size_t i=0;i<29;++i)
  {
    auto changed=k;changed.input[i]=std::nextafter(changed.input[i],1.);
    EXPECT_EQ(c.find(changed,true),nullptr)<<i;
  }
  auto other=k;other.point.first++;EXPECT_EQ(c.find(other,true),nullptr);
  other=k;other.point.second++;EXPECT_EQ(c.find(other,true),nullptr);
  other=k;other.point.side=0;EXPECT_EQ(c.find(other,true),nullptr);
}
TEST(AbaqusCDPExactResultCache, SignedZeroNaNAndInfinityCannotAlias)
{
  Cache c(4);Cache::Result r{};auto k=key();c.store(k,r);
  auto neg=k;neg.input[0]=-0.;EXPECT_EQ(c.find(neg,true),nullptr);
  for(double v:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})
  {auto bad=k;bad.input[0]=v;c.store(bad,r);EXPECT_EQ(c.find(bad,true),nullptr);}
  EXPECT_EQ(c.entries(),1u);
}
TEST(AbaqusCDPExactResultCache, EpochInvalidationAndCapacityAreBounded)
{
  Cache c(1);Cache::Result r{};auto k=key();c.store(k,r);
  auto next=k;next.point.first++;c.store(next,r);EXPECT_EQ(c.entries(),1u);
  EXPECT_EQ(c.counters().capacity_skips,1u);EXPECT_EQ(c.find(next,false),nullptr);
  c.invalidate();EXPECT_EQ(c.find(k,true),nullptr);EXPECT_EQ(c.entries(),0u);
  c.store(next,r);EXPECT_NE(c.find(next,false),nullptr);EXPECT_EQ(c.counters().reused_hits,0u);
}
TEST(AbaqusCDPExactResultCache, CompleteIntegratedResultIsCopiedWithoutChangingOriginal)
{
  const std::string dir="test/tests/cdp_material_table/data/";
  CDPMaterialTable table(dir+"compression_hardening.csv",dir+"compression_damage.csv",
                         dir+"tension_stiffening.csv",dir+"tension_damage.csv",3.04e10);
  AbaqusCDPLocalIntegrator local(table,{3.04e10,.2,36.31,.1,1.16,.667,40,1e-9,1e-7,1e-6});
  AbaqusCDPStateIntegrator state(local,{1.,0.,5e-4,1e-12});
  AbaqusCDPSubstepIntegrator integrator(state,{16,0.,1e-8});
  const double tension=table.stressValues(CDPMaterialTable::Branch::TENSION).front();
  const double compression=table.stressValues(CDPMaterialTable::Branch::COMPRESSION).front();
  for(double stress:{1e6,1.05*tension,-1.05*compression})
  {
    AbaqusCDPSubstepIntegrator::SymmetricTensor strain={stress/3.04e10,-.2*stress/3.04e10,
                                                       -.2*stress/3.04e10,0.,0.,0.};
    const auto result=integrator.integrateLinearized({},strain,.001,{});
    Cache cache(4);const auto k=Cache::key(1,0,{},strain,.001,{});
    cache.store(k,result);const auto * hit=cache.find(k,true);ASSERT_NE(hit,nullptr);
    EXPECT_EQ(fingerprint(result),fingerprint(*hit));
    Cache separate_material(4);EXPECT_EQ(separate_material.find(k,true),nullptr);
  }
}
TEST(AbaqusCDPExactResultCache, PassiveTraceAndTraceLimitDoNotChangeSolution)
{
  const std::string dir="test/tests/cdp_material_table/data/";
  CDPMaterialTable table(dir+"compression_hardening.csv",dir+"compression_damage.csv",
                         dir+"tension_stiffening.csv",dir+"tension_damage.csv",3.04e10);
  AbaqusCDPLocalIntegrator::Parameters p{3.04e10,.2,36.31,.1,1.16,.667,40,1e-9,1e-7,1e-6};
  std::ostringstream log;AbaqusCDPLocalIntegrator::DiagnosticTrace trace;trace.stream=&log;trace.limit=1;
  AbaqusCDPLocalIntegrator base(table,p),observed(table,p,&trace);
  AbaqusCDPStateIntegrator bs(base,{1.,0.,5e-4,1e-12}),os(observed,{1.,0.,5e-4,1e-12});
  AbaqusCDPSubstepIntegrator bi(bs,{16,0.,1e-8}),oi(os,{16,0.,1e-8});
  const double stress=1.05*table.stressValues(CDPMaterialTable::Branch::TENSION).front();
  AbaqusCDPSubstepIntegrator::SymmetricTensor e={stress/3.04e10,-.2*stress/3.04e10,-.2*stress/3.04e10,0.,0.,0.};
  EXPECT_EQ(fingerprint(bi.integrateLinearized({},e,.001,{})),fingerprint(oi.integrateLinearized({},e,.001,{})));
  EXPECT_EQ(trace.written,1u);EXPECT_GT(trace.dropped,0u);EXPECT_FALSE(trace.io_failed);

  std::ostringstream broken;broken.setstate(std::ios::badbit);
  AbaqusCDPLocalIntegrator::DiagnosticTrace failed_io;failed_io.stream=&broken;
  AbaqusCDPLocalIntegrator io_local(table,p,&failed_io);
  AbaqusCDPStateIntegrator io_state(io_local,{1.,0.,5e-4,1e-12});
  AbaqusCDPSubstepIntegrator io_integrator(io_state,{16,0.,1e-8});
  EXPECT_EQ(fingerprint(bi.integrateLinearized({},e,.001,{})),
            fingerprint(io_integrator.integrateLinearized({},e,.001,{})));
  EXPECT_TRUE(failed_io.io_failed);
}

TEST(AbaqusCDPExactResultCache, TraceDoesNotChangeExceptionAndFailedResultIsNotStored)
{
  const std::string dir="test/tests/cdp_material_table/data/";
  CDPMaterialTable table(dir+"compression_hardening.csv",dir+"compression_damage.csv",
                         dir+"tension_stiffening.csv",dir+"tension_damage.csv",3.04e10);
  auto p=AbaqusCDPLocalIntegrator::Parameters{3.04e10,.2,36.31,.1,1.16,.667,0,1e-9,1e-7,1e-6};
  std::ostringstream log;AbaqusCDPLocalIntegrator::DiagnosticTrace trace;trace.stream=&log;
  std::string errors[2];Cache cache(4);
  const double stress=1.05*table.stressValues(CDPMaterialTable::Branch::TENSION).front();
  AbaqusCDPSubstepIntegrator::SymmetricTensor e={stress/3.04e10,-.2*stress/3.04e10,-.2*stress/3.04e10,0.,0.,0.};
  const auto k=Cache::key(1,0,{},e,.001,{});
  for(unsigned int i=0;i<2;++i)
  {
    AbaqusCDPLocalIntegrator local(table,p,i?&trace:nullptr);
    AbaqusCDPStateIntegrator state(local,{1.,0.,5e-4,1e-12});
    AbaqusCDPSubstepIntegrator integration(state,{4,0.,1e-8});
    try{cache.store(k,integration.integrateLinearized({},e,.001,{}));}
    catch(const std::exception & error){errors[i]=error.what();cache.invalidate();}
  }
  ASSERT_FALSE(errors[0].empty());EXPECT_EQ(errors[0],errors[1]);
  EXPECT_EQ(cache.entries(),0u);EXPECT_EQ(cache.find(k,true),nullptr);
}
