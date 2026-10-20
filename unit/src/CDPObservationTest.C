#include "CDPObservation.h"
#include "AbaqusCDPSubstepIntegrator.h"
#include "gtest/gtest.h"
#include <stdexcept>
namespace
{
struct Observer : CDPObservation::Sink
{
  unsigned int records=0, matrices=0, failures=0, candidates=0;
  bool enterPoint(std::uint64_t e,unsigned int qp,bool) override { return (e==1881 || e==767) && qp<4; }
  void write(const char * stream,const char * event,const char *,const std::vector<double> & values,const std::string &) noexcept override
  {
    ++records;
    if(std::string(stream)=="local-jacobian" && values.size()==81) ++matrices;
    if(std::string(event)=="exception" || std::string(event)=="partition_failed") ++failures;
    if(std::string(event)=="candidate") ++candidates;
  }
};
struct Install
{
  Install(CDPObservation::Sink & s) { CDPObservation::configured=&s; }
  ~Install() { CDPObservation::configured=nullptr; }
};
CDPMaterialTable table()
{
  const std::string d="test/tests/cdp_material_table/data/";
  return CDPMaterialTable(d+"compression_hardening.csv",d+"compression_damage.csv",d+"tension_stiffening.csv",d+"tension_damage.csv",3.04e10);
}
}
TEST(CDPObservation, ScopedSelectionRestoresDisabledState)
{
  Observer o; Install installed(o);
  { CDPObservation::PointScope p(9999,0); EXPECT_FALSE(CDPObservation::active()); }
  { CDPObservation::PointScope p(1881,3); EXPECT_TRUE(CDPObservation::active()); }
  EXPECT_FALSE(CDPObservation::active());
}
TEST(CDPObservation, EightPointLabelsDoNotAlterMaterialMapOrTangent)
{
  const auto t=table();
  AbaqusCDPLocalIntegrator local(t,{3.04e10,.2,36.31,.1,1.16,.667});
  AbaqusCDPStateIntegrator state(local,{1.,0.,.0005});
  AbaqusCDPSubstepIntegrator integrator(state,{256,2.5e-5,1e-8});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor strain={1.5e-4,-3e-5,-3e-5,2e-5,0.,0.};
  const auto baseline=integrator.integrateLinearized({},strain,.002,{});
  Observer o; Install installed(o);
  for(const auto label : {1881,767})
  for(unsigned int qp=0;qp<4;++qp)
  {
    CDPObservation::PointScope p(label,qp);
    const auto observed=integrator.integrateLinearized({},strain,.002,{});
    EXPECT_EQ(observed.algorithmic_tangent,baseline.algorithmic_tangent);
    EXPECT_EQ(observed.result.final_result.cauchy_stress,baseline.result.final_result.cauchy_stress);
    EXPECT_EQ(observed.result.final_result.state.backbone.plastic_strain,baseline.result.final_result.state.backbone.plastic_strain);
    EXPECT_EQ(observed.result.final_result.state.viscous_plastic_strain,baseline.result.final_result.state.viscous_plastic_strain);
    EXPECT_EQ(observed.result.final_result.state.backbone.tensile_equivalent_plastic_strain,baseline.result.final_result.state.backbone.tensile_equivalent_plastic_strain);
    EXPECT_EQ(observed.result.final_result.state.viscous_tension_damage,baseline.result.final_result.state.viscous_tension_damage);
    EXPECT_EQ(observed.result.total_local_iterations,baseline.result.total_local_iterations);
    EXPECT_EQ(observed.result.total_automatic_jacobian_evaluations,baseline.result.total_automatic_jacobian_evaluations);
    EXPECT_EQ(observed.result.total_local_factorizations,baseline.result.total_local_factorizations);
    EXPECT_EQ(observed.result.accepted_substeps,baseline.result.accepted_substeps);
  }
  EXPECT_GT(o.matrices,0u); EXPECT_GT(o.candidates,0u);
}
TEST(CDPObservation, FailureRecordsSurviveWithoutChangingOldState)
{
  const auto t=table();
  AbaqusCDPLocalIntegrator::Parameters parameters={3.04e10,.2,36.31,.1,1.16,.667};
  parameters.maximum_iterations=0;
  AbaqusCDPLocalIntegrator local(t,parameters);
  AbaqusCDPStateIntegrator state(local,{1.,0.,.0005});
  AbaqusCDPSubstepIntegrator integrator(state,{4,0.,1e-8});
  AbaqusCDPStateIntegrator::State old;
  Observer o; Install installed(o); CDPObservation::PointScope p(1881,0);
  EXPECT_THROW(integrator.integrateLinearized({},{1e-3,0.,0.,0.,0.,0.},.002,old),std::runtime_error);
  EXPECT_GT(o.failures,0u);
  EXPECT_EQ(old.backbone.plastic_strain,(AbaqusCDPSubstepIntegrator::SymmetricTensor{}));
  EXPECT_DOUBLE_EQ(old.viscous_tension_damage,0.);
}
