#include "AbaqusCDPSubstepIntegrator.h"
#include "CDPObservation.h"
#include "gtest/gtest.h"
#include <map>
#include <string>
#include <vector>

namespace
{
struct Sink : CDPObservation::Sink
{
  std::map<std::string, std::vector<double>> values;
  std::vector<std::uint64_t> trial_ids;
  bool enterPoint(std::uint64_t e, unsigned int qp, bool) override { return e == 2597 && qp == 0; }
  void write(const char *, const char * event, const char * field,
             const std::vector<double> & v, const std::string &) noexcept override
  {
    try {
      values[std::string(event) + ":" + field] = v;
      if (std::string(event) == "trial" && std::string(field) == "trial_stress")
        trial_ids.push_back(CDPObservation::context.local_solve);
    } catch (...) {}
  }
};
}

TEST(DS04G2597Observation, DisabledAndSelectedPathsHaveIdenticalMaterialResult)
{
  const std::string d = "unit/data/ds04g-sp1/";
  CDPMaterialTable table(d+"compression_hardening.csv", d+"compression_damage.csv",
                         d+"tension_stiffening.csv", d+"tension_damage.csv", 2.97915e10);
  AbaqusCDPLocalIntegrator local(table, {2.97915e10, .2, 30., .1, 1.16, .6667});
  AbaqusCDPStateIntegrator state(local, {0, 1, .0005});
  AbaqusCDPSubstepIntegrator integrator(state, {256, 2.5e-5, 1e-8});
  Sink sink;
  for (const double strain : {1e-6, 6e-5})
  {
    const AbaqusCDPSubstepIntegrator::SymmetricTensor target = {strain, 0., 0., 0., 0., 0.};
    const auto expected = integrator.integrateLinearized({}, target, .002, {});
    CDPObservation::configured = &sink;
    {
      CDPObservation::PointScope scope(2597, 0);
      const auto actual = integrator.integrateLinearized({}, target, .002, {});
      EXPECT_EQ(actual.algorithmic_tangent, expected.algorithmic_tangent);
      EXPECT_EQ(actual.result.final_result.cauchy_stress, expected.result.final_result.cauchy_stress);
      EXPECT_EQ(actual.result.final_result.state.backbone.plastic_strain, expected.result.final_result.state.backbone.plastic_strain);
      EXPECT_EQ(actual.result.final_result.state.viscous_plastic_strain, expected.result.final_result.state.viscous_plastic_strain);
      EXPECT_EQ(actual.result.final_result.state.backbone.tensile_equivalent_plastic_strain, expected.result.final_result.state.backbone.tensile_equivalent_plastic_strain);
      EXPECT_EQ(actual.result.final_result.state.backbone.compressive_equivalent_plastic_strain, expected.result.final_result.state.backbone.compressive_equivalent_plastic_strain);
      EXPECT_EQ(actual.result.final_result.backbone.active_branch, expected.result.final_result.backbone.active_branch);
      EXPECT_EQ(actual.result.final_result.state.viscous_tension_damage, expected.result.final_result.state.viscous_tension_damage);
      EXPECT_EQ(actual.result.final_result.state.viscous_compression_damage, expected.result.final_result.state.viscous_compression_damage);
      EXPECT_EQ(actual.result.accepted_substeps, expected.result.accepted_substeps);
      EXPECT_EQ(actual.result.cutback_count, expected.result.cutback_count);
    }
    CDPObservation::configured = nullptr;
    EXPECT_FALSE(CDPObservation::active());
    ASSERT_EQ(sink.values["trial:trial_stress"].size(), 6u);
    ASSERT_EQ(sink.values["trial:ft_fc_scale_f_threshold"].size(), 5u);
  }
  ASSERT_GE(sink.trial_ids.size(), 4u);
  for (std::size_t i = 1; i < sink.trial_ids.size(); ++i)
    EXPECT_GT(sink.trial_ids[i], sink.trial_ids[i - 1]);
  CDPObservation::configured = &sink;
  { CDPObservation::PointScope other(2597, 1); EXPECT_FALSE(CDPObservation::active()); }
  CDPObservation::configured = nullptr;
}
