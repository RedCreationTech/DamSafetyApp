#include "AbaqusCDPSubstepIntegrator.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <fstream>
namespace
{
using State = AbaqusCDPStateIntegrator::State;
const std::string data = "unit/data/ds04g-sp1/";
std::array<double, 16> pack(const State &s)
{
  std::array<double, 16> a;
  std::copy(s.backbone.plastic_strain.begin(), s.backbone.plastic_strain.end(), a.begin());
  a[6] = s.backbone.tensile_equivalent_plastic_strain;
  a[7] = s.backbone.compressive_equivalent_plastic_strain;
  std::copy(s.viscous_plastic_strain.begin(), s.viscous_plastic_strain.end(), a.begin() + 8);
  a[14] = s.viscous_tension_damage;
  a[15] = s.viscous_compression_damage;
  return a;
}
State unpack(const std::array<double, 16> &a)
{
  State s;
  std::copy_n(a.begin(), 6, s.backbone.plastic_strain.begin());
  s.backbone.tensile_equivalent_plastic_strain = a[6];
  s.backbone.compressive_equivalent_plastic_strain = a[7];
  std::copy_n(a.begin() + 8, 6, s.viscous_plastic_strain.begin());
  s.viscous_tension_damage = a[14];
  s.viscous_compression_damage = a[15];
  return s;
}
template <class A> void read(std::istream &s, A &a)
{
  for (auto &x : a)
    s >> x;
}
class DS04GSP1MaterialReplay : public testing::TestWithParam<const char *>
{
};
} // namespace
TEST_P(DS04GSP1MaterialReplay, FrozenCompleteUpdateAndOldStateReplay)
{
  const std::string id = GetParam();
  std::ifstream in(data + id + ".txt");
  std::string got;
  std::array<double, 6> strain0, strain1, stress;
  std::array<double, 16> old, expected;
  std::array<double, 36> tangent;
  std::array<double, 2> bone;
  double dt;
  int branch;
  in >> got;
  read(in, strain0);
  read(in, strain1);
  read(in, old);
  in >> dt;
  read(in, stress);
  read(in, expected);
  read(in, tangent);
  read(in, bone);
  in >> branch;
  ASSERT_TRUE(bool(in));
  ASSERT_EQ(id, got);
  CDPMaterialTable table(data + "compression_hardening.csv", data + "compression_damage.csv",
                         data + "tension_stiffening.csv", data + "tension_damage.csv", 2.97915e10);
  AbaqusCDPLocalIntegrator local(table, {2.97915e10, .2, 30., .1, 1.16, .6667});
  AbaqusCDPStateIntegrator core(local, {0, 1, .0005});
  AbaqusCDPSubstepIntegrator sub(core, {256, 2.5e-5, 1e-8});
  const auto s = unpack(old);
  const auto a = sub.integrateLinearized(strain0, strain1, dt, s);
  double peak = 0, tpeak = 0;
  for (double x : stress)
    peak = std::max(peak, std::abs(x));
  for (double x : tangent)
    tpeak = std::max(tpeak, std::abs(x));
  const auto &r = a.result.final_result;
  for (unsigned i = 0; i < 6; ++i)
    EXPECT_NEAR(r.cauchy_stress[i], stress[i], std::max(.01, 1e-10 * peak));
  const auto state = pack(r.state);
  for (unsigned i = 0; i < 16; ++i)
    EXPECT_NEAR(state[i], expected[i], i < 14 ? 1e-12 : 1e-10);
  for (unsigned c = 0; c < 6; ++c)
    for (unsigned j = 0; j < 6; ++j)
      EXPECT_NEAR(a.algorithmic_tangent[c][j], tangent[6 * c + j], std::max(.01, 1e-10 * tpeak));
  EXPECT_EQ(int(r.backbone.active_branch), branch);
  EXPECT_EQ(pack(s), old);
  // A discarded longer-step trial must not alter a repeat from the same accepted old state.
  try
  {
    (void)sub.integrate(strain0, strain1, dt * 2, s);
  }
  catch (const std::exception &)
  {
  }
  const auto b = sub.integrateLinearized(strain0, strain1, dt, s);
  EXPECT_EQ(pack(b.result.final_result.state), state);
  EXPECT_EQ(b.algorithmic_tangent, a.algorithmic_tangent);
  EXPECT_EQ(pack(s), old);
}
INSTANTIATE_TEST_SUITE_P(RealInputs, DS04GSP1MaterialReplay,
                         testing::Values("P01", "P02", "P03", "P04", "P05", "P06"));
