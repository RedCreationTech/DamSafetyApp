#include "AbaqusCDPSubstepIntegrator.h"
#include "CDPObservation.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

class AbaqusCDPRealStateTestAccess
{
public:
  using Local = AbaqusCDPLocalIntegrator;
  static Local::LocalVector residual(const Local & core, const Local::LocalVector & z,
                                    const Local::SymmetricTensor & strain,
                                    const Local::State & old, double scale)
  { return core.evaluate(z, strain, old, scale).residual; }
  static Local::LocalMatrix ad(const Local & core, const Local::LocalVector & z,
                              const Local::SymmetricTensor & strain,
                              const Local::State & old, double scale)
  { return core.automaticDifferentiationJacobian(z, strain, old, scale); }
  static Local::LocalMatrix productionFD(const Local & core, const Local::LocalVector & z,
                                        const Local::SymmetricTensor & strain,
                                        const Local::State & old, double scale)
  { return core.numericalJacobian(z, strain, old, scale); }
  static double strainScale(const Local & core) { return core._strain_scale; }
  static std::vector<double> diagnostics(const Local & core, const Local::LocalVector & z,
                                        const Local::SymmetricTensor & strain,
                                        const Local::State & old, double scale)
  {
    const auto e = core.evaluate(z, strain, old, scale);
    const auto stress = AbaqusCDPFormula::stressInvariants(e.stress);
    const auto plastic = AbaqusCDPFormula::stressInvariants(e.plastic_increment);
    std::vector<double> values{e.tensile_equivalent_plastic_strain, e.compressive_equivalent_plastic_strain,
                               e.plastic_multiplier, e.tensile_increment, e.compressive_increment};
    values.insert(values.end(), stress.principal_stress.begin(), stress.principal_stress.end());
    values.insert(values.end(), plastic.principal_stress.begin(), plastic.principal_stress.end());
    values.push_back(stress.tension_weight);
    return values;
  }
};

namespace
{
using Local = AbaqusCDPLocalIntegrator;
using StateCore = AbaqusCDPStateIntegrator;
using Sub = AbaqusCDPSubstepIntegrator;
using Tensor = Local::SymmetricTensor;
using State = StateCore::State;
const std::string data = "unit/data/ds04g-realstate-tangent/";
constexpr double E = 2.97915e10;

template <class Array> void read(std::istream & in, Array & values)
{ for (auto & x : values) in >> x; }
template <class Array> std::vector<double> vec(const Array & a)
{ return {a.begin(), a.end()}; }

std::array<double, 16> packed(const State & s)
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
State unpack(const std::array<double, 16> & a)
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
struct Fixture
{
  std::string id;
  Tensor old_strain, target, stress;
  State old;
  double dt;
  unsigned int N;
  std::array<double, 16> state;
  std::array<double, 36> tangent;
  bool has_local;
  Tensor local_target;
  Local::State local_old;
  double stress_scale, strain_scale;
  Local::LocalVector z, residual;
  std::array<double, 81> ad, fd;
};
Fixture load(const std::string & id)
{
  std::ifstream in(data + id + ".txt");
  Fixture f{};
  std::array<double, 16> old;
  in >> f.id;
  read(in, f.old_strain); read(in, f.target); read(in, old);
  f.old = unpack(old);
  in >> f.dt >> f.N;
  read(in, f.stress); read(in, f.state); read(in, f.tangent);
  in >> f.has_local;
  if (f.has_local)
  {
    read(in, f.local_target); read(in, f.local_old.plastic_strain);
    in >> f.local_old.tensile_equivalent_plastic_strain
       >> f.local_old.compressive_equivalent_plastic_strain >> f.stress_scale >> f.strain_scale;
    read(in, f.z); read(in, f.residual); read(in, f.ad); read(in, f.fd);
  }
  if (!in || f.id != id) throw std::runtime_error("invalid real-state fixture " + id);
  return f;
}
struct Output
{
  std::ofstream file;
  Output(const std::string & id)
  {
    const char * directory = std::getenv("DS04G_TANGENT_OUTPUT");
    if (directory)
    {
      file.open(std::string(directory) + "/" + id + ".csv");
      if (!file) throw std::runtime_error("cannot open output directory");
      file << "layer,mode,direction,h,side,N,signature,status,values,note\n" << std::setprecision(17);
    }
  }
  void row(const std::string & layer, const std::string & mode, int direction, double h,
           int side, unsigned int n, const std::string & signature,
           const std::vector<double> & values, std::string note = {}, const char * status = "ok")
  {
    if (!file.is_open()) return;
    for (char & c : note) if (c == ',' || c == '\n' || c == '\r' || c == '"') c = ';';
    file << layer << ',' << mode << ',' << direction << ',' << h << ',' << side << ',' << n
         << ',' << signature << ',' << status << ',';
    for (const auto x : values) file << x << ' ';
    file << ',' << note << '\n';
  }
};

// Observe existing evaluations only; retain baseline roots and actual partition failures.
struct Trace : CDPObservation::Sink
{
  Output & out;
  explicit Trace(Output & value) : out(value) {}
  bool enterPoint(std::uint64_t, unsigned int, bool) override { return true; }
  void write(const char * stream, const char * event, const char * field,
             const std::vector<double> & values, const std::string & note) noexcept override
  {
    try
    {
      const auto & c = CDPObservation::context;
      out.row(std::string("trace-") + stream, std::string(event) + ":" + field,
              c.substep, c.iteration, c.partition, c.subdivisions,
              std::to_string(c.local_solve), values, note);
    }
    catch (...) { failed = true; }
  }
  bool failed = false;
};

int sign(double x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
std::string tableSegment(const CDPMaterialTable & table, CDPMaterialTable::Branch branch, double k)
{
  const auto & a = table.equivalentPlasticStrain(branch);
  const auto lower = std::lower_bound(a.begin(), a.end(), k);
  const auto upper = std::upper_bound(a.begin(), a.end(), k);
  return std::to_string(lower - a.begin()) + ":" + std::to_string(upper - a.begin());
}
std::string signature(const CDPMaterialTable & table, const StateCore::Result & r,
                      const State & old)
{
  std::ostringstream s;
  s << Local::branchName(r.backbone.active_branch) << ':'
    << tableSegment(table, CDPMaterialTable::Branch::TENSION,
                    r.state.backbone.tensile_equivalent_plastic_strain) << ':'
    << tableSegment(table, CDPMaterialTable::Branch::COMPRESSION,
                    r.state.backbone.compressive_equivalent_plastic_strain);
  Tensor increment;
  for (unsigned int j = 0; j < 6; ++j)
    increment[j] = r.state.backbone.plastic_strain[j] - old.backbone.plastic_strain[j];
  for (const auto & tensor : {r.backbone.effective_stress, r.viscous_effective_stress, increment})
  {
    const auto v = AbaqusCDPFormula::stressInvariants(tensor).principal_stress;
    for (const auto x : v) s << ':' << sign(x);
    s << ':' << (std::abs(v[1] - v[0]) <= 1e-10 * std::max(1.0, std::abs(v[2])))
      << ':' << (std::abs(v[2] - v[1]) <= 1e-10 * std::max(1.0, std::abs(v[2])));
  }
  return s.str();
}
std::vector<double> resultValues(const StateCore::Result & r)
{
  auto values = vec(r.cauchy_stress);
  const auto state = packed(r.state);
  values.insert(values.end(), state.begin(), state.end());
  return values;
}
struct Map
{
  StateCore::Result result;
  unsigned int n;
  std::string signature;
};
Map fixed(const StateCore & core, const CDPMaterialTable & table, const Fixture & f,
          const Tensor & target, unsigned int n, Output * out = nullptr,
          const std::string & mode = {}, int direction = -1, double h = 0, int side = 0)
{
  State working = f.old;
  std::optional<StateCore::Result> final;
  std::string sig;
  for (unsigned int i = 1; i <= n; ++i)
  {
    Tensor strain;
    for (unsigned int j = 0; j < 6; ++j)
      strain[j] = f.old_strain[j] + (double(i) / n) * (target[j] - f.old_strain[j]);
    StateCore::Result r;
    try { r = core.integrate(strain, f.dt / n, working); }
    catch (const std::exception & e)
    { throw std::runtime_error("fixed_N=" + std::to_string(n) + ";failed_substep=" + std::to_string(i) + ";" + e.what()); }
    const auto step_sig = signature(table, r, working);
    sig += (i == 1 ? "" : "|") + step_sig;
    if (out) out->row("substate", mode, direction, h, side, n, step_sig, resultValues(r), std::to_string(i));
    working = r.state;
    final = std::move(r);
  }
  return {*final, n, sig};
}
Map adaptive(const Sub & sub, const StateCore & core, const CDPMaterialTable & table,
             const Fixture & f, const Tensor & target, Output * out, int direction, double h, int side)
{
  const auto r = sub.integrate(f.old_strain, target, f.dt, f.old);
  unsigned int starting_n = 1;
  double increment = 0;
  for (unsigned int j = 0; j < 6; ++j) increment = std::max(increment, std::abs(target[j] - f.old_strain[j]));
  while (increment / starting_n > 2.5e-5) starting_n *= 2;
  for (unsigned int n = starting_n; n < r.accepted_substeps; n *= 2)
  {
    try
    {
      fixed(core, table, f, target, n);
      ADD_FAILURE() << "primal adaptive selector skipped a successful fixed partition";
    }
    catch (const std::exception & e)
    { if (out) out->row("replayed-failed-partition", "adaptive", direction, h, side, n, "", {}, e.what(), "failed"); }
  }
  const auto replay = fixed(core, table, f, target, r.accepted_substeps, out, "adaptive", direction, h, side);
  EXPECT_EQ(resultValues(r.final_result), resultValues(replay.result));
  if (out) out->row("path", "adaptive", direction, h, side, r.accepted_substeps, replay.signature,
                    {double(r.attempted_partitions), double(r.cutback_count),
                     double(r.total_local_iterations), double(r.total_jacobian_fallbacks)});
  return {r.final_result, r.accepted_substeps, replay.signature};
}

void localCheck(const Local & core, const CDPMaterialTable & table, const Fixture & f, Output & out)
{
  const auto r = AbaqusCDPRealStateTestAccess::residual(core, f.z, f.local_target, f.local_old, f.stress_scale);
  const auto ad = AbaqusCDPRealStateTestAccess::ad(core, f.z, f.local_target, f.local_old, f.stress_scale);
  const auto fd = AbaqusCDPRealStateTestAccess::productionFD(core, f.z, f.local_target, f.local_old, f.stress_scale);
  out.row("local-replay", "residual", -1, 0, 0, 0, "", vec(r));
  out.row("local-expected", "residual", -1, 0, 0, 0, "", vec(f.residual));
  for (unsigned int i = 0; i < 9; ++i)
  {
    EXPECT_NEAR(r[i], f.residual[i], 1e-12);
    out.row("local-ad", "row", i, 0, 0, 0, "", vec(ad[i]));
    out.row("local-production-fd", "row", i, 0, 0, 0, "", vec(fd[i]));
    for (unsigned int j = 0; j < 9; ++j)
    {
      EXPECT_NEAR(ad[i][j], f.ad[9 * i + j], 1e-10 * std::max(1.0, std::abs(f.ad[9*i+j])));
      EXPECT_NEAR(fd[i][j], f.fd[9 * i + j], 1e-10 * std::max(1.0, std::abs(f.fd[9*i+j])));
    }
  }
  const auto evaluate = [&](const Local::LocalVector & z) {
    return AbaqusCDPRealStateTestAccess::residual(core, z, f.local_target, f.local_old, f.stress_scale);
  };
  const auto local_signature = [&](const Local::LocalVector & z, int column, double h, int side) {
    const auto d = AbaqusCDPRealStateTestAccess::diagnostics(core, z, f.local_target, f.local_old, f.stress_scale);
    std::string s = tableSegment(table, CDPMaterialTable::Branch::TENSION, d[0]) + ":" +
                    tableSegment(table, CDPMaterialTable::Branch::COMPRESSION, d[1]);
    for (unsigned int i = 2; i < 11; ++i) s += ":" + std::to_string(sign(d[i]));
    out.row("local-diagnostics", "raw", column, h, side, 0, s, d);
    return s;
  };
  const auto center_signature = local_signature(f.z, -1, 0, 0);
  for (unsigned int j = 0; j < 9; ++j)
    for (const auto factor : {1.0, .3, .1, .03, .01})
      for (int restricted = 0; restricted < 2; ++restricted)
      {
        double h = 1e-7 * std::max(1.0, std::abs(f.z[j])) * factor;
        const double kc = f.local_old.compressive_equivalent_plastic_strain + f.z[8] * f.strain_scale;
        if (restricted && j == 8 && kc > 0)
          h = std::min(h, kc / (4 * f.strain_scale) * factor);
        auto plus = f.z, minus = f.z, twice = f.z;
        plus[j] += h; minus[j] -= h; twice[j] += 2*h;
        const auto a = evaluate(plus), b = evaluate(minus), c = evaluate(twice);
        const auto sa = local_signature(plus, j, h, 1);
        const auto sb = local_signature(minus, j, h, -1);
        const auto sc = local_signature(twice, j, h, 2);
        Local::LocalVector central, forward;
        for (unsigned int i = 0; i < 9; ++i)
        {
          central[i] = (a[i] - b[i]) / (2*h);
          forward[i] = (-3*r[i] + 4*a[i] - c[i]) / (2*h);
        }
        const std::string mode = restricted ? "same-side" : "original-scale";
        out.row("local-fd-central", mode, j, h, 0, 0, center_signature+"|"+sa+"|"+sb, vec(central));
        out.row("local-fd-forward", mode, j, h, 1, 0, center_signature+"|"+sa+"|"+sc, vec(forward));
        out.row("local-sample", mode, j, h, 1, 0, "", vec(a));
        out.row("local-sample", mode, j, h, -1, 0, "", vec(b));
      }
}

class AbaqusCDPRealState : public testing::TestWithParam<const char *> {};
TEST_P(AbaqusCDPRealState, FrozenInputReplayAndIndependentDifferences)
{
  const auto f = load(GetParam());
  Output out(f.id);
  CDPMaterialTable table(data + "compression_hardening.csv", data + "compression_damage.csv",
                         data + "tension_stiffening.csv", data + "tension_damage.csv", E);
  Local local(table, {E, .2, 30., .1, 1.16, .6667});
  StateCore core(local, {0., 1., .0005});
  Sub sub(core, {256, 2.5e-5, 1e-8});
  const auto old_before = packed(f.old);
  if (f.has_local)
  {
    EXPECT_EQ(AbaqusCDPRealStateTestAccess::strainScale(local), f.strain_scale);
    localCheck(local, table, f, out);
  }
  Trace trace(out);
  auto * previous = CDPObservation::configured;
  CDPObservation::configured = &trace;
  std::optional<Sub::LinearizedResult> linear;
  {
    CDPObservation::PointScope scope(0, 0);
    linear = sub.integrateLinearized(f.old_strain, f.target, f.dt, f.old);
  }
  CDPObservation::configured = previous;
  ASSERT_FALSE(trace.failed);
  const auto baseline = adaptive(sub, core, table, f, f.target, &out, -1, 0, 0);
  out.row("replay", "primal", -1, 0, 0, baseline.n, baseline.signature, resultValues(baseline.result));
  out.row("replay", "linearized", -1, 0, 0, linear->result.accepted_substeps, "", resultValues(linear->result.final_result));
  std::vector<double> expected = vec(f.stress);
  expected.insert(expected.end(), f.state.begin(), f.state.end());
  out.row("expected", "recorded", -1, 0, 0, f.N, "", expected);
  std::vector<double> tangent;
  for (const auto & column : linear->algorithmic_tangent)
    tangent.insert(tangent.end(), column.begin(), column.end());
  out.row("tangent", "candidate", -1, 0, 0, linear->result.accepted_substeps, "", tangent);
  out.row("tangent", "recorded", -1, 0, 0, f.N, "", vec(f.tangent));
  bool g0 = baseline.n == f.N && linear->result.accepted_substeps == f.N;
  for (unsigned int i = 0; i < 6; ++i)
  {
    const double peak = std::abs(*std::max_element(f.stress.begin(), f.stress.end(),
                                      [](double a, double b){return std::abs(a)<std::abs(b);}));
    const double tol = std::max(.01, 1e-10 * peak);
    g0 &= std::abs(baseline.result.cauchy_stress[i] - f.stress[i]) <= tol;
    g0 &= std::abs(linear->result.final_result.cauchy_stress[i] - f.stress[i]) <= tol;
  }
  for (unsigned int i = 0; i < 16; ++i)
  {
    const double tol = i < 14 ? 1e-12 : 1e-10;
    g0 &= std::abs(packed(baseline.result.state)[i] - f.state[i]) <= tol;
    g0 &= std::abs(packed(linear->result.final_result.state)[i] - f.state[i]) <= tol;
  }
  ASSERT_TRUE(g0) << "Stop differences: frozen input replay mismatch";
  const auto repeated = sub.integrate(f.old_strain, f.target, f.dt, f.old);
  EXPECT_EQ(resultValues(repeated.final_result), resultValues(baseline.result));
  EXPECT_EQ(repeated.accepted_substeps, baseline.n);
  // Only the four accepted points are full-map derivative targets.
  if (!f.has_local)
  {
    std::optional<std::pair<Tensor, Tensor>> bracket;
    for (const unsigned int n : {1u, 8u})
      if (f.id == "F1" || f.id == "F2")
        try
        {
          const auto r = fixed(core, table, f, f.target, n);
          out.row("counterfactual", "fixed", -1, 0, 0, n, r.signature, resultValues(r.result));
        }
        catch (const std::exception & e)
        { out.row("counterfactual", "fixed", -1, 0, 0, n, "", {}, e.what(), "failed"); }
    for (int d = 0; d < 8; ++d)
    {
      Tensor direction{};
      if (d < 6) direction[d] = 1;
      else if (d == 6) direction[0] = direction[1] = 1/std::sqrt(2.0);
      else { direction[0] = direction[3] = 1/std::sqrt(3.0); direction[1] = -direction[0]; }
      out.row("direction", "unit", d, 0, 0, f.N, "", vec(direction));
      for (const double h : {1e-6, 1e-7, 1e-8, 1e-9, 1e-10})
        for (const auto & mode : {std::string("fixed"), std::string("adaptive")})
        {
          unsigned int positive_n = 0, negative_n = 0;
          for (const int side : {1, -1, 2, -2})
          {
            auto target = f.target;
            for (unsigned int i = 0; i < 6; ++i) target[i] += side*h*direction[i];
            try
            {
              const auto r = mode == "fixed" ? fixed(core, table, f, target, f.N, &out, mode, d, h, side)
                                             : adaptive(sub, core, table, f, target, &out, d, h, side);
              out.row("sample", mode, d, h, side, r.n, r.signature, resultValues(r.result));
              if (mode == "adaptive" && side == 1) positive_n = r.n;
              if (mode == "adaptive" && side == -1) negative_n = r.n;
            }
            catch (const std::exception & e)
            { out.row("sample", mode, d, h, side, 0, "", {}, e.what(), "failed"); }
          }
          if (!bracket && positive_n && negative_n && positive_n != negative_n)
          {
            auto plus = f.target, minus = f.target;
            for (unsigned int j = 0; j < 6; ++j)
            { plus[j] += h*direction[j]; minus[j] -= h*direction[j]; }
            bracket = std::make_pair(minus, plus);
          }
        }
    }
    if (bracket)
    {
      auto left = bracket->first, right = bracket->second;
      auto l = adaptive(sub, core, table, f, left, nullptr, -2, 0, 0);
      auto r = adaptive(sub, core, table, f, right, nullptr, -2, 0, 0);
      for (unsigned int i = 0; i < 12; ++i)
      {
        Tensor mid;
        for (unsigned int j = 0; j < 6; ++j) mid[j] = .5*(left[j]+right[j]);
        const auto m = adaptive(sub, core, table, f, mid, nullptr, -2, 0, 0);
        if (m.n == l.n) { left = mid; l = m; }
        else { right = mid; r = m; }
      }
      for (const int side : {-1, 1})
      {
        const auto & target = side == -1 ? left : right;
        const auto & map = side == -1 ? l : r;
        auto values = vec(target);
        const auto ending = resultValues(map.result);
        values.insert(values.end(), ending.begin(), ending.end());
        out.row("selector-bracket", "12-bisections", -2, 0, side, map.n, map.signature, values);
      }
    }
    auto plus = f.target, minus = f.target;
    plus[0] += 1e-8; minus[0] -= 1e-8;
    const auto a = sub.integrate(f.old_strain, plus, f.dt, f.old);
    const auto b = sub.integrate(f.old_strain, minus, f.dt, f.old);
    const auto reverse_b = sub.integrate(f.old_strain, minus, f.dt, f.old);
    const auto reverse_a = sub.integrate(f.old_strain, plus, f.dt, f.old);
    EXPECT_EQ(resultValues(a.final_result), resultValues(reverse_a.final_result));
    EXPECT_EQ(resultValues(b.final_result), resultValues(reverse_b.final_result));
  }
  EXPECT_EQ(packed(f.old), old_before);
}
INSTANTIATE_TEST_CASE_P(DS04G, AbaqusCDPRealState,
    testing::Values("F1", "F2", "F3", "F4", "L1", "L2", "L3"));
}
