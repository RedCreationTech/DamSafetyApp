#include "AbaqusCDPSubstepIntegrator.h"
#include "CDPObservation.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
using Local = AbaqusCDPLocalIntegrator;
using Core = AbaqusCDPStateIntegrator;
using Sub = AbaqusCDPSubstepIntegrator;
using Tensor = Local::SymmetricTensor;
using State = Core::State;
const std::string data = "unit/data/ds04g-sixstate-update/";
constexpr double E = 2.97915e10;
template <class A> std::vector<double> vec(const A & a) { return {a.begin(), a.end()}; }
template <class A> void read(std::istream & in, A & a) { for (auto & v : a) in >> v; }
std::array<double,16> pack(const State & s)
{
  std::array<double,16> a;
  std::copy(s.backbone.plastic_strain.begin(), s.backbone.plastic_strain.end(), a.begin());
  a[6]=s.backbone.tensile_equivalent_plastic_strain; a[7]=s.backbone.compressive_equivalent_plastic_strain;
  std::copy(s.viscous_plastic_strain.begin(), s.viscous_plastic_strain.end(), a.begin()+8);
  a[14]=s.viscous_tension_damage; a[15]=s.viscous_compression_damage; return a;
}
State unpack(const std::array<double,16> & a)
{
  State s; std::copy_n(a.begin(),6,s.backbone.plastic_strain.begin());
  s.backbone.tensile_equivalent_plastic_strain=a[6]; s.backbone.compressive_equivalent_plastic_strain=a[7];
  std::copy_n(a.begin()+8,6,s.viscous_plastic_strain.begin());
  s.viscous_tension_damage=a[14]; s.viscous_compression_damage=a[15]; return s;
}
std::vector<double> resultValues(const Core::Result & r)
{
  auto a=vec(r.cauchy_stress); const auto s=pack(r.state); a.insert(a.end(),s.begin(),s.end()); return a;
}
struct Fixture
{
  std::string id; Tensor old_strain,target,stress; State old; double dt;
  std::array<double,16> expected; std::array<double,36> tangent; std::array<double,2> bone_damage; int branch;
};
Fixture load(const std::string & id)
{
  Fixture f{}; std::array<double,16> old; std::ifstream in(data+id+".txt");
  in>>f.id; read(in,f.old_strain); read(in,f.target); read(in,old); f.old=unpack(old);
  in>>f.dt; read(in,f.stress); read(in,f.expected); read(in,f.tangent); read(in,f.bone_damage); in>>f.branch;
  if (!in || f.id!=id) throw std::runtime_error("Invalid fixture "+id); return f;
}
std::string directory()
{
  const char * p=std::getenv("DS04G_SIXSTATE_OUTPUT");
  if (!p) throw std::runtime_error("DS04G_SIXSTATE_OUTPUT required"); return p;
}
struct Output
{
  std::ofstream file; std::string arm; unsigned int n=0,j=0;
  Output(const std::string & id, const std::string & phase) : file(directory()+"/"+id+"-"+phase+".csv")
  {
    if(!file) throw std::runtime_error("Cannot open output");
    file<<"arm,N,substep,local_solve,local_iteration,field,status,values,note\n"<<std::setprecision(17);
  }
  void row(const std::string & field,const std::vector<double> & v,const std::string & status="ok",std::string note={})
  {
    for(char & c:note) if(c==',' || c=='\n' || c=='\r' || c=='"') c=';';
    file<<arm<<','<<n<<','<<j<<','<<CDPObservation::context.local_solve<<','<<CDPObservation::context.iteration
        <<','<<field<<','<<status<<',';
    for(double x:v) file<<x<<' ';
    file<<','<<note<<'\n'; if(!file) throw std::runtime_error("Output failure");
  }
};
struct Trace : CDPObservation::Sink
{
  Output & out; bool failed=false; CDPObservation::Sink * previous;
  explicit Trace(Output & value) : out(value),previous(CDPObservation::configured) { CDPObservation::configured=this; }
  ~Trace() { CDPObservation::configured=previous; }
  bool enterPoint(std::uint64_t,unsigned int,bool) override { return true; }
  void write(const char * stream,const char * event,const char * field,const std::vector<double> & values,const std::string & note) noexcept override
  { try { out.row(std::string("trace/")+stream+"/"+event+"/"+field,values,"ok",note); } catch(...) { failed=true; } }
};
std::string segment(const CDPMaterialTable & t,CDPMaterialTable::Branch b,double k)
{
  const auto & x=t.equivalentPlasticStrain(b);
  return std::to_string(std::lower_bound(x.begin(),x.end(),k)-x.begin())+":"+
         std::to_string(std::upper_bound(x.begin(),x.end(),k)-x.begin());
}
void emit(Output & out,const Core::Result & r,const State & before,const CDPMaterialTable & table)
{
  out.row("result22",resultValues(r)); out.row("old16",vec(pack(before)));
  out.row("backbone_plastic",vec(r.state.backbone.plastic_strain));
  out.row("backbone_effective_stress",vec(r.backbone.effective_stress));
  out.row("viscous_effective_stress",vec(r.viscous_effective_stress));
  out.row("bone_damage",{r.backbone.backbone_tension_damage,r.backbone.backbone_compression_damage});
  out.row("stiffness_weight",{r.damage.stiffness_factor});
  out.row("solve",{double(r.backbone.active_branch),r.backbone.residual_norm,double(r.backbone.iterations),double(r.backbone.jacobian_fallbacks)});
  out.row("table_segments",{},"ok",segment(table,CDPMaterialTable::Branch::TENSION,r.state.backbone.tensile_equivalent_plastic_strain)+"|"+segment(table,CDPMaterialTable::Branch::COMPRESSION,r.state.backbone.compressive_equivalent_plastic_strain));
  if (r.backbone.plastic) EXPECT_LE(r.backbone.residual_norm,1e-9);
}
Core::Result fixed(const Core & core,const CDPMaterialTable & table,const Fixture & f,unsigned int n,Output & out,const std::string & arm)
{
  State working=f.old; std::optional<Core::Result> final; out.arm=arm; out.n=n;
  for(unsigned int i=1;i<=n;++i)
  {
    out.j=i; Tensor target;
    for(unsigned int k=0;k<6;++k) target[k]=f.old_strain[k]+(double(i)/n)*(f.target[k]-f.old_strain[k]);
    const auto before=pack(working); out.row("target_dt",[&](){auto v=vec(target);v.push_back(f.dt/n);return v;}());
    CDPObservation::PointScope point(0,0); CDPObservation::context.subdivisions=n; CDPObservation::context.substep=i;
    auto result=core.integrate(target,f.dt/n,working);
    EXPECT_EQ(pack(working),before); emit(out,result,working,table); working=result.state; final=std::move(result);
  }
  return *final;
}
bool matches(const Fixture & f,const Core::Result & r)
{
  bool ok=true; double peak=0; for(double x:f.stress) peak=std::max(peak,std::abs(x));
  for(unsigned int i=0;i<6;++i) ok &= std::isfinite(r.cauchy_stress[i]) && std::abs(r.cauchy_stress[i]-f.stress[i])<=std::max(.01,1e-10*peak);
  const auto state=pack(r.state);
  for(unsigned int i=0;i<16;++i) ok &= std::isfinite(state[i]) && std::abs(state[i]-f.expected[i])<=(i<14?1e-12:1e-10);
  ok &= int(r.backbone.active_branch)==f.branch;
  ok &= (state[7]==0)==(f.expected[7]==0);
  return ok;
}
void saveG0(const std::string & id,const Core::Result & r)
{
  std::ofstream out(directory()+"/"+id+".g0");out<<id<<'\n'<<std::setprecision(17);
  for(double x:resultValues(r)) out<<x<<' ';out<<'\n';
  if(!out) throw std::runtime_error("Cannot save G0 marker");
}
std::vector<double> readG0(const std::string & id)
{
  std::ifstream in(directory()+"/"+id+".g0");std::string got;in>>got;std::vector<double> v(22);read(in,v);
  if(!in || got!=id) throw std::runtime_error("Missing matching G0 result "+id);return v;
}
class DS04GSixStateUpdate : public testing::TestWithParam<const char *> {};
TEST_P(DS04GSixStateUpdate,G0FrozenReplay)
{
  const auto f=load(GetParam());Output out(f.id,"g0");Trace trace(out);
  CDPMaterialTable table(data+"compression_hardening.csv",data+"compression_damage.csv",data+"tension_stiffening.csv",data+"tension_damage.csv",E);
  Local local(table,{E,.2,30.,.1,1.16,.6667});Core core(local,{0.,1.,.0005});Sub sub(core,{256,2.5e-5,1e-8});const auto original=pack(f.old);
  out.arm="G0-L";out.n=1;out.j=1;const auto start=std::chrono::steady_clock::now();
  std::optional<Sub::LinearizedResult> linear;
  { CDPObservation::PointScope point(0,0);linear=sub.integrateLinearized(f.old_strain,f.target,f.dt,f.old); }
  ASSERT_EQ(linear->result.accepted_substeps,1u);ASSERT_EQ(linear->result.cutback_count,0u);ASSERT_TRUE(matches(f,linear->result.final_result));emit(out,linear->result.final_result,f.old,table);
  std::vector<double> tangent;for(const auto & c:linear->algorithmic_tangent)tangent.insert(tangent.end(),c.begin(),c.end());out.row("C6_column_major",tangent);
  double peak=0;for(double x:f.tangent)peak=std::max(peak,std::abs(x));
  for(unsigned int i=0;i<36;++i)ASSERT_LE(std::abs(tangent[i]-f.tangent[i]),std::max(.01,1e-10*peak));
  out.arm="G0-P";std::optional<Sub::Result> primal;
  { CDPObservation::PointScope point(0,0);primal=sub.integrate(f.old_strain,f.target,f.dt,f.old); }
  ASSERT_EQ(primal->accepted_substeps,1u);ASSERT_EQ(primal->cutback_count,0u);ASSERT_TRUE(matches(f,primal->final_result));emit(out,primal->final_result,f.old,table);
  const auto one=fixed(core,table,f,1,out,"F1-before");ASSERT_TRUE(matches(f,one));ASSERT_EQ(pack(f.old),original);
  const auto onevalues=resultValues(one),linearvalues=resultValues(linear->result.final_result),primalvalues=resultValues(primal->final_result);
  for(unsigned int i=0;i<22;++i)
  {
    const double tol=i<6?std::max(.01,1e-10*std::abs(onevalues[i])):i<20?1e-12:1e-10;
    ASSERT_LE(std::abs(onevalues[i]-linearvalues[i]),tol);ASSERT_LE(std::abs(onevalues[i]-primalvalues[i]),tol);
  }
  ASSERT_FALSE(trace.failed);ASSERT_FALSE(HasFailure());out.arm="gate";out.row("G0_passed",{1});
  out.row("wall_seconds",{std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()});saveG0(f.id,one);
}
TEST_P(DS04GSixStateUpdate,FixedN1N2CompleteUpdate)
{
  const char * gate=std::getenv("DS04G_SIXSTATE_GATE");ASSERT_TRUE(gate && std::string(gate)=="G0_ALL_VERIFIED");
  for(const auto & id:{"P01","P02","P03","P04","P05","P06"})ASSERT_EQ(readG0(id).size(),22u);
  const auto f=load(GetParam());const auto expected1=readG0(f.id);Output out(f.id,"compare");Trace trace(out);
  CDPMaterialTable table(data+"compression_hardening.csv",data+"compression_damage.csv",data+"tension_stiffening.csv",data+"tension_damage.csv",E);
  Local local(table,{E,.2,30.,.1,1.16,.6667});Core core(local,{0.,1.,.0005});const auto original=pack(f.old);bool failed=false;
  const auto start=std::chrono::steady_clock::now();
  try { fixed(core,table,f,2,out,"F2"); }
  catch(const std::exception & e) { failed=true;out.row("exception",{},"failed",e.what()); }
  const auto after=fixed(core,table,f,1,out,"F1-after");EXPECT_EQ(resultValues(after),expected1);EXPECT_EQ(pack(f.old),original);
  out.arm="gate";out.row("old_unchanged_repeat_equal",{double(pack(f.old)==original),double(resultValues(after)==expected1)});
  out.row("wall_seconds",{std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()});EXPECT_FALSE(trace.failed);EXPECT_FALSE(failed);
}
INSTANTIATE_TEST_CASE_P(DS04G6,DS04GSixStateUpdate,testing::Values("P01","P02","P03","P04","P05","P06"));
}
