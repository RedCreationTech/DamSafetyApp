#include "CDPObservationProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"
#include "MaterialPropertyStorage.h"
#include "RankTwoTensor.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include <cmath>
#include <iomanip>
#include <sstream>
registerMooseObject("DamSafetyApp", CDPObservationProblem);
namespace
{
std::string quoted(const std::string & s)
{
  std::string r = "\"";
  for (char c : s) { if (c == '"') r += '"'; r += c; }
  return r + '"';
}
}
InputParameters CDPObservationProblem::validParams()
{
  auto p = FEProblem::validParams();
  p.addParam<bool>("observation_enabled", false, "Enable passive recorders, never extra assembly");
  p.addParam<Real>("observation_start", 2.698, "First observed cumulative target time");
  p.addParam<Real>("observation_end", 2.720, "Last observed cumulative target time");
  p.addParam<Real>("earthquake_offset", 1., "Cumulative time origin of EQ");
  p.addParam<std::string>("observation_prefix", "ds04g-observation", "Rank-file prefix");
  p.addParam<unsigned long long>("observation_byte_limit", 65536000, "Per-rank detailed byte cap; mark incomplete, keep solve running");
  p.addParam<std::vector<unsigned int>>("observation_labels", {}, "Source element labels");
  p.addParam<std::vector<Real>>("observation_vertices", {}, "Four x/y vertices per label in SI; compare unordered physical sets");
  return p;
}
CDPObservationProblem::CDPObservationProblem(const InputParameters & p)
  : FEProblem(p), _enabled(getParam<bool>("observation_enabled")),
    _start(getParam<Real>("observation_start")), _end(getParam<Real>("observation_end")),
    _offset(getParam<Real>("earthquake_offset")), _prefix(getParam<std::string>("observation_prefix")),
    _limit(getParam<unsigned long long>("observation_byte_limit"))
{
  if (_enabled && libMesh::n_threads() != 1) paramError("observation_enabled", "Requires one thread per MPI rank");
  if (_start > _end) paramError("observation_start", "Must not exceed end");
}
CDPObservationProblem::~CDPObservationProblem()
{
  if (_enabled) status();
  if (CDPObservation::configured == this) CDPObservation::configured = nullptr;
}
void CDPObservationProblem::initialSetup()
{
  FEProblem::initialSetup();
  if (!_enabled) return;
  const auto & labels = getParam<std::vector<unsigned int>>("observation_labels");
  const auto & xy = getParam<std::vector<Real>>("observation_vertices");
  if (labels.empty() || xy.size() != labels.size()*8) paramError("observation_vertices", "Expected eight coordinates per label");
  for (unsigned int k=0; k<labels.size(); ++k)
  {
    unsigned int matches=0;
    for (const auto * e : mesh().getMesh().active_local_element_ptr_range())
    {
      if (e->dim()!=2 || e->n_nodes()!=4) continue;
      unsigned int hit=0;
      for (unsigned int i=0;i<4;++i)
        for (unsigned int j=0;j<4;++j)
          if (std::abs(e->point(i)(0)-xy[k*8+j*2])<1e-8 &&
              std::abs(e->point(i)(1)-xy[k*8+j*2+1])<1e-8) ++hit;
      if (hit==4) { ++matches; _targets.emplace(e->id(),labels[k]); }
    }
    comm().sum(matches);
    if (matches != 1) paramError("observation_vertices", "Source element mapping must be unique: ",labels[k]);
  }
  if (CDPObservation::configured) mooseError("Only one observation problem per process is supported");
  CDPObservation::configured = this;
  for (const auto & target : _targets)
  {
    _element=target.first;
    const auto * e=mesh().getMesh().elem_ptr(target.first);
    for (unsigned int i=0;i<4;++i)
      write("probe-map", "vertex", "node_id_x_y_z", {static_cast<double>(e->node_id(i)),e->point(i)(0),e->point(i)(1),e->point(i)(2)}, {});
  }
  _element=0;
  status();
}
bool CDPObservationProblem::window() const { return _enabled && time()>=_start-1e-9 && time()<=_end+1e-9; }
bool CDPObservationProblem::enterPoint(std::uint64_t e, unsigned int qp, bool new_call)
{
  if (!window() || !_targets.count(e) || _truncated || _io_failed) return false;
  _element=e; _qp=qp;
  if (new_call) { _call=++_next_call; _last_call[{e,qp}]={_evaluation,_call}; }
  else
  {
    const auto it=_last_call.find({e,qp});
    _call=it!=_last_call.end() && it->second.first==_evaluation ? it->second.second : 0;
  }
  return true;
}
void CDPObservationProblem::write(const char * stream, const char * event, const char * name,
                                 const std::vector<double> & v, const std::string & note) noexcept
{
  try
  {
    if (!_enabled || _io_failed) return;
    const bool light=std::string(stream)=="attempt-events";
    if (!light && (_truncated || _bytes>=_limit)) { _truncated=true; return; }
    auto & f=_streams[stream];
    if (!f.is_open())
    {
      f.open(_prefix+"-"+stream+"-rank"+std::to_string(processor_id())+".csv",std::ios::out);
      f << "record,time,eq_time,dt,dt_old,step,attempt,evaluation,kind,element,source_element,qp,call,partition,N,substep,local_solve,local_iteration,event,field,values,note\n";
    }
    std::ostringstream a;
    a << std::setprecision(17);
    for (unsigned int i=0;i<v.size();++i) { if(i) a << ' '; a << v[i]; }
    const auto & c=CDPObservation::context;
    std::ostringstream row;
    row << std::setprecision(17) << ++_record << ',' << time() << ',' << time()-_offset << ',' << dt() << ',' << dtOld() << ',' << timeStep() << ',' << _attempt << ',' << _evaluation << ',' << _kind << ',' << _element << ',' << (_targets.count(_element)?_targets.at(_element):0) << ',' << _qp << ',' << _call << ',' << c.partition << ',' << c.subdivisions << ',' << c.substep << ',' << c.local_solve << ',' << c.iteration << ',' << event << ',' << name << ',' << quoted(a.str()) << ',' << quoted(note) << '\n';
    const auto text=row.str(); f << text; _bytes+=text.size();
    if (!f) _io_failed=true;
    if (std::string(event).find("exception")!=std::string::npos) f.flush();
  }
  catch (...) { _io_failed=true; }
}
void CDPObservationProblem::status() noexcept
{
  try
  {
    for(auto & f:_streams) f.second.flush();
    std::ofstream f(_prefix+"-status-rank"+std::to_string(processor_id())+".json");
    f << "{\"schema\":\"ds04g-observation-v1\",\"rank\":" << processor_id() << ",\"records\":" << _record << ",\"bytes\":" << _bytes << ",\"capture_truncated\":" << (_truncated?"true":"false") << ",\"io_failed\":" << (_io_failed?"true":"false") << ",\"history_incomplete\":" << (_history_incomplete?"true":"false") << ",\"target_elements_on_rank\":" << _targets.size() << "}\n";
  } catch (...) { _io_failed=true; }
}
void CDPObservationProblem::history(const char * event)
{
  const auto & storage=getMaterialPropertyStorage();
  for (auto state:storage.stateIndexRange())
    for (const auto & entry:storage.props(state))
    {
      if (!_targets.count(entry.first->id())) continue;
      _element=entry.first->id();
      for (const auto & side:entry.second)
        for(unsigned int id=0;id<side.second.size();++id)
        {
          const auto * value=side.second.queryValue(id);
          if(!value) continue;
          const auto name=storage.queryStatefulPropName(storage.statefulProps().at(id));
          const auto * scalar=dynamic_cast<const MaterialProperty<Real> *>(value);
          const auto * tensor=dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
          if (!name || (!scalar && !tensor)) { _history_incomplete=true; continue; }
          for(unsigned int qp=0;qp<value->size();++qp)
          {
            _qp=qp; std::vector<double> values;
            if(scalar) values.push_back((*scalar)[qp]);
            else for(unsigned int c=0;c<9;++c) values.push_back((*tensor)[qp](c/3,c%3));
            const auto field="state"+std::to_string(state)+"_side"+std::to_string(side.first)+"_"+*name;
            write("history-events",event,field.c_str(),values,{});
          }
        }
    }
  for (auto * system : {static_cast<SystemBase *>(&getNonlinearSystem(0)),static_cast<SystemBase *>(&getAuxiliarySystem())})
  {
    const auto & sys=system->system();
    const auto * u=system->currentSolution();
    for(const auto * node:mesh().getMesh().local_node_ptr_range())
    {
      bool selected=false;
      const auto & xy=getParam<std::vector<Real>>("observation_vertices");
      for(unsigned int j=0;j<xy.size();j+=2)
        if(std::abs((*node)(0)-xy[j])<1e-8 && std::abs((*node)(1)-xy[j+1])<1e-8) selected=true;
      if(!selected) continue;
      for(unsigned int var=0;var<sys.n_vars();++var)
      {
        const auto name=sys.variable_name(var);
        if(name!="disp_x" && name!="disp_y" && name!="vel_x" && name!="vel_y" && name!="accel_x" && name!="accel_y") continue;
        if(!node->n_dofs(sys.number(),var)) continue;
        const auto dof=node->dof_number(sys.number(),var,0);
        if(dof<u->first_local_index() || dof>=u->last_local_index()) continue;
        write("history-events",event,name.c_str(),{static_cast<double>(node->id()),static_cast<double>(dof),(*u)(dof),system->solutionOld()(dof),system->solutionOlder()(dof)},"node_dof_current_old_older");
      }
    }
  }
  _element=0; _qp=0;
}
void CDPObservationProblem::attemptEvent(const char * event)
{
  if(!_enabled) return;
  if(std::string(event)=="attempt_enter") ++_attempt;
  _kind="event"; _element=0; _qp=0; CDPObservation::context={};
  write("attempt-events",event,"time",{time()},{});
  if(window() || std::abs(time()-(_offset+1.680))<1e-9 || std::abs(time()-_offset)<1e-9) history(event);
  if (std::string(event).find("solve_")==0)
  {
    auto * solver=dynamic_cast<libMesh::PetscNonlinearSolver<Number> *>(getNonlinearSystem(0).nonlinearSolver());
    if(solver)
    {
      SNESConvergedReason reason; KSPConvergedReason linear_reason; KSP ksp;
      PetscInt iterations; PetscReal norm;
      const auto snes=solver->snes();
      const bool bad=SNESGetConvergedReason(snes,&reason) || SNESGetIterationNumber(snes,&iterations) || SNESGetFunctionNorm(snes,&norm) || SNESGetKSP(snes,&ksp) || KSPGetConvergedReason(ksp,&linear_reason);
      if(bad) _history_incomplete=true;
      else write("attempt-events",event,"SNES_reason_iterations_norm_KSP_reason",{static_cast<double>(reason),static_cast<double>(iterations),norm,static_cast<double>(linear_reason)},{});
    }
  }
  if(window()) status();
}
void CDPObservationProblem::computeResidual(const NumericVector<Number> & x, NumericVector<Number> & r,unsigned int n)
{
  if(!_enabled) { FEProblem::computeResidual(x,r,n); return; }
  ++_evaluation; _kind="residual";
  if(window())
  {
    auto * solver=dynamic_cast<libMesh::PetscNonlinearSolver<Number> *>(getNonlinearSystem(n).nonlinearSolver());
    if(!solver || SNESMonitorSet(solver->snes(),monitor,this,nullptr)) mooseError("Cannot attach observation SNES monitor");
  }
  FEProblem::computeResidual(x,r,n);
  if(window()) { _element=0; write("attempt-events","residual_complete","l2_norm",{r.l2_norm()},{}); }
}
void CDPObservationProblem::computeJacobian(const NumericVector<Number> & x,SparseMatrix<Number> & j,unsigned int n)
{
  if(!_enabled) { FEProblem::computeJacobian(x,j,n); return; }
  ++_evaluation; _kind="jacobian";
  FEProblem::computeJacobian(x,j,n);
}
PetscErrorCode CDPObservationProblem::monitor(SNES snes,PetscInt iteration,PetscReal norm,void * context)
{
  auto & p=*static_cast<CDPObservationProblem *>(context);
  if(!p.window()) return PETSC_SUCCESS;
  SNESConvergedReason reason; PetscCall(SNESGetConvergedReason(snes,&reason));
  p._element=0;
  p.write("attempt-events","SNES_accepted_iterate","iteration_norm_reason",{static_cast<double>(iteration),norm,static_cast<double>(reason)},{});
  return PETSC_SUCCESS;
}
