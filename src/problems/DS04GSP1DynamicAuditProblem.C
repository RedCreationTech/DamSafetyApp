#include "DS04GSP1DynamicAuditProblem.h"
#include "DS04GSP1.h"
#include "CDPAssemblyProbe.h"
#include "AuxiliarySystem.h"
#include "KernelBase.h"
#include "MaterialPropertyStorage.h"
#include "MooseException.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "RankTwoTensor.h"
#include "libmesh/elem.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>

registerMooseObject("DamSafetyApp", DS04GSP1DynamicAuditProblem);

InputParameters DS04GSP1DynamicAuditProblem::validParams()
{
  auto p = FEProblem::validParams();
  p.addClassDescription("TEST ONLY: four prescribed dynamic histories on one QUAD4; no solve");
  p.addRequiredParam<std::vector<Real>>("snapshot_u", "Eight interleaved real snapshot U values");
  p.addRequiredParam<std::vector<Real>>("snapshot_v", "Eight interleaved real snapshot V values");
  p.addRequiredParam<std::vector<Real>>("snapshot_a", "Eight interleaved real snapshot A values");
  p.addRequiredParam<std::vector<Real>>("snapshot_stress", "Real snapshot S11 S22 S33 S12");
  return p;
}

DS04GSP1DynamicAuditProblem::DS04GSP1DynamicAuditProblem(const InputParameters & p) : FEProblem(p) {}

void DS04GSP1DynamicAuditProblem::initialSetup()
{
  FEProblem::initialSetup();
  auto & nl = getNonlinearSystem(0);
  auto & aux = getAuxiliarySystem();
  auto & sys = dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  const auto & mesh = this->mesh().getMesh();
  if (comm().size() != 1 || libMesh::n_threads() != 1 || mesh.n_elem() != 1 ||
      mesh.n_nodes() != 4 || sys.n_dofs() != 8 || mesh.elem_ptr(0)->type() != QUAD4)
    mooseError("Dynamic audit requires serial one-QUAD4/eight DOFs; never use on a structure");
  for (const auto & name : {"snapshot_u", "snapshot_v", "snapshot_a", "snapshot_stress"})
    if (getParam<std::vector<Real>>(name).size() != (std::string(name) == "snapshot_stress" ? 4 : 8))
      paramError(name, "Wrong prescribed snapshot size");

  const auto * elem = mesh.elem_ptr(0);
  const auto * original = nl.currentSolution();
  auto x = original->clone();
  auto residual = x->zero_clone();
  auto J = SparseMatrix<Number>::build(comm());
  J->init(8, 8, 8, 8, 8, 0);
  std::array<dof_id_type, 8> dofs, vdofs, adofs;
  DS04GSP1::Nodes nodes;
  std::ofstream mapping("dynamic_dof_map.csv");
  mapping << std::setprecision(17) << "local_node,node_id,x,y,dof_x,dof_y\n";
  for (unsigned i = 0; i < 4; ++i)
  {
    const auto * node = elem->node_ptr(i);
    for (unsigned c = 0; c < 2; ++c)
    {
      const auto axis = c ? "y" : "x";
      nodes[i][c] = (*node)(c);
      dofs[2*i+c] = node->dof_number(sys.number(), sys.variable_number(c ? "disp_y" : "disp_x"), 0);
      vdofs[2*i+c] = node->dof_number(aux.number(), aux.system().variable_number(std::string("vel_")+axis), 0);
      adofs[2*i+c] = node->dof_number(aux.number(), aux.system().variable_number(std::string("accel_")+axis), 0);
    }
    mapping << i << ',' << node->id() << ',' << (*node)(0) << ',' << (*node)(1) << ','
            << dofs[2*i] << ',' << dofs[2*i+1] << '\n';
  }
  const auto geom = DS04GSP1::geometry(nodes);
  std::array<std::unique_ptr<NumericVector<Number>>, 3> original_nl, original_aux;
  for (unsigned state = 0; state < 3; ++state)
  {
    original_nl[state] = nl.solutionState(state).clone();
    original_aux[state] = aux.solutionState(state).clone();
  }
  const Real original_time = time(), original_oldtime = timeOld(), original_dt = dt(), original_dtold = dtOld();
  const int original_step = timeStep();
  auto & storage = const_cast<MaterialPropertyStorage &>(getMaterialPropertyStorage());
  // Test injection is done once per case, never during a residual/Jacobian evaluation.
  auto visit = [&](auto action)
  {
    for (auto state : storage.stateIndexRange())
      for (const auto & entry : storage.props(state))
        for (const auto & side : entry.second)
        {
          auto & props = storage.setProps(entry.first, side.first, state);
          for (unsigned id = 0; id < props.size(); ++id)
          {
            auto * value = props.queryValue(id);
            if (!value) continue;
            const auto prop = storage.queryStatefulPropName(storage.statefulProps().at(id));
            auto * scalar = dynamic_cast<MaterialProperty<Real> *>(value);
            auto * tensor = dynamic_cast<MaterialProperty<RankTwoTensor> *>(value);
            if (!prop || (!scalar && !tensor)) mooseError("Unsupported diagnostic history property");
            for (unsigned qp = 0; qp < value->size(); ++qp)
              for (unsigned c = 0; c < (scalar ? 1u : 9u); ++c)
                action(state, *prop, qp, c, scalar ? (*scalar)[qp] : (*tensor)[qp](c/3,c%3));
          }
        }
  };
  auto material_values = [&](bool old_only)
  {
    std::vector<Real> values;
    visit([&](unsigned state, const std::string &, unsigned, unsigned, Real & value)
          { if (!old_only || state) values.push_back(value); });
    return values;
  };
  const auto original_material = material_values(false);
  auto accepted_values = [&]()
  {
    auto values = material_values(true);
    for (unsigned state = 1; state < 3; ++state)
      for (auto * system : std::array<SystemBase *, 2>{{&nl, &aux}})
      {
        std::vector<Number> data;
        system->solutionState(state).localize(data);
        values.insert(values.end(), data.begin(), data.end());
      }
    values.push_back(dtOld());
    return values;
  };
  std::map<std::string, bool> enabled_before;
  for (const auto & k : nl.getKernelWarehouse().getObjects()) enabled_before[k->name()] = k->enabled();
  std::map<std::string, std::set<std::string>> groups{
    {"entity", {"TM_all0", "TM_all1"}},
    {"mass_base", {"TM_all_inertia_x", "TM_all_inertia_y"}},
    {"mass_correction", {"ds04g_sp1_mass_0", "ds04g_sp1_mass_1"}},
    {"mass", {"TM_all_inertia_x", "TM_all_inertia_y", "ds04g_sp1_mass_0", "ds04g_sp1_mass_1"}},
    {"gravity_base", {"gravity_y"}}, {"gravity_correction", {"ds04g_sp1_gravity_1"}},
    {"gravity", {"gravity_y", "ds04g_sp1_gravity_1"}},
    {"hg", {"ds04g_sp1_hourglass_0", "ds04g_sp1_hourglass_1"}},
    {"all", {"TM_all0", "TM_all1", "TM_all_inertia_x", "TM_all_inertia_y",
             "ds04g_sp1_mass_0", "ds04g_sp1_mass_1", "gravity_y", "ds04g_sp1_gravity_1",
             "ds04g_sp1_hourglass_0", "ds04g_sp1_hourglass_1"}}};
  if (enabled_before.count("ds04g_all_velocity_damping_x"))
  {
    groups["entity_rate"] = {"ds04g_all_velocity_damping_x", "ds04g_all_velocity_damping_y"};
    groups["entity_total"] = groups["entity"];
    groups["entity_total"].insert(groups["entity_rate"].begin(), groups["entity_rate"].end());
    groups["all"].insert(groups["entity_rate"].begin(), groups["entity_rate"].end());
  }
  std::ofstream vectors("dynamic_residuals.csv"), matrices("dynamic_matrices.csv"),
    histories("dynamic_histories.csv"), cases("dynamic_cases.csv"), activity("dynamic_active_kernels.csv");
  vectors << std::setprecision(17) << "case,term,evaluation,h_index,perturb_dof,sign,row,value\n";
  matrices << std::setprecision(17) << "case,term,row,col,value\n";
  histories << std::setprecision(17) << "case,stage,state,property,qp,component,value\n";
  cases << std::setprecision(17) << "case,h,h_previous,old_before_count,old_unchanged_after_reject,replay_exact,replay_jacobian_exact,status\n";
  activity << "case,term,kernel,enabled,cached_active,use_displaced_mesh\n";
  std::ofstream intervals("dynamic_interval_weights.csv");
  intervals << std::setprecision(17) << "case,stage,kernel,h,h_previous,two_interval,A,B,C,zeta_for_weight_check\n";
  auto record_intervals = [&](unsigned hcase, const std::string & stage)
  {
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      if (k->name() == "TM_all0" || k->name() == "ds04g_sp1_hourglass_0")
      {
        const bool on = k->getParam<bool>("use_previous_accepted_interval");
        const Real zeta = k->name() == "ds04g_sp1_hourglass_0" ? k->getParam<Real>("zeta")
            : std::stod(k->getParam<MaterialPropertyName>("zeta"));
        const auto w = on ? DS04GSP1::twoIntervalHistory(-.05,zeta,dt(),dtOld())
                          : DS04GSP1::history(-.05,zeta,dt());
        intervals << hcase << ',' << stage << ',' << k->name() << ',' << dt() << ','
                  << dtOld() << ',' << on << ',' << w.current << ',' << w.old << ',' << w.older
                  << ',' << zeta << '\n';
      }
  };
  auto select = [&](const std::string & name, unsigned history)
  {
    const auto & wanted = groups.at(name);
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      const_cast<bool &>(k->getParam<bool>("enable")) = wanted.count(k->name());
    updateActiveObjects();
    std::set<std::string> active;
    for (const auto & k : nl.getKernelWarehouse().getActiveObjects()) active.insert(k->name());
    if (active != wanted) mooseError("Actual dynamic audit kernel cache differs from requested group");
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      activity << history << ',' << name << ',' << k->name() << ',' << k->enabled() << ','
               << active.count(k->name()) << ',' << k->getParam<bool>("use_displaced_mesh") << '\n';
  };
  auto evaluate = [&](const NumericVector<Number> & u, unsigned history, const std::string & term,
                      const std::string & label, int hi, int column, int sign)
  {
    FEProblem::computeResidual(u, *residual, 0);
    std::array<Real, 8> result;
    for (unsigned i = 0; i < 8; ++i)
    {
      result[i] = (*residual)(dofs[i]);
      vectors << history << ',' << term << ',' << label << ',' << hi << ',' << column << ','
              << sign << ',' << i << ',' << result[i] << '\n';
    }
    return result;
  };
  const auto & U = getParam<std::vector<Real>>("snapshot_u");
  const auto & V = getParam<std::vector<Real>>("snapshot_v");
  const auto & Q = getParam<std::vector<Real>>("snapshot_a");
  const auto & S = getParam<std::vector<Real>>("snapshot_stress");
  constexpr Real motion_scale = 1e-3, stress_scale = .05, E = 2.97915e10, nu = .2;
  const std::array<Real, 4> hs{{.002,.001,.002,.002}}, hps{{.002,.002,.001,.00025}};
  for (unsigned history = 0; history < 4; ++history)
  {
    const Real h = hs[history], hp = hps[history];
    time() = 2+h; timeOld() = 2; timeStep() = 10; dt() = h; dtOld() = hp;
    std::array<Real, 8> uold, uolder, vold, aold;
    for (unsigned j = 0; j < 8; ++j)
    {
      aold[j] = motion_scale*Q[j]; vold[j] = motion_scale*V[j]-h*aold[j];
      uold[j] = U[j]-h*motion_scale*V[j]+.5*h*h*aold[j];
      uolder[j] = uold[j]-hp*vold[j]+.5*hp*hp*aold[j];
      x->set(dofs[j], U[j]);
      nl.solutionOld().set(dofs[j], uold[j]); nl.solutionOlder().set(dofs[j], uolder[j]);
      aux.solutionOld().set(vdofs[j], vold[j]); aux.solutionOld().set(adofs[j], aold[j]);
      aux.solutionOlder().set(vdofs[j], vold[j]-hp*aold[j]);
      aux.solutionOlder().set(adofs[j], aold[j]);
    }
    x->close(); nl.solutionOld().close(); nl.solutionOlder().close();
    aux.solutionOld().close(); aux.solutionOlder().close();
    RankTwoTensor sold, slope, Fold(RankTwoTensor::initIdentity), Folder(RankTwoTensor::initIdentity);
    sold(0,0)=stress_scale*S[0]; sold(1,1)=stress_scale*S[1]; sold(2,2)=nu*(sold(0,0)+sold(1,1));
    sold(0,1)=sold(1,0)=stress_scale*S[3];
    slope(0,0)=2e5; slope(1,1)=-5e5; slope(2,2)=nu*(slope(0,0)+slope(1,1));
    slope(0,1)=slope(1,0)=3e5;
    for (unsigned i=0;i<4;++i)
      for (unsigned c=0;c<2;++c)
        for (unsigned b=0;b<2;++b)
        {
          Fold(c,b)+=uold[2*i+c]*geom.gradient[i][b];
          Folder(c,b)+=uolder[2*i+c]*geom.gradient[i][b];
        }
    auto seed = [&](unsigned state, const std::string & property, unsigned, unsigned c, Real & value)
    {
      const auto stress = state==2 ? sold-hp*slope : sold;
      const auto F = state==2 ? Folder : Fold;
      value = 0;
      if (property=="stress") value=stress(c/3,c%3);
      else if (property=="deformation_gradient") value=F(c/3,c%3);
      else if (property.find("rotation")!=std::string::npos) value=c/3==c%3;
      else if (property=="elastic_strain" || property=="total_strain" || property=="mechanical_strain")
        value=((1+nu)*stress(c/3,c%3)-(c/3==c%3 ? nu*stress.trace() : 0.))/E;
    };
    visit(seed);
    auto record = [&](const std::string & stage)
    {
      visit([&](unsigned state,const std::string & p,unsigned qp,unsigned c,Real & value)
            { histories << history << ',' << stage << ',' << state << ',' << p << ',' << qp
                        << ',' << c << ',' << value << '\n'; });
    };
    record("seeded");
    const auto before = accepted_values();
    record_intervals(history,"clean");
    for (const auto & group : groups)
    {
      select(group.first, history);
      const auto stem="dynamic_H"+std::to_string(history)+"_"+group.first;
      CDPAssemblyProbe::beginCapture(stem+"_qp.csv",true,true,true);
      try { evaluate(*x,history,group.first,"base",-1,-1,0); }
      catch (...) { CDPAssemblyProbe::endCapture(); throw; }
      CDPAssemblyProbe::endCapture();
      FEProblem::computeJacobian(*x,*J,0);
      for (unsigned i=0;i<8;++i)
        for (unsigned j=0;j<8;++j)
          matrices << history << ',' << group.first << ',' << i << ',' << j << ',' << (*J)(dofs[i],dofs[j]) << '\n';
      for (unsigned hi=0;hi<3;++hi)
        for (unsigned j=0;j<8;++j)
          for (int sign : {-1,1})
          {
            const Real eps=hi==0 ? 1e-6 : hi==1 ? 3e-7 : 1e-7;
            auto trial=x->clone(); trial->add(dofs[j],sign*eps); trial->close();
            evaluate(*trial,history,group.first,"FD",hi,j,sign);
            FEProblem::computeResidual(*x,*residual,0);
          }
      if (accepted_values()!=before) mooseError("Accepted state mutated by frozen FD evaluation");
    }
    select("all",history);
    const auto clean=evaluate(*x,history,"all","clean_before_reject",-1,-1,0);
    FEProblem::computeJacobian(*x,*J,0);
    std::array<std::array<Real,8>,8> clean_matrix;
    for (unsigned i=0;i<8;++i)
      for (unsigned j=0;j<8;++j) clean_matrix[i][j]=(*J)(dofs[i],dofs[j]);
    const auto clean_current=material_values(false);
    // Explicit, bounded use of the same restoreSolutions API as a rejected solve.
    // This does not claim an executioner-driven failed CDP trajectory.
    dt()=2*h; time()=2+2*h;
    record_intervals(history,"rejected_trial");
    auto trial=x->clone();
    for (unsigned j=0;j<8;++j)
      trial->set(dofs[j],uold[j]+2*h*vold[j]+.5*4*h*h*aold[j]);
    trial->close();
    evaluate(*trial,history,"all","intentional_rejected_trial",-1,-1,0);
    try
    {
      throw MooseException("TEST ONLY: reject assembled trial before any state advance");
    }
    catch (const MooseException &) { restoreSolutions(); }
    const bool old_kept=accepted_values()==before;
    record("after_reject_restore");
    dt()=h; time()=2+h;
    record_intervals(history,"replay");
    const auto replay=evaluate(*x,history,"all","clean_after_reject",-1,-1,0);
    const bool exact=clean==replay && clean_current==material_values(false);
    FEProblem::computeJacobian(*x,*J,0);
    bool jexact=true;
    for (unsigned i=0;i<8;++i)
      for (unsigned j=0;j<8;++j) jexact=jexact && clean_matrix[i][j]==(*J)(dofs[i],dofs[j]);
    cases << history << ',' << h << ',' << hp << ',' << before.size() << ',' << old_kept
          << ',' << exact << ',' << jexact << ",complete\n";
    if (!old_kept || !exact || !jexact) mooseError("Bounded reject/restore replay changed frozen state/response");
    record("replayed");
  }
  // Restore the initial diagnostic environment. No test history is promoted.
  for (unsigned state=0;state<3;++state)
  {
    nl.solutionState(state)=*original_nl[state]; aux.solutionState(state)=*original_aux[state];
    nl.solutionState(state).close(); aux.solutionState(state).close();
  }
  std::size_t index=0;
  visit([&](unsigned,const std::string &,unsigned,unsigned,Real & value)
        { value=original_material.at(index++); });
  for (const auto & k : nl.getKernelWarehouse().getObjects())
    const_cast<bool &>(k->getParam<bool>("enable"))=enabled_before.at(k->name());
  updateActiveObjects();
  time()=original_time; timeOld()=original_oldtime; timeStep()=original_step; dt()=original_dt; dtOld()=original_dtold;
  FEProblem::computeResidual(*original,*residual,0);
  if (!mapping || !vectors || !matrices || !histories || !cases || !activity || !intervals)
    mooseError("Could not write complete dynamic audit evidence");
  mooseInfo("DS04G dynamic audit complete: four prescribed histories, zero accepted steps");
}
