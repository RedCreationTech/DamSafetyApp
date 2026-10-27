#include "CDPAssemblyProbe.h"
#include "KernelBase.h"
#include "DS04GSP1AuditProblem.h"
#include "MaterialPropertyStorage.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "RankTwoTensor.h"
#include "libmesh/elem.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", DS04GSP1AuditProblem);
InputParameters DS04GSP1AuditProblem::validParams()
{
  auto p = FEProblem::validParams();
  p.addClassDescription("TEST ONLY: frozen serial one-QUAD4 actual assembly; no time advance");
  p.addParam<bool>("audit_kinematics", false, "Capture same-call F/Fold/rotation and true material tangent");
  p.addParam<bool>("refine_compression", false, "Add smaller compression FD scales within N64");
  p.addParam<bool>("refine_nonlinear", false, "Add smaller affine/tension/compression FD scales");
  return p;
}
DS04GSP1AuditProblem::DS04GSP1AuditProblem(const InputParameters &p) : FEProblem(p) {}

void DS04GSP1AuditProblem::initialSetup()
{
  FEProblem::initialSetup();
  auto &nl = getNonlinearSystem(0);
  auto &sys = dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  const auto &mesh = this->mesh().getMesh();
  if (comm().size() != 1 || libMesh::n_threads() != 1 || mesh.n_elem() != 1 ||
      mesh.n_nodes() != 4 || sys.n_dofs() != 8 || mesh.elem_ptr(0)->type() != QUAD4)
    mooseError("DS04GSP1AuditProblem is restricted to serial one-QUAD4/eight-DOF diagnostics");
  const auto *elem = mesh.elem_ptr(0);
  const auto *original = nl.currentSolution();
  auto x = original->clone();
  auto residual = x->zero_clone();
  // At INITIAL the live matrix has no populated sparsity to clone. Allocate
  // the complete eight-DOF audit pattern without changing the solver matrix.
  auto J = SparseMatrix<Number>::build(comm());
  J->init(8, 8, 8, 8, 8, 0);
  time() = .002;
  dt() = .002;
  dtOld() = .002;
  executeControls(EXEC_TIMESTEP_BEGIN);
  // Normal Transient initialization refreshes this warehouse after INITIAL.
  // This audit runs inside initialSetup, so perform that lifecycle step here.
  std::ofstream active("audit_active_kernels.csv");
  active << "stage,kernel,enabled,cached_active\n";
  auto record_active = [&](const std::string & stage)
  {
    std::set<std::string> names;
    for (const auto & k : nl.getKernelWarehouse().getActiveObjects()) names.insert(k->name());
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      active << stage << ',' << k->name() << ',' << k->enabled() << ',' << names.count(k->name()) << '\n';
    return names;
  };
  record_active("before_refresh");
  updateActiveObjects();
  const auto active_names = record_active("after_refresh");
  for (const auto & name : active_names)
    if (name.find("TM_all") != std::string::npos || name.find("gravity") != std::string::npos)
      mooseError("Frozen internal-stiffness audit must exclude dynamic/gravity kernels");
  active.close();
  std::array<std::array<dof_id_type, 2>, 4> dofs;
  std::array<double, 2> center{};
  std::ofstream mapping("audit_dof_map.csv");
  mapping << std::setprecision(17) << "local_node,node_id,x,y,dof_x,dof_y\n";
  for (unsigned i = 0; i < 4; ++i)
  {
    const auto *node = elem->node_ptr(i);
    for (unsigned c = 0; c < 2; ++c)
    {
      dofs[i][c] = node->dof_number(sys.number(), sys.variable_number(c ? "disp_y" : "disp_x"), 0);
      center[c] += (*node)(c) / 4.;
    }
    mapping << i << ',' << node->id() << ',' << (*node)(0) << ',' << (*node)(1) << ',' << dofs[i][0]
            << ',' << dofs[i][1] << '\n';
  }
  mapping.close();
  auto history = [&](const std::string &name)
  {
    std::ofstream f(name);
    f << std::setprecision(17) << "element,side,state,property,qp,component,value\n";
    const auto &storage = getMaterialPropertyStorage();
    for (auto state : storage.stateIndexRange())
      for (const auto &entry : storage.props(state))
        for (const auto &side : entry.second)
          for (unsigned id = 0; id < side.second.size(); ++id)
          {
            const auto *value = side.second.queryValue(id);
            if (!value)
              continue;
            const auto prop = storage.queryStatefulPropName(storage.statefulProps().at(id));
            const auto *scalar = dynamic_cast<const MaterialProperty<Real> *>(value);
            const auto *tensor = dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
            if (!prop || (!scalar && !tensor))
              mooseError("Unsupported audit state property");
            for (unsigned qp = 0; qp < value->size(); ++qp)
              for (unsigned c = 0; c < (scalar ? 1u : 9u); ++c)
                f << entry.first->id() << ',' << side.first << ',' << state << ',' << *prop << ','
                  << qp << ',' << c << ',' << (scalar ? (*scalar)[qp] : (*tensor)[qp](c / 3, c % 3))
                  << '\n';
          }
    if (!f)
      mooseError("Cannot write audit state");
  };
  auto evaluate = [&](const NumericVector<Number> &u, const std::string &name)
  {
    CDPAssemblyProbe::beginCapture(name + "_qp.csv", false, true,
                                   getParam<bool>("audit_kinematics"));
    try
    {
      FEProblem::computeResidual(u, *residual, 0);
    }
    catch (...)
    {
      CDPAssemblyProbe::endCapture();
      throw;
    }
    CDPAssemblyProbe::endCapture();
    residual->print_matlab(name + "_r.m");
  };
  history("audit_history_before.csv");
  std::ofstream cases("audit_cases.csv");
  cases << "case,status\n";
  cases.flush();
  for (const auto &name : {"zero", "micro", "affine", "rigid", "tension", "compression"})
  {
    const std::string stem = std::string("audit_") + name;
    x->zero();
    for (unsigned i = 0; i < 4; ++i)
    {
      const Real X = elem->point(i)(0) - center[0], Y = elem->point(i)(1) - center[1];
      Real u = 0, v = 0;
      if (std::string(name) == "micro")
      {
        u = (i % 2 ? -1 : 1) * 1e-7;
        v = u * .3;
      }
      if (std::string(name) == "affine")
      {
        u = 1e-4 * X + 2e-4 * Y;
        v = -1e-4 * Y;
      }
      if (std::string(name) == "rigid")
      {
        u = (std::cos(.3) - 1) * X - std::sin(.3) * Y + .02;
        v = std::sin(.3) * X + (std::cos(.3) - 1) * Y - .03;
      }
      if (std::string(name) == "tension")
      {
        u = 2e-4 * X;
        v = -4e-5 * Y;
      }
      if (std::string(name) == "compression")
      {
        u = -8e-4 * X;
        v = 1.6e-4 * Y;
      }
      x->set(dofs[i][0], u);
      x->set(dofs[i][1], v);
    }
    x->close();
    try
    {
      x->print_matlab(stem + "_u.m");
      evaluate(*x, stem + "_base");
      if (getParam<bool>("audit_kinematics"))
        CDPAssemblyProbe::beginCapture(stem + "_jac_material.csv", true, true, true);
      try
      {
        FEProblem::computeJacobian(*x, *J, 0);
      }
      catch (...)
      {
        if (getParam<bool>("audit_kinematics"))
          CDPAssemblyProbe::endCapture();
        throw;
      }
      if (getParam<bool>("audit_kinematics"))
        CDPAssemblyProbe::endCapture();
      J->print_matlab(stem + "_J.m");
      for (unsigned i = 0; i < 4; ++i)
        for (unsigned c = 0; c < 2; ++c)
          for (unsigned k = 0; k < ((getParam<bool>("refine_compression") &&
                                      std::string(name) == "compression") ||
                                     (getParam<bool>("refine_nonlinear") &&
                                      (std::string(name) == "affine" || std::string(name) == "tension" ||
                                       std::string(name) == "compression")) ? 5u : 3u); ++k)
          {
            const Real h = k == 3 ? 3e-8 : k == 4 ? 1e-8 :
                                                     std::pow(10., -5 - static_cast<int>(k));
            for (const Real sign : {1., -1.})
            {
              auto trial = x->clone();
              trial->add(dofs[i][c], sign * h);
              trial->close();
              evaluate(*trial, stem + "_d" + std::to_string(2 * i + c) + "_h" + std::to_string(k) +
                                   (sign > 0 ? "_plus" : "_minus"));
              // Restore the current-solution pointer before destroying the scratch vector.
              FEProblem::computeResidual(*x, *residual, 0);
            }
          }
      evaluate(*x, stem + "_replay");
      cases << name << ",complete\n";
      cases.flush();
    }
    catch (const std::exception &e)
    {
      cases << name << ",failed\n";
      std::ofstream error(stem + "_failure.txt");
      error << e.what() << '\n';
    }
    FEProblem::computeResidual(*original, *residual, 0);
  }
  history("audit_history_after.csv");
  if (!mapping || !cases)
    mooseError("Cannot write element audit index");
  time() = 0.;
}
