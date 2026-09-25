#!/usr/bin/env python3
"""Apply the dam-2b weak-plane-stress x B-bar Jacobian fix to a pristine MOOSE 4bce02d tree.

Target: MOOSE modules/solid_mechanics kernels
  - StressDivergenceTensors.{h,C}: add volume-averaged shape-function values and the missing
    B-bar correction terms in the off-diagonal Jacobian w.r.t. the out-of-plane strain variable.
  - WeakPlaneStress.{h,C}: new `volumetric_locking_correction` flag (forwarded automatically by
    the SolidMechanics Physics action), volume-averaged shape-function helpers, and the missing
    B-bar correction terms in the diagonal and displacement off-diagonal Jacobians.

Root cause (D2EQ 009 stall / Q2): Compute2DSmallStrain replaces the strain trace with its
element volume average when volumetric_locking_correction is active. The true derivative of
stress w.r.t. the out-of-plane strain variable therefore contains (avg - local)/3 type terms.
StressDivergenceTensors' main diagonal block carries the matching 4-term B-bar expansion, but
the strain_zz-coupled Jacobian entries (StressDivergence off-diagonal, WeakPlaneStress diagonal
and off-diagonals) did not, making the assembled Jacobian inconsistent with the residual
(measured ||J-Jfd||_F/||J||_F ~ 3.5% unscaled, ~20% with automatic_scaling; PJFNK unaffected).

Usage: python3 apply_moose_zz_bbar_fix.py /path/to/moose
Idempotent: skips files already patched.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(sys.argv[1]) / "modules/solid_mechanics"
K = ROOT / "src/kernels"
H = ROOT / "include/kernels"


def patch_file(path: Path, replacements: list[tuple[str, str]]) -> None:
    text = path.read_text()
    if all(new in text for _, new in replacements):
        print(f"already patched: {path}")
        return
    for old, new in replacements:
        if old not in text:
            raise SystemExit(f"anchor not found in {path}:\n{old[:200]}")
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique in {path}:\n{old[:200]}")
        text = text.replace(old, new)
    path.write_text(text)
    print(f"patched: {path}")


# ---------------- WeakPlaneStress.h ----------------
patch_file(H / "WeakPlaneStress.h", [
    ("""  /// d(strain)/d(temperature), if computed by ComputeThermalExpansionEigenstrain
  std::vector<const MaterialProperty<RankTwoTensor> *> _deigenstrain_dT;
};""",
     """  /// d(strain)/d(temperature), if computed by ComputeThermalExpansionEigenstrain
  std::vector<const MaterialProperty<RankTwoTensor> *> _deigenstrain_dT;

  /// Volumetric locking (B-bar) correction flag; must match the strain calculator material
  const bool _volumetric_locking_correction;

  /// Volume-averaged shape function values (per trial function) for the B-bar Jacobian terms
  std::vector<Real> _avg_phi;

  /// Volume-averaged shape function gradients for the B-bar Jacobian terms
  std::vector<std::vector<Real>> _avg_grad_phi;

  /// Compute volume-averaged shape function values over the current element
  void computeAveragePhi();

  /// Compute volume-averaged shape function gradients over the current element
  void computeAverageGradientPhi();

  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;
};"""),
])

# ---------------- WeakPlaneStress.C ----------------
patch_file(K / "WeakPlaneStress.C", [
    ("""  params.set<bool>("use_displaced_mesh") = false;

  return params;""",
     """  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>("volumetric_locking_correction",
                        false,
                        "Set to true when the strain calculator applies the volumetric locking "
                        "correction, so that the Jacobian includes the matching B-bar terms");

  return params;"""),
    ("""    _disp_coupled(isCoupled("displacements")),
    _ndisp(_disp_coupled ? coupledComponents("displacements") : 0),""",
     """    _disp_coupled(isCoupled("displacements")),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _ndisp(_disp_coupled ? coupledComponents("displacements") : 0),"""),
    ("""Real
WeakPlaneStress::computeQpJacobian()
{
  return _Jacobian_mult[_qp](_direction, _direction, _direction, _direction) * _test[_i][_qp] *
         _phi[_j][_qp];
}""",
     """void
WeakPlaneStress::computeAveragePhi()
{
  _avg_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_phi[_i] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_phi[_i] += _phi[_i][_qp] * _JxW[_qp] * _coord[_qp];
    _avg_phi[_i] /= _current_elem_volume;
  }
}

void
WeakPlaneStress::computeAverageGradientPhi()
{
  _avg_grad_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_grad_phi[_i].resize(3);
    for (unsigned int component = 0; component < 3; ++component)
    {
      _avg_grad_phi[_i][component] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_phi[_i][component] += _grad_phi[_i][_qp](component) * _JxW[_qp] * _coord[_qp];
      _avg_grad_phi[_i][component] /= _current_elem_volume;
    }
  }
}

void
WeakPlaneStress::computeJacobian()
{
  if (_volumetric_locking_correction)
    computeAveragePhi();
  Kernel::computeJacobian();
}

void
WeakPlaneStress::computeOffDiagJacobian(unsigned int jvar)
{
  if (_volumetric_locking_correction)
  {
    computeAveragePhi();
    computeAverageGradientPhi();
  }
  Kernel::computeOffDiagJacobian(jvar);
}

Real
WeakPlaneStress::computeQpJacobian()
{
  Real jacobian = _Jacobian_mult[_qp](_direction, _direction, _direction, _direction) *
                  _test[_i][_qp] * _phi[_j][_qp];

  // B-bar correction: the strain calculator replaces the strain trace with its element volume
  // average, so d(sigma_dd)/d(strain_zz) also contains the tangent trace-sum times the
  // difference between the volume-averaged and the local shape function value.
  if (_volumetric_locking_correction)
  {
    Real sum_C3x1 = 0.0;
    for (unsigned int k = 0; k < 3; ++k)
      sum_C3x1 += _Jacobian_mult[_qp](_direction, _direction, k, k);
    jacobian += sum_C3x1 * _test[_i][_qp] * (_avg_phi[_j] - _phi[_j][_qp]) / 3.0;
  }

  return jacobian;
}"""),
    ("""        val = _Jacobian_mult[_qp](
                  _direction, _direction, coupled_direction_index, coupled_direction_index) *
              _test[_i][_qp] * _grad_phi[_j][_qp](coupled_direction_index);
      }""",
     """        val = _Jacobian_mult[_qp](
                  _direction, _direction, coupled_direction_index, coupled_direction_index) *
              _test[_i][_qp] * _grad_phi[_j][_qp](coupled_direction_index);

        // B-bar correction for the displacement off-diagonal (strain trace volume-averaged by
        // the strain calculator).
        if (_volumetric_locking_correction)
        {
          Real sum_C3x1 = 0.0;
          for (unsigned int k = 0; k < 3; ++k)
            sum_C3x1 += _Jacobian_mult[_qp](_direction, _direction, k, k);
          val += sum_C3x1 * _test[_i][_qp] *
                 (_avg_grad_phi[_j][coupled_direction_index] -
                  _grad_phi[_j][_qp](coupled_direction_index)) / 3.0;
        }
      }"""),
])

# ---------------- StressDivergenceTensors.h: add _avg_phi + helper ----------------
sdt_h = H / "StressDivergenceTensors.h"
text = sdt_h.read_text()
if "_avg_phi" in text:
    print(f"already patched: {sdt_h}")
else:
    anchor = "_avg_grad_phi;"
    if text.count(anchor) != 1:
        raise SystemExit("StressDivergenceTensors.h anchor not unique")
    text = text.replace(anchor, anchor + """

  /// Volume-averaged shape function values (per trial function) for the B-bar Jacobian terms
  std::vector<Real> _avg_phi;""")
    anchor2 = "void computeAverageGradientPhi();"
    if text.count(anchor2) != 1:
        raise SystemExit("computeAverageGradientPhi declaration anchor not unique")
    text = text.replace(anchor2, anchor2 + """

  /// Compute volume-averaged shape function values over the current element
  void computeAveragePhi();""")
    sdt_h.write_text(text)
    print(f"patched: {sdt_h}")

# ---------------- StressDivergenceTensors.C ----------------
patch_file(K / "StressDivergenceTensors.C", [
    ("""void
StressDivergenceTensors::computeAverageGradientPhi()
{""",
     """void
StressDivergenceTensors::computeAveragePhi()
{
  _avg_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_phi[_i] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_phi[_i] += _phi[_i][_qp] * _JxW[_qp] * _coord[_qp];
    _avg_phi[_i] /= _current_elem_volume;
  }
}

void
StressDivergenceTensors::computeAverageGradientPhi()
{"""),
    ("""  // off-diagonal Jacobian with respect to a coupled out_of_plane_strain variable
  if (_out_of_plane_strain_coupled && jvar == _out_of_plane_strain_var)
    return _Jacobian_mult[_qp](
               _component, _component, _out_of_plane_direction, _out_of_plane_direction) *
           _grad_test[_i][_qp](_component) * _phi[_j][_qp];""",
     """  // off-diagonal Jacobian with respect to a coupled out_of_plane_strain variable
  if (_out_of_plane_strain_coupled && jvar == _out_of_plane_strain_var)
  {
    Real jacobian = _Jacobian_mult[_qp](
                        _component, _component, _out_of_plane_direction, _out_of_plane_direction) *
                    _grad_test[_i][_qp](_component) * _phi[_j][_qp];

    // B-bar correction: the strain calculator volume-averages the strain trace, so the
    // derivative w.r.t. the out-of-plane strain variable contains (avg - local)/3 type terms on
    // both the trial side and the B-bar test-correction side.
    if (_volumetric_locking_correction)
    {
      Real sum_C3x3 = _Jacobian_mult[_qp].sum3x3();
      RealGradient sum_C3x1 = _Jacobian_mult[_qp].sum3x1();

      // trial-side trace correction
      jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                  (_avg_phi[_j] - _phi[_j][_qp]) / 3.0;

      // derivative of the B-bar test correction (tr(sigma)/3 * (avg_grad_test - grad_test))
      const Real d_tr_sigma = sum_C3x1(_out_of_plane_direction) * _phi[_j][_qp] +
                              sum_C3x3 * (_avg_phi[_j] - _phi[_j][_qp]) / 3.0;
      jacobian += d_tr_sigma * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) /
                  3.0;
    }

    return jacobian;
  }"""),
    ("""void
StressDivergenceTensors::computeOffDiagJacobian(const unsigned int jvar)
{
  if (_volumetric_locking_correction)
  {
    computeAverageGradientPhi();
    computeAverageGradientTest();
  }""",
     """void
StressDivergenceTensors::computeOffDiagJacobian(const unsigned int jvar)
{
  if (_volumetric_locking_correction)
  {
    computeAverageGradientPhi();
    computeAverageGradientTest();
    if (jvar == _out_of_plane_strain_var)
      computeAveragePhi();
  }"""),
])

print("MOOSE weak-plane-stress x B-bar Jacobian patch applied.")
