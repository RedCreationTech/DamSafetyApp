// Build this unchanged driver against the frozen and candidate headers/libraries.
// It performs local unit-level calls only, never an FE solve or a LIMS Job.
#include "AbaqusCDPStateIntegrator.h"
#include "AbaqusCDPSubstepIntegrator.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
  constexpr double E = 3.04e10, nu = .2;
  const std::string p = "test/tests/cdp_material_table/data/";
  CDPMaterialTable table(p + "compression_hardening.csv", p + "compression_damage.csv",
                        p + "tension_stiffening.csv", p + "tension_damage.csv", E);
  const AbaqusCDPLocalIntegrator local(table, {E, nu, 36.31, .1, 1.16, .667, 40, 1e-9, 1e-7, 1e-6});
  std::cout << std::setprecision(17);
  const double ft = table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.).stress.value;
  const double fc = table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.).stress.value;
  for (double mu : {0., .0005})
  {
    const AbaqusCDPStateIntegrator core(local, {1., 1., mu, 1e-12});
    AbaqusCDPStateIntegrator::State old;
    for (double s : {.5*ft, 1.05*ft, .4*ft, -.5*fc, -1.05*fc, -.2*fc, .3*ft})
    {
      const auto result = core.integrateLinearized({s/E, -nu*s/E, -nu*s/E, 0., 0., 0.}, .00025, old);
      for (double value : result.result.cauchy_stress) std::cout << value << ',';
      for (double value : result.result.state.backbone.plastic_strain) std::cout << value << ',';
      std::cout << result.result.state.backbone.tensile_equivalent_plastic_strain << ','
                << result.result.state.backbone.compressive_equivalent_plastic_strain << ',';
      for (double value : result.result.state.viscous_plastic_strain) std::cout << value << ',';
      std::cout << result.result.state.viscous_tension_damage << ','
                << result.result.state.viscous_compression_damage << ',';
      for (unsigned int col = 0; col < 22; ++col)
        for (unsigned int row = 0; row < 22; ++row)
          std::cout << result.derivative[col][row] << ',';
      std::cout << result.result.damage.stiffness_factor << '\n';
      old = result.result.state;
    }
    // Preserve the original failing centered-direction case as a diagnostic.
    for (double s : {1.05*ft, -1.05*fc})
    {
      const AbaqusCDPStateIntegrator::SymmetricTensor target = {s/E, -nu*s/E, -nu*s/E, 0., 0., 0.};
      const AbaqusCDPSubstepIntegrator steps(core, {16, std::abs(target[0])/3., 1e-8});
      const auto linear = steps.integrateLinearized({}, target, .001, {});
      const AbaqusCDPStateIntegrator::SymmetricTensor direction = {.2,.3,-.4,0.,.05,.03};
      const auto a = steps.applyTangent(linear.algorithmic_tangent, direction);
      const auto b = steps.directionalDerivative({}, target, .001, {}, direction, 1e-9);
      std::cerr << "uniaxial-diagnostic mu=" << mu << " stress=" << s;
      for (unsigned int k = 0; k < 6; ++k) std::cerr << ',' << a[k] << ',' << b[k];
      std::cerr << '\n';
    }
  }
}
