#pragma once

#include "ComputeFiniteStrainElasticStress.h"

/**
 * Explicit opt-in, isotropic elastic diagnostic for Hughes-Winget kinematics.
 *
 * The Hughes-Winget increment is expressed in the spatial midpoint basis. This
 * research material rotates the old stress and then adds C:D; it is NOT a CDP
 * model, does not transport plastic history and does not supply an exact finite
 * strain Jacobian. It must not be used as a production constitutive model.
 */
class DS04GHWSpatialElasticStress : public ComputeFiniteStrainElasticStress
{
public:
  static InputParameters validParams();
  DS04GHWSpatialElasticStress(const InputParameters & parameters);
  void initialSetup() override;

protected:
  void computeQpStress() override;

private:
  // Requiring these named properties catches the standard Taylor-only configuration;
  // their presence alone cannot establish their producer or configuration.
  const MaterialProperty<RankTwoTensor> & _hw_midpoint_gradient;
  const MaterialProperty<RankTwoTensor> & _hw_displacement_increment_gradient;
};
