#include "CDPQuasiStaticPhysics.h"
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "meta_action");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_variable");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_material");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_master_action_material");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_scalar_kernel");
registerMooseAction("DamSafetyApp", CDPQuasiStaticPhysics, "add_user_object");
