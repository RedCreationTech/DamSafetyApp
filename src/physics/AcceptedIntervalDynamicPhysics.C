#include "AcceptedIntervalDynamicPhysics.h"
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "meta_action");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_variable");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_material");
registerMooseAction("DamSafetyApp", AcceptedIntervalDynamicPhysics, "add_master_action_material");
