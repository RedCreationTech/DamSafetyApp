import hashlib
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('prepare', ROOT / 'tools/diagnostics/prepare_ds04g_2597_observer.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class ObserverPreparationTest(unittest.TestCase):
    def test_only_inserts_observer_and_subclassed_stepper(self):
        source = b'[Executioner]\n  [TimeStepper]\n    type = TimeSequenceStepper\n    time_sequence = 1\n  []\n[]\n'
        with patch.object(prepare, 'Q0_SHA', hashlib.sha256(source).hexdigest()):
            observed = prepare.prepare(source).decode()
        self.assertEqual(observed.replace(prepare.PROBLEM, '').replace('type = CDPObservationTimeStepper', 'type = TimeSequenceStepper'), source.decode())
        self.assertIn('observation_start = 3.018', observed)
        self.assertIn('observation_end = 3.202', observed)
        self.assertIn("observation_labels = '2597'", observed)

    def test_reject_unlocked_input(self):
        with self.assertRaises(ValueError):
            prepare.prepare(b'untrusted input')

    def test_hooks_are_guarded_and_no_monitor_is_installed(self):
        problem = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        self.assertNotIn('SNESMonitorSet', problem)
        self.assertIn('qp != 0', problem)
        local = (ROOT / 'src/utils/AbaqusCDPLocalIntegrator.C').read_text()
        self.assertIn('"trial_stress", trial_stress', local)
        self.assertIn('"ft_fc_scale_f_threshold"', local)
        material = (ROOT / 'src/materials/AbaqusCDPStressUpdate.C').read_text()
        self.assertLess(material.index('"raw_increment", increment'), material.index('strain_increment -= inelastic_strain_increment'))
        self.assertIn('if (CDPObservation::active())', material)

    def test_trial_hook_is_on_the_actual_linearized_call_chain(self):
        local = (ROOT / 'src/utils/AbaqusCDPLocalIntegrator.C').read_text()
        state = (ROOT / 'src/utils/AbaqusCDPStateIntegrator.C').read_text()
        sub = (ROOT / 'src/utils/AbaqusCDPSubstepIntegrator.C').read_text()
        material = (ROOT / 'src/materials/AbaqusCDPStressUpdate.C').read_text()
        self.assertIn('_substep_integrator.integrateLinearized(old_total_strain, new_total_strain, _dt, old_state)', material)
        self.assertIn('_state_integrator.integrateLinearized(', sub)
        self.assertIn('_backbone_integrator.integrateLinearized(total_strain, old_state.backbone)', state)
        self.assertIn('LinearizedResult linearized{integrate(total_strain, old_state), {}};', local)
        trial = local.index('"trial_stress", trial_stress')
        self.assertLess(local.index('const auto trial_stress = elasticStress', local.index('AbaqusCDPLocalIntegrator::integrate(')), trial)
        self.assertLess(trial, local.index('if (trial_yield <= _parameters.residual_tolerance * stress_scale)', trial))
        self.assertIn('CDPObservation::context.local_solve = ++CDPObservation::sequence', local)
        self.assertIn('"partition_failed"', sub)
        self.assertIn('"partition_succeeded"', sub)

    def test_residual_observer_does_not_add_norm_or_monitor(self):
        problem = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        residual = problem.split('void CDPObservationProblem::computeResidual(', 1)[1].split('void CDPObservationProblem::computeJacobian(', 1)[0]
        self.assertNotIn('l2_norm(', residual)
        self.assertNotIn('SNESMonitorSet', problem)
        self.assertNotIn('comm().', residual)
        self.assertIn('if(!_enabled) { FEProblem::computeResidual(x,r,n); return; }', residual)
        self.assertIn('vector identity is not alone proof of line-search acceptance', problem)
        stepper = (ROOT / 'src/timesteppers/CDPObservationTimeStepper.C').read_text()
        self.assertLess(stepper.index('attemptEvent("accepted_end")'), stepper.index('TimeSequenceStepper::acceptStep()'))


if __name__ == '__main__':
    unittest.main()
