from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


class AcceptedOutputObserverSourceTest(unittest.TestCase):
    def test_default_off_paths_call_only_the_base_implementation(self):
        source = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        self.assertIn('"observation_enabled", false', source)
        self.assertIn('"trace_state_swap", false', source)
        for name, args in (
            ('initialSetup', ''),
            ('advanceState', ''),
            ('execute', 'exec_type'),
            ('outputStep', 'type'),
        ):
            body = source.split('CDPObservationProblem::' + name + '(', 1)[1].split('\n}\n', 1)[0]
            self.assertIn('if (!_enabled)', body)
            self.assertIn('FEProblem::' + name + '(' + args + ');', body)
        for name in ('computeResidual', 'computeJacobian'):
            body = source.split('CDPObservationProblem::' + name + '(', 1)[1].split('\n}\n', 1)[0]
            self.assertIn('if (!_enabled)', body)
            self.assertIn('FEProblem::' + name + '(soln,', body)

    def test_evaluation_context_is_raii_restored_on_all_paths(self):
        source = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        for name in ('computeResidual', 'computeJacobian'):
            body = source.split('CDPObservationProblem::' + name + '(', 1)[1].split('\n}\n', 1)[0]
            self.assertIn('CDPObservationFEContextGuard restore', body)
            self.assertIn('catch (...)', body)
            self.assertIn('throw;', body)
        guard = (ROOT / 'include/utils/CDPObservationContextGuard.h').read_text()
        self.assertIn('~CDPObservationFEContextGuard() noexcept', guard)
        self.assertIn('_evaluation = _old_evaluation;', guard)
        self.assertIn('_context = _old_context;', guard)

    def test_constant_dt_forwards_attempt_accept_and_reject(self):
        header = (ROOT / 'include/timesteppers/CDPObservationConstantDT.h').read_text()
        source = (ROOT / 'src/timesteppers/CDPObservationConstantDT.C').read_text()
        self.assertIn('public ConstantDT', header)
        self.assertIn('attemptEvent("attempt_enter")', source)
        self.assertIn('ConstantDT::step()', source)
        self.assertIn('attemptEvent("reject_before")', source)
        self.assertIn('ConstantDT::rejectStep()', source)
        self.assertIn('attemptEvent("reject_after")', source)
        self.assertIn('attemptEvent("accepted_end")', source)
        self.assertIn('ConstantDT::acceptStep()', source)

    def test_material_records_full_state_rotation_and_exception_without_retry(self):
        source = (ROOT / 'src/materials/AbaqusCDPStressUpdate.C').read_text()
        update = source.split('AbaqusCDPStressUpdate::updateState(', 1)[1]
        self.assertIn('"old16", old16', update)
        self.assertIn('"new16", new16', update)
        self.assertIn('rotation[row][column] = rotation_increment(row, column)', update)
        self.assertIn('"rotation_row_major", rotation', update)
        self.assertIn('final.backbone.residual_norm', update)
        self.assertEqual(update.count('_substep_integrator.integrateLinearized('), 1)
        catch = update.split('catch (const std::exception & error)', 1)[1]
        self.assertLess(catch.index('"material_exception"'), catch.index('throw MooseException(error.what())'))
        self.assertNotIn('integrateLinearized', catch)

    def test_observer_does_not_add_solver_or_material_evaluations(self):
        source = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        self.assertNotIn('SNESMonitorSet', source)
        self.assertNotIn('computeResidual(soln, residual, sys_num);\n    FEProblem::computeResidual', source)
        stepper = (ROOT / 'src/timesteppers/CDPObservationConstantDT.C').read_text()
        self.assertNotIn('computeResidual', stepper)
        self.assertNotIn('computeJacobian', stepper)
        self.assertNotIn('integrateLinearized', stepper)

    def test_residual_preparation_failure_cannot_skip_or_repeat_base(self):
        source = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        body = source.split('CDPObservationProblem::computeResidual(', 1)[1].split('\n}\n', 1)[0]
        preparation = body.split('if (window())', 1)[1].split('FEProblem::computeResidual', 1)[0]
        self.assertIn('try', preparation)
        self.assertIn('catch (...)', preparation)
        self.assertIn('_history_incomplete = true;', preparation)
        self.assertNotIn('FEProblem::computeResidual', preparation)
        self.assertEqual(body.count('FEProblem::computeResidual(soln, residual, sys_num);'), 2)

    def test_csv_records_quote_fields_and_consumer_handles_escaped_quotes(self):
        problem = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        self.assertIn('quoted(name)', problem)
        self.assertIn('quoted(encoded_values.str())', problem)
        self.assertIn('quoted(note)', problem)
        consumer = (ROOT / 'src/timesteppers/DS04GAcceptedIntervalTestStepper.C').read_text()
        self.assertIn('readCsvRow', consumer)
        self.assertIn("input.peek() == '\"'", consumer)
        self.assertIn("field += '\"'", consumer)
        self.assertNotIn('std::getline(stream, column', consumer)
        self.assertNotIn('const auto row = columns(', consumer)

    def test_byte_cap_includes_headers_and_all_event_rows(self):
        source = (ROOT / 'src/problems/CDPObservationProblem.C').read_text()
        writer = source.split('CDPObservationProblem::write(', 1)[1].split('\n}\n', 1)[0]
        self.assertNotIn('const bool light', writer)
        self.assertIn('header.size() > _limit - _bytes', writer)
        self.assertIn('text.size() > _limit - _bytes - header.size()', writer)
        self.assertIn('file << header << text;', writer)
        self.assertIn('_bytes += header.size() + text.size();', writer)
        self.assertLess(writer.index('const auto text = row.str()'), writer.index('_truncated = true;'))
        self.assertLess(writer.index('_truncated = true;'), writer.index('file << header << text;'))

    def test_real_cdp_fixture_covers_success_reject_swap_accept_and_output(self):
        fixture = (ROOT / 'test/tests/ds04g_sp1_dynamic/cdp_accepted_output_observer.i').read_text()
        self.assertIn('type = AbaqusCDPStressUpdate', fixture)
        self.assertIn('trace_state_swap = true', fixture)
        self.assertIn('check_material_observation = true', fixture)
        self.assertIn('type = MaterialRealAux', fixture)
        self.assertIn('maximum_strain_increment = 2.5e-5', fixture)
        checker = (ROOT / 'src/timesteppers/DS04GAcceptedIntervalTestStepper.C').read_text()
        for event in ('material_enter', 'material_success', 'trial', 'partition_succeeded',
                      'swapback_after', 'swapback_candidate', 'material_swapback_after'):
            self.assertIn(event, checker)
        readme = (ROOT / 'test/tests/ds04g_sp1_dynamic/README.md').read_text()
        self.assertIn('未可靠注入真实材料异常', readme)
        self.assertIn('完整OFF/OFF与OFF/ON数值对照也未在此声称覆盖', readme)


if __name__ == '__main__':
    unittest.main()
