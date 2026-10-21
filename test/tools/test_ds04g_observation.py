"""Synthetic tooling checks; these fixtures are never physical/recovery evidence."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

TOOLS = Path(__file__).resolve().parents[2] / 'tools/diagnostics'


def load(name):
    spec = importlib.util.spec_from_file_location(name, TOOLS / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = load('prepare_ds04g_observation')
auditor = load('audit_ds04g_observation')
LABELS = [833, 1312, 959]


class PreparationTests(unittest.TestCase):
    def test_exact_source_times_and_physics_survive_observation_only_edit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'parent.i'
            original = ("[Executioner]\n  end_time = 2.720000\n  time_sequence = '0 1 1.000500 2.720000'\n"
                        "  type = CDPObservationTimeStepper\n[]\n"
                        "[Physics]\n  enable_finite_newmark_damping = false\n[]\n"
                        "[Problem]\n  type = CDPObservationProblem\n  observation_enabled = false\n"
                        "  observation_start = 2.698\n  observation_end = 2.720\n"
                        "  observation_labels = '1881 767'\n  observation_vertices = '0 0'\n[]\n")
            source.write_text(original)
            nodes = root / 'nodes.csv'; elements = root / 'elements.csv'
            nodes.write_text('node_label,X_mm,Y_mm\n1,0,0\n2,1000,0\n3,1000,1000\n4,0,1000\n')
            with elements.open('w', newline='') as f:
                w = csv.writer(f);w.writerow(['element_label', 'connectivity'])
                for label in LABELS:w.writerow([label, '[1,2,3,4]'])
            out = root / 'out'
            result = prepare.prepare(source, nodes, elements, out, labels=LABELS,
                source_sha=hashlib.sha256(source.read_bytes()).hexdigest(), preserve_source=True,
                observation_start=2.68, observation_end=2.72, byte_limit=4*1024**3)
            self.assertEqual((out/'ds04g_prefix_control.i').read_bytes(), source.read_bytes())
            observed = (out/'ds04g_observe_eq1p720.i').read_text()
            self.assertIn("observation_labels = '833 1312 959'", observed)
            self.assertIn('observation_byte_limit = 4294967296', observed)
            self.assertEqual(observed.split('[Problem]')[0], original.split('[Problem]')[0])
            self.assertTrue(result['outside_observation_block_byte_identical'])
            self.assertEqual(result['expected_native_QPs'], 12)
            with self.assertRaises(AssertionError):
                prepare.prepare(source, nodes, elements, root/'bad', source_sha='0'*64)


class IntegrityTests(unittest.TestCase):
    def fixture(self, directory, *, missing_point=False, bad_matrix=False, truncated=False,
                duplicate_owner=False, missing_target=False, corrupt_sequence=False):
        by_rank = {r: [] for r in range(4)}

        def add(rank, label, qp, event, field, values, time=2.68, call=1):
            by_rank[rank].append(dict(record=len(by_rank[rank])+1, time=time,
                source_element=label, qp=qp, evaluation=1, call=call,
                event=event, field=field, values=' '.join(map(str, values))))

        for i, label in enumerate(LABELS):
            for qp in range(4):
                if missing_point and (i, qp) == (2, 3):continue
                for t in [2.68, 2.6805]:
                    if missing_target and t > 2.68:continue
                    add(i, label, qp, 'qp', 'x_y_z_JxW_coord', [1,2,0,0.25,1], t)
                    add(i, label, qp, 'qp', 'reference_x_y_z', [qp%2,qp//2,0], t)
            add(i, label, 0, 'return', 'column_major_6x6', [1]*36)
            add(i, label, 0, 'assembly', 'RankFour_ijkl', [1]*(80 if bad_matrix else 81))
        if duplicate_owner:
            add(3, LABELS[0], 0, 'qp', 'x_y_z_JxW_coord', [1,2,0,0.25,1])
        add(0, 0, 0, 'root', 'row_major_9x9', [1]*81)
        for t in [2.68, 2.6805]:
            if not (missing_target and t > 2.68):add(0,0,0,'accepted_end','time',[t],t)
        for rank, rows in by_rank.items():
            if rows:
                if corrupt_sequence and rank == 0:rows[-1]['record'] = rows[-2]['record']
                with (directory/f'trace-test-rank{rank}.csv').open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            status=dict(rank=rank,records=len(rows),capture_truncated=truncated and rank==0,
                        io_failed=False,history_incomplete=False,streams={'test':len(rows)} if rows else {})
            (directory/f'trace-status-rank{rank}.json').write_text(json.dumps(status))

    def review(self, **options):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root,**options)
            return auditor.audit(root,'trace',job_id='synthetic-only',labels=LABELS,
                require_reject=False,window=[2.68,2.6805],nominal_dt=.0005)

    def test_no_global_reject_does_not_claim_rollback_coverage(self):
        r=self.review();self.assertTrue(r['complete_for_declared_scope'],r['errors'])
        self.assertFalse(r['rollback_observed'])
        self.assertEqual(r['expected_native_QPs'],12)

    def test_missing_QP_rejected(self):self.assertFalse(self.review(missing_point=True)['complete_for_declared_scope'])
    def test_bad_matrix_rejected(self):self.assertIn('bad_RankFour_shape',self.review(bad_matrix=True)['errors'])
    def test_truncation_rejected(self):self.assertFalse(self.review(truncated=True)['complete_for_declared_scope'])
    def test_duplicate_owner_rejected(self):self.assertIn('duplicate_owner',self.review(duplicate_owner=True)['errors'])
    def test_missing_target_rejected(self):self.assertIn('missing_accepted_target_in_window',self.review(missing_target=True)['errors'])
    def test_record_gap_or_duplicate_rejected(self):self.assertIn('duplicate_or_invalid_record',self.review(corrupt_sequence=True)['errors'])
    def test_default_eight_QP_scope_cannot_accept_twelve_QP(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root)
            self.assertFalse(auditor.audit(root,'trace',development=True)['complete_for_declared_scope'])


if __name__ == '__main__':unittest.main()
