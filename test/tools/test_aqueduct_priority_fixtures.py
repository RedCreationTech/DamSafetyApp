import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import netCDF4
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
import build_aqueduct_priority_fixtures as fixtures


class PriorityFixtureTest(unittest.TestCase):
    def test_three_orthogonal_internal_bars_have_unique_runtime_hosts(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'fixture';fixtures.embedded(out)
            with netCDF4.Dataset(out/'fixture.e') as nc:
                xyz=np.column_stack([nc[v][:] for v in ('coordx','coordy','coordz')])
                vectors=[]
                for j in (2,3,4):
                    conn=nc[f'connect{j}'][0]-1;vectors.append(xyz[conn[1]]-xyz[conn[0]])
                np.testing.assert_allclose(vectors,np.eye(3)*.5)
            rows=list(csv.DictReader((out/'embedded-map.csv').open()))
            self.assertEqual(len({r['secondary_node'] for r in rows}),6)
            self.assertEqual({r['primary_element'] for r in rows},{'0'})
            text=(out/'embedded_three_axis.i').read_text()
            self.assertEqual(text.count('type = AqueductEmbeddedConstraint'),9)

    def test_shell_solid_rows_reproduce_rigid_motion_and_keep_functors_distinct(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'fixture';fixtures.shell_solid(out)
            data=json.loads((out/'coupling-rows.json').read_text());xyz=np.array(data['coordinates'])
            theta=np.array([.01,-.02,.03]);u=np.cross(theta,xyz)+[.1,.2,.3]
            values=np.column_stack([u,np.tile(theta,(len(xyz),1))])
            for row in data['rows']:
                residual=sum(t['coefficient']*values[t['node'],fixtures.FIELDS.index(t['variable'])] for t in row['terms'])
                self.assertAlmostEqual(residual,0,places=12)
            text=(out/'shell_solid_transfer.i').read_text()
            self.assertEqual(text.count("displacements = ''"),6)
            for v in fixtures.FIELDS:self.assertIn('functor = load_'+v,text)

    def test_plastic_oracle_keeps_cumulative_hardening_after_reverse_unloading(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'fixture';fixtures.shell_plastic(out)
            d=json.loads((out/'fixture-identity.json').read_text());a,b,c,z=d['expected_history']
            self.assertGreater(b['right_u'],0)
            self.assertEqual(z['right_u'],0)
            self.assertAlmostEqual(z['peq_1'],2*a['peq_1'])
            self.assertEqual(a['right_moment'],-450000.)


if __name__=='__main__':unittest.main()
