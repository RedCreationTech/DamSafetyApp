import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import netCDF4

TOOLS = Path(__file__).resolve().parents[2] / 'tools' / 'abaqus'
sys.path.insert(0, str(TOOLS))
spec = importlib.util.spec_from_file_location('prepare_aqueduct_case', TOOLS / 'prepare_aqueduct_case.py')
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)

FIXTURE = '''*Part, name=混凝土
*Node
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
*Element, type=C3D8R
1,1,2,3,4,5,6,7,8
*Elset, elset=ALL, generate
1,1,1
*Solid Section, elset=ALL, material=CON
,
*End Part
*Part, name=BAR
*Node
1,0.2,0.2,0.2
2,0.8,0.2,0.2
*Element, type=T3D2
1,1,2
*Elset, elset=ALL, generate
1,1,1
*Solid Section, elset=ALL, material=STEEL
3.14,
*End Part
*Part, name=PLATE
*Node
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
*Element, type=S4R
1,1,2,3,4
*Elset, elset=ALL, generate
1,1,1
*Shell Section, elset=ALL, material=STEEL, offset=SPOS
100,3
*End Part
*Assembly, name=Assembly
*Instance, name=A, part=混凝土
*End Instance
*Instance, name=B, part=混凝土
*End Instance
*Instance, name=R, part=BAR
*End Instance
*Instance, name=P, part=PLATE
*End Instance
*Element, type=CONN3D2
1, P.1, P.2
*Elset, elset=SPRING, generate
1,1,1
*Connector Section, elset=SPRING, behavior=CB
Cartesian, Align
*Elset, elset=HOST, instance=A
1
*Elset, elset=REBAR, instance=R
1
*Elset, elset=SOLIDFACE, instance=A
1
*Elset, elset=SHELLFACE, instance=P
1
*Surface, type=ELEMENT, name=S1
SOLIDFACE,S1
*Surface, type=ELEMENT, name=S2
SHELLFACE,E1
*Tie, name=TIE, adjust=yes
S1,S2
*Embedded Element, host elset=HOST
REBAR
*Shell to Solid Coupling, constraint name=SS
S2,S1
*End Assembly
*Connector Behavior, name=CB
*Connector Elasticity, component=1
605,
*Material, name=CON
*Density
2.64e-9
*Elastic
30000,0.2
*Material, name=STEEL
*Density
7.85e-9
*Elastic
200000,0.3
*Plastic
400,0
696,0.1
*Surface Interaction, name=FRIC
1,
*Friction, slip tolerance=0.005
0.55,
*Surface Behavior, pressure-overclosure=HARD
*Contact Pair, interaction=FRIC, type=SURFACE TO SURFACE
S1,S2
*Amplitude, name=X
0,0,1,1
*Baseline Correction
*Step, name=GRA, nlgeom=NO
*Static, stabilize=0.0002
0.01,1,1e-15,1
*End Step
*Step, name=ELC, nlgeom=NO
*Dynamic
0.001,27.5,1e-8
*End Step
'''


class PrepareAqueductTest(unittest.TestCase):
    def test_confirmed_status_cannot_override_conflicting_units_or_source(self):
        approval = dict(case_id='AQ1023-FULL-ELC', source_sha256=tool.SOURCE_SHA,
                        A01=dict(status='confirmed', authoritative_complete_input=True,
                                 external_overrides=False, add_gravity_or_water=False),
                        A02=dict(status='confirmed',
                                 source_units=dict(length='mm', force='N', time='s', stress='MPa', mass='tonne'),
                                 preserve_materials_sections_and_connector_parameters=True,
                                 cdp_viscosity_seconds=0.05, output_units='SI'))
        tool.validate_confirmation(approval, tool.SOURCE_SHA)
        with self.assertRaisesRegex(ValueError, 'source'):
            tool.validate_confirmation(approval, 'changed')
        approval['A02']['source_units']['length'] = 'm'
        with self.assertRaisesRegex(ValueError, 'units'):
            tool.validate_confirmation(approval, tool.SOURCE_SHA)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.source = self.work / 'model.inp'
        self.source.write_bytes(FIXTURE.encode('gbk'))

    def test_coincident_instances_and_shell_keep_their_own_nodes(self):
        out = self.work / 'out'
        tool.prepare(self.source, out)
        with netCDF4.Dataset(out / 'assembly-geometry.e') as nc:
            self.assertEqual(len(nc.dimensions['num_nodes']), 22)
            self.assertEqual(len(nc.dimensions['num_elem']), 4)
            self.assertEqual(float(nc['coordx'][1]), 0.001)
            self.assertEqual(float(nc['coordx'][9]), 0.001)
            self.assertNotEqual(int(nc['node_num_map'][1]), int(nc['node_num_map'][9]))
            types = {nc[name].elem_type for name in nc.variables if name.startswith('connect')}
            self.assertEqual(types, {'HEX8', 'TRUSS', 'SHELL4'})
        self.assertIn('混凝土', (out / 'source-utf8.inp').read_text())
        self.assertEqual(self.source.read_bytes(), FIXTURE.encode('gbk'))

    def test_interfaces_keep_host_identity_shell_offset_and_connector_endpoints(self):
        out = self.work / 'out'
        tool.prepare(self.source, out)
        d = json.loads((out / 'assembly-contract.json').read_text())
        self.assertEqual(d['interfaces']['embedded'][0]['host'], [['A', 1]])
        self.assertEqual(d['interfaces']['embedded'][0]['embedded'], [['R', 1]])
        self.assertEqual(d['interfaces']['connectors'][0]['endpoints'][0]['instance'], 'P')
        self.assertEqual(len(d['interfaces']['contacts']), 1)
        shell = next(s for s in d['sections'] if s['keyword'] == 'shell section')
        self.assertEqual(shell['offset'], 'SPOS')
        self.assertEqual(shell['thickness_m'], 0.1)
        self.assertEqual(shell['thickness_integration_points'], 3)
        self.assertEqual(len(d['surfaces']['S2'][0]['nodes']), 2)
        self.assertFalse(d['ready_for_solve'])
        self.assertFalse(d['raw_amplitudes_are_effective_motion'])

    def test_s6_surface_order_is_cyclic_for_numerical_integration(self):
        self.source.write_bytes(FIXTURE.replace('SOLIDFACE,S1','SOLIDFACE,S6').encode('gbk'))
        out=self.work/'s6';tool.prepare(self.source,out)
        d=json.loads((out/'assembly-contract.json').read_text())
        self.assertEqual(d['surfaces']['S1'][0]['nodes'],[4,1,5,8])

    def test_units_scale_stress_density_area_but_not_plastic_strain(self):
        out = self.work / 'out'
        tool.prepare(self.source, out)
        d = json.loads((out / 'assembly-contract.json').read_text())
        self.assertAlmostEqual(d['materials_si']['CON']['density'][0][0], 2640)
        self.assertEqual(d['materials_si']['CON']['elastic'][0], [3e10, 0.2])
        self.assertEqual(d['materials_si']['STEEL']['plastic'][1], [696e6, 0.1])
        truss = next(s for s in d['sections'] if 'area_m2' in s)
        self.assertAlmostEqual(truss['area_m2'], 3.14e-6)

    def test_unknown_physics_is_rejected_before_output(self):
        self.source.write_bytes((FIXTURE + '*Visco\n1\n').encode('gbk'))
        with self.assertRaisesRegex(ValueError, 'unknown keyword'):
            tool.prepare(self.source, self.work / 'out')
        self.assertFalse((self.work / 'out').exists())

    def test_invalid_numeric_or_duplicate_geometry_is_not_silently_skipped(self):
        for replacement in ['1,nan,0,0', '1,0,0,0\n1,0,0,0']:
            with self.subTest(replacement=replacement):
                self.source.write_bytes(FIXTURE.replace('1,0,0,0', replacement, 1).encode('gbk'))
                with self.assertRaises(ValueError):
                    tool.prepare(self.source, self.work / 'out')

    def test_missing_surface_member_fails(self):
        self.source.write_bytes(FIXTURE.replace('SOLIDFACE,S1', 'UNKNOWN,S1').encode('gbk'))
        with self.assertRaisesRegex(ValueError, 'missing assembly element set'):
            tool.prepare(self.source, self.work / 'out')

    def test_nonmonotonic_motion_is_rejected(self):
        self.source.write_bytes(FIXTURE.replace('0,0,1,1', '0,0,0,1').encode('gbk'))
        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            tool.prepare(self.source, self.work / 'out')

    def test_second_run_does_not_overwrite_outputs(self):
        out = self.work / 'out'
        tool.prepare(self.source, out)
        before = (out / 'assembly-contract.json').read_bytes()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            tool.prepare(self.source, out)
        self.assertEqual((out / 'assembly-contract.json').read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
