import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path

import netCDF4


PROJECT_DIR = Path(__file__).resolve().parents[2]
TOOL_PATH = PROJECT_DIR / 'tools' / 'abaqus' / 'abaqus2exodus.py'
SPEC = importlib.util.spec_from_file_location('abaqus2exodus', TOOL_PATH)
CONVERTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONVERTER)


class Abaqus2Exodus2DTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        (self.work / 'mass.inc').write_text(textwrap.dedent('''\
            *Element, type=MASS, elset=emass
            101, 2
            *Mass, elset=emass, TYPE=ANISOTROPIC
            10, 20, 0
            *Element, type=MASS, elset=emass
            102, 3
            *Mass, elset=emass, TYPE=ANISOTROPIC
            30, 40, 0
        '''), encoding='utf-8')
        (self.work / 'mixed.inp').write_text(textwrap.dedent('''\
            *Part, name=DAM
            *Node
            1, 0, 0
            2, 1, 0
            3, 1, 1
            4, 0, 1
            5, 2, 0
            *Element, type=CPS4R
            1, 1, 2, 3, 4
            *Element, type=CPS3
            2, 2, 5, 3
            *Elset, elset=ALL
            1, 2
            *Solid Section, elset=ALL, material=CONCRETE
            1,
            *Include, input=mass.inc
            *End Part
            *Assembly, name=ASSEMBLY
            *Instance, name=DAM-1, part=DAM
            *End Instance
            *Nset, nset=MASS_NODES, instance=DAM-1
            2, 3
            *Elset, elset=LOAD_EDGE, instance=DAM-1
            1
            *Surface, type=ELEMENT, name=LOAD
            LOAD_EDGE, S2
            *End Assembly
            *Material, name=CONCRETE
            *Density
            2260.74
            *Elastic
            3.04e10, 0.2
            *Boundary
            MASS_NODES, 1, 2
            *Amplitude, name=EQ, time=TOTAL TIME
            0.0, -0.1, 0.01, 0.2
            *Step, name=GRA, nlgeom=NO, inc=100
            *Static
            0.1, 1.0, 1e-15, 1.0
            *Dload
            , GRAV, 9.8, 0.0, -1.0
            *Dsload
            LOAD, HP, 592116.0, 60.42, 0.0
            *End Step
            *Step, name=EQ, nlgeom=NO, inc=1000
            *Dynamic
            0.01, 50.0, 5e-7
            *Boundary, op=NEW, amplitude=EQ, type=ACCELERATION
            MASS_NODES, 1, 1, 3.44
            *End Step
        '''), encoding='utf-8')

    def tearDown(self):
        self.temp.cleanup()

    def test_mixed_2d_include_sideset_and_anisotropic_mass(self):
        model = CONVERTER.parse_inp(self.work / 'mixed.inp')
        self.assertEqual(len(model.source_files), 2)
        part = model.parts['DAM']
        self.assertEqual(part.elem_types, {1: 'CPS4R', 2: 'CPS3'})
        self.assertEqual(len(part.point_elems), 2)
        self.assertEqual(part.point_mass[101], [10.0, 20.0, 0.0])
        self.assertEqual(part.point_mass[102], [30.0, 40.0, 0.0])
        self.assertEqual(model.amplitude_options['EQ']['time'], 'TOTAL TIME')
        self.assertEqual(model.initial_boundaries[0]['dof2'], 2)
        self.assertEqual(model.steps[0]['options']['inc'], '100')
        self.assertEqual(model.steps[0]['loads'], [
            {'region': '', 'type': 'GRAV', 'value': 9.8,
             'parameters': [0.0, -1.0]},
            {'surface': 'LOAD', 'type': 'HP', 'value': 592116.0,
             'parameters': [60.42, 0.0]},
        ])
        self.assertEqual(model.steps[1]['boundaries'][0]['type'],
                         'ACCELERATION')
        self.assertEqual(model.steps[1]['boundaries'][0]['op'], 'NEW')

        result = CONVERTER.build_global_mesh(model, 1e-9)
        gm, blocks, block_types, block_meta, nodesets, sidesets = result
        self.assertEqual(len(gm.coords), 5)
        self.assertEqual(sum(len(v) for v in blocks.values()), 2)
        self.assertEqual(sorted(block_types.values()), ['QUAD4', 'TRI3'])
        self.assertEqual(len(sidesets['LOAD']), 1)
        self.assertEqual(sidesets['LOAD'][0][1], 2)

        props = CONVERTER._resolve_point_props(model, gm)
        self.assertEqual(len(props), 2)
        self.assertEqual(
            [sum(row['mass'][i] for row in props) for i in range(3)],
            [40.0, 60.0, 0.0])

        mass_files = CONVERTER.write_nodal_mass_csv(
            props, self.work / 'added_mass')
        self.assertEqual(len(Path(mass_files['x']).read_text().splitlines()), 2)
        self.assertTrue(Path(mass_files['y']).read_text().endswith(',40\n'))

        output = self.work / 'mixed.e'
        CONVERTER.write_exodus(
            output, gm, blocks, block_types, block_meta,
            nodesets, sidesets, 'mixed test')
        with netCDF4.Dataset(output) as exodus:
            self.assertEqual(len(exodus.dimensions['num_nodes']), 5)
            self.assertEqual(len(exodus.dimensions['num_elem']), 2)
            self.assertEqual(len(exodus.dimensions['num_el_blk']), 2)
            self.assertEqual(len(exodus.dimensions['num_side_sets']), 1)
            element_types = sorted(
                getattr(exodus.variables[name], 'elem_type')
                for name in exodus.variables if name.startswith('connect'))
            self.assertEqual(element_types, ['QUAD4', 'TRI3'])
            self.assertEqual(list(exodus.variables['side_ss1'][:]), [2])

    def test_recursive_include_cycle_is_rejected(self):
        (self.work / 'cycle-a.inp').write_text(
            '*Include, input=cycle-b.inp\n', encoding='utf-8')
        (self.work / 'cycle-b.inp').write_text(
            '*Include, input=cycle-a.inp\n', encoding='utf-8')

        with self.assertRaisesRegex(ValueError, '循环引用'):
            CONVERTER.parse_inp(self.work / 'cycle-a.inp')

    def test_mass_without_matching_property_is_rejected(self):
        (self.work / 'missing-mass.inp').write_text(textwrap.dedent('''\
            *Part, name=DAM
            *Node
            1, 0, 0
            *Element, type=MASS, elset=emass
            1, 1
            *End Part
            *Assembly, name=ASSEMBLY
            *Instance, name=DAM-1, part=DAM
            *End Instance
            *End Assembly
        '''), encoding='utf-8')

        model = CONVERTER.parse_inp(self.work / 'missing-mass.inp')
        gm, *_ = CONVERTER.build_global_mesh(model, 1e-9)
        with self.assertRaisesRegex(ValueError, '缺少对应 \\*Mass'):
            CONVERTER._resolve_point_props(model, gm)


if __name__ == '__main__':
    unittest.main()
