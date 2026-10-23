import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

path=Path(__file__).resolve().parents[2]/'tools/abaqus/build_aqueduct_interface_cases.py'
spec=importlib.util.spec_from_file_location('aqueduct_interfaces',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class InterfaceCasesTest(unittest.TestCase):
    def test_connector_retains_three_springs_and_three_align_components(self):
        text=module.connector_case([605000,605000,21965000])
        self.assertEqual(text.count('type = LinearNodalConstraint'),6)
        self.assertEqual(text.count('formulation = kinematic'),3)
        self.assertEqual(text.count('formulation = penalty'),3)
        self.assertIn('penalty = 21965000',text)
        self.assertEqual(text.count('weights = 1'),6)
        self.assertEqual(text.count('primary = 0'),6)
        self.assertEqual(text.count('secondary_node_ids = 1'),6)

    def test_contact_fixture_does_not_close_or_tie_initial_gap(self):
        text=module.contact_case()
        self.assertIn('zmin = 0.16',text)
        self.assertIn('zmax = 0.1',text)
        self.assertIn('friction_coefficient = .55',text)
        self.assertIn('formulation = mortar',text)
        self.assertNotIn('glued',text)
        self.assertNotIn('tangential_tolerance',text)

    def test_contract_and_overwrite_are_checked(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);contract=root/'contract.json';out=root/'out'
            d=dict(case_id='AQ1023-FULL-ELC',interfaces=dict(connector_behaviors_si={str(i):dict({'1':605000.,'2':605000.,'3':21965000.}) for i in range(4)}))
            contract.write_text(json.dumps(d));summary=module.build(contract,out)
            self.assertEqual(len(summary['artifacts']),3)
            self.assertFalse((out/'full-aqueduct.i').exists())
            self.assertIn('not yet reproduced',summary['contact_limitation'])
            with self.assertRaises(ValueError):module.build(contract,out)
            d['interfaces']['connector_behaviors_si']['0']['1']=1
            contract.write_text(json.dumps(d))
            with self.assertRaises(ValueError):module.build(contract,root/'out2')

if __name__=='__main__':unittest.main()
