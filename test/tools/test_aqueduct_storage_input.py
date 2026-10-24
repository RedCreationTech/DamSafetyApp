import importlib.util
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/abaqus'))
from build_aqueduct_full_input import configure_mpc_storage

class AqueductStorageInputTest(unittest.TestCase):
    def test_only_explicit_mpc_types_change_and_default_bytes_preserved(self):
        source='[Constraints]\n  [row]\n    type = AqueductLinearMPC\n  []\n  [embedded]\n    type = AqueductEmbeddedMPC\n  []\n  [contact]\n    type = AqueductElasticSlipConstraint\n  []\n[]\n'
        self.assertEqual(configure_mpc_storage(source),source)
        changed=configure_mpc_storage(source,True)
        self.assertEqual(changed.count('reuse_assembly_storage = true'),2)
        self.assertEqual(changed.replace('    reuse_assembly_storage = true\n',''),source)
    def test_physical_coefficients_and_existing_couplings_are_untouched(self):
        source='type = AqueductLinearMPC\ncoefficients = \'1 -0.125 -0.875\'\nmultiplier_compliance = 1.652892561983471e-6\n'
        changed=configure_mpc_storage(source,True)
        self.assertEqual(changed.replace('reuse_assembly_storage = true\n',''),source)

if __name__=='__main__':unittest.main()
