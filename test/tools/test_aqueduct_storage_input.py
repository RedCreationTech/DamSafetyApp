import importlib.util
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/abaqus'))
from build_aqueduct_full_input import configure_mpc_storage, configure_mpc_cross_blocks, build

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

    def test_cross_blocks_default_is_byte_preserving(self):
        source='type = AqueductLinearMPC\nreuse_assembly_storage = true\ncoefficients = \'1 -1\'\n'
        self.assertEqual(configure_mpc_cross_blocks(source),source)
        changed=configure_mpc_cross_blocks(source,True)
        self.assertEqual(changed.replace('assemble_cross_blocks = true\n',''),source)

    def test_cross_blocks_only_changes_the_two_fixed_mpc_types(self):
        source='  type = AqueductEmbeddedMPC\n  reuse_assembly_storage = true\n  type = AqueductLinearMPC\n  reuse_assembly_storage = true\n  type = AqueductElasticSlipConstraint\n'
        changed=configure_mpc_cross_blocks(source,True)
        self.assertEqual(changed.count('assemble_cross_blocks = true'),2)
        self.assertEqual(changed.replace('  assemble_cross_blocks = true\n',''),source)

    def test_cross_blocks_requires_storage_before_creating_output(self):
        with self.assertRaisesRegex(ValueError,'require --reuse-mpc-storage'):
            build('/not-read','/not-read','/not-created',mpc_cross_blocks=True)

if __name__=='__main__':unittest.main()
