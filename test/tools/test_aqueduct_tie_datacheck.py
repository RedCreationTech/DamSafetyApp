import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
from audit_aqueduct_tie_datacheck import parse_dat


class TiePrintTest(unittest.TestCase):
    text = '''LOCAL TO GLOBAL NODE AND ELEMENT MAPS
    P-1  1  101
    Q-1  2  202
                 element    element
    BAD-1 99 999
T I E  C O N S T R A I N T S
         101  TRANSLATION         202
                              203 204
E M B E D D E D   E L E M E N T
         777    8  123  0.5
*tie, name=ASSEMBLY_PJZ1, adjust=YES
 101  1.1E-12 ADJUSTED TO ( 1.0, 2.0, 3.0 ).
C O N T A C T   P A I R (S)
 -273- -274- -272- 50.662 FINITE SLIDING
'''

    def test_maps_continuations_and_embedded_separation(self):
        nodes, ties, adjustments, contacts = parse_dat(self.text)
        self.assertEqual(nodes, {101: ('P-1', 1), 202: ('Q-1', 2)})
        self.assertEqual(ties[101]['main_global_labels'], [202, 203, 204])
        self.assertNotIn(777, ties)
        self.assertEqual(adjustments['ASSEMBLY_PJZ1'][0]['distance_mm_printed'], 1.1e-12)
        self.assertEqual(contacts[0]['characteristic_length_mm_printed'], 50.662)

    def test_duplicate_labels_fail_closed(self):
        with self.assertRaisesRegex(ValueError, 'duplicate DAT'):
            parse_dat(self.text.replace('Q-1  2  202', 'Q-1  2  101'))
        with self.assertRaisesRegex(ValueError, 'duplicate secondary'):
            parse_dat(self.text.replace('                              203 204', '         101  TRANSLATION         203'))


if __name__ == '__main__':
    unittest.main()
