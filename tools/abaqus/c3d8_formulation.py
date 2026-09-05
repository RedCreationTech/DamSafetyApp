"""Apply the small-strain Abaqus C3D8 volumetric formulation to a MOOSE input.

Only call after verifying that the source uses C3D8 (not C3D8R/C3D8I).
The physics-level parameter updates both strain and stress-divergence terms.
"""

def apply_c3d8_bbar(text: str) -> str:
    anchor = '[Physics/SolidMechanics/QuasiStatic/concrete]\n'
    if text.count(anchor) != 1:
        raise ValueError('Expected one concrete quasi-static physics block')
    start = text.index(anchor) + len(anchor)
    end = text.index('\n[]', start)
    block = text[start:end]
    if 'strain = SMALL' not in block or 'incremental = true' not in block:
        raise ValueError('Only incremental small-strain C3D8 is supported')
    if 'volumetric_locking_correction' in text:
        raise ValueError('Existing volumetric formulation requires explicit review')
    return text[:start] + ('  # Abaqus C3D8: element-averaged volumetric strain and virtual strain.\n'
                          '  volumetric_locking_correction = true\n') + text[start:]
