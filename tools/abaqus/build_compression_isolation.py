#!/usr/bin/env python3
"""Prepare a free-lateral HEX8 diagnostic from the sealed 2026-08-31 Job.

No submission or solver execution. Tables are copied byte-for-byte. Optional
zero viscosity is a labelled isolation experiment, never the original baseline.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

BASELINE_INPUT_SHA = 'f8f558b93ad674192ca51f4e7fae0ce85eeb774d9127f00fdeb095e478034e9a'
SOLVER_SHA = 'fecadd3e13cdb7003d8c7fc7a0365bd588254730'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(baseline, output, zero_viscosity=False):
    manifest = json.loads((baseline/'manifest.json').read_text())
    snap = manifest['input_snapshot']
    source = baseline/'input'/snap['input_file']
    assert sha(source) == BASELINE_INPUT_SHA == snap['input_sha256']
    assert manifest['solver_versions']['dam_safety_app'] == SOLVER_SHA
    original = source.read_text()
    common_path = Path(__file__).resolve().parents[2]/'test/tests/abaqus_cdp_stress_update/single_hex8_common.i'
    common = common_path.read_text()
    prefix = common.split('[Materials]', 1)[0]
    prefix = prefix.replace('  elem_type = HEX8', '  xmin = 0\n  xmax = 0.15\n  ymin = 0\n  ymax = 0.15\n  zmin = 0\n  zmax = 0.15\n  elem_type = HEX8')
    prefix = re.sub(r"generate_output = '[^']+'", "generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_xx strain_yy strain_zz strain_xy strain_yz strain_xz vonmises_stress'", prefix)
    material = '[Materials]'+original.split('[Materials]', 1)[1].split('[Postprocessors]', 1)[0]
    material = material.replace('    block = concrete_cube__concrete\n', '')
    if zero_viscosity:
        material, count = re.subn(r'(?m)^    viscosity = .*$', '    viscosity = 0', material)
        assert count == 1
    post = '[Postprocessors]'+common.split('[Postprocessors]', 1)[1].split('[Preconditioning]', 1)[0]
    extra = ''
    for field in ('stress_zz','stress_yz','stress_xz','strain_zz','strain_yz','strain_xz','combined_damage'):
        extra += f'  [average_{field}]\n    type = ElementAverageValue\n    variable = {field}\n  []\n'
    end = post.rfind('[]')
    post = post[:end]+extra+post[end:]
    mode = 'mu0' if zero_viscosity else 'mu0005'
    case = 'tpl-cdpc-uc-free-hex8-'+mode+'-v1'
    driver = original[original.index('[Preconditioning/smp]'):].replace('file_base = uniaxial_compression', 'file_base = compression_isolation')
    boundary = '''[Functions/load_path]
  type = ParsedFunction
  expression = '-0.0025*t'
[]
[BCs]
  [x_symmetry]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [y_symmetry]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [z_base]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [z_compression]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = load_path
  []
[]
'''
    output.mkdir(parents=True, exist_ok=False)
    inputs = output/'input'
    inputs.mkdir()
    filename = case+'.i'
    (inputs/filename).write_text('# Free-lateral single HEX8; diagnostic, not original constrained cube.\n'+prefix+material+post+boundary+driver)
    copied = {}
    for item in snap['extra_files']:
        if item['name'] not in ('compression_hardening.csv','compression_damage.csv','tension_stiffening.csv','tension_damage.csv','uniaxial_compression.inp','release-audit.json'):
            continue
        src = baseline/'input'/item['name']
        assert sha(src) == item['sha256']
        (inputs/item['name']).write_bytes(src.read_bytes())
        copied[item['name']] = item['sha256']
    assert len(copied) == 6
    diagnostic = {'test_id': 'TEST-CDPC-UC-D03' if zero_viscosity else 'TEST-CDPC-UC-D02',
        'classification': 'mock/prototype diagnostic; no engineering acceptance',
        'baseline_job': manifest['job_id'], 'baseline_input_sha256': BASELINE_INPUT_SHA,
        'expected_solver_sha': SOLVER_SHA, 'builder_sha256': sha(Path(__file__)),
        'common_input_sha256': sha(common_path), 'copied_files': copied,
        'changed_from_full_case': ['GeneratedMesh 1 HEX8, 0.15m cube', 'free lateral / symmetry faces, original top lateral clamp removed', 'diagnostic state outputs']+(['viscosity=0 isolation'] if zero_viscosity else []),
        'preserved': ['four tables','E/nu','dilation/eccentricity/strength ratios','wt=0,wc=1','loading -0.0025*t through 1s','global/local tolerances','material maximum strain increment 2.5e-5, maximum substeps 256','adaptive time control and 0.01s output'],
        'limitations': ['single-element output is still projected, not an instrumented raw-IP dump','no Abaqus result CSV is solver input','source INP retained for material provenance, not executed as this diagnostic geometry']}
    (inputs/'diagnostic-manifest.json').write_text(json.dumps(diagnostic, ensure_ascii=False, indent=2)+'\n')
    submission = {'project_id':'gmp-ise','case_name':case,'input_file':filename,'input_sha256':sha(inputs/filename),
        'mesh_files':[], 'extra_files':[{'name':p.name,'sha256':sha(p)} for p in sorted(inputs.iterdir()) if p.name != filename],
        'command':'DamSafetyApp-opt -i '+filename}
    (output/'submission.json').write_text(json.dumps(submission, ensure_ascii=False, indent=2)+'\n')
    return submission


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--zero-viscosity', action='store_true')
    args = parser.parse_args()
    print(json.dumps(build(args.baseline, args.output, args.zero_viscosity), indent=2))
