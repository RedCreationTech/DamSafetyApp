#!/usr/bin/env python3
"""Create an opt-in input from an unmodified Q0 snapshot; never run a solver."""
import argparse
import hashlib
from pathlib import Path

Q0_SHA = '54117c1e85bd8d0f92cfb63f45458dad4908c9be6f2a32cbfad85b0fb8535447'
PROBLEM = '''[Problem]
  type = CDPObservationProblem
  observation_enabled = true
  observation_start = 3.018
  observation_end = 3.202
  earthquake_offset = 1
  observation_prefix = ds04g-2597-observation
  observation_byte_limit = 65536000
  observation_labels = '2597'
  observation_vertices = '39.92139453125 25.22817578125 40.18651953125 24.391552734375 41.13575390625 24.682125 40.742546875 25.558146484375'
[]

'''


def prepare(source):
    if hashlib.sha256(source).hexdigest() != Q0_SHA:
        raise ValueError('input is not the frozen Q0 snapshot')
    text = source.decode('utf-8')
    if '[Problem]' in text or text.count('type = TimeSequenceStepper') != 1 or text.count('[Executioner]') != 1:
        raise ValueError('unexpected Q0 input structure')
    return text.replace('[Executioner]', PROBLEM + '[Executioner]', 1).replace(
        'type = TimeSequenceStepper', 'type = CDPObservationTimeStepper', 1).encode('utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    output = prepare(args.input.read_bytes())
    with args.output.open('xb') as f:
        f.write(output)
    print('output_sha256=' + hashlib.sha256(output).hexdigest())


if __name__ == '__main__':
    main()
