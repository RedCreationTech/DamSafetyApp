#!/usr/bin/env python3
"""Check actual passive small-model outputs; never launch a solver."""
import csv
import json
import math
import re
from pathlib import Path


def numeric_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        value = line.strip().rstrip(';')
        if re.fullmatch(r'[-+0-9.eE\s]+', value or '!'):
            rows.append([float(x) for x in value.split()])
    assert rows and all(math.isfinite(x) for row in rows for x in row), path
    return rows


def vector(path, n):
    values = numeric_rows(path)
    assert len(values) == n and all(len(row) == 1 for row in values), path
    return [row[0] for row in values]


def matrix(path, n):
    result = [[0.0] * n for _ in range(n)]
    for i, j, value in numeric_rows(path):
        assert int(i) == i and int(j) == j and 1 <= i <= n and 1 <= j <= n
        result[int(i)-1][int(j)-1] += value
    return result


def norm(values):
    return math.sqrt(math.fsum(x*x for x in values))


def matvec(a, x):
    return [math.fsum(v*w for v, w in zip(row, x)) for row in a]


def verify(root):
    evaluations = list(csv.reader((root/'passive_capture_evaluations.csv').open()))
    accepted = list(csv.reader((root/'passive_capture_accepted.csv').open()))
    linear = list(csv.reader((root/'passive_capture_linear_solves.csv').open()))
    assert evaluations and accepted
    assert {round(float(row[1]), 8) for row in evaluations} == {0.04, 0.08}
    outputs = []
    for last, time, iteration, reported_norm in accepted:
        last, iteration = int(last), int(iteration)
        stem = 'passive_capture_eval'+str(last)
        rows = list(csv.DictReader((root/(stem+'_map_rank0.csv')).open()))
        rows.sort(key=lambda row: int(row['dof']))
        n = len(rows)
        assert n == 12 and [int(row['dof']) for row in rows] == list(range(n))
        assert {row['variable'] for row in rows} == {'q_x', 'q_y', 'strain_zz'}
        scales = [float(row['scaling']) for row in rows]
        assert all(math.isfinite(x) and x > 0 for x in scales)
        state = vector(root/('passive_capture_accepted_eval'+str(last)+'_it'+str(iteration)+'_u.m'), n)
        # Match real evaluations by their complete state, not just their index.
        candidates = [int(row[0]) for row in evaluations
                      if abs(float(row[1])-float(time)) < 1e-9 and int(row[0]) <= last]
        matches = [eid for eid in candidates
                   if max(abs(a-b) for a, b in zip(state, vector(root/('passive_capture_eval'+str(eid)+'_u.m'), n)))
                   <= max(1e-25, max(map(abs, state))*1e-12)]
        assert matches, ('No actual residual identity', last, iteration)
        r = vector(root/('passive_capture_eval'+str(max(matches))+'_r.m'), n)
        assert abs(norm(r)-float(reported_norm)) <= max(1e-25, abs(float(reported_norm))*1e-9)
        if iteration == 0:
            continue  # No completed linear solve exists yet.
        completed = [int(row[0]) for row in linear if int(row[3]) == 1
                     and abs(float(row[1])-float(time)) < 1e-9 and int(row[0]) <= last]
        assert completed, ('No completed linear capture', last, iteration)
        linear_id = max(completed)
        linear_stem = 'passive_capture_linear_eval'+str(linear_id)
        k = matrix(root/(linear_stem+'_J.m'), n)
        pc = matrix(root/('passive_capture_actual_jac_eval'+str(last)+'_J.m'), n)
        assert max(abs(a-b) for ka, pa in zip(k, pc) for a, b in zip(ka, pa)) <= max(
            1e-25, max(abs(x) for row in k for x in row)*1e-12)
        rhs = vector(root/(linear_stem+'_rhs.m'), n)
        solution = vector(root/(linear_stem+'_solution.m'), n)
        actual_rhs_state = vector(root/('passive_capture_eval'+str(linear_id)+'_r.m'), n)
        sign_error = min(norm([a-sign*b for a, b in zip(rhs, actual_rhs_state)]) for sign in (-1, 1))
        assert sign_error <= max(1e-25, max(norm(rhs), norm(actual_rhs_state))*1e-9)
        error = norm([a-b for a, b in zip(matvec(k, solution), rhs)])
        k_inf = max(math.fsum(map(abs, row)) for row in k)
        denominator = math.sqrt(n)*k_inf*max(map(abs, solution))+norm(rhs)
        backward_error = error/max(denominator, 1e-300)
        assert backward_error <= 1e-9, ('Stored KSP solution changed or mismatched', backward_error)
        raw = [[v/scales[i] for v in row] for i, row in enumerate(k)]
        raw_inf = max(math.fsum(map(abs, row)) for row in raw)
        z = [i for i, row in enumerate(rows) if row['variable'] == 'strain_zz']
        assert len(z) == 4
        right = left = 0.0
        for dof in z[1:]:
            direction = [0.0] * n
            direction[dof], direction[z[0]] = 1.0, -1.0
            right = max(right, max(map(abs, matvec(raw, direction)))/max(raw_inf, 1e-300))
            left = max(left, max(abs(raw[dof][j]-raw[z[0]][j]) for j in range(n))/max(raw_inf, 1e-300))
        assert max(right, left) <= 1e-9, (right, left)
        # Record raw residuals by variable; constrained displacement rows are
        # algebraic BC equations and are not direct physical reaction forces.
        raw_residual = {variable: max(abs(r[i]/scales[i]) for i, row in enumerate(rows)
                                     if row['variable'] == variable)
                        for variable in ('q_x', 'q_y', 'strain_zz')}
        outputs.append(dict(time=float(time), iteration=iteration,
                            linear_evaluation=linear_id,
                            matching_actual_evaluations=matches,
                            linear_backward_error=backward_error,
                            right_null_normalized=right, left_null_normalized=left,
                            unscaled_residual_by_variable=raw_residual))
    assert {round(row['time'], 8) for row in outputs} == {0.04, 0.08}
    control = list(csv.DictReader((root/'control-history.csv').open()))
    capture = list(csv.DictReader((root/'capture-history.csv').open()))
    assert len(control) == len(capture) == 21
    maximum = 0.0
    for a, b in zip(control, capture):
        assert a.keys() == b.keys()
        for key in a:
            if not key.startswith('strain_zz_'):
                maximum = max(maximum, abs(float(a[key])-float(b[key])))
        center_a = math.fsum(float(a['strain_zz_'+str(i)]) for i in range(4))/4
        center_b = math.fsum(float(b['strain_zz_'+str(i)]) for i in range(4))/4
        maximum = max(maximum, abs(center_a-center_b))
    assert maximum <= 1e-9, ('Passive observation changed motion', maximum)
    summary = dict(passed=True, actual_capture_times=[0.04, 0.08], systems=outputs,
                   maximum_observable_control_difference=maximum,
                   limitations=['Small one-point QUAD4 only, no Abaqus reference',
                                'Row unscaling does not recover constrained reactions',
                                'No full-dam gate or deficient-pivot exception waived'])
    (root/'capture-verification.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(verify(Path('.')), indent=2))
