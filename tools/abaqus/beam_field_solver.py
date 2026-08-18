#!/usr/bin/env python3
# beam_field_solver.py — 梁模型地震时程与三场结果求解器
#
# 由 2026-08-10 PR-RG-400gal-X 算例创建。动机: hongchuang-opt (MOOSE) 在
# i5-1038NG7 上每残差评估 ~2s, 6500 步不可行; 且发现 PresetAcceleration
# scale_factor 不被使用 (MOOSE bug)。本求解器对线性体系精确组装 K/M/C,
# MPC 用约束变换精确消元, Newmark 常加速度, 一次性稀疏 LU 分解后逐步回代
# — 全 65s (6500 步) 分钟级完成。
#
# DamSafetyApp 版本源自 demo-process@0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8
# 的 tools/beam_direct_solver.py，并增加节点绝对加速度恢复、梁截面应力包络
# 恢复，以及位移/加速度/应力三份独立 Exodus 输出。
# 2026-08-17: Newmark γ/β 改为 CLI 可配 (--newmark-gamma/--newmark-beta,
# 默认 0.5/0.25 不变), 并新增 --allfields-out 单场合并 Exodus:
# 节点位移/转角、绝对加速度、基底约束反力, 单元应力包络、主应力与应变。
#
# 模型语义与 tools/gen_beam_case.py 一致 (同一 report.json):
#   - 3D Timoshenko 梁 (φ 剪切因子), 截面特性来自 Abaqus I/CIRC 公式
#   - y_orientation = report n1 (已按 Abaqus 投影语义处理)
#   - *Damping → C = α·M + β·K
#   - *Mass / *Rotary Inertia → 主节点集中质量/惯量
#   - *Nonstructural Mass → 附加密度 Δρ = m'/A
#   - *MPC BEAM → 从节点 6 自由度刚体约束消元 (u_s = u_m + θ_m × r)
#   - *RELEASE s1/s2 allm → 端部铰接静力凝聚 (精确)
#   - 基底: dofs 2-6 固定, x 向 prescribed 加速度 = scale × amplitude(t)
#
# 应力语义:
#   - axial_stress: 梁端轴力绝对值包络 / A
#   - bending_stress: 两个主弯矩在截面最外缘的绝对应力包络
#   - torsional_shear: T*c/J 的截面外缘近似值
#   - vonmises_stress: sqrt((axial+bending)^2 + 3*torsion^2)
#   - s_max/mid/min_principal: 以最外缘 σ=axial+bending 与 τ=torsional_shear
#     按平面应力主值公式 σ/2±sqrt(σ²/4+τ²)、σ2=0 (与 vonmises 口径自洽:
#     sqrt(σ1²-σ1σ3+σ3²) == sqrt(σ²+3τ²))
#   - axial_strain/bending_strain: σ/E 的线弹性外缘应变
#   这是线性梁截面的保守可视化包络，不是实体积分点连续应力。

import argparse
import json
import time

import numpy as np
import netCDF4
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import splu, eigsh

NEWMARK_BETA = 0.25
NEWMARK_GAMMA = 0.5


def newmark_coefficients(gamma, beta, dt):
    """返回 Newmark 更新式系数 (a0..a5)。γ=0.5/β=0.25 为常加速度法;
    γ>0.5 引入算法阻尼 (等价 HHT-α, α=0.5-γ), 精度降为一阶。"""
    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt * (gamma / (2 * beta) - 1.0)
    return a0, a1, a2, a3, a4, a5


# ---------------------------------------------------------------- 截面特性
def circ_props(r):
    A = np.pi * r * r
    I = np.pi * r ** 4 / 4.0
    return A, I, I, 2 * I          # A, Iy, Iz, J


def isect_props(dims):
    l, h, b1, b2, t1, t2, t3 = dims
    hw = h - t1 - t2
    Af1, Af2, Aw = b1 * t1, b2 * t2, hw * t3
    A = Af1 + Af2 + Aw
    ybar = (Af1 * t1 / 2 + Aw * (t1 + hw / 2) + Af2 * (h - t2 / 2)) / A
    In1 = (b1 * t1 ** 3 / 12 + Af1 * (ybar - t1 / 2) ** 2 +
           t3 * hw ** 3 / 12 + Aw * (t1 + hw / 2 - ybar) ** 2 +
           b2 * t2 ** 3 / 12 + Af2 * (h - t2 / 2 - ybar) ** 2)
    In2 = t1 * b1 ** 3 / 12 + t2 * b2 ** 3 / 12 + hw * t3 ** 3 / 12
    J = (b1 * t1 ** 3 + b2 * t2 ** 3 + hw * t3 ** 3) / 3.0
    return A, In1, In2, J


def section_fiber_limits(section, dims):
    """返回与 Iy/Iz 对应的最外缘距离 (cy, cz)。"""
    if section == 'CIRC':
        return dims[0], dims[0]
    _, h, b1, b2, t1, t2, t3 = dims
    hw = h - t1 - t2
    af1, af2, aw = b1 * t1, b2 * t2, hw * t3
    area = af1 + af2 + aw
    ybar = (af1 * t1 / 2 + aw * (t1 + hw / 2) +
            af2 * (h - t2 / 2)) / area
    return max(b1, b2) / 2.0, max(ybar, h - ybar)


# ---------------------------------------------------------------- 梁单元
def beam_ke_local(E, G, A, Iy, Iz, J, L, Ay, Az):
    """3D Timoshenko 梁局部刚度 (12×12), 局部 x=轴向, y=n1"""
    k = np.zeros((12, 12))

    def put(idofs, mat):
        for a in range(len(idofs)):
            for b in range(len(idofs)):
                k[idofs[a], idofs[b]] += mat[a][b]

    put([0, 6], (E * A / L) * np.array([[1, -1], [-1, 1]]))
    put([3, 9], (G * J / L) * np.array([[1, -1], [-1, 1]]))
    phi_z = 12.0 * E * Iz / (G * Ay * L * L) if Ay > 0 else 0.0
    c = E * Iz / (L ** 3 * (1 + phi_z))
    put([1, 5, 7, 11], c * np.array([
        [12, 6 * L, -12, 6 * L],
        [6 * L, (4 + phi_z) * L * L, -6 * L, (2 - phi_z) * L * L],
        [-12, -6 * L, 12, -6 * L],
        [6 * L, (2 - phi_z) * L * L, -6 * L, (4 + phi_z) * L * L]]))
    phi_y = 12.0 * E * Iy / (G * Az * L * L) if Az > 0 else 0.0
    c = E * Iy / (L ** 3 * (1 + phi_y))
    put([2, 4, 8, 10], c * np.array([
        [12, -6 * L, -12, -6 * L],
        [-6 * L, (4 + phi_y) * L * L, 6 * L, (2 - phi_y) * L * L],
        [-12, 6 * L, 12, 6 * L],
        [-6 * L, (2 - phi_y) * L * L, 6 * L, (4 + phi_y) * L * L]]))
    return k


def beam_me_local(rho, A, Iy, Iz, L):
    """一致质量 (局部坐标, MOOSE InertialForceBeam 同约定):
    平动 ρAL (1/3,1/6); 转动 ρL·diag(Iy+Iz, Iz, Iy) (1/3,1/6)"""
    m = np.zeros((12, 12))
    mt = rho * A * L
    for d in range(3):
        m[d, d] += mt / 3
        m[d + 6, d + 6] += mt / 3
        m[d, d + 6] += mt / 6
        m[d + 6, d] += mt / 6
    for a, Ia in enumerate([Iy + Iz, Iz, Iy]):
        w = rho * L * Ia
        m[a + 3, a + 3] += w / 3
        m[a + 9, a + 9] += w / 3
        m[a + 3, a + 9] += w / 6
        m[a + 9, a + 3] += w / 6
    return m


def condense_release(k, m, end):
    """端部弯矩释放 (allm): 释放端转动 DOF 静力凝聚。
    返回 (kc, mc, keep) — 凝聚后矩阵 + 保留的局部 DOF"""
    rel = [9, 10, 11] if end == 's2' else [3, 4, 5]
    keep = [i for i in range(12) if i not in rel]
    k_aa = k[np.ix_(keep, keep)]
    k_ab = k[np.ix_(keep, rel)]
    k_bb = k[np.ix_(rel, rel)]
    kc = k_aa - k_ab @ np.linalg.solve(k_bb, k_ab.T)
    mc = m[np.ix_(keep, keep)]
    # 被释放端转动惯量并入该端平动对角 (保持总转动惯量量级, 防零质量 DOF)
    for r in rel:
        tnode = [i for i in keep if (i // 6) == (r // 6) and i % 6 < 3]
        for t in tnode:
            mc[t, t] += m[r, r] / max(1, len(tnode))
    return kc, mc, keep


# ---------------------------------------------------------------- Exodus 写出
def _put_names(variable, names):
    for index, name in enumerate(names):
        encoded = name.encode()
        variable[index, :len(encoded)] = list(
            np.frombuffer(encoded, dtype='S1'))


def _prepare_result(mesh_path, out_path, times):
    import shutil
    shutil.copy(mesh_path, out_path)
    nc = netCDF4.Dataset(out_path, 'r+')
    vt = nc.variables['time_whole']
    vt[0] = 0.0
    vt[1:len(times) + 1] = np.array(times)
    return nc


def _write_nodal_fields(nc, fields, times):
    """在打开的 nc 上写节点场，fields 为 name -> (frame,node)。"""
    names = list(fields)
    nc.createDimension('num_nod_var', len(names))
    variable = nc.createVariable(
        'name_nod_var', 'S1', ('num_nod_var', 'len_name'))
    _put_names(variable, names)
    for index, name in enumerate(names, start=1):
        values = nc.createVariable(
            f'vals_nod_var{index}', 'f8', ('time_step', 'num_nodes'))
        values[0, :] = 0.0
        values[1:len(times) + 1, :] = fields[name]


def _write_element_fields(nc, fields, times, block_names):
    """在打开的 nc 上写单元场，fields 为 name -> block -> (frame,elem)。"""
    names = list(fields)
    nc.createDimension('num_elem_var', len(names))
    variable = nc.createVariable(
        'name_elem_var', 'S1', ('num_elem_var', 'len_name'))
    _put_names(variable, names)
    truth = nc.createVariable(
        'elem_var_tab', 'i4', ('num_el_blk', 'num_elem_var'))
    truth[:, :] = 0
    for field_index, name in enumerate(names, start=1):
        for block_index, block_name in enumerate(block_names, start=1):
            block_values = fields[name].get(block_name)
            if block_values is None:
                continue
            truth[block_index - 1, field_index - 1] = 1
            values = nc.createVariable(
                f'vals_elem_var{field_index}eb{block_index}', 'f8',
                ('time_step', f'num_el_in_blk{block_index}'))
            values[0, :] = 0.0
            values[1:len(times) + 1, :] = block_values


def write_nodal_result(mesh_path, out_path, fields, times):
    """克隆网格并写节点场，fields 为 name -> (frame,node)。"""
    nc = _prepare_result(mesh_path, out_path, times)
    _write_nodal_fields(nc, fields, times)
    nc.close()


def write_element_result(mesh_path, out_path, fields, times, block_names):
    """克隆网格并写单元场，fields 为 name -> block -> (frame,elem)。"""
    nc = _prepare_result(mesh_path, out_path, times)
    _write_element_fields(nc, fields, times, block_names)
    nc.close()


def write_all_result(mesh_path, out_path, nodal_fields, element_fields,
                     times, block_names):
    """克隆网格并把节点场与单元场写入同一份合并 Exodus。"""
    nc = _prepare_result(mesh_path, out_path, times)
    _write_nodal_fields(nc, nodal_fields, times)
    _write_element_fields(nc, element_fields, times, block_names)
    nc.close()


STRESS_ENVELOPE_NAMES = (
    'axial_stress', 'bending_stress', 'torsional_shear', 'vonmises_stress')


def recover_beam_stress(displacements, element_records):
    """由节点位移恢复每个物理梁 block 的截面应力包络 (MPa)、
    平面应力主值与线弹性外缘应变。"""
    fields = {name: {} for name in STRESS_ENVELOPE_NAMES + (
        's_max_principal', 's_mid_principal', 's_min_principal',
        'axial_strain', 'bending_strain')}
    for block_name, records in element_records.items():
        shape = (displacements.shape[0], len(records))
        values = {name: np.zeros(shape) for name in fields}
        for element_index, record in enumerate(records):
            n0, n1 = record['nodes']
            global_displacement = np.concatenate(
                (displacements[:, n0, :], displacements[:, n1, :]), axis=1)
            local_displacement = global_displacement @ record['transform'].T
            keep = record['keep']
            local_force_kept = local_displacement[:, keep] @ record['ke'].T
            local_force = np.zeros((displacements.shape[0], 12))
            local_force[:, keep] = local_force_kept
            axial_force = np.maximum(
                np.abs(local_force[:, 0]), np.abs(local_force[:, 6]))
            torsion = np.maximum(
                np.abs(local_force[:, 3]), np.abs(local_force[:, 9]))
            moment_y = np.maximum(
                np.abs(local_force[:, 4]), np.abs(local_force[:, 10]))
            moment_z = np.maximum(
                np.abs(local_force[:, 5]), np.abs(local_force[:, 11]))
            axial = axial_force / record['area']
            bending = (moment_y * record['cz'] / record['iy'] +
                       moment_z * record['cy'] / record['iz'])
            torsional = torsion * max(record['cy'], record['cz']) / record['j']
            sigma = axial + bending
            radius = np.sqrt(sigma ** 2 / 4.0 + torsional ** 2)
            values['axial_stress'][:, element_index] = axial
            values['bending_stress'][:, element_index] = bending
            values['torsional_shear'][:, element_index] = torsional
            values['vonmises_stress'][:, element_index] = np.sqrt(
                sigma ** 2 + 3.0 * torsional ** 2)
            values['s_max_principal'][:, element_index] = sigma / 2.0 + radius
            values['s_mid_principal'][:, element_index] = 0.0
            values['s_min_principal'][:, element_index] = sigma / 2.0 - radius
            values['axial_strain'][:, element_index] = axial / record['e']
            values['bending_strain'][:, element_index] = bending / record['e']
        for name in fields:
            fields[name][block_name] = values[name]
    return fields


def recover_base_reactions(stiffness, mass, damping_alpha, damping_beta,
                           displacements, velocities, accelerations,
                           constrained_dofs):
    """在约束平动自由度上恢复反力 R = M·a + C·v + K·u (外载为运动激励,
    约束处残差即反力), 非约束节点置 0。

    displacements/velocities/accelerations 为 (frame, node, 6) 全自由度重构。
    返回 (frame, node, 3) 的 rf_x/y/z。"""
    nfr, nnode = displacements.shape[0], displacements.shape[1]
    damping = damping_alpha * mass + damping_beta * stiffness
    u_flat = displacements.reshape(nfr, -1).T
    v_flat = velocities.reshape(nfr, -1).T
    a_flat = accelerations.reshape(nfr, -1).T
    residual = mass @ a_flat + damping @ v_flat + stiffness @ u_flat
    reactions = np.zeros((nfr, nnode, 3))
    for dof in constrained_dofs:
        if dof % 6 >= 3:
            continue
        reactions[:, dof // 6, dof % 6] = residual[dof, :]
    return reactions


# ---------------------------------------------------------------- 主流程
def build_parser():
    ap = argparse.ArgumentParser()
    ap.add_argument('--report', required=True)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--displacement-out', required=True)
    ap.add_argument('--acceleration-out', required=True)
    ap.add_argument('--stress-out', required=True)
    ap.add_argument('--allfields-out',
                    help='合并全部节点/单元场的单场 Exodus 输出')
    ap.add_argument('--csv', help='输出帧的顶点与全场峰值摘要 CSV')
    ap.add_argument('--newmark-gamma', type=float, default=NEWMARK_GAMMA,
                    help='Newmark γ, 默认 0.5 (常加速度, 无数值阻尼)')
    ap.add_argument('--newmark-beta', type=float, default=NEWMARK_BETA,
                    help='Newmark β, 默认 0.25 (常加速度, 无数值阻尼)')
    ap.add_argument('--dt', type=float, default=None)
    ap.add_argument('--end-time', type=float, default=None)
    ap.add_argument('--output-interval', type=int, default=5)
    ap.add_argument('--save-matrices', help='导出 Kf/Mf/Cf npz (验证用)')
    ap.add_argument('--no-releases', action='store_true',
                    help='忽略 *RELEASE 端部释放 (与 MOOSE 模型对照用)')
    return ap


def main():
    args = build_parser().parse_args()

    t_start = time.time()
    r = json.load(open(args.report))

    # ---- 读网格 ----
    nc = netCDF4.Dataset(args.mesh)
    coords = np.column_stack([nc.variables['coordx'][:],
                              nc.variables['coordy'][:],
                              nc.variables['coordz'][:]])
    nnode = coords.shape[0]
    nblk = nc.dimensions['num_el_blk'].size
    eb_raw = nc.variables['eb_names']
    eb_names = [eb_raw[i].tobytes().decode('ascii', 'replace').strip('\x00 ')
                for i in range(eb_raw.shape[0])]
    blocks = {}
    for bi, nm in enumerate(eb_names, start=1):
        blocks[nm] = nc.variables[f'connect{bi}'][:] - 1
    nsets = {}
    if 'ns_names' in nc.variables:
        ns_raw = nc.variables['ns_names']
        ns_names = [ns_raw[i].tobytes().decode('ascii', 'replace').strip('\x00 ')
                    for i in range(ns_raw.shape[0])]
        for ni, nm in enumerate(ns_names, start=1):
            nsets[nm] = nc.variables[f'node_ns{ni}'][:]
    nc.close()
    print(f"[load] 节点={nnode} 块={len(blocks)} "
          f"单元={sum(len(v) for v in blocks.values())}")

    ndof = nnode * 6
    step = r['steps'][0]
    dyn = step.get('dynamic') or [0.01, 1.0]
    dt = args.dt or dyn[0]
    T = args.end_time or (dyn[1] if len(dyn) > 1 else 1.0)
    nsteps = int(round(T / dt))

    mats = r['materials']
    eta = zeta = 0.0
    for mt in mats.values():
        if mt.get('damping'):
            eta = mt['damping'].get('alpha', 0.0)
            zeta = mt['damping'].get('beta', 0.0)
            break
    print(f"[mat] Rayleigh: α={eta} β={zeta}")

    extra_rho = {}
    for ns in r.get('nonstructural_mass', []):
        for blk in ns['per_block']:
            bs = r['beam_sections'].get(blk)
            if not bs:
                continue
            A = circ_props(bs['dims'][0])[0] if bs['section'] == 'CIRC' \
                else isect_props(bs['dims'])[0]
            extra_rho[blk] = extra_rho.get(blk, 0.0) + ns['value'] / A

    # ---- 单元组装 (全局 K, M; MPC 连杆跳过) ----
    K = lil_matrix((ndof, ndof))
    M = lil_matrix((ndof, ndof))
    release_map = {rel['geid']: rel['end']
                   for rel in r.get('releases_global', [])}
    if args.no_releases:
        release_map = {}
    geid = 0
    n_hinge = 0
    element_records = {}
    for bname, conn in sorted(blocks.items()):
        if bname.startswith('mpc_beam_links'):
            geid += len(conn)
            continue
        bs = r['beam_sections'][bname]
        mat = mats[bs['material']]
        E, nu = mat['elastic'][0][0], mat['elastic'][0][1]
        G = E / (2 * (1 + nu))
        rho = mat['density'][0][0] + extra_rho.get(bname, 0.0)
        if bs['section'] == 'CIRC':
            A, Iy, Iz, J = circ_props(bs['dims'][0])
            Ay = Az = 0.9 * A
        else:
            A, Iy, Iz, J = isect_props(bs['dims'])
            Ay = Az = A / 1.2
        cy, cz = section_fiber_limits(bs['section'], bs['dims'])
        element_records[bname] = []
        n1 = np.array(bs['n1'], dtype=float)
        for e in conn:
            geid += 1
            p0, p1 = coords[e[0]], coords[e[1]]
            t = p1 - p0
            L = np.linalg.norm(t)
            t = t / L
            y = n1 - np.dot(n1, t) * t
            ny = np.linalg.norm(y)
            if ny < 1e-8:
                y = np.array([-t[1], t[0], 0.0])
                ny = np.linalg.norm(y)
                if ny < 1e-8:
                    y = np.array([0.0, -t[2], t[1]])
                    ny = np.linalg.norm(y)
            y = y / ny
            z = np.cross(t, y)
            R = np.vstack([t, y, z])
            T12 = np.zeros((12, 12))
            for i in range(4):
                T12[3 * i:3 * i + 3, 3 * i:3 * i + 3] = R
            ke = beam_ke_local(E, G, A, Iy, Iz, J, L, Ay, Az)
            me = beam_me_local(rho, A, Iy, Iz, L)
            end = release_map.get(geid)
            if end:
                ke, me, keep = condense_release(ke, me, end)
                n_hinge += 1
            else:
                keep = list(range(12))
            element_records[bname].append({
                'nodes': tuple(e), 'transform': T12, 'ke': ke.copy(),
                'keep': list(keep), 'area': A, 'iy': Iy, 'iz': Iz,
                'j': J, 'cy': cy, 'cz': cz, 'e': E})
            Tk = T12[np.ix_(keep, keep)]
            kg = Tk.T @ ke @ Tk
            mg = Tk.T @ me @ Tk
            gd = [(e[0] if i < 6 else e[1]) * 6 + i % 6 for i in keep]
            for a in range(len(gd)):
                K[gd[a], gd[a]] += kg[a, a]
                M[gd[a], gd[a]] += mg[a, a]
                for b in range(a + 1, len(gd)):
                    K[gd[a], gd[b]] += kg[a, b]
                    K[gd[b], gd[a]] += kg[a, b]
                    M[gd[a], gd[b]] += mg[a, b]
                    M[gd[b], gd[a]] += mg[a, b]
    print(f"[assemble] 单元组装完成 (铰接凝聚 {n_hinge}), "
          f"{time.time() - t_start:.1f}s")

    # ---- 点质量/转动惯量 ----
    for pm in r.get('point_mass', []):
        n0 = (pm['gid'] - 1) * 6
        if pm['kind'] == 'mass':
            mass = pm['mass']
            directional_mass = ([mass, mass, mass]
                                if np.isscalar(mass) else mass[:3])
            for d in range(3):
                M[n0 + d, n0 + d] += directional_mass[d]
            print(f"[mass] 节点 {pm['gid']}: m={directional_mass} t")
        else:
            for d in range(3):
                M[n0 + 3 + d, n0 + 3 + d] += pm['inertia'][d]
            print(f"[mass] 节点 {pm['gid']}: I={pm['inertia'][:3]}")

    K = K.tocsr()
    M = M.tocsr()

    # ---- MPC BEAM 约束变换矩阵 G (ndof × nfree) ----
    master_of = {}
    for link in r.get('mpc_links', []):
        master_of[link['slave']] = link['master']

    def exo_set(name):
        for cand in nsets:
            if cand == name or cand.startswith(f"{name}__"):
                return cand
        return name

    fixed_dofs = []
    accel_dofs = []
    amp_name = scale = None
    for b in step['boundaries']:
        nodes = nsets[exo_set(b['set'])]
        if b.get('amplitude'):
            amp_name, scale = b['amplitude'], b['value']
            for nid in nodes:
                accel_dofs.append((nid - 1) * 6 + b['dof1'] - 1)
        else:
            for nid in nodes:
                for dof in range(b['dof1'], b['dof2'] + 1):
                    fixed_dofs.append((nid - 1) * 6 + dof - 1)
    accel_dofs = sorted(set(accel_dofs))
    fixed_dofs = sorted(set(fixed_dofs))
    constrained = set(fixed_dofs) | set(accel_dofs)
    elim_dofs = set()
    elim_map = {}                        # elim dof -> [(free dof, w)]
    for s, m0 in master_of.items():
        rvec = coords[s - 1] - coords[m0 - 1]
        sd, md = (s - 1) * 6, (m0 - 1) * 6
        # u_s = u_m + θ_m × rvec;  θ×r = [θy rz−θz ry, θz rx−θx rz, θx ry−θy rx]
        W = np.array([[0, rvec[2], -rvec[1]],
                      [-rvec[2], 0, rvec[0]],
                      [rvec[1], -rvec[0], 0]])    # W @ θ = θ × r
        for d in range(3):
            elim_map[sd + d] = [(md + d, 1.0)] + [(md + 3 + k, W[d, k])
                                                  for k in range(3)]
            elim_map[sd + 3 + d] = [(md + 3 + d, 1.0)]
    elim_dofs = set(elim_map)

    free = [d for d in range(ndof)
            if d not in constrained and d not in elim_dofs]
    Kd = K.diagonal()
    Md = M.diagonal()
    dangling = [d for d in free if abs(Kd[d]) < 1e-12 and abs(Md[d]) < 1e-12]
    if dangling:
        print(f"[bc] {len(dangling)} 个悬空 dof → 固定")
        constrained.update(dangling)
        free = [d for d in free if d not in constrained]
    fidx = {d: i for i, d in enumerate(free)}
    nfree = len(free)
    print(f"[dof] 总={ndof} 自由={nfree} MPC消元={len(elim_dofs)} "
          f"约束={len(constrained)} (其中加速度 {len(accel_dofs)})")

    # G 的行: free→单位; elim→master 组合; constrained→空
    rows, cols, vals = [], [], []
    for d in free:
        rows.append(d)
        cols.append(fidx[d])
        vals.append(1.0)
    for d, combo in elim_map.items():
        for md, w in combo:
            if md in fidx:
                rows.append(d)
                cols.append(fidx[md])
                vals.append(w)
    G = csr_matrix((vals, (rows, cols)), shape=(ndof, nfree))

    Kf = (G.T @ K @ G).tocsr()
    Mf = (G.T @ M @ G).tocsr()
    acol = np.array(accel_dofs)
    Kfc = (G.T @ K[:, acol]).tocsr()
    Mfc = (G.T @ M[:, acol]).tocsr()
    print(f"[reduce] {time.time() - t_start:.1f}s")

    # ---- 阻尼 & Newmark ----
    gamma = args.newmark_gamma
    beta = args.newmark_beta
    print(f"[int] Newmark γ={gamma} β={beta}"
          + (" (常加速度)" if (gamma, beta) == (0.5, 0.25) else ""))
    Cf = eta * Mf + zeta * Kf
    Cfc = eta * Mfc + zeta * Kfc
    a0, a1, a2, a3, a4, a5 = newmark_coefficients(gamma, beta, dt)
    A = a0 * Mf + a1 * Cf + Kf
    print("[factor] 稀疏 LU ...")
    lu = splu(A.tocsc())
    print(f"[factor] 完成 {time.time() - t_start:.1f}s")

    # ---- 特征频率校核 ----
    try:
        w2 = eigsh(Kf, k=6, M=Mf, sigma=0.01, which='LM', maxiter=5000,
                   return_eigenvectors=False)
        freqs = np.sqrt(np.maximum(w2, 0)) / (2 * np.pi)
        print(f"[eig] 前 6 阶频率 (Hz): {np.sort(freqs).round(3)}")
    except Exception as ex:
        print(f"[eig] 跳过: {ex}")
    if args.save_matrices:
        from scipy.sparse import save_npz
        save_npz(args.save_matrices.replace('.npz', '_K.npz'), Kf)
        save_npz(args.save_matrices.replace('.npz', '_M.npz'), Mf)
        print(f"[save] 矩阵已导出 {args.save_matrices}_{{K,M}}.npz")

    # ---- 地面加速度插值 ----
    amp = np.array(r['amplitudes'][amp_name])
    at, av = amp[:, 0], amp[:, 1]

    def ground_accel(t):
        return scale * float(np.interp(t, at, av))

    # ---- 时程积分 ----
    u = np.zeros(nfree)
    v = np.zeros(nfree)
    acc = np.zeros(nfree)
    ug = vg = 0.0
    out_every = args.output_interval
    frames = []
    veloc_frames = []
    accel_frames = []
    times = []
    ug_frames = []
    vg_frames = []
    ag_frames = []
    top_dof = None
    for pm in r.get('point_mass', []):
        top_dof = (pm['gid'] - 1) * 6
        break
    top_hist = []
    print(f"[solve] {nsteps} 步 dt={dt} T={T} ...")
    ts0 = time.time()
    for n in range(nsteps):
        t_new = (n + 1) * dt
        ag0 = ground_accel(n * dt)
        ag1 = ground_accel(t_new)
        ug1 = ug + dt * vg + dt * dt * ((0.5 - beta) * ag0
                                        + beta * ag1)
        vg1 = vg + dt * ((1 - gamma) * ag0 + gamma * ag1)
        b = (Mf @ (a0 * u + a2 * v + a3 * acc) +
             Cf @ (a1 * u + a4 * v + a5 * acc) -
             Kfc @ (np.full(len(accel_dofs), ug1)) -
             Cfc @ (np.full(len(accel_dofs), vg1)) -
             Mfc @ (np.full(len(accel_dofs), ag1)))
        u1 = lu.solve(b)
        acc1 = a0 * (u1 - u) - a2 * v - a3 * acc
        v1 = v + dt * ((1 - gamma) * acc + gamma * acc1)
        u, v, acc = u1, v1, acc1
        ug, vg = ug1, vg1
        top_hist.append(u[fidx[top_dof]] if top_dof in fidx else np.nan)
        if (n + 1) % out_every == 0 or n == nsteps - 1:
            frames.append(u.copy())
            veloc_frames.append(v.copy())
            accel_frames.append(acc.copy())
            ug_frames.append(ug)
            vg_frames.append(vg)
            ag_frames.append(ag1)
            times.append(t_new)
        if (n + 1) % 500 == 0:
            print(f"  {n + 1}/{nsteps} t={t_new:.1f}s "
                  f"({time.time() - ts0:.0f}s) max|u|={np.abs(u).max():.3g}")

    peak = np.nanmax(np.abs(top_hist))
    print(f"[solve] 完成 {time.time() - ts0:.0f}s, "
          f"顶点 |disp_x| 峰值 = {peak:.4g} mm")

    # ---- 重构全自由度并写 Exodus ----
    print("[recover] 位移、绝对加速度与梁截面应力 ...")
    nfr = len(frames)
    Fr = np.array(frames)                # (nfr, nfree)
    Vr = np.array(veloc_frames)          # (nfr, nfree)
    Ar = np.array(accel_frames)          # (nfr, nfree)

    def reconstruct(free_frames, ground_frames):
        """把 (nfr, nfree) 自由自由度帧重构为 (nfr, nnode, 6) 全场。"""
        full = np.zeros((nfr, nnode, 6))
        for d in free:
            full[:, d // 6, d % 6] = free_frames[:, fidx[d]]
        for d in accel_dofs:
            full[:, d // 6, d % 6] = np.array(ground_frames)
        for d, combo in elim_map.items():
            vec = np.zeros(nfr)
            for md, w in combo:
                if md in fidx:
                    vec += w * free_frames[:, fidx[md]]
                elif md in accel_dofs:
                    vec += w * np.array(ground_frames)
            full[:, d // 6, d % 6] = vec
        return full

    U6 = reconstruct(Fr, ug_frames)
    V6 = reconstruct(Vr, vg_frames)
    A6 = reconstruct(Ar, ag_frames)

    all_fields = recover_beam_stress(U6, element_records)
    stress_fields = {name: all_fields[name] for name in STRESS_ENVELOPE_NAMES}
    displacement_fields = {
        name: U6[:, :, index] for index, name in enumerate(
            ['disp_x', 'disp_y', 'disp_z', 'rot_x', 'rot_y', 'rot_z'])}
    acceleration_fields = {
        name: A6[:, :, index] for index, name in enumerate([
            'accel_x', 'accel_y', 'accel_z',
            'rot_accel_x', 'rot_accel_y', 'rot_accel_z'])}

    print("[write] 三份 Exodus ...")
    write_nodal_result(
        args.mesh, args.displacement_out, displacement_fields, times)
    write_nodal_result(
        args.mesh, args.acceleration_out, acceleration_fields, times)
    write_element_result(
        args.mesh, args.stress_out, stress_fields, times, eb_names)

    if args.allfields_out:
        print("[write] 合并全场 Exodus ...")
        rf = recover_base_reactions(
            K, M, eta, zeta, U6, V6, A6, sorted(constrained))
        reaction_fields = {
            name: rf[:, :, index]
            for index, name in enumerate(['rf_x', 'rf_y', 'rf_z'])}
        write_all_result(
            args.mesh, args.allfields_out,
            {**displacement_fields, **acceleration_fields,
             **reaction_fields},
            all_fields, times, eb_names)

    max_accel = np.sqrt(np.sum(A6[:, :, :3] ** 2, axis=2)).max(axis=1)
    max_stress = np.zeros(nfr)
    for values in stress_fields['vonmises_stress'].values():
        max_stress = np.maximum(max_stress, values.max(axis=1))
    top_node = top_dof // 6 if top_dof is not None else 0
    if args.csv:
        with open(args.csv, 'w') as summary:
            summary.write(
                'time,top_disp_x,top_accel_x,max_accel,max_vonmises_stress\n')
            for index, time_value in enumerate(times):
                summary.write(
                    f'{time_value},{U6[index, top_node, 0]},'
                    f'{A6[index, top_node, 0]},{max_accel[index]},'
                    f'{max_stress[index]}\n')
    print(f"[peak] max|accel|={max_accel.max():.6g} mm/s^2, "
          f"max von Mises envelope={max_stress.max():.6g} MPa")
    print(f"完成 ✓ 总耗时 {time.time() - t_start:.0f}s")


if __name__ == '__main__':
    main()
