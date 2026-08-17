#!/usr/bin/env python3
"""abaqus2exodus.py — Abaqus .inp → Exodus II (MOOSE 兼容) 转换工具

红创科技多物理场仿真平台 / 数据接入工具

DamSafetyApp 版本源自 demo-process@0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8，
在保留既有三维/梁转换能力的基础上增加二维混合单元、二维 sideset、递归
*Include 和各向异性节点附加质量支持。

将 Abaqus 输入卡（.inp）解析并转换为 MOOSE 可直接读取的 Exodus II 网格文件：
  - Part/Instance 装配（含平移 + 旋转）
  - C3D8R → HEX8 实体单元块（按 part × material 分块）
  - CPS4/CPS4R/CPS3 → QUAD4/TRI3 二维实体单元块
  - T3D2  → TRUSS 钢筋单元块
  - Nset / Elset / Surface → Exodus nodesets（供 BC / 荷载 / 后处理引用）
  - 跨 instance 重合节点合并（等效 Abaqus *Tie 绑定约束）
  - 材料 / 幅值 / 分析步导出为 JSON 报告，供生成 MOOSE .i 使用

用法:
  python abaqus2exodus.py --inp /path/to/job.inp --out mesh.e \
      --report report.json --merge-tol 1e-9

依赖: netCDF4 (moose conda env)
"""

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Abaqus 关键字解析
# ---------------------------------------------------------------------------

def _parse_kv(line):
    """解析 'name=value, name=value' 形式的参数串"""
    kv = {}
    for tok in line.split(','):
        tok = tok.strip()
        if not tok:
            continue
        if '=' in tok:
            k, v = tok.split('=', 1)
            kv[k.strip().lower()] = v.strip()
        else:
            kv[tok.lower()] = True
    return kv


def _nums(line):
    """解析逗号分隔的数值行"""
    out = []
    for tok in line.rstrip(',').split(','):
        tok = tok.strip()
        if tok:
            out.append(float(tok))
    return out


def _ids(line):
    return [int(float(t)) for t in line.rstrip(',').split(',') if t.strip()]


class Part:
    def __init__(self, name):
        self.name = name
        self.nodes = {}          # id -> (x, y, z)
        self.elems = {}          # id -> [node ids]
        self.elem_types = {}     # id -> Abaqus element type
        self.etype = None        # backward-compatible last structural type
        self.point_elems = {}    # id -> {type, elset, conn}
        self.point_mass = {}     # point element id -> [mx, my, mz]
        self.nsets = {}
        self.elsets = {}
        self.sections = []       # [(elset, material, area)]
        self.beam_sections = []  # [{elset, material, section, dims, n1}]
        self.releases = []       # [(eid, end, code)]  *RELEASE 端部释放


class Instance:
    def __init__(self, name, part):
        self.name = name
        self.part = part
        self.trans = (0.0, 0.0, 0.0)
        self.rot = None          # (ax,ay,az, bx,by,bz, angle_deg)
        self.data_lines = []     # 原始变换行


class InpModel:
    def __init__(self):
        self.parts = {}
        self.instances = []
        self.asm_nodes = {}      # assembly 级节点 (参考点)
        self.asm_nsets = {}
        self.asm_elsets = {}     # name -> [(instance, ids)] (instance 限定)
        self.asm_node_nsets = {} # nset name -> [(instance, ids)]
        self.surfaces = {}       # name -> [(instance, elset, face)]
        self.materials = defaultdict(dict)   # name -> {density, elastic, plastic, cdp...}
        self.amplitudes = {}
        self.amplitude_options = {}  # name -> *Amplitude keyword parameters
        self.initial_boundaries = []  # model-level *Boundary definitions
        self.steps = []
        self.ties = []
        self.couplings = []
        self.embedded = []
        self.mpcs = []           # [{type, slave, master}]  *MPC 数据行
        self.asm_elems = {}      # assembly 级元素: elset -> {eid: [nids]} (MASS/ROTARYI)
        self.asm_elem_types = {} # elset -> 元素类型 (MASS/ROTARYI)
        self.asm_mass = {}       # elset -> 质量值 *Mass
        self.asm_rotary = {}     # elset -> [I11,I22,I33,I12,I13,I23] *Rotary Inertia
        self.nonstruct_mass = [] # [{elset, units, value}] *Nonstructural Mass
        self.source_files = []   # 主输入及递归 *Include 文件（绝对路径）


# Abaqus C3D8 各面的局部节点 (1-based, 仅用于取节点集合, 不关心绕向)
C3D8_FACES = {
    'S1': (1, 2, 3, 4),
    'S2': (5, 6, 7, 8),
    'S3': (1, 2, 6, 5),
    'S4': (2, 3, 7, 6),
    'S5': (3, 4, 8, 7),
    'S6': (1, 4, 5, 8),
}

ELEMENT_TYPE_MAP = {
    'B31': ('BEAM2', 2),
    'C3D8R': ('HEX8', 8),
    'T3D2': ('TRUSS', 2),
    'CPS4': ('QUAD4', 4),
    'CPS4R': ('QUAD4', 4),
    'CPS3': ('TRI3', 3),
}

# Abaqus 与 libMesh/Exodus 的一阶二维边编号在这些单元上同为 1-based。
ELEMENT_FACES = {
    'CPS4': {'S1': (1, 2), 'S2': (2, 3), 'S3': (3, 4), 'S4': (4, 1)},
    'CPS4R': {'S1': (1, 2), 'S2': (2, 3), 'S3': (3, 4), 'S4': (4, 1)},
    'CPS3': {'S1': (1, 2), 'S2': (2, 3), 'S3': (3, 1)},
    'C3D8R': C3D8_FACES,
}

EXODUS_SIDE_MAP = {
    'CPS4': {'S1': 1, 'S2': 2, 'S3': 3, 'S4': 4},
    'CPS4R': {'S1': 1, 'S2': 2, 'S3': 3, 'S4': 4},
    'CPS3': {'S1': 1, 'S2': 2, 'S3': 3},
}

# 材料子关键字 -> 报告字段名
MAT_KEYS = {
    'density': 'density',
    'elastic': 'elastic',
    'plastic': 'plastic',
    'concrete damaged plasticity': 'cdp',
    'concrete compression hardening': 'cdp_compression_hardening',
    'concrete tension stiffening': 'cdp_tension_stiffening',
    'concrete compression damage': 'cdp_compression_damage',
    'concrete tension damage': 'cdp_tension_damage',
}


def _read_with_includes(path, stack=()):
    """读取 Abaqus 输入并在原位置递归展开 *Include。"""
    source = Path(path).expanduser().resolve()
    if source in stack:
        chain = ' -> '.join(str(p) for p in (*stack, source))
        raise ValueError(f'*Include 循环引用: {chain}')
    if not source.is_file():
        raise FileNotFoundError(f'Abaqus 输入文件不存在: {source}')

    lines = []
    sources = [str(source)]
    with source.open('r', encoding='utf-8', errors='ignore') as f:
        raw = [ln.rstrip('\r\n') for ln in f]
    for line in raw:
        stripped = line.strip()
        if stripped.lower().startswith('*include'):
            kw_line = stripped[1:]
            head = kw_line.split(',')[0]
            kv = _parse_kv(kw_line[len(head):].lstrip(','))
            include_name = kv.get('input') or kv.get('file')
            if not include_name:
                raise ValueError(f'*Include 缺少 input/file: {source}: {line}')
            include_name = str(include_name).strip('"\'')
            child_lines, child_sources = _read_with_includes(
                source.parent / include_name, (*stack, source))
            lines.extend(child_lines)
            sources.extend(child_sources)
        else:
            lines.append(line)
    return lines, sources


def parse_inp(path):
    m = InpModel()
    raw, m.source_files = _read_with_includes(path)

    i = 0
    n = len(raw)
    cur_part = None          # 当前 part (Part)
    cur_inst = None          # 当前 instance
    in_assembly = False
    cur_mat = None
    cur_mat_key = None
    cur_amp = None
    cur_step = None
    cur_nset = None          # (dict, name, instance, is_generate, is_node)
    cur_surf = None
    pending_constraint = None  # '** Constraint: <name>' 注释暂存
    cur_constraint = None      # 当前关键字所属约束名
    last_point_part = None     # 最近 *Element,type=MASS 所属 Part
    last_point_eids = []       # 最近点单元块，用于紧随其后的 *Mass

    def close_data_block():
        nonlocal cur_nset, cur_surf, cur_mat_key, cur_amp
        cur_nset = None
        cur_surf = None
        cur_mat_key = None
        cur_amp = None

    while i < n:
        line = raw[i]
        i += 1
        if not line.strip():
            continue
        if line.startswith('**'):
            cm = re.match(r'^\*\*\s*Constraint:\s*(.+?)\s*$', line)
            if cm:
                pending_constraint = cm.group(1)
            continue
        if line.startswith('*'):
            kw_line = line[1:].strip()
            kw_name = kw_line.split(',')[0].strip().lower()
            kv = _parse_kv(kw_line[len(kw_line.split(',')[0]):].lstrip(','))
            close_data_block()
            # '** Constraint:' 注释归属于紧随其后的关键字块
            cur_constraint, pending_constraint = pending_constraint, None

            if kw_name == 'part':
                cur_part = Part(kv['name'])
                m.parts[cur_part.name] = cur_part
            elif kw_name == 'end part':
                cur_part = None
            elif kw_name == 'assembly':
                in_assembly = True
            elif kw_name == 'end assembly':
                in_assembly = False
            elif kw_name == 'instance':
                cur_inst = Instance(kv['name'], kv['part'])
                m.instances.append(cur_inst)
            elif kw_name == 'end instance':
                cur_inst = None
            elif kw_name == 'node':
                cur_nset = ('nodes', None, kv.get('instance'), False)
            elif kw_name == 'element':
                elem_type = kv.get('type', '').upper()
                if cur_part is not None:
                    if elem_type in ('MASS', 'ROTARYI'):
                        last_point_part = cur_part
                        last_point_eids = []
                        cur_nset = ('point_elems', elem_type,
                                    kv.get('elset'), False)
                    else:
                        cur_part.etype = elem_type
                        last_point_part = None
                        last_point_eids = []
                        cur_nset = ('elems', elem_type,
                                    kv.get('instance'), False)
                else:
                    # assembly 级元素 (MASS / ROTARYI 等): elset 名存入 name 槽
                    es = kv.get('elset')
                    if es:
                        m.asm_elem_types[es] = elem_type
                    cur_nset = ('elems', es, None, False)
            elif kw_name in ('nset', 'elset'):
                name = kv.get('nset') or kv.get('elset')
                gen = 'generate' in kv
                inst = kv.get('instance')
                kind = 'nset' if kw_name == 'nset' else 'elset'
                cur_nset = (kind, name, inst, gen)
            elif kw_name == 'beam section':
                # 数据行 1: 截面尺寸; 数据行 2: n1 方向矢量
                dl = []
                while i < n and len(dl) < 2:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    dl.append(ln.strip())
                    i += 1
                dims = _nums(dl[0]) if dl else []
                n1v = _nums(dl[1]) if len(dl) > 1 else []
                if cur_part is not None:
                    cur_part.beam_sections.append({
                        'elset': kv.get('elset'),
                        'material': kv.get('material'),
                        'section': (kv.get('section') or '').upper(),
                        'dims': dims, 'n1': n1v})
            elif kw_name == 'release':
                # 数据行: eid,端点(s1/s2),释放代码 (如 allm)
                while i < n:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    toks = [t.strip() for t in ln.rstrip(',').split(',')]
                    try:
                        eid = int(float(toks[0]))
                        if cur_part is not None:
                            cur_part.releases.append(
                                (eid, toks[1] if len(toks) > 1 else '',
                                 toks[2] if len(toks) > 2 else ''))
                    except (ValueError, IndexError):
                        pass
                    i += 1
            elif kw_name == 'mpc':
                # 数据行: 类型,从集,主集 (如 BEAM, _PickedSet171, _PickedSet172)
                while i < n:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    toks = [t.strip() for t in ln.rstrip(',').split(',')
                            if t.strip()]
                    if len(toks) >= 3:
                        m.mpcs.append({'type': toks[0].upper(),
                                       'slave': toks[1], 'master': toks[2],
                                       'constraint': cur_constraint})
                    i += 1
            elif kw_name == 'mass':
                # 数据行: 标量或 TYPE=ANISOTROPIC 的 [mx,my,mz]。
                while i < n:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    try:
                        vals = _nums(ln)
                        if str(kv.get('type', '')).upper() == 'ANISOTROPIC':
                            mass = (vals + [0.0, 0.0, 0.0])[:3]
                        else:
                            mass = [vals[0], vals[0], vals[0]]
                        if last_point_part is not None and last_point_eids:
                            for eid in last_point_eids:
                                last_point_part.point_mass[eid] = list(mass)
                        else:
                            # 兼容原 assembly 级标量质量结构。
                            m.asm_mass[kv.get('elset')] = mass[0]
                    except (ValueError, IndexError):
                        pass
                    i += 1
                last_point_part = None
                last_point_eids = []
            elif kw_name == 'rotary inertia':
                # 数据行: I11 I22 I33 I12 I13 I23
                while i < n:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    try:
                        m.asm_rotary[kv.get('elset')] = _nums(ln)
                    except (ValueError, IndexError):
                        pass
                    i += 1
            elif kw_name == 'nonstructural mass':
                # 数据行: 每单位长度/面积质量
                while i < n:
                    ln = raw[i]
                    if not ln.strip() or ln.startswith('**'):
                        i += 1
                        continue
                    if ln.strip().startswith('*'):
                        break
                    try:
                        m.nonstruct_mass.append({
                            'elset': kv.get('elset'),
                            'units': kv.get('units'),
                            'value': _nums(ln)[0]})
                    except (ValueError, IndexError):
                        pass
                    i += 1
            elif kw_name == 'dynamic':
                if cur_step is not None:
                    while i < n and (not raw[i].strip()
                                     or raw[i].startswith('**')):
                        i += 1
                    if i < n and not raw[i].strip().startswith('*'):
                        cur_step['dynamic'] = _nums(raw[i])
                        i += 1
            elif kw_name == 'damping' and cur_mat is not None:
                # Rayleigh 阻尼在 kwargs 中 (alpha/beta), 无数据行
                m.materials[cur_mat]['damping'] = {
                    k: float(v) for k, v in kv.items()
                    if k in ('alpha', 'beta')}
            elif kw_name == 'solid section':
                # 数据行在下一行 (桁架为截面积)
                data = ''
                while i < n and not raw[i].strip().startswith('*'):
                    if raw[i].strip() and not raw[i].startswith('**'):
                        data = raw[i].strip()
                    i += 1
                if cur_part is not None:
                    area = None
                    try:
                        vals = _nums(data)
                        area = vals[0] if vals else None
                    except ValueError:
                        pass
                    cur_part.sections.append(
                        (kv.get('elset'), kv.get('material'), area))
            elif kw_name == 'surface':
                name = kv.get('name')
                m.surfaces.setdefault(name, [])
                cur_surf = name
            elif kw_name == 'material':
                cur_mat = kv.get('name')
            elif kw_name == 'amplitude':
                cur_amp = kv.get('name')
                m.amplitudes[cur_amp] = []
                m.amplitude_options[cur_amp] = dict(kv)
            elif kw_name == 'step':
                cur_step = {'name': kv.get('name'), 'nlgeom': kv.get('nlgeom'),
                            'options': dict(kv), 'boundaries': [],
                            'loads': [], 'outputs': []}
                m.steps.append(cur_step)
            elif kw_name == 'end step':
                cur_step = None
            elif kw_name == 'static':
                if cur_step is not None:
                    # 下一行是步长控制
                    while i < n and (not raw[i].strip() or raw[i].startswith('**')):
                        i += 1
                    if i < n:
                        cur_step['static'] = _nums(raw[i])
                        i += 1
            elif kw_name == 'boundary':
                cur_nset = ('boundary', kv.get('amplitude'), {
                    'type': kv.get('type'), 'op': kv.get('op')}, False)
            elif kw_name == 'dload':
                cur_nset = ('dload', None, None, False)
            elif kw_name == 'dsload':
                cur_nset = ('dsload', None, None, False)
            elif kw_name == 'tie':
                m.ties.append(dict(kv))
                cur_nset = ('tie_surfs', len(m.ties) - 1, None, False)
            elif kw_name == 'coupling':
                m.couplings.append(dict(kv))
            elif kw_name == 'embedded element':
                m.embedded.append(dict(kv))
                cur_nset = ('embedded_sets', len(m.embedded) - 1, None, False)
            elif kw_name in MAT_KEYS and cur_mat is not None:
                cur_mat_key = MAT_KEYS[kw_name]
                m.materials[cur_mat].setdefault(cur_mat_key, [])
            elif kw_name in ('output', 'node output', 'element output',
                             'contact output', 'restart', 'preprint',
                             'heading', 'end instance', 'surface interaction',
                             'friction', 'surface behavior', 'cohesive behavior',
                             'damage initiation', 'damage evolution',
                             'damage stabilization', 'contact pair',
                             'el print', 'node print', 'monitor', 'controls',
                             'cload', 'dsflux', 'kinematic'):
                # 已知但不需要数据行的关键字 (cload/dload 等如需可扩展)
                if kw_name == 'kinematic':
                    pass
                cur_nset = ('skip', None, None, False) if kw_name in (
                    'cload', 'dload') else cur_nset
            else:
                # 未识别关键字: 安全跳过其数据行由下一轮 * 处理
                pass
            continue

        # ---------------- 数据行 ----------------
        if cur_nset is not None:
            kind, name, inst, gen = cur_nset
            try:
                if kind == 'nodes':
                    vals = _nums(line)
                    nid = int(vals[0])
                    xyz = tuple(vals[1:4]) + (0.0,) * (3 - len(vals[1:4]))
                    if cur_inst is not None:
                        pass  # instance 内不会出现 *Node 数据 (part 内才有)
                    elif cur_part is not None and not in_assembly:
                        cur_part.nodes[nid] = xyz
                    else:
                        m.asm_nodes[nid] = xyz
                elif kind == 'elems':
                    vals = _ids(line)
                    if cur_part is not None:
                        cur_part.elems[vals[0]] = vals[1:]
                        cur_part.elem_types[vals[0]] = name
                    elif name:
                        m.asm_elems.setdefault(name, {})[vals[0]] = vals[1:]
                elif kind == 'point_elems':
                    vals = _ids(line)
                    if cur_part is not None:
                        eid = vals[0]
                        cur_part.point_elems[eid] = {
                            'type': name,
                            'elset': inst,
                            'conn': vals[1:],
                        }
                        last_point_eids.append(eid)
                elif kind in ('nset', 'elset'):
                    ids = _ids(line)
                    if gen:
                        # generate: start, end, step
                        s, e = ids[0], ids[1]
                        st = ids[2] if len(ids) > 2 else 1
                        ids = list(range(s, e + 1, st))
                    tgt = None
                    if in_assembly and inst:
                        tgt = m.asm_node_nsets if kind == 'nset' else m.asm_elsets
                        tgt.setdefault(name, {})
                        tgt[name].setdefault(inst, [])
                        tgt[name][inst].extend(ids)
                    elif in_assembly:
                        tgt = m.asm_nsets if kind == 'nset' else m.asm_elsets
                        tgt.setdefault(name, []).extend(ids)
                    elif cur_part is not None:
                        tgt = cur_part.nsets if kind == 'nset' else cur_part.elsets
                        tgt.setdefault(name, []).extend(ids)
                elif kind == 'boundary':
                    vals = line.rstrip(',').split(',')
                    if len(vals) >= 3:
                        boundary = {
                            'set': vals[0].strip(),
                            'dof1': int(float(vals[1])),
                            'dof2': int(float(vals[2])),
                            'value': float(vals[3]) if len(vals) > 3 and vals[3].strip() else 0.0,
                            'amplitude': name,
                            'type': inst.get('type') if isinstance(inst, dict) else None,
                            'op': inst.get('op') if isinstance(inst, dict) else None,
                            'encastre': any('encastre' in v.lower() for v in vals[1:]),
                        }
                        if cur_step is not None:
                            cur_step['boundaries'].append(boundary)
                        else:
                            m.initial_boundaries.append(boundary)
                elif kind == 'dload':
                    vals = [value.strip() for value in line.rstrip(',').split(',')]
                    if cur_step is not None and len(vals) >= 3:
                        cur_step['loads'].append({
                            'region': vals[0],
                            'type': vals[1].upper(),
                            'value': float(vals[2]),
                            'parameters': [float(value) for value in vals[3:]
                                           if value],
                        })
                elif kind == 'dsload':
                    vals = line.rstrip(',').split(',')
                    if cur_step is not None and len(vals) >= 3:
                        cur_step['loads'].append({
                            'surface': vals[0].strip(),
                            'type': vals[1].strip(),
                            'value': float(vals[2]),
                            'parameters': [float(value) for value in vals[3:]
                                           if value.strip()],
                        })
                elif kind == 'tie_surfs':
                    # *Tie 数据行: slave 面, master 面
                    surfs = [t.strip() for t in line.rstrip(',').split(',') if t.strip()]
                    if len(surfs) >= 2:
                        m.ties[name]['slave'] = surfs[0]
                        m.ties[name]['master'] = surfs[1]
                elif kind == 'embedded_sets':
                    # *Embedded Element 数据行: 被嵌入的 elset (钢筋)
                    sets_ = [t.strip() for t in line.rstrip(',').split(',') if t.strip()]
                    m.embedded[name].setdefault('embedded', []).extend(sets_)
            except (ValueError, IndexError):
                pass
            continue

        if cur_surf is not None:
            parts = [t.strip() for t in line.rstrip(',').split(',') if t.strip()]
            if len(parts) == 2:
                m.surfaces[cur_surf].append((parts[0], parts[1].upper()))
            continue

        if cur_amp is not None:
            try:
                vals = _nums(line)
                m.amplitudes[cur_amp].extend(
                    [(vals[j], vals[j + 1]) for j in range(0, len(vals) - 1, 2)])
            except ValueError:
                pass
            continue

        if cur_mat_key is not None and cur_mat is not None:
            try:
                m.materials[cur_mat][cur_mat_key].append(_nums(line))
            except ValueError:
                pass
            continue

        # instance 变换行
        if cur_inst is not None:
            try:
                vals = _nums(line)
                if len(vals) == 3:
                    cur_inst.trans = tuple(vals)
                elif len(vals) == 7:
                    cur_inst.rot = tuple(vals)
            except ValueError:
                pass
            continue

    return m


# ---------------------------------------------------------------------------
# 装配: 变换 + 节点合并
# ---------------------------------------------------------------------------

def _rot_matrix(axis, angle_deg):
    ux, uy, uz = axis
    n = math.sqrt(ux * ux + uy * uy + uz * uz)
    if n == 0:
        return None
    ux, uy, uz = ux / n, uy / n, uz / n
    th = math.radians(angle_deg)
    c, s, t = math.cos(th), math.sin(th), 1 - math.cos(th)
    return [
        [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
        [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
        [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
    ]


def rotate_vector(v, inst):
    """将方向矢量按 instance 旋转 (不平移) 到全局坐标系"""
    if not inst.rot:
        return tuple(v)
    ax, ay, az, bx, by, bz, ang = inst.rot
    R = _rot_matrix((bx - ax, by - ay, bz - az), ang)
    if not R:
        return tuple(v)
    x, y, z = v[0], v[1], v[2]
    return (R[0][0] * x + R[0][1] * y + R[0][2] * z,
            R[1][0] * x + R[1][1] * y + R[1][2] * z,
            R[2][0] * x + R[2][1] * y + R[2][2] * z)


def transform_point(xyz, inst):
    # Abaqus *Instance 语义: 先平移, 再绕 (平移后坐标系中的) 轴 a→b 旋转
    # (旋转轴两点 a, b 以装配/全局坐标给出, 已用 6-15.inp 实测验证)
    tx, ty, tz = inst.trans
    x, y, z = xyz[0] + tx, xyz[1] + ty, xyz[2] + tz
    if inst.rot:
        ax, ay, az, bx, by, bz, ang = inst.rot
        R = _rot_matrix((bx - ax, by - ay, bz - az), ang)
        if R:
            px, py, pz = x - ax, y - ay, z - az
            x = R[0][0] * px + R[0][1] * py + R[0][2] * pz + ax
            y = R[1][0] * px + R[1][1] * py + R[1][2] * pz + ay
            z = R[2][0] * px + R[2][1] * py + R[2][2] * pz + az
    return (x, y, z)


MAX_NAME = 28  # Exodus/MOOSE 名称安全长度 (MOOSE 按 32 截断, 留余量防冲突)


def _shorten(s, budget=MAX_NAME):
    """超长名称智能缩短: 丢弃中间 token, 保留首 token + 末两个 token 保证区分度
    如 'AA_dinglaing_zongjin_D12' → 'AA_zongjin_D12'"""
    if len(s) <= budget:
        return s
    toks = s.split('_')
    if len(toks) > 2:
        s2 = toks[0] + '_' + '_'.join(toks[-2:])
        if len(s2) <= budget:
            return s2
    return s[:budget]


class NameSanitizer:
    """清洗名称并保证 ≤28 字符且全局唯一"""

    def __init__(self):
        self._seen = {}

    def __call__(self, name):
        s = re.sub(r'[^A-Za-z0-9_]', '_', str(name))
        s = _shorten(s)
        if s in self._seen and self._seen[s] != name:
            i = 2
            while f"{s[:MAX_NAME - 2]}_{i}" in self._seen:
                i += 1
            s = f"{s[:MAX_NAME - 2]}_{i}"
        self._seen[s] = name
        return s


_sanitize_global = NameSanitizer()


def sanitize(name):
    return _sanitize_global(name)


class GlobalMesh:
    """全局合并网格"""

    def __init__(self, tol):
        self.tol = tol
        self.coords = []          # global_id-1 -> (x,y,z)
        self._hash = {}
        self.node_map = {}        # (instance, local_id) -> global_id
        self.asm_node_map = {}    # assembly node id -> global_id
        self.merged_count = 0

    def add(self, key, xyz):
        h = (round(xyz[0] / self.tol), round(xyz[1] / self.tol),
             round(xyz[2] / self.tol))
        # 允许邻近 8 格碰撞检查
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    gid = self._hash.get((h[0] + dx, h[1] + dy, h[2] + dz))
                    if gid is not None:
                        gx, gy, gz = self.coords[gid - 1]
                        if (abs(gx - xyz[0]) <= self.tol and
                                abs(gy - xyz[1]) <= self.tol and
                                abs(gz - xyz[2]) <= self.tol):
                            self.node_map[key] = gid
                            self.merged_count += 1
                            return gid
        gid = len(self.coords) + 1
        self.coords.append(xyz)
        self._hash[h] = gid
        self.node_map[key] = gid
        return gid


def build_global_mesh(model, tol):
    gm = GlobalMesh(tol)
    gm.elem_origin = {}          # 全局单元 id -> (instance, part 单元 id)
    gm.elem_block = {}           # 全局单元 id -> block 名
    gm.block_beam = {}           # block -> {material, section, dims, n1}
    blocks = defaultdict(list)    # block_name -> [(elem_global_id, [gids])]
    block_etype = {}
    block_meta = {}               # block -> (part, material)
    elem_counter = 0

    # 元素 -> 材料: part.sections (elset, material, area)
    for inst in model.instances:
        part = model.parts.get(inst.part)
        if part is None:
            print(f"  WARN: instance {inst.name} 引用不存在的 part {inst.part}",
                  file=sys.stderr)
            continue
        # 局部节点 -> 全局
        local2g = {}
        for nid, xyz in part.nodes.items():
            g = gm.add((inst.name, nid), transform_point(xyz, inst))
            local2g[nid] = g
        # 元素分块
        el2mat = {}
        for elset, material, area in part.sections:
            ids = part.elsets.get(elset, [])
            for eid in ids:
                el2mat[eid] = (material, area)
        # 同一 part 多个截面积时需按面积分块 (如 zgjl: D8+D12)
        part_areas = {a for _, _, a in part.sections if a is not None}
        multi_area = len(part_areas) > 1
        # 梁截面: eid -> 截面索引
        el2sec = {}
        for si, bs in enumerate(part.beam_sections):
            for eid in part.elsets.get(bs['elset'] or '', []):
                el2sec[eid] = si
        sec_key2name = {}        # 截面 key -> block 名 (同 part 内)
        part_types = set(part.elem_types.values())
        for eid, conn in part.elems.items():
            elem_counter += 1
            mat, area = el2mat.get(eid, ('UNASSIGNED', None))
            beam_meta = None
            abaqus_type = part.elem_types.get(eid, part.etype)
            if abaqus_type not in ELEMENT_TYPE_MAP:
                raise ValueError(
                    f'不支持的 Abaqus 单元类型 {abaqus_type}: '
                    f'{inst.name}/{eid}')
            exo_type, expected_nodes = ELEMENT_TYPE_MAP[abaqus_type]
            if len(conn) != expected_nodes:
                raise ValueError(
                    f'{inst.name}/{eid} ({abaqus_type}) 需要 '
                    f'{expected_nodes} 个节点，实际 {len(conn)} 个')
            if part.beam_sections and abaqus_type == 'B31':
                si = el2sec.get(eid)
                if si is not None:
                    bs = part.beam_sections[si]
                    n1g = rotate_vector(bs['n1'], inst)
                    # 单元轴向 (全局)
                    p0 = gm.coords[local2g[conn[0]] - 1]
                    p1 = gm.coords[local2g[conn[1]] - 1]
                    t = [p1[k] - p0[k] for k in range(3)]
                    tl = math.sqrt(sum(c * c for c in t)) or 1.0
                    t = [c / tl for c in t]
                    # n1 投影到垂直于轴的平面 (同 Abaqus 投影语义)
                    d = sum(n1g[k] * t[k] for k in range(3))
                    proj = [n1g[k] - d * t[k] for k in range(3)]
                    pl = math.sqrt(sum(c * c for c in proj))
                    if pl < 1e-6:
                        # n1 ∥ 轴: 任取垂直矢量
                        proj = [-t[1], t[0], 0.0]
                        pl = math.sqrt(sum(c * c for c in proj))
                        if pl < 1e-6:
                            proj = [0.0, -t[2], t[1]]
                            pl = math.sqrt(sum(c * c for c in proj))
                    n1_eff = tuple(round(c / pl, 6) for c in proj)
                    key = (bs['material'], bs['section'],
                           tuple(bs['dims']), n1_eff,
                           tuple(round(c, 5) for c in t))
                    if key not in sec_key2name:
                        sect = (f"{bs['section']}"
                                f"_g{len(sec_key2name) + 1}")
                        clean = re.sub(r'[^A-Za-z0-9_]', '_', inst.part)
                        pshort = _shorten(
                            clean, MAX_NAME - len(str(bs['material']))
                            - len(sect) - 3)
                        sec_key2name[key] = sanitize(
                            f"{pshort}__{bs['material']}_{sect}")
                        gm.block_beam[sec_key2name[key]] = {
                            'material': bs['material'],
                            'section': bs['section'],
                            'dims': bs['dims'], 'n1': list(n1_eff)}
                    bname = sec_key2name[key]
                    beam_meta = gm.block_beam[bname]
                else:
                    clean = re.sub(r'[^A-Za-z0-9_]', '_', inst.part)
                    bname = sanitize(
                        f"{_shorten(clean, MAX_NAME - 13)}__UNASSIGNED")
            else:
                area_tag = (f"_A{int(round(area))}"
                            if (multi_area and area) else '')
                bname = sanitize(f"{_shorten(re.sub(r'[^A-Za-z0-9_]', '_', inst.part), MAX_NAME - len(str(mat)) - 2 - len(area_tag))}__{mat}{area_tag}")
            if len(part_types) > 1:
                bname = sanitize(
                    f"{_shorten(bname, MAX_NAME - len(exo_type) - 1)}_"
                    f"{exo_type}")
            blocks[bname].append((elem_counter, [local2g[n] for n in conn]))
            block_etype[bname] = exo_type
            block_meta[bname] = (inst.part, mat if beam_meta is None
                                 else beam_meta['material'])
            gm.elem_origin[elem_counter] = (inst.name, eid)
            gm.elem_block[elem_counter] = bname

    # assembly 级节点 (参考点)
    for nid, xyz in model.asm_nodes.items():
        g = gm.add(('__assembly__', nid), xyz)
        gm.asm_node_map[nid] = g

    # ---------------- nodesets ----------------
    nodesets = {}   # name -> sorted list of global ids

    # 1) assembly 级 nset (如 rp-1)
    for name, ids in model.asm_nsets.items():
        gids = [gm.asm_node_map.get(i) or gm.node_map.get(('__assembly__', i))
                for i in ids]
        nodesets[sanitize(name)] = sorted(g for g in gids if g)

    # 2) instance 限定的 nset/elset (instance 名按预算缩短, 保留集合名)
    for name, per_inst in model.asm_node_nsets.items():
        for inst, ids in per_inst.items():
            gids = [gm.node_map.get((inst, i)) for i in ids]
            gids = [g for g in gids if g]
            if gids:
                iclean = re.sub(r'[^A-Za-z0-9_]', '_', inst)
                nodesets[sanitize(
                    f"{name}__{_shorten(iclean, MAX_NAME - len(name) - 2)}")] = sorted(set(gids))
    for name, per_inst in model.asm_elsets.items():
        nodes = set()
        for inst, ids in per_inst.items():
            part = model.parts.get(model_inst_part(model, inst))
            if part is None:
                continue
            for eid in ids:
                for nid in part.elems.get(eid, []):
                    g = gm.node_map.get((inst, nid))
                    if g:
                        nodes.add(g)
        if nodes:
            nodesets[sanitize(f"ELSET_{name}")] = sorted(nodes)

    # 3) surfaces -> nodesets + 二维真实 sidesets
    sidesets = defaultdict(list)  # name -> [(global element id, Exodus side)]
    origin2geid = {origin: geid for geid, origin in gm.elem_origin.items()}
    for name, entries in model.surfaces.items():
        nodes = set()
        for elset_name, face in entries:
            # elset 是 instance 限定的 (elset_name 形如 __PickedSurf837_S2,
            # 在 asm_elsets 里有 instance 映射)
            per_inst = model.asm_elsets.get(elset_name, {})
            for inst, eids in per_inst.items():
                part_name = model_inst_part(model, inst)
                part = model.parts.get(part_name)
                if part is None:
                    continue
                for eid in eids:
                    conn = part.elems.get(eid)
                    if not conn:
                        continue
                    abaqus_type = part.elem_types.get(eid, part.etype)
                    face_nodes = ELEMENT_FACES.get(abaqus_type, {}).get(face)
                    sel = ([conn[f - 1] for f in face_nodes]
                           if face_nodes else conn)
                    for nid in sel:
                        g = gm.node_map.get((inst, nid))
                        if g:
                            nodes.add(g)
                    side = EXODUS_SIDE_MAP.get(abaqus_type, {}).get(face)
                    geid = origin2geid.get((inst, eid))
                    if side and geid:
                        sidesets[sanitize(name)].append((geid, side))
        if nodes:
            nodesets[sanitize(f"SURF_{name}")] = sorted(nodes)

    return (gm, blocks, block_etype, block_meta, nodesets,
            {name: sorted(set(entries))
             for name, entries in sidesets.items()})


def model_inst_part(model, inst_name):
    for inst in model.instances:
        if inst.name == inst_name:
            return inst.part
    return None


# ---------------------------------------------------------------------------
# *Tie 约束缝合 (等效 Abaqus adjust=yes): slave 面节点并入最近 master 节点
# ---------------------------------------------------------------------------

def apply_tie_stitch(model, gm, blocks, block_etype, nodesets, tie_tol,
                     blacklist=None):
    """返回 (合并对数, 跳过对数, merged_slave_roots)。

    Abaqus *Tie, adjust=yes 会把 slave 面节点吸附到 master 面。
    这里用并查集把距离 ≤ tie_tol 的 slave/master 节点合并，
    然后压实节点编号并重映射全部连接关系。
    blacklist: 禁止合并的 slave 节点 (防单元翻转回滚用)。"""
    blacklist = blacklist or set()
    ties = [t for t in model.ties if t.get('slave') and t.get('master')]
    if not ties:
        return 0, 0, set()

    nn = len(gm.coords)
    parent = list(range(nn + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    # 节点 → 单元索引 (防退化: 同单元内节点不得互并)
    # node_elems[根] 累积整个 union 的关联单元; members[根] 累积成员
    node_elems = defaultdict(set)
    elem_conn = {}
    for bname, elems in blocks.items():
        for ei, (_, conn) in enumerate(elems):
            elem_conn[(bname, ei)] = conn
            for g in conn:
                node_elems[g].add((bname, ei))
    members = {g: {g} for g in range(1, nn + 1)}

    # 节点 → 关联单元最小边长 (防止节点移动过远导致薄层单元翻转)
    HEX_EDGES = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                 (0,4),(1,5),(2,6),(3,7)]
    node_min_edge = {}
    hex_elems = []           # (bname, ei, conn) 便于精确翻转检查
    node_hexes = defaultdict(list)
    for bname, elems in blocks.items():
        if block_etype[bname] != 'HEX8':
            continue
        for ei, (_, conn) in enumerate(elems):
            hex_elems.append(conn)
            for g in conn:
                node_hexes[g].append(conn)
            for a, b in HEX_EDGES:
                pa, pb = gm.coords[conn[a]-1], gm.coords[conn[b]-1]
                L = math.dist(pa, pb)
                for g in (conn[a], conn[b]):
                    if L < node_min_edge.get(g, 1e30):
                        node_min_edge[g] = L

    import numpy as np
    _SIG = np.array(_HEX_SIGNS)

    def hex_min_corner_det(conn, override=None):
        """8 角点 det(J) 最小值; override={节点id: 新坐标} 用于试算"""
        p = []
        for g in conn:
            if override and g in override:
                p.append(override[g])
            else:
                p.append(gm.coords[g - 1])
        p = np.array(p)
        worst = 1e30
        for xi, eta, zeta in _HEX_CORNERS:
            ds = _SIG[:, 0] * (1 + _SIG[:, 1] * eta) * (1 + _SIG[:, 2] * zeta) / 8.0
            dt = (1 + _SIG[:, 0] * xi) * _SIG[:, 1] * (1 + _SIG[:, 2] * zeta) / 8.0
            du = (1 + _SIG[:, 0] * xi) * (1 + _SIG[:, 1] * eta) * _SIG[:, 2] / 8.0
            J = np.array([ds @ p, dt @ p, du @ p])
            worst = min(worst, np.linalg.det(J))
        return worst

    # master 节点空间哈希
    cell = max(tie_tol, 1.0)

    def hkey(p):
        return (int(p[0] // cell), int(p[1] // cell), int(p[2] // cell))

    merged, skipped = 0, 0
    merged_slaves = []
    for tie in ties:
        slave_key = 'SURF_' + sanitize(tie['slave'])
        master_key = 'SURF_' + sanitize(tie['master'])
        slaves = nodesets.get(slave_key, [])
        masters = nodesets.get(master_key, [])
        if not slaves or not masters:
            print(f"  WARN: tie {tie.get('name')} 面节点为空 "
                  f"({slave_key}:{len(slaves)}, {master_key}:{len(masters)})",
                  file=sys.stderr)
            continue
        grid = defaultdict(list)
        for g in masters:
            grid[hkey(gm.coords[g - 1])].append(g)
        for sg in slaves:
            if sg in blacklist:
                skipped += 1
                continue
            sp = gm.coords[sg - 1]
            k = hkey(sp)
            # 节点允许的最大移动距离: 关联单元最小边长的 1/4
            max_move = 0.25 * node_min_edge.get(sg, tie_tol)
            best, bestd = None, min(tie_tol, max_move)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for mg in grid.get((k[0] + dx, k[1] + dy, k[2] + dz), []):
                            mp = gm.coords[mg - 1]
                            d = math.dist(sp, mp)
                            if d <= bestd:
                                # 防退化检查
                                rs, rm = find(sg), find(mg)
                                if rs == rm:
                                    continue
                                if node_elems[rs] & node_elems[rm]:
                                    continue
                                best, bestd = rm, d
            if best is not None:
                # 精确翻转检查: 整个 union (共享一个位置) 试移到 best 坐标,
                # 检查 union 全部关联 hex 的角点 det
                # (已接受的合并即时更新 gm.coords, 后续检查看到最新几何)
                rs = find(sg)
                bc = gm.coords[best - 1]
                ok = True
                for key in node_elems.get(rs, set()):
                    conn = elem_conn[key]
                    if hex_min_corner_det(
                            conn, {g: bc for g in conn if find(g) == rs}) <= 1e-9:
                        ok = False
                        break
                if ok:
                    parent[rs] = best
                    node_elems[best] |= node_elems.get(rs, set())
                    node_elems[rs] = set()
                    mb = members.setdefault(best, {best})
                    mb |= members.get(rs, {rs})
                    members[rs] = set()
                    for g in mb:
                        gm.coords[g - 1] = bc  # union 即时生效
                    merged += 1
                    merged_slaves.append(sg)
                else:
                    skipped += 1
            else:
                skipped += 1

    if merged == 0:
        return merged, skipped, {}

    # 压实: 保留各连通分量的根节点; 坐标以实体(HEX8)成员为准,
    # 防止实体节点被合到钢筋节点坐标上导致单元畸变/翻转
    merged_roots = _compact_with_solid_coords(gm, blocks, block_etype,
                                              nodesets, parent)
    slave_new_ids = {sg: merged_roots[sg] for sg in merged_slaves
                     if sg in merged_roots}
    return merged, skipped, slave_new_ids


def _compact_with_solid_coords(gm, blocks, block_etype, nodesets, parent):
    """并查集压实 + 重编号, 坐标优先取实体单元成员节点。
    返回 {被并掉的节点原 id: 压实后的新 id}"""
    nn = len(gm.coords)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    solid_members = set()
    for b, elems in blocks.items():
        if block_etype[b] == 'HEX8':
            for _, conn in elems:
                solid_members.update(conn)

    canonical = {}
    for g in range(1, nn + 1):
        if g in solid_members:
            canonical.setdefault(find(g), g)
    # 根本身是实体成员时, 一律以根的坐标为准 (与精确翻转检查一致)
    for g in range(1, nn + 1):
        r = find(g)
        if r == g and r in solid_members:
            canonical[r] = r

    keep = sorted({find(g) for g in range(1, nn + 1)})
    newid = {g: i + 1 for i, g in enumerate(keep)}

    def remap(g):
        return newid[find(g)]

    merged_map = {g: remap(g) for g in range(1, nn + 1) if find(g) != g}

    gm.coords = [gm.coords[canonical.get(g, g) - 1] for g in keep]
    for bname, elems in blocks.items():
        blocks[bname] = [(eid, [remap(g) for g in conn]) for eid, conn in elems]
    for name, ids in nodesets.items():
        nodesets[name] = sorted({remap(g) for g in ids})
    return merged_map


# ---------------------------------------------------------------------------
# *Embedded Element 等效: 钢筋节点缝合至最近实体节点 (共享自由度, 完美粘结)
# ---------------------------------------------------------------------------

def apply_rebar_stitch(gm, blocks, block_etype, nodesets):
    """TRUSS 块节点并入最近的 HEX8 实体节点。

    等效 Abaqus *Embedded Element: 钢筋与混凝土共享位移 (完美粘结),
    钢筋 truss 刚度自然叠加到实体节点上。缝合后零长度 truss 单元被剔除。
    为避免 truss 节点并入 BC/荷载面节点后, MOOSE 由 nodeset 生成 sideset 时
    把 truss 的 0D 侧面混入 (Pressure BC 会因 0D 元素崩溃),
    属于 SURF_/固定 nodesets 的实体节点不作为缝合目标。
    返回 (缝合节点数, 剔除单元数, 最大缝合距离)。"""
    from scipy.spatial import cKDTree

    protected = set()
    for name, ids in nodesets.items():
        if name.startswith('SURF_'):
            protected.update(ids)

    # 全部实体节点 (判断 union 是否已连到结构)
    solid_nodes_all = sorted({g for b, elems in blocks.items()
                              if block_etype[b] == 'HEX8'
                              for _, conn in elems for g in conn})
    # 可作为缝合目标的实体节点 (排除 SURF_ 面节点)
    solid_nodes = [g for g in solid_nodes_all if g not in protected]
    truss_blocks = [b for b in blocks if block_etype[b] == 'TRUSS']
    if not truss_blocks or not solid_nodes:
        return 0, 0, 0.0

    import numpy as np
    solid_arr = np.array(solid_nodes)
    coords = np.array(gm.coords)
    tree = cKDTree(coords[solid_arr - 1])

    nn = len(gm.coords)
    parent = list(range(nn + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    truss_nodes = sorted({g for b in truss_blocks
                          for _, conn in blocks[b] for g in conn})
    # 初始合并可能把实体节点并入钢筋 union (钢筋 instance 先处理)。
    # 含实体成员的 union 已与结构相连 — 跳过, 避免实体节点被拉走翻转 hex
    solid_set = set(solid_nodes_all)
    solid_root = set()
    for g in range(1, nn + 1):
        if g in solid_set:
            solid_root.add(find(g))

    tp = coords[np.array(truss_nodes) - 1]
    dists, idxs = tree.query(tp)
    maxd = 0.0
    stitched = 0
    for g, d, ii in zip(truss_nodes, dists, idxs):
        rg = find(g)
        if rg in solid_root:
            continue  # 已含实体成员, 天然相连
        sg = int(solid_arr[ii])
        if rg == find(sg):
            continue
        parent[rg] = find(sg)
        solid_root.add(find(sg))
        stitched += 1
        maxd = max(maxd, float(d))

    # 压实重编号 (坐标以实体成员为准)
    _compact_with_solid_coords(gm, blocks, block_etype, nodesets, parent)

    # 删除零长度 / 近零长度 truss 单元
    dropped = 0
    coords_new = gm.coords
    for bname, elems in list(blocks.items()):
        new_elems = []
        for eid, conn in elems:
            if len(set(conn)) < 2:
                dropped += 1
                continue
            if block_etype[bname] == 'TRUSS':
                p, q = coords_new[conn[0] - 1], coords_new[conn[1] - 1]
                if math.dist(p, q) < 5.0:
                    dropped += 1
                    continue
            new_elems.append((eid, conn))
        blocks[bname] = new_elems
    return stitched, dropped, maxd


# ---------------------------------------------------------------------------
# *Embedded Element → MPC (LinearNodalConstraint) 导出
# 每个钢筋节点定位宿主 hex, 形函数插值: u_truss = Σ N_i · u_solid_i
# ---------------------------------------------------------------------------

_HEX_NODE_SIGNS = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                   (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]


def _hex_shape(xi, eta, zeta):
    return [0.125 * (1 + a * xi) * (1 + b * eta) * (1 + c * zeta)
            for a, b, c in _HEX_NODE_SIGNS]


def _hex_map(conn_pts, xi, eta, zeta):
    """三线性映射 + 雅可比。返回 (物理点, J 3x3)"""
    N = _hex_shape(xi, eta, zeta)
    p = [sum(N[i] * conn_pts[i][d] for i in range(8)) for d in range(3)]
    dN = [0.125 * a * (1 + b * eta) * (1 + c * zeta) for a, b, c in _HEX_NODE_SIGNS]
    eN = [0.125 * (1 + a * xi) * b * (1 + c * zeta) for a, b, c in _HEX_NODE_SIGNS]
    zN = [0.125 * (1 + a * xi) * (1 + b * eta) * c for a, b, c in _HEX_NODE_SIGNS]
    J = [[sum(dN[i] * conn_pts[i][d] for i in range(8)) for d in range(3)],
         [sum(eN[i] * conn_pts[i][d] for i in range(8)) for d in range(3)],
         [sum(zN[i] * conn_pts[i][d] for i in range(8)) for d in range(3)]]
    return p, J


def _solve3(J, b):
    """解 3x3 线性系统 J x = b (高斯消元)"""
    import copy as _c
    A = [row[:] + [b[i]] for i, row in enumerate(J)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(A[r][col]))
        if abs(A[piv][col]) < 1e-14:
            return None
        A[col], A[piv] = A[piv], A[col]
        for r in range(col + 1, 3):
            f = A[r][col] / A[col][col]
            for c in range(col, 4):
                A[r][c] -= f * A[col][c]
    x = [0.0] * 3
    for r in (2, 1, 0):
        x[r] = (A[r][3] - sum(A[r][c] * x[c] for c in range(r + 1, 3))) / A[r][r]
    return x


def point_in_hex(pt, conn_pts, tol=0.05, max_it=50, best_effort=False):
    """Newton 反解自然坐标。
    best_effort=False: 仅内部点命中返回 (xi,eta,zeta), 否则 None。
    best_effort=True: 不收敛也返回截断到 [-1,1] 的最近自然坐标 (投影)。"""
    center = [sum(p[d] for p in conn_pts) / 8.0 for d in range(3)]
    _, J0 = _hex_map(conn_pts, 0.0, 0.0, 0.0)
    dx = _solve3(J0, [pt[d] - center[d] for d in range(3)])
    if dx is None:
        xi = eta = zeta = 0.0
    else:
        xi, eta, zeta = dx
    converged = False
    for _ in range(max_it):
        p, J = _hex_map(conn_pts, xi, eta, zeta)
        r = [pt[d] - p[d] for d in range(3)]
        if max(abs(r[0]), abs(r[1]), abs(r[2])) < 1e-8:
            converged = True
            break
        dx = _solve3(J, r)
        if dx is None:
            break
        scale = 1.0
        m = max(abs(dx[0]), abs(dx[1]), abs(dx[2]))
        if m > 1.0:
            scale = 1.0 / m
        xi, eta, zeta = xi + dx[0] * scale, eta + dx[1] * scale, zeta + dx[2] * scale
        if best_effort:
            # 迭代中约束在略大范围内, 防发散
            xi = max(-1.5, min(1.5, xi))
            eta = max(-1.5, min(1.5, eta))
            zeta = max(-1.5, min(1.5, zeta))
        elif max(abs(xi), abs(eta), abs(zeta)) > 5.0:
            return None
    if converged and max(abs(xi), abs(eta), abs(zeta)) <= 1.0 + tol:
        return (xi, eta, zeta)
    if best_effort:
        return (max(-1.0, min(1.0, xi)),
                max(-1.0, min(1.0, eta)),
                max(-1.0, min(1.0, zeta)))
    return None


def emit_rebar_mpc(gm, blocks, block_etype, out_path, penalty=1e7,
                   vars_=('disp_x', 'disp_y', 'disp_z')):
    """为每个钢筋节点生成 LinearNodalConstraint (宿主 hex 形函数插值)。
    找不到宿主的节点退化为最近实体节点 w=1。返回 (总数, 插值数, 退化数)。"""
    import numpy as np
    from scipy.spatial import cKDTree

    solid_elems = []
    for b, elems in blocks.items():
        if block_etype[b] == 'HEX8':
            solid_elems.extend(elems)
    truss_nodes = sorted({g for b, e in blocks.items()
                          if block_etype[b] == 'TRUSS'
                          for _, conn in e for g in conn})
    C = np.array(gm.coords)

    # hex 质心 KD 树粗筛候选
    cents = np.array([C[np.array(conn) - 1].mean(axis=0)
                      for _, conn in solid_elems])
    sizes = np.array([C[np.array(conn) - 1].max(axis=0) -
                      C[np.array(conn) - 1].min(axis=0)
                      for _, conn in solid_elems])
    radii = np.linalg.norm(sizes, axis=1) / 2.0
    tree = cKDTree(cents)

    n_interp = n_degen = 0
    lines = []
    for g in truss_nodes:
        pt = C[g - 1]
        cand = tree.query_ball_point(pt, r=200.0)
        cand.sort(key=lambda i: np.linalg.norm(cents[i] - pt) - radii[i])
        weights = None
        for i in cand[:24]:
            conn = solid_elems[i][1]
            pts = [tuple(C[n - 1]) for n in conn]
            res = point_in_hex(pt, pts)
            if res is not None:
                # 截断到 [-1,1]: 外插形函数可能为负 → penalty 对角线为负 →
                # 矩阵不定 → MUMPS PC_FAILED
                xi, eta, zeta = (max(-1.0, min(1.0, c)) for c in res)
                weights = list(zip([int(n) for n in conn],
                                   _hex_shape(xi, eta, zeta)))
                break
        if weights is None:
            # 退化: 最近实体节点 w=1
            all_solid = sorted({n for _, conn in solid_elems for n in conn})
            d, ii = cKDTree(C[np.array(all_solid) - 1]).query(pt)
            weights = [(int(all_solid[ii]), 1.0)]
            n_degen += 1
        else:
            n_interp += 1
        prim = ' '.join(str(n) for n, _ in weights)
        wts = ' '.join(f'{w:.8g}' for _, w in weights)
        for v in vars_:
            lines.append(f"  [mpc_n{g}_{v[-1]}]\n"
                         f"    type = LinearNodalConstraint\n"
                         f"    variable = {v}\n"
                         f"    primary = '{prim}'\n"
                         f"    secondary_node_ids = '{g}'\n"
                         f"    weights = '{wts}'\n"
                         f"    penalty = {penalty:g}\n"
                         f"  []\n")
    with open(out_path, 'w') as f:
        f.write('[Constraints]\n')
        f.writelines(lines)
        f.write('[]\n')
    return len(truss_nodes), n_interp, n_degen


# ---------------------------------------------------------------------------
# HEX8 角点 Jacobian 检查 (检测单元翻转)
# ---------------------------------------------------------------------------

_HEX_SIGNS = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
              (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
_HEX_CORNERS = [(xi, eta, zeta) for xi in (-1, 1) for eta in (-1, 1)
                for zeta in (-1, 1)]


def find_inverted_hexes(gm, blocks, block_etype):
    """返回 {(bname, elem_index)}: 角点 det(J) <= 0 的单元"""
    import numpy as np
    C = np.array(gm.coords)
    bad = set()
    for bname, elems in blocks.items():
        if block_etype[bname] != 'HEX8':
            continue
        for ei, (_, conn) in enumerate(elems):
            p = C[np.array(conn) - 1]
            for xi, eta, zeta in _HEX_CORNERS:
                ds = np.array([a*(1+b*eta)*(1+c*zeta)
                               for a, b, c in _HEX_SIGNS]) / 8.0
                dt = np.array([(1+a*xi)*b*(1+c*zeta)
                               for a, b, c in _HEX_SIGNS]) / 8.0
                du = np.array([(1+a*xi)*(1+b*eta)*c
                               for a, b, c in _HEX_SIGNS]) / 8.0
                J = np.array([ds @ p, dt @ p, du @ p])
                if np.linalg.det(J) <= 0:
                    bad.add((bname, ei))
                    break
    return bad


# ---------------------------------------------------------------------------
# Exodus II 写出
# ---------------------------------------------------------------------------

def write_exodus(path, gm, blocks, block_etype, block_meta, nodesets,
                  sidesets, title):
    import netCDF4
    import numpy as np

    nc = netCDF4.Dataset(path, 'w', format='NETCDF4')

    nblk = len(blocks)
    nns = len(nodesets)
    nss = len(sidesets)
    nelem = sum(len(v) for v in blocks.values())
    block_names = sorted(blocks.keys())
    element_order = [eid for b in block_names for eid, _ in blocks[b]]
    exodus_elem_index = {eid: i for i, eid in enumerate(element_order, 1)}

    nc.createDimension('len_name', 256)
    nc.createDimension('len_line', 81)
    nc.createDimension('four', 4)
    nc.createDimension('time_step', None)
    nc.createDimension('num_dim', 3)
    nc.createDimension('num_nodes', len(gm.coords))
    nc.createDimension('num_elem', nelem)
    nc.createDimension('num_el_blk', nblk)
    nc.createDimension('num_node_sets', nns)
    nc.createDimension('num_side_sets', nss)

    nc.setncattr('api_version', np.float32(8.11))
    nc.setncattr('version', np.float32(8.11))
    nc.setncattr('floating_point_word_size', np.int32(8))
    nc.setncattr('file_size', np.int32(1))
    nc.setncattr('maximum_name_length', np.int32(32))
    nc.setncattr('int64_status', np.int32(0))
    nc.setncattr('title', title)

    def _names_var(vname, dimname, count, names):
        v = nc.createVariable(vname, 'S1', (dimname, 'len_name'))
        for i, nm in enumerate(names):
            b = nm.encode('ascii', 'replace')[:255]
            v[i, :len(b)] = list(np.frombuffer(b, dtype='S1'))
        return v

    # 时间
    vt = nc.createVariable('time_whole', 'f8', ('time_step',))
    vt[0] = 0.0

    # 坐标
    coords = np.array(gm.coords, dtype='f8')
    for d, vn in enumerate(('coordx', 'coordy', 'coordz')):
        v = nc.createVariable(vn, 'f8', ('num_nodes',))
        v[:] = coords[:, d]

    v = nc.createVariable('coor_names', 'S1', ('num_dim', 'len_name'))
    for d, nm in enumerate(('x', 'y', 'z')):
        v[d, 0] = np.frombuffer(nm.encode(), dtype='S1')[0]

    # node/elem num map (1:1)
    v = nc.createVariable('node_num_map', 'i4', ('num_nodes',))
    v[:] = np.arange(1, len(gm.coords) + 1)
    v = nc.createVariable('elem_num_map', 'i4', ('num_elem',))
    v[:] = np.array(element_order, dtype='i4')

    # 单元块
    nc.createVariable('eb_status', 'i4', ('num_el_blk',))[:] = 1
    vebp = nc.createVariable('eb_prop1', 'i4', ('num_el_blk',))
    vebp.setncattr('name', 'ID')
    vebp[:] = np.arange(1, nblk + 1)
    _names_var('eb_names', 'num_el_blk', nblk, block_names)

    eoff = 0
    for bi, bname in enumerate(block_names, start=1):
        elems = blocks[bname]
        nc.createDimension(f'num_el_in_blk{bi}', len(elems))
        npe = len(elems[0][1])
        nc.createDimension(f'num_nod_per_el{bi}', npe)
        vc = nc.createVariable(f'connect{bi}', 'i4',
                               (f'num_el_in_blk{bi}', f'num_nod_per_el{bi}'))
        vc.setncattr('elem_type', block_etype[bname])
        vc[:] = np.array([c for _, c in elems], dtype='i4')
        eoff += len(elems)

    # nodesets
    if nns:
        nc.createVariable('ns_status', 'i4', ('num_node_sets',))[:] = 1
        vnsp = nc.createVariable('ns_prop1', 'i4', ('num_node_sets',))
        vnsp.setncattr('name', 'ID')
        vnsp[:] = np.arange(1, nns + 1)
        ns_names = sorted(nodesets.keys())
        _names_var('ns_names', 'num_node_sets', nns, ns_names)
        for ni, name in enumerate(ns_names, start=1):
            ids = nodesets[name]
            nc.createDimension(f'num_nod_ns{ni}', len(ids))
            v = nc.createVariable(f'node_ns{ni}', 'i4', (f'num_nod_ns{ni}',))
            v[:] = np.array(ids, dtype='i4')

    # sidesets: elem_ss 使用 Exodus 文件内部单元序号，side_ss 为 1-based。
    if nss:
        nc.createVariable('ss_status', 'i4', ('num_side_sets',))[:] = 1
        vssp = nc.createVariable('ss_prop1', 'i4', ('num_side_sets',))
        vssp.setncattr('name', 'ID')
        vssp[:] = np.arange(1, nss + 1)
        ss_names = sorted(sidesets.keys())
        _names_var('ss_names', 'num_side_sets', nss, ss_names)
        for si, name in enumerate(ss_names, start=1):
            entries = sidesets[name]
            nc.createDimension(f'num_side_ss{si}', len(entries))
            ve = nc.createVariable(
                f'elem_ss{si}', 'i4', (f'num_side_ss{si}',))
            vs = nc.createVariable(
                f'side_ss{si}', 'i4', (f'num_side_ss{si}',))
            ve[:] = np.array(
                [exodus_elem_index[eid] for eid, _ in entries], dtype='i4')
            vs[:] = np.array([side for _, side in entries], dtype='i4')

    nc.close()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def _resolve_point_props(model, gm):
    """*Mass / *Rotary Inertia → 全局节点 (供 MOOSE NodalKernel)"""
    out = []
    for inst in model.instances:
        part = model.parts.get(inst.part)
        if part is None:
            continue
        for eid, elem in part.point_elems.items():
            if elem['type'] != 'MASS':
                continue
            if eid not in part.point_mass:
                raise ValueError(
                    f'{inst.name}/{eid} MASS 单元缺少对应 *Mass 数据')
            conn = elem['conn']
            if len(conn) != 1:
                raise ValueError(
                    f'{inst.name}/{eid} MASS 单元应连接 1 个节点，'
                    f'实际 {len(conn)} 个')
            gid = gm.node_map.get((inst.name, conn[0]))
            if not gid:
                raise ValueError(
                    f'{inst.name}/{eid} MASS 节点 {conn[0]} 无法映射')
            out.append({
                'kind': 'mass',
                'instance': inst.name,
                'element': eid,
                'source_node': conn[0],
                'gid': gid,
                'mass': list(part.point_mass[eid]),
                'xyz': list(gm.coords[gid - 1]),
            })
    for elset, val in model.asm_mass.items():
        for eid, conn in model.asm_elems.get(elset, {}).items():
            gid = (gm.asm_node_map.get(conn[0])
                   or gm.node_map.get(('__assembly__', conn[0])))
            if gid:
                out.append({'kind': 'mass', 'gid': gid,
                            'mass': [val, val, val],
                            'xyz': list(gm.coords[gid - 1])})
    for elset, vals in model.asm_rotary.items():
        for eid, conn in model.asm_elems.get(elset, {}).items():
            gid = (gm.asm_node_map.get(conn[0])
                   or gm.node_map.get(('__assembly__', conn[0])))
            if gid:
                out.append({'kind': 'rotary', 'gid': gid, 'inertia': vals,
                            'xyz': list(gm.coords[gid - 1])})
    return out


def write_nodal_mass_csv(point_props, prefix):
    """按 x/y/z 方向导出 NodalTranslationalInertia 所需的 4 列 CSV。"""
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    aggregated = {}
    for entry in point_props:
        if entry.get('kind') != 'mass':
            continue
        gid = entry['gid']
        if gid not in aggregated:
            aggregated[gid] = {
                'xyz': entry['xyz'],
                'mass': [0.0, 0.0, 0.0],
            }
        for direction in range(3):
            aggregated[gid]['mass'][direction] += entry['mass'][direction]

    outputs = {}
    for direction, suffix in enumerate(('x', 'y', 'z')):
        path = prefix.parent / f'{prefix.name}_{suffix}.csv'
        with path.open('w', encoding='utf-8', newline='') as f:
            for gid in sorted(aggregated):
                row = aggregated[gid]
                x, y, z = row['xyz']
                f.write(f'{x:.17g},{y:.17g},{z:.17g},'
                        f'{row["mass"][direction]:.17g}\n')
        outputs[suffix] = str(path)
    return outputs


def _resolve_releases(model, gm):
    """*RELEASE → 全局单元 id 列表 [{geid, end, code}]"""
    origin2geid = {v: k for k, v in gm.elem_origin.items()}
    out = []
    for p in model.parts.values():
        for eid, end, code in p.releases:
            for inst in model.instances:
                if inst.part != p.name:
                    continue
                geid = origin2geid.get((inst.name, eid))
                if geid:
                    out.append({'geid': geid, 'end': end, 'code': code})
    return out


def _resolve_nonstruct(model, gm):
    """*Nonstructural Mass → 每 block 元素计数 (供附加密度换算)"""
    origin2geid = {v: k for k, v in gm.elem_origin.items()}
    out = []
    for entry in model.nonstruct_mass:
        per_block = {}
        missing = 0
        for inst, eids in model.asm_elsets.get(entry['elset'], {}).items():
            for eid in eids:
                geid = origin2geid.get((inst, eid))
                if geid is None:
                    missing += 1
                    continue
                b = gm.elem_block.get(geid)
                per_block[b] = per_block.get(b, 0) + 1
        out.append({'elset': entry['elset'], 'units': entry['units'],
                    'value': entry['value'], 'per_block': per_block,
                    'missing': missing})
    return out


def main():
    ap = argparse.ArgumentParser(description='Abaqus .inp → Exodus II 转换器')
    ap.add_argument('--inp', required=True, help='Abaqus 输入文件')
    ap.add_argument('--out', required=True, help='输出 Exodus .e 文件')
    ap.add_argument('--report', help='JSON 报告输出路径')
    ap.add_argument('--merge-tol', type=float, default=1e-9,
                    help='跨 instance 节点合并容差 (默认 1e-9, 与模型单位一致)')
    ap.add_argument('--tie-tol', type=float, default=20.0,
                    help='*Tie 面节点缝合容差 (默认 20, 等效 adjust=yes)')
    ap.add_argument('--no-rebar-stitch', action='store_true',
                    help='跳过钢筋节点缝合 (保持原始几何, 配合 MOOSE embedded 约束)')
    ap.add_argument('--mpc', metavar='OUT_I',
                    help='导出钢筋 embedded MPC 约束片段 (LinearNodalConstraint)')
    ap.add_argument('--render-map', metavar='OUT_JSON',
                    help='导出原始钢筋几何→缝合节点渲染映射 (配合钢筋缝合)')
    ap.add_argument('--mass-csv-prefix', metavar='PREFIX',
                    help='导出 PREFIX_x/y/z.csv，供 NodalTranslationalInertia 使用')
    ap.add_argument('--add-nodeset', action='append', default=[],
                    metavar='NAME:BLOCK:COND',
                    help="按几何条件追加 nodeset, 如 'CAOSHEN_TOP:caoshen__con:z>=3940' "
                         "(COND 支持 x/y/z 与 >= <= == > <; 可重复; 用于接触/绑定界面)")
    ap.add_argument('--allow-unresolved-constraints', action='store_true',
                    help='约束 (*MPC 等) 无法解析时仅警告 '
                         '(默认: 硬失败退出 — 约束丢失会导致模型连接不完整)')
    args = ap.parse_args()

    print(f"[1/3] 解析 {args.inp} ...")
    model = parse_inp(args.inp)
    print(f"      parts={len(model.parts)} instances={len(model.instances)} "
          f"materials={len(model.materials)} surfaces={len(model.surfaces)} "
          f"ties={len(model.ties)} steps={len(model.steps)}")

    print(f"[2/3] 装配全局网格 (merge_tol={args.merge_tol}) ...")
    (gm, blocks, block_etype, block_meta,
     nodesets, sidesets) = build_global_mesh(model, args.merge_tol)
    print(f"      全局节点={len(gm.coords)} (合并 {gm.merged_count}) "
          f"单元={sum(len(v) for v in blocks.values())} "
          f"块={len(blocks)} nodesets={len(nodesets)} "
          f"sidesets={len(sidesets)}")

    # *Tie 绑定约束缝合 (等效 Abaqus adjust=yes), 带防翻转回滚
    import copy
    snapshot = (list(gm.coords), copy.deepcopy(blocks),
                copy.deepcopy(nodesets))
    blacklist = set()
    tie_merged = tie_skipped = 0
    for attempt in range(4):
        tie_merged, tie_skipped, slave_ids = apply_tie_stitch(
            model, gm, blocks, block_etype, nodesets, args.tie_tol,
            blacklist)
        bad = find_inverted_hexes(gm, blocks, block_etype)
        if not bad:
            break
        # 回滚: 把导致翻转的合并节点加入黑名单, 恢复快照重缝
        root_to_slave = {v: k for k, v in slave_ids.items()}
        n_before = len(blacklist)
        for bname, ei in bad:
            for g in blocks[bname][ei][1]:
                if g in root_to_slave:
                    blacklist.add(root_to_slave[g])
        print(f"      *Tie 缝合导致 {len(bad)} 个单元翻转, "
              f"回滚 {len(blacklist) - n_before} 处合并并重试 "
              f"(attempt {attempt + 1})")
        gm.coords, blocks, nodesets = (list(snapshot[0]),
                                       copy.deepcopy(snapshot[1]),
                                       copy.deepcopy(snapshot[2]))
    if tie_merged or tie_skipped:
        print(f"      *Tie 缝合: 合并 {tie_merged} 对, 跳过 {tie_skipped} "
              f"(tol={args.tie_tol}), 压实后节点={len(gm.coords)}")

    # *Embedded Element 等效: 钢筋节点缝合至最近实体节点
    has_truss = any(block_etype[b] == 'TRUSS' for b in blocks)
    if args.no_rebar_stitch or not has_truss:
        rebar_stitched, rebar_dropped, rebar_maxd = 0, 0, 0.0
        # 不缝合也快照全部钢筋几何 (v4: 求解不含钢筋, 渲染全量保留)
        truss_orig = {}
        truss_orig_coords = {}
        for b, elems in blocks.items():
            if block_etype[b] == 'TRUSS':
                truss_orig[b] = list(elems)
                truss_orig_coords[b] = {
                    eid: (gm.coords[conn[0] - 1], gm.coords[conn[1] - 1])
                    for eid, conn in elems}
    else:
        # 快照原始钢筋几何 (渲染映射用)
        truss_orig = {}
        truss_orig_coords = {}
        for b, elems in blocks.items():
            if block_etype[b] == 'TRUSS':
                truss_orig[b] = list(elems)
                truss_orig_coords[b] = {
                    eid: (gm.coords[conn[0] - 1], gm.coords[conn[1] - 1])
                    for eid, conn in elems}
        rebar_stitched, rebar_dropped, rebar_maxd = apply_rebar_stitch(
            gm, blocks, block_etype, nodesets)
    if rebar_stitched:
        print(f"      *Embedded 缝合: {rebar_stitched} 个钢筋节点并入实体节点 "
              f"(最大距离 {rebar_maxd:.1f}), 剔除零长度单元 {rebar_dropped}, "
              f"压实后节点={len(gm.coords)}")

    # 渲染映射: 原始直线钢筋单元 → 求解网格节点 (渲染时按位移场回弹)
    if args.render_map and (rebar_stitched or args.no_rebar_stitch):
        rmap = {'blocks': []}
        for b, orig_elems in truss_orig.items():
            final_by_eid = {eid: conn for eid, conn in blocks.get(b, [])}
            entry = {'name': b, 'elements': []}
            for eid, oconn in orig_elems:
                if eid in final_by_eid:
                    p0, p1 = truss_orig_coords[b][eid]
                    entry['elements'].append({
                        'eid': eid, 'p0': list(p0), 'p1': list(p1),
                        'n': final_by_eid[eid]})
            if entry['elements']:
                rmap['blocks'].append(entry)
        import json as _json
        with open(args.render_map, 'w') as f:
            _json.dump(rmap, f)
        n_el = sum(len(e['elements']) for e in rmap['blocks'])
        print(f"      渲染映射: {n_el} 个钢筋单元 → {args.render_map}")

    # MPC embedded 约束片段导出 (LinearNodalConstraint, 形函数插值)
    if args.mpc:
        n_tot, n_interp, n_degen = emit_rebar_mpc(
            gm, blocks, block_etype, args.mpc)
        print(f"      MPC 约束: {n_tot} 个钢筋节点 "
              f"(插值 {n_interp}, 最近点退化 {n_degen}) → {args.mpc}")

    # 追加几何条件 nodeset (接触/绑定界面, 供 TiedValueConstraint 等使用)
    for spec in args.add_nodeset:
        import re as _re
        m = _re.match(r'([^:]+):([^:]+):(x|y|z)(>=|<=|==|>|<)([-\d.eE]+)', spec)
        if not m:
            raise SystemExit(f'✗ 无法解析 --add-nodeset: {spec}')
        name, blk, axis, op, val = m.groups()
        val = float(val)
        if blk not in blocks:
            raise SystemExit(f'✗ 块不存在: {blk} (可选: {sorted(blocks)})')
        nodes = sorted({n for _, conn in blocks[blk] for n in conn})
        a = 'xyz'.index(axis)
        ops = {'>=': lambda c: c >= val, '<=': lambda c: c <= val,
               '==': lambda c: abs(c - val) < 1e-6,
               '>': lambda c: c > val, '<': lambda c: c < val}
        sel = [n for n in nodes if ops[op](gm.coords[n - 1][a])]
        nodesets[sanitize(name)] = sel
        print(f'      追加 nodeset {sanitize(name)}: block={blk} '
              f'{axis}{op}{val} → {len(sel)} 节点')

    # *MPC BEAM: 生成主-从刚性连杆 (spider) 单元, 等效刚体运动约束
    mpc_links = []
    constraints_audit = []       # 约束提取审核记录 → report['constraints']
    if model.mpcs:
        def _resolve_nset(name):
            if name in model.asm_node_nsets:
                g = []
                for inst, ids in model.asm_node_nsets[name].items():
                    for nid in ids:
                        gid = gm.node_map.get((inst, nid))
                        if gid:
                            g.append(gid)
                return sorted(set(g))
            if name in model.asm_nsets:
                g = [gm.asm_node_map.get(i)
                     or gm.node_map.get(('__assembly__', i))
                     for i in model.asm_nsets[name]]
                return sorted(set(x for x in g if x))
            return []
        eid_max = max((e for elems in blocks.values() for e, _ in elems),
                      default=0)
        mpc_group = {}           # (轴向, n1) -> block 名
        mpc_failures = []        # 约束解析失败记录 (默认硬失败)
        mpc_link_count = {}      # id(mpc) -> 生成连杆数
        for mpc in model.mpcs:
            cname = mpc.get('constraint') or '(无名约束)'
            if mpc['type'] != 'BEAM':
                mpc_failures.append(
                    f"{cname}: 不支持的 MPC 类型 {mpc['type']} "
                    f"({mpc['slave']} -> {mpc['master']})")
                continue
            slaves = _resolve_nset(mpc['slave'])
            masters = _resolve_nset(mpc['master'])
            if not slaves or not masters:
                mpc_failures.append(
                    f"{cname}: MPC {mpc['slave']}/{mpc['master']} "
                    f"节点集解析为空 (slave={len(slaves)} master={len(masters)})")
                continue
            mpc_link_count[id(mpc)] = 0
            m0 = masters[0]
            for s in slaves:
                if s == m0:
                    continue
                # 按连杆轴向分组, 取垂直 y_orientation (MOOSE 要求)
                p0 = gm.coords[s - 1]
                p1 = gm.coords[m0 - 1]
                t = [p1[k] - p0[k] for k in range(3)]
                tl = math.sqrt(sum(c * c for c in t)) or 1.0
                t = tuple(round(c / tl, 6) for c in t)
                ref = [1.0, 0.0, 0.0]
                d = sum(ref[k] * t[k] for k in range(3))
                proj = [ref[k] - d * t[k] for k in range(3)]
                pl = math.sqrt(sum(c * c for c in proj))
                if pl < 1e-6:
                    proj = [0.0, 1.0, 0.0]
                    d = sum(proj[k] * t[k] for k in range(3))
                    proj = [proj[k] - d * t[k] for k in range(3)]
                    pl = math.sqrt(sum(c * c for c in proj))
                n1 = [round(c / pl, 6) for c in proj]
                nb = mpc_group.setdefault(
                    (t, tuple(n1)), f"mpc_beam_links_g{len(mpc_group) + 1}")
                eid_max += 1
                blocks[nb].append((eid_max, [s, m0]))
                gm.elem_block[eid_max] = nb
                gm.block_beam[nb] = {'material': 'reactor',
                                     'section': 'RIGID', 'dims': [],
                                     'n1': n1}
                mpc_links.append({'slave': s, 'master': m0, 'block': nb,
                                  'constraint': mpc.get('constraint')})
                mpc_link_count[id(mpc)] += 1
        if mpc_links:
            for nb in set(gm.elem_block[e] for e in gm.elem_block
                          if gm.elem_block[e].startswith('mpc_beam_links')):
                block_etype[nb] = 'BEAM2'
                block_meta[nb] = ('__mpc__', 'RIGID')
            print(f"      MPC BEAM: 生成 {len(mpc_links)} 根刚性连杆 "
                  f"({len(mpc_group)} 个方向组)")
        # 约束完整性: 任何约束丢失都是致命错误 (除非显式允许)
        for mpc in model.mpcs:
            cst = {'name': mpc.get('constraint'), 'type': mpc['type'],
                   'slave': mpc['slave'], 'master': mpc['master'],
                   'links': mpc_link_count.get(id(mpc), 0)}
            cst['status'] = 'ok' if cst['links'] > 0 else 'failed'
            constraints_audit.append(cst)
        if mpc_failures:
            for f in mpc_failures:
                print(f"      \u2717 约束转换失败: {f}", file=sys.stderr)
            if not args.allow_unresolved_constraints:
                print("      约束是模型连接完整性必需字段, 拒绝静默丢弃。"
                      " 如确认可忽略, 加 --allow-unresolved-constraints",
                      file=sys.stderr)
                sys.exit(2)
            print("      WARN: --allow-unresolved-constraints 已启用, "
                  "以上约束被丢弃", file=sys.stderr)

    point_props = _resolve_point_props(model, gm)
    point_mass_nodes = sorted({entry['gid'] for entry in point_props
                               if entry.get('kind') == 'mass'})
    if point_mass_nodes:
        nodesets[sanitize('POINT_MASS')] = point_mass_nodes
    mass_csv = (write_nodal_mass_csv(point_props, args.mass_csv_prefix)
                if args.mass_csv_prefix else {})

    print(f"[3/3] 写出 Exodus: {args.out}")
    write_exodus(args.out, gm, blocks, block_etype, block_meta, nodesets,
                 sidesets, args.inp)

    # 报告
    xs = [c[0] for c in gm.coords]
    ys = [c[1] for c in gm.coords]
    zs = [c[2] for c in gm.coords]
    report = {
        'source': args.inp,
        'source_files': model.source_files,
        'merge_tol': args.merge_tol,
        'tie_tol': args.tie_tol,
        'tie_merged': tie_merged,
        'rebar_stitched': rebar_stitched,
        'rebar_dropped': rebar_dropped,
        'rebar_max_stitch_dist': rebar_maxd,
        'num_nodes': len(gm.coords),
        'num_merged': gm.merged_count,
        'num_elems': sum(len(v) for v in blocks.values()),
        'bbox': {'x': [min(xs), max(xs)], 'y': [min(ys), max(ys)],
                 'z': [min(zs), max(zs)]},
        'blocks': {b: {'count': len(blocks[b]), 'type': block_etype[b],
                       'part': block_meta[b][0], 'material': block_meta[b][1]}
                   for b in sorted(blocks)},
        'nodesets': {n: len(v) for n, v in sorted(nodesets.items())},
        'sidesets': {n: len(v) for n, v in sorted(sidesets.items())},
        'materials': {k: v for k, v in model.materials.items()},
        'amplitudes': model.amplitudes,
        'amplitude_options': model.amplitude_options,
        'initial_boundaries': model.initial_boundaries,
        'steps': model.steps,
        'ties': model.ties,
        'couplings': model.couplings,
        'embedded': model.embedded,
        'beam_sections': getattr(gm, 'block_beam', {}),
        'mpcs': model.mpcs,
        'mpc_links': mpc_links,
        'constraints': constraints_audit,
        'releases': {p.name: p.releases
                     for p in model.parts.values() if p.releases},
        'releases_global': _resolve_releases(model, gm),
        'point_mass': point_props,
        'point_mass_csv': mass_csv,
        'nonstructural_mass': _resolve_nonstruct(model, gm),
    }
    if args.report:
        with open(args.report, 'w', encoding='utf-8') as f:
            json.dump(report, f, ensure_ascii=False, indent=2, default=str)
        print(f"      报告: {args.report}")
    print("完成 ✓")


if __name__ == '__main__':
    main()
