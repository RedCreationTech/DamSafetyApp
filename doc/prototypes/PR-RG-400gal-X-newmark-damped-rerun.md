# PR-RG-400gal-X Newmark 数值阻尼重算与全场输出记录

> 复算日期：2026-08-17
>
> 状态：`prototype`
>
> 适用范围：响应"峰值应力比 Abaqus 隐式动力结果偏大"的反馈，用带算法阻尼的 Newmark 参数重算并补齐全场输出；不是 Abaqus 等价验证或工程安全结论

## 1. 背景与动机

客户反馈：既有复算（见 [PR-RG-400gal-X-three-field-rerun.md](PR-RG-400gal-X-three-field-rerun.md)）应力总体趋势正确，但峰值比 Abaqus 大。Abaqus 隐式动力分析默认使用 HHT-α 积分（含数值阻尼），而既有复算使用 Newmark 常加速度法（γ=0.5、β=0.25），完全没有数值阻尼，高频人为振荡不被衰减，推高了加速度与应力的包络峰值。

本次把求解器 Newmark γ/β 改为 CLI 可配参数，用 γ=0.55、β=0.2756 重算完整 65 s，并把该线弹性梁模型能真实产出的全部场合并输出到单份 Exodus。

γ=0.55、β=0.2756 与 HHT-α 的 α=-0.05 等价：γ=0.5-α=0.55，β=(1-α)²/4=0.275625≈0.2756。

## 2. 输入与来源

| 对象 | 版本或 SHA-256 |
|---|---|
| Abaqus 输入 `/home/kevin/Abaqus/PR-RG-400gal-X/PR-RG-400gal-X.inp` | `c746e1766ee39de50c2999223c65cbc7eef3c9e2805f11dca087af9e2c9ccf35`（与上次一致，未变） |
| DamSafetyApp 基线 | `b371dcea62332c9e772c3d4e41aac63d428db5c3` + 本次工作区改动 |
| 改造前 `beam_field_solver.py` | `e8a0ccee9dcea9bfb43d407658e58ec3ecf9929ca0a698a639c5f8c89173703c` |
| 改造后 `beam_field_solver.py` | `0aa9884ce4cfd7f7151766f4adf2b439ece6eac81dc3540a31daa2c16a120022` |
| 改造后 `test_beam_field_solver.py` | `728adb787053989497025d5a97d6f07660a410c595cb3c4cbb2d163ee5e7295b` |

原输入只读。模型、网格、Rayleigh 阻尼（α=0.0702、β=0.005697）、MPC 消元、基底加速度施加、`--no-releases`（已知连接近似）均与上次完全一致，唯一变化是 Newmark γ/β 与输出内容。

## 3. 求解器改造

`tools/abaqus/beam_field_solver.py`：

- 新增 `--newmark-gamma`/`--newmark-beta`（默认 0.5/0.25，常加速度法不变）；积分系数抽取为 `newmark_coefficients()`；
- 新增 `--allfields-out`：单份合并 Exodus，节点场与单元场写进同一文件；
- 节点场新增 `rf_x/y/z`：由 R = M·a + C·v + K·u 在约束平动自由度上恢复（外载为基底运动激励，约束处残差即反力），非约束节点置 0；为此时程循环增加速度帧记录与全自由度速度重构；
- 单元场新增 `s_max/mid/min_principal`、`axial_strain`、`bending_strain`（见第 4 节口径）；
- Exodus 写出重构为 `_write_nodal_fields`/`_write_element_fields` 共用内部函数，三份独立输出行为不变。

`test/tools/test_beam_field_solver.py`：新增 5 项测试——默认参数与旧内联公式逐项一致、默认参数高频自由振动幅值守恒、γ=0.55/β=0.2756 高频幅值衰减、主应力/应变恢复（含 σ-τ 组合与 vonmises 自洽）、约束反力恢复。`tools/abaqus/run-tests.sh` 全部 11 项通过。

## 4. 积分参数对照

| 项目 | 上次复算 | 本次重算 |
|---|---|---|
| Newmark γ | 0.5 | 0.55 |
| Newmark β | 0.25 | 0.2756 |
| 算法阻尼 | 无 | 有（等价 HHT-α α=-0.05） |
| 精度阶 | 二阶 | 一阶（γ>0.5 的固有损失，预期行为） |
| Rayleigh α/β | 0.0702 / 0.005697 | 相同 |
| dt / 步数 / 时长 | 0.01 s / 6,500 / 65 s | 相同 |
| `--no-releases` | 是 | 是 |

单元场应力口径（与既有恢复自洽，σ 为截面最外缘正应力包络 axial+bending，τ 为扭转剪应力 T·c/J）：

- `vonmises_stress` = sqrt(σ² + 3τ²)（不变）；
- `s_max/min_principal` = σ/2 ± sqrt(σ²/4 + τ²)，`s_mid_principal` = 0（平面应力主值；恒有 σ1≥0≥σ3，故 0 为中值）。恒等式 sqrt(σ1²-σ1σ3+σ3²) = sqrt(σ²+3τ²) 保证与 vonmises 口径严格自洽，合并文件实测最大偏差 1.1e-13 MPa；
- `axial_strain` = axial_stress/E，`bending_strain` = bending_stress/E（线弹性外缘应变，与应力同口径）。

## 5. 新旧峰值对比

同口径（顶点 = 点质量节点 701；应力/加速度为全场包络）：

| 指标 | 上次 γ=0.5/β=0.25 | 本次 γ=0.55/β=0.2756 | 变化 |
|---|---:|---:|---:|
| 顶点绝对水平位移 disp_x 最小值 | -1590.8877 mm @ 10.45 s | -1588.3281 mm @ 10.45 s | -0.16% |
| 顶点 \|accel_x\| 峰值 | 10916.4549 mm/s² @ 25.75 s | 10567.3151 mm/s² @ 25.75 s | -3.2% |
| 全场绝对加速度幅值峰值 | 19733.8561 mm/s² @ 25.80 s（节点 161） | 19180.9228 mm/s² @ 25.80 s（节点 161） | -2.8% |
| 梁 von Mises 应力包络峰值 | 591.8902 MPa @ 23.45 s | 578.1414 MPa @ 23.45 s | -2.3% |

应力峰值位置不变（block `insulator___aluminium_I_g166`）。以低频为主的位移几乎不变，高频含量更高的加速度与应力峰值被压低 2%~3%，与算法阻尼只衰减高频分量的预期一致。峰值时刻均未移动。

## 6. 交付文件

新结果位于独立子目录，既有交付物未覆盖：

```text
/home/kevin/DamSafetyApp/.build/cases/pr-rg-400gal-x/deliverables/newmark-g0.55-b0.2756/
```

| 文件 | 内容 | SHA-256 |
|---|---|---|
| `PR-RG-400gal-X_allfields.e` | **合并全场**：15 个节点变量 + 10 个单元变量（见下） | `c56ec088eed93e3d973c319a8a6b670d64a0b5574a8ab872380f5f9936ef0608` |
| `PR-RG-400gal-X_displacement.e` | 6 个节点位移/转角变量 | `abd83ece0df5d5e80e3d96343d7ec585ea48af873ec7819635c2dcd11540af09` |
| `PR-RG-400gal-X_stress.e` | 4 个梁截面应力包络变量 | `45ae815dc0bfb1cac14079afeb8a0efbe2867a3f4d9317ed274f855c73d4a137` |
| `PR-RG-400gal-X_acceleration.e` | 6 个节点平动/转动加速度变量 | `fd406aca9758040c8eba5484967f80bd300ee9b43f8099ec55fdf086bb163223` |
| `PR-RG-400gal-X_summary.csv` | 峰值时程摘要 | `278626d094b115f952a8828c2eafadc188692b5e0965ea94c65699a7b3653595` |

计算日志：`.build/cases/pr-rg-400gal-x/logs/full-solve-newmark-g0.55-b0.2756.log`。

`PR-RG-400gal-X_allfields.e` 内容（701 节点、768 单元、1,301 个时间状态，含 t=0，间隔 0.05 s）：

- 节点场：`disp_x/y/z`、`rot_x/y/z`、`accel_x/y/z`、`rot_accel_x/y/z`、`rf_x/y/z`（仅 12 个基底约束节点非零，峰值 \|rf\| ≈ 7.19e5 N）；
- 单元场：`axial_stress`、`bending_stress`、`torsional_shear`、`vonmises_stress`、`s_max_principal`、`s_mid_principal`、`s_min_principal`、`axial_strain`、`bending_strain`。

**未产出的客户清单场**：`DAMAGEC`、`DAMAGET`、`PE`、`PEEQ`、`PEMAG`、`AC_YIELD` 是塑性/损伤相关场，本模型为线弹性，这些量在物理上不存在，不以零场充数。若需要，应引入塑性/损伤本构后另行计算。

## 7. 验证

- `tools/abaqus/run-tests.sh`：11 项测试全部通过（含新增 5 项）；
- 回归：改造后以默认参数完整重算 65 s，位移/加速度/应力三份 Exodus 与既有交付物逐数组比较，最大绝对差均为 0（复现能力不回归）；
- 合并 Exodus 用 Python 读回校验：701 节点、768 单元、1,301 状态、0.05 s 间隔；15 个节点变量、10 个单元变量命名齐全（snake_case）；全部数值有限；
- 主应力与 vonmises 自洽性：峰值 block 上 sqrt(σ1²-σ1σ3+σ3²) 与 `vonmises_stress` 最大偏差 1.1e-13 MPa；`axial_stress/axial_strain` 比值中位数 = 72000 MPa，与铝材 E 一致；
- 基底反力仅出现在 12 个基底约束节点，非约束节点全为 0；
- 前六阶频率与上次相同（1.076、1.076、1.141、3.667、3.667、4.261 Hz），模型未变。

## 8. 已知限制

- γ=0.55>0.5 时 Newmark 降为一阶精度，这是引入算法阻尼的固有代价；低频响应（位移）受影响极小（-0.16%），但严格收敛性需以 dt 减半复核；
- 本次峰值仍比 Abaqus 高还是低未做数值对标——没有 Abaqus `.odb` 同工况结果，差异方向和幅度需要正式对照确认；
- 仍为线弹性模型：`vonmises_stress` 峰值 578 MPa 已超常见铝材/结构钢线弹性适用范围，只能作为非线性复核触发项；塑性/损伤场（DAMAGEC/DAMAGET/PE/PEEQ/PEMAG/AC_YIELD）未产出（见第 6 节）；
- 未模拟 Abaqus `*RELEASE`（沿用 `--no-releases`），基底反力为运动激励下的约束残差，未与 Abaqus RF 输出对标；
- 求解器仍是 DamSafetyApp 内的 SciPy 稀疏线性梁积分器，不是 `blackbear-opt`；
- 未生成新视频；如需渲染可对 `deliverables/newmark-g0.55-b0.2756/` 下 Exodus 复用 `render_beam_fields.py`。
