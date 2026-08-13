# PR-RG-400gal-X 位移、应力和加速度三场复算记录

> 复算日期：2026-08-13
>
> 状态：`prototype`
>
> 适用范围：复现 `demo-process` 既有线性梁计算并补齐结果场，不是 Abaqus 等价验证或工程安全结论

## 1. 输入与来源

| 对象 | 版本或 SHA-256 |
|---|---|
| Abaqus 输入 `/home/kevin/Abaqus/PR-RG-400gal-X/PR-RG-400gal-X.inp` | `c746e1766ee39de50c2999223c65cbc7eef3c9e2805f11dca087af9e2c9ccf35` |
| DamSafetyApp 基线 | `2188251fc1c53929a3da2c74a185d5533e469bf0` |
| `demo-process` 来源提交 | `0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8` |
| 原 `beam_direct_solver.py` | `bb61a35440399705a5d6c97d6a43225e0f007bc91da94cc66fb7e6931fa6ee12` |
| DamSafetyApp `beam_field_solver.py` | `e8a0ccee9dcea9bfb43d407658e58ec3ecf9929ca0a698a639c5f8c89173703c` |

原输入只读，副本、网格、日志和大结果全部位于 `DamSafetyApp/.build/cases/pr-rg-400gal-x/`。`demo-process` 的源码和既有 `outputs/`、`renders/` 均未修改。

## 2. 模型与计算定义

| 项目 | 数值 |
|---|---:|
| 节点 | 701 |
| 梁/MPC 单元 | 768 |
| Exodus block | 188 |
| 总自由度 | 4,206 |
| 消元后自由度 | 4,062 |
| MPC 消元自由度 | 72 |
| 基底加速度自由度 | 12 |
| 时间步长 | 0.01 s |
| 计算时长 | 65 s，6,500 步 |
| 输出间隔 | 0.05 s，共 1,301 个状态（含 t=0） |

复算沿用既有 `demo-process` 算法和参数：三维 Timoshenko 梁、一致质量、Rayleigh 阻尼 `α=0.0702`、`β=0.005697`、集中质量/转动惯量、MPC BEAM 刚体消元和 Newmark 常加速度法。前六阶频率为 `1.076、1.076、1.141、3.667、3.667、4.261 Hz`。

既有结果日志明确显示铰接凝聚数为 0，原 MOOSE 输入也注明 `*RELEASE` 未模拟。为保持复算可比性，本次显式使用 `--no-releases`；这是一项已知连接近似，不表示 Abaqus `*RELEASE` 已迁移。

## 3. 三场恢复语义

- 位移场：节点 `disp_x/y/z` 和 `rot_x/y/z`；
- 加速度场：Newmark 状态恢复的节点绝对 `accel_x/y/z` 和 `rot_accel_x/y/z`；
- 应力场：由梁端局部内力恢复 `axial_stress`、`bending_stress`、`torsional_shear`，并形成截面最外缘 `vonmises_stress` 包络。

应力包络采用 `sqrt((|N|/A + |My|cz/Iy + |Mz|cy/Iz)^2 + 3(Tc/J)^2)`。它适合梁模型热点筛查和可视化，但不是 Abaqus 梁截面点或实体积分点应力。

## 4. 结果

| 指标 | 峰值 | 时间/位置 |
|---|---:|---|
| 顶点绝对水平位移 | `-1590.8877 mm` | 10.45 s |
| 可视化相对位移幅值 | `368.430 mm` | 全时程最大值 |
| 绝对加速度幅值 | `19733.8561 mm/s²` | 25.80 s，节点 161，坐标 `[2806.16284, 0.000124, 9787.70898] mm` |
| 梁 von Mises 等效应力包络 | `591.8902 MPa` | 23.45 s，block `insulator___aluminium_I_g166`，单元 696 |

`591.8902 MPa` 已超过常见铝材/结构钢的线弹性适用范围。由于当前模型是线弹性的，该数值只能作为非线性复核触发项，不能直接作为真实应力或安全裕度结论。

## 5. 交付文件

服务器目录：

```text
/home/kevin/DamSafetyApp/.build/cases/pr-rg-400gal-x/deliverables/
```

| 文件 | 内容 | SHA-256 |
|---|---|---|
| `PR-RG-400gal-X_displacement.e` | 6 个节点位移/转角变量 | `02e930095d4b3d7aa9c3478489199cde6fd7bcf0c12186267ac2081749c70afa` |
| `PR-RG-400gal-X_displacement.mp4` | 相对位移云图与放大变形 | `4c1f9a7a4a40e9068091a73afeaa9e64442e459673362c67c66a1fdbb8646a67` |
| `PR-RG-400gal-X_stress.e` | 4 个梁截面应力包络变量 | `7c9f96a1dcb502025cd9b8ba7e0860134ba8b7c171922f7e078a20e7543336d3` |
| `PR-RG-400gal-X_stress.mp4` | 梁 von Mises 等效应力包络云图 | `1b3187af6e4889c866ed030cddafc2914c7188a65e43076de93b72c7e8b6426c` |
| `PR-RG-400gal-X_acceleration.e` | 6 个节点平动/转动加速度变量 | `0729ebcba0ed19d304c053e6267def2066394a92d54a04251ae3caa750a3bd1d` |
| `PR-RG-400gal-X_acceleration.mp4` | 绝对加速度幅值云图 | `6366e1a3e58996a0492475d8d6ec0b747a8b04abb539a4f451b3c6caad0f9383` |
| `PR-RG-400gal-X_summary.csv` | 三场峰值时程摘要 | `0ab93ec6d2691478a20078370b16734bf7c8bfa09b4b66b904aafef310701e5f` |

每个 `.e` 都包含 701 节点、768 单元和 1,301 个相同时间状态；所有数值均为有限值。三份 MP4 均为 1600×900、H.264、yuv420p、20 fps、40 帧、2.0 s。

## 6. 验证

- DamSafetyApp 位移结果与 `demo-process/outputs/pr_rg_400gal_x/pr_rg_400gal_x_out.e` 的时间、坐标和六个节点变量逐数组比较，最大绝对差均为 0；
- 12 个基底节点的输出 `accel_x` 与 `RG-X` 幅值乘缩放系数逐状态比较，最大绝对差为 0；
- 应力恢复最小轴拉测试得到解析值，工具测试共 5 项通过；
- ParaView 正确识别三份 Exodus 的变量，三份视频分别在独立 pvpython 进程生成；
- 位移、应力和加速度代表帧均已人工检查，模型、标题、时间和单一对应色标显示正常。

## 7. 已知限制

- 完整计算由 DamSafetyApp 中迁移并扩展的 SciPy 稀疏线性梁积分器完成，不是 `blackbear-opt` 求解；这是为了严格复现既有 `demo-process` 结果；
- 未模拟 Abaqus `*RELEASE`，未处理材料屈服、塑性、接触或几何非线性；
- 没有 Abaqus `.odb` 或同工况结果可做应力、加速度数值等价对照；
- MP4 每 33 个 Exodus 状态抽取一帧，只用于浏览；完整数据以 `.e` 和摘要 CSV 为准；
- 工程交付前应确认材料强度、`*RELEASE` 连接语义和所需应力口径，并执行 Abaqus/MOOSE 或试验对照。
