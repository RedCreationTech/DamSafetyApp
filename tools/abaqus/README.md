# Abaqus 转换工具

本目录的 `abaqus2exodus.py` 源自 `demo-process@0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8`，原脚本 SHA-256 为 `7d01cdc2560b5fc6d3eb089c8fece6a3d066df3f1d9192844c754d0bebda1c19`。DamSafetyApp 在保留原有 `B31/C3D8R/T3D2` 能力的基础上增加：

- 同一 Part 内逐单元保存类型，支持混合 `CPS4/CPS4R/CPS3`；
- `CPS4/CPS4R → QUAD4`、`CPS3 → TRI3`；
- 二维 `S1～S4` 转换为 Exodus sideset；
- 相对主输入路径递归展开 `*Include`，循环引用硬失败；
- Part 内 `MASS` 与 `TYPE=ANISOTROPIC` 三方向质量解析；
- 生成 `NodalTranslationalInertia` 可读取的 `x/y/z` 三份节点质量 CSV；
- 报告模型级初始边界、`*Dload, GRAV`、`*Dsload, HP` 参数、`*Boundary` 选项、分析步参数和 `*Amplitude, time=...`；
- 不支持的结构单元类型、错误连接数或缺失质量数据硬失败。

## 环境

转换与梁时程求解工具使用独立环境，不修改 P0 锁定的 `.build/env`。环境包含 Python、NumPy、netCDF4 和 SciPy：

```bash
./tools/abaqus/create-env.sh
```

## 转换

```bash
.build/abaqus-converter-env/bin/python \
  tools/abaqus/abaqus2exodus.py \
  --inp /path/to/model.inp \
  --out /path/to/model.e \
  --report /path/to/model-conversion.json \
  --mass-csv-prefix /path/to/model-added-mass
```

默认 `merge_tol=1e-9`，单位与模型一致。需要跨 instance 合并时必须依据模型尺度显式指定，不能沿用原工具的 `0.5` 默认值。

## 附加质量语义

`*_x.csv`、`*_y.csv`、`*_z.csv` 每行均为：

```text
x,y,z,mass
```

三份文件分别连接到 `disp_x/disp_y/disp_z` 的 `NodalTranslationalInertia`。附加质量进入惯性项，影响模态和动力响应；它不产生静水压力，不能替代 Abaqus `*Dsload, HP`。静水压力与附加质量必须分别建模。

当前 MOOSE `NodalTranslationalInertia` 可以在瞬态动力计算中读取这些 CSV。将方向性点质量纳入广义特征值质量矩阵仍需专门的模态实现和回归，不应直接宣称与 Abaqus Lanczos 结果等价。

## 验证

```bash
./tools/abaqus/run-tests.sh
```

测试覆盖递归 include、混合 `QUAD4/TRI3`、二维 sideset、重复 `elset` 的各向异性质量、Exodus 写出和分方向 CSV。

## 二维坝体静力—动力输入生成

`gen_dam_two_step_case.py` 将增强后的转换报告生成为两个受控 MOOSE 输入：先求解重力与静水压力静力平衡，再从静力 Exodus 读取位移和面外应变作为动力初值。动力步保留静力荷载，同时施加节点附加质量、Rayleigh 阻尼与基底加速度。

```bash
.build/abaqus-converter-env/bin/python \
  tools/abaqus/gen_dam_two_step_case.py \
  --report /path/to/model-conversion.json \
  --mesh /path/to/model.e \
  --mass-x /path/to/model-added-mass_x.csv \
  --mass-y /path/to/model-added-mass_y.csv \
  --output-dir /path/to/generated

cd /path/to/generated
/path/to/blackbear-opt -i dam_2d_full_static.i
mpiexec -n 4 /path/to/blackbear-opt -i dam_2d_full_dynamic.i
```

求解命令必须在生成目录中执行，因为 `file_base=results/...` 相对当前工作目录解析；输入内的网格、CSV 和初值文件路径则相对 `.i` 文件解析。

当前生成器只支持已校验的二维、单材料、X 向基底加速度场景，采用小变形线弹性平面应力、Newmark 平均加速度法和 Abaqus 材料中的 Rayleigh `alpha/beta`。它不转换 CDP 本构，不执行 Abaqus 的 Lanczos 频率步，也不能替代与 Abaqus ODB 的正式数值对标。

## 线性梁三场时程

`beam_field_solver.py` 源自锁定的 `demo-process` 线性梁直接积分器，支持把同一次 Newmark 时程计算拆分为位移、梁截面应力包络和绝对加速度三份 Exodus：

```bash
.build/abaqus-converter-env/bin/python \
  tools/abaqus/beam_field_solver.py \
  --report /path/to/report.json \
  --mesh /path/to/mesh.e \
  --displacement-out /path/to/displacement.e \
  --stress-out /path/to/stress.e \
  --acceleration-out /path/to/acceleration.e \
  --csv /path/to/summary.csv
```

输出语义：

- 位移：`disp_x/y/z` 与 `rot_x/y/z`；
- 绝对加速度：`accel_x/y/z` 与 `rot_accel_x/y/z`；
- 梁应力：`axial_stress`、`bending_stress`、`torsional_shear` 和 `vonmises_stress`，单位随 mm-t-N-s 模型为 MPa。

`vonmises_stress` 是由两端内力恢复的截面最外缘保守包络，不是实体积分点应力。是否启用 `--no-releases` 必须依据目标基线决定，不能静默改变连接语义。

Newmark 积分参数可用 `--newmark-gamma`/`--newmark-beta` 调整，默认 `0.5`/`0.25`（常加速度法，无数值阻尼），与既有基线复现结果逐数组一致；`γ>0.5` 引入算法阻尼（等价 HHT-α，`α=0.5-γ`）但精度降为一阶。`--allfields-out` 额外输出单份合并 Exodus，含全部节点场（位移/转角、绝对加速度、基底约束反力 `rf_x/y/z`）与全部单元场（应力包络、主应力 `s_max/mid/min_principal`、外缘应变 `axial/bending_strain`）。

三场视频使用：

```bash
/home/kevin/miniforge3/envs/moose/bin/pvpython \
  tools/abaqus/render_beam_fields.py \
  displacement.e stress.e acceleration.e /path/to/render-output <field>
```

`<field>` 分别为 `displacement`、`stress` 或 `acceleration`。建议每个场使用独立 pvpython 进程，避免 ParaView 保留前一场的色标。

## ParaView 动画

服务器安装 ParaView 的环境中，可参考 `demo-process` 的渲染链生成变形与 von Mises 应力动画：

```bash
/home/kevin/miniforge3/envs/moose/bin/pvpython \
  tools/abaqus/render_dam_dynamic.py \
  /path/to/result.e \
  /path/to/render-output \
  --max-frames 200 \
  --seconds-per-state 0.04 \
  --video-name dam_2d_full_static_dynamic.mp4
```

脚本只使用 Exodus 中真实存在的时间状态；时间状态多于 `--max-frames` 时在全时程内等间隔抽样。脚本按模型尺度自动计算位移显示放大倍数，输出 1600×900、H.264/yuv420p MP4。显示放大只用于可视化，不改变 `.e` 中的计算值。

## 已知边界

- 转换器输出网格、集合、二维 sideset、附加质量清单及受控分析语义报告；MOOSE `.i` 仍由场景生成器生成，不能把 JSON 报告本身当作求解输入；
- Abaqus CDP 参数不能机械映射为 BlackBear 损伤塑性参数；
- 三维 surface 到 Exodus side 编号仍沿用旧工具的 nodeset 行为，未在本次扩展中声明支持；
- 转换成功不等于物理等价，必须继续做拓扑、质量合计、边界、反力和关键响应对比。
- 2D 坝体已完成线弹性静力预载和完整 50 秒动力复算，但状态仍为 `prototype`；证据与限制见 [`../../doc/prototypes/abaqus-2d-static-dynamic-rerun.md`](../../doc/prototypes/abaqus-2d-static-dynamic-rerun.md)。
- 线性梁直接积分器用于复算已审核的 B31 时程模型，不替代 BlackBear 非线性材料或正式 Abaqus 对照；截面应力超过材料线性范围时必须转入非线性复核。
