# Abaqus 二维混合网格与附加质量转换原型

> 验证日期：2026-08-13
>
> 状态：`prototype`
>
> 关联任务：TASK-MOOSE-017 前置技术验证，任务本身仍为 `planned`

## 1. 目标与范围

本原型验证以下技术问题：

- REQ-CONV-2D-001：同一 Abaqus Part 内混合 `CPS4/CPS4R/CPS3` 时，保持逐单元拓扑并输出可由 MOOSE 读取的 Exodus；
- REQ-CONV-SURF-001：将二维 `S1～S4` element surface 映射为 Exodus sideset；
- REQ-CONV-MASS-001：解析 `*Include` 内逐节点 `MASS` 与 `*Mass, TYPE=ANISOTROPIC`，保留三个方向的附加质量；
- REQ-CONV-AUDIT-001：转换报告包含源文件、块、集合、质量合计和输出文件，未知或不完整输入硬失败。

转换器源自 `demo-process@0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8` 的 `tools/abaqus2exodus.py`，原脚本 SHA-256 为 `7d01cdc2560b5fc6d3eb089c8fece6a3d066df3f1d9192844c754d0bebda1c19`。原目录未被修改。

## 2. 实现结果

项目内工具位于 `tools/abaqus/`，主要变化如下：

1. 结构单元类型改为逐 element 保存，`CPS4/CPS4R → QUAD4`、`CPS3 → TRI3`，不同拓扑写入不同 Exodus block；
2. 递归展开相对路径 `*Include`，缺失文件和循环引用直接失败；
3. 二维 Abaqus 边号映射为 Exodus side，并生成 sideset；
4. `MASS` 与结构单元分离，重复的同名 `elset` 按其紧邻 `*Mass` 正确关联；
5. 各向异性质量导出为 `x/y/z` 三份 `x,y,z,mass` CSV，供 MOOSE `NodalTranslationalInertia` 使用；
6. 默认节点合并容差从原工具的 `0.5` 收紧为 `1e-9`，避免按本模型尺度错误合并节点。

转换器隔离环境为 Python 3.11.15、NumPy 2.3.5、netCDF4 1.7.4，位于被忽略的 `.build/abaqus-converter-env`，未改变 P0 上游构建环境。

## 3. 自动化测试

执行：

```bash
./tools/abaqus/run-tests.sh
```

测试 ID 与覆盖范围：

| 测试 | 结果 | 覆盖内容 |
|---|---|---|
| TEST-CONV-2D-001 | 通过 | 单 Part 混合 QUAD4/TRI3、递归 include、sideset、Exodus 写出 |
| TEST-CONV-MASS-001 | 通过 | 重复 `elset`、各向异性质量合计、三方向 CSV |
| TEST-CONV-NEG-001 | 通过 | include 循环、缺失质量属性硬失败 |

## 4. 用户坝体模型转换

输入副本位于 `.build/cases/abaqus-2d-dam-p0/source/`，运行产物位于同级 `converted-v2/`，均不进入 Git。

| 项目 | 转换结果 |
|---|---:|
| 节点 | 11,190，合并 0 |
| 结构单元 | 10,956 |
| QUAD4 | 10,955 |
| TRI3 | 1 |
| `_PickedSurf7` | 110 条边 |
| `POINT_MASS` / `UP_STREAM__DAM_1` | 均为同一组 113 节点 |
| 附加质量 x 向合计 | 2,129,652.257 |
| 附加质量 y 向合计 | 123,419.370 |
| 附加质量 z 向合计 | 0 |

MOOSE 对转换后 Exodus 的 `--mesh-only` 检查成功，识别两个子域、`QUAD4/TRI3`、6 个节点集和 1 个 sideset。将既有静力原型切换到转换后 Exodus 后，CSV 与直接读取 Abaqus 网格的结果逐字节一致，SHA-256 均为 `7fbdd017ec24585e15f50d5cd5aeb1c1f71fed023be49a2a2747de704ba97756`。

## 5. 附加质量动力冒烟

通过 `NodalTranslationalInertia` 将 x/y 两方向 CSV 接入 0.02 s 的线弹性动力冒烟计算；`--check-input` 和两个时间步均成功，结果均为有限值。关闭全部 `NodalKernels` 后运行同一对照算例：

| t (s) | 有附加质量 `max_abs_accel_x` | 无附加质量 `max_abs_accel_x` |
|---:|---:|---:|
| 0.01 | -4.35994888e-3 | -4.35994888e-3 |
| 0.02 | -6.11800839e-3 | -5.85375499e-3 |

第二时间步响应发生变化，证明附加质量已参与动力方程。此输入采用恒定小幅基底加速度，仅用于连接性冒烟；没有迁移完整 50 s 地震时程、阻尼、CDP、静力预加载或 Abaqus Lanczos 模态过程。

### 5.1 Exodus 与 MP4 交付副本

参考 `demo-process` 的 ParaView 渲染链，使用 `WarpByVector` 放大显示位移，并按单元 `vonmises_stress` 着色。视频只包含 Exodus 中 `0、0.01、0.02 s` 三个真实状态，每个状态保持 1 秒，没有插值或伪造计算帧。

交付目录：`.build/cases/abaqus-2d-dam-p0/deliverables/`。

| 文件 | 内容 | SHA-256 |
|---|---|---|
| `dam_dynamic_added_mass_prototype.e` | 11,190 节点、10,956 单元、3 个状态的 Exodus 结果 | `263f6e39331912e43a1dce1e586220b6d792d72ac1c95b6a5d4472bb94aac0bf` |
| `dam_dynamic_added_mass_prototype.mp4` | 1600×900、H.264/yuv420p、25 fps、3.0 s | `45bc31491fe931149af97579c890abf6a30beb4634d0964d6cc809d3951d3e41` |

渲染脚本 `tools/abaqus/render_dam_dynamic.py` 的 SHA-256 为 `b2b40047a1d2c19d6fea2bbb1c6173ec1aee51bff199063e1f5dabb7336c223d`。最终帧已人工检查，网格、云图、标题、时间状态和图例可见。

## 6. 关键证据哈希

| 证据 | SHA-256 |
|---|---|
| 转换器 | `c72d629dafcd56f14057dcb4e6752e9c0c84dc3b094b918965ad30a1fc767dfc` |
| `dam-2b.e` | `21a1e7a333b995e33e4bd0bae216b8ab8bbf2530799abc124e560e4d5bce462b` |
| `conversion.json` | `bea86adeb53118f37c6812ee6b57aa2eabf9f974fac85d8a63747de4a25b7a75` |
| `added-mass_x.csv` | `b50db5df3bcf01bda10ac09ad60b30d21c1c54f72f5bb6bd757b74bc369bf42a` |
| `added-mass_y.csv` | `99d140392fcef85315169bb4364b057b4acb5a9762b92c54d8ee888f8a235719` |
| 转换网格 MOOSE 读取日志 | `baa948a1eaea4fa87206e49d4bfed7dabee21c9f9b005f31dd930f502580fe87` |
| 转换网格静力日志 | `8619f75875924025023fa466483d483a9be16911a558bf0bcaac35f6cc85a727` |
| 附加质量动力日志 | `8c39e77f873788588cdab3653eb1b66b2812ef52e938a250644da2a7507702ad` |
| 无附加质量对照日志 | `8a21c4dabe1d7b74d07418445e2adc74adb3304796f480c4f07b63a8565f2437` |

## 7. 物理含义与剩余门禁

`add.dat` 表示上游节点的方向性水体附加质量。它模拟库水与坝体共同振动带来的水动力惯性，影响模态与地震响应；它本身不会产生静态水头压力，不能替代 Abaqus `*Dsload, HP`。工程模型需要分别表达静水压力和附加质量。

本原型只证明转换与瞬态动力接入可行。TASK-MOOSE-017 仍须补齐版本化转换 Manifest、更多负向/金标准测试、节点与边界差异报告、单位约束和 Abaqus 对照；TASK-MOOSE-012/013 仍须完成结构动力模板、材料、阻尼、模态及地震时程迁移。上述门禁未通过前，结果不得标记为 `verified` 或用于工程安全结论。
