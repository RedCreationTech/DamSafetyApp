# Abaqus 转换工具

本目录的 `abaqus2exodus.py` 源自 `demo-process@0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8`，原脚本 SHA-256 为 `7d01cdc2560b5fc6d3eb089c8fece6a3d066df3f1d9192844c754d0bebda1c19`。DamSafetyApp 在保留原有 `B31/C3D8R/T3D2` 能力的基础上增加：

- 同一 Part 内逐单元保存类型，支持混合 `CPS4/CPS4R/CPS3`；
- `CPS4/CPS4R → QUAD4`、`CPS3 → TRI3`；
- 二维 `S1～S4` 转换为 Exodus sideset；
- 相对主输入路径递归展开 `*Include`，循环引用硬失败；
- Part 内 `MASS` 与 `TYPE=ANISOTROPIC` 三方向质量解析；
- 生成 `NodalTranslationalInertia` 可读取的 `x/y/z` 三份节点质量 CSV；
- 不支持的结构单元类型、错误连接数或缺失质量数据硬失败。

## 环境

转换器使用独立环境，不修改 P0 锁定的 `.build/env`：

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

## ParaView 动画

服务器安装 ParaView 的环境中，可参考 `demo-process` 的渲染链生成变形与 von Mises 应力动画：

```bash
/home/kevin/miniforge3/envs/moose/bin/pvpython \
  tools/abaqus/render_dam_dynamic.py \
  /path/to/result.e \
  /path/to/render-output
```

脚本只使用 Exodus 中真实存在的时间状态，按模型尺度自动计算位移显示放大倍数，输出 1600×900、H.264/yuv420p MP4。显示放大只用于可视化，不改变 `.e` 中的计算值。

## 已知边界

- 只转换网格、集合、二维 sideset 和附加质量清单；材料与分析步主要进入 JSON 报告，MOOSE `.i` 仍需受控模板生成；
- Abaqus CDP 参数不能机械映射为 BlackBear 损伤塑性参数；
- 三维 surface 到 Exodus side 编号仍沿用旧工具的 nodeset 行为，未在本次扩展中声明支持；
- 转换成功不等于物理等价，必须继续做拓扑、质量合计、边界、反力和关键响应对比。
- 当前动力冒烟结果只有三个真实时间状态；其 MP4 是产出链路示例，不是完整地震响应动画。
