# DamSafetyApp

水工混凝土结构安全评估与寿命预测领域求解器。项目基于固定提交的 BlackBear/MOOSE 构建，当前处于 P0 可复现构建与官方门禁阶段。

## P0 快速入口

```bash
./scripts/p0/fetch-upstream.sh
./scripts/p0/create-env.sh
./scripts/p0/build-blackbear.sh
./scripts/p0/probe-environment.sh
./scripts/p0/smoke-runtime.sh
./scripts/p0/run-official-gate.sh
```

默认约束：

- 上游源码：`.upstream/blackbear`；
- 隔离依赖环境：`.build/env`；
- 构建与测试日志：`.build/logs`；
- 编译方法：`METHOD=opt`；
- 默认编译并行度：4；
- 官方测试默认并行度：2；
- 不读取或修改 `/home/kevin/gt/demo/mayor/rig/build/moose`。

具体版本、命令和完成证据见 [`doc/p0/`](doc/p0/)。当前状态不代表求解器已通过工程验收。

## Abaqus 二维转换原型

`tools/abaqus/abaqus2exodus.py` 是从锁定的 `demo-process` 转换器派生的项目内原型，现已支持同一 Part 中混合 `CPS4/CPS4R/CPS3`、二维 sideset，以及 `*Include` 中各向异性节点附加质量的分方向 CSV 导出。

```bash
./tools/abaqus/create-env.sh
./tools/abaqus/run-tests.sh
```

使用方法、来源和已知边界见 [`tools/abaqus/README.md`](tools/abaqus/README.md)。该能力当前状态为 `prototype`，是 TASK-MOOSE-017 的前置技术验证，不表示完整 Abaqus 材料、分析步或动力响应已经等价迁移。

同目录的 `beam_field_solver.py` 可复算已审核的线性 B31 梁时程，并分别输出位移、梁截面应力包络和绝对加速度 Exodus；`render_beam_fields.py` 生成对应三场 MP4。PR-RG-400gal-X 实例证据见 [`doc/prototypes/PR-RG-400gal-X-three-field-rerun.md`](doc/prototypes/PR-RG-400gal-X-three-field-rerun.md)。
