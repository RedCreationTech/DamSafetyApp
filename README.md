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

## 在其他电脑部署运行

在 Windows（WSL2）或其他已有 MOOSE 开发环境的机器上，从源码编译、验证并运行
DamSafetyApp，请参阅 [`部署运行.md`](部署运行.md)。该流程不依赖当前
`192.168.0.138` 计算机。

## Abaqus 二维转换原型

`tools/abaqus/abaqus2exodus.py` 是从锁定的 `demo-process` 转换器派生的项目内原型，现已支持同一 Part 中混合 `CPS4/CPS4R/CPS3`、二维 sideset，以及 `*Include` 中各向异性节点附加质量的分方向 CSV 导出。

```bash
./tools/abaqus/create-env.sh
./tools/abaqus/run-tests.sh
```

使用方法、来源和已知边界见 [`tools/abaqus/README.md`](tools/abaqus/README.md)。该能力当前状态为 `prototype`，是 TASK-MOOSE-017 的前置技术验证，不表示完整 Abaqus 材料、分析步或动力响应已经等价迁移。

`gen_dam_two_step_case.py` 可从增强后的转换报告生成线弹性二维坝体“静力预载 → 完整动力”两阶段 MOOSE 输入，保留重力、静水压力、节点附加质量、Rayleigh 阻尼和 Abaqus `TOTAL TIME` 加速度幅值。2D 坝 50 秒复算证据见 [`doc/prototypes/abaqus-2d-static-dynamic-rerun.md`](doc/prototypes/abaqus-2d-static-dynamic-rerun.md)。该能力仍为 `prototype`，尚未迁移混凝土损伤塑性和 Lanczos 模态步。

同目录的 `beam_field_solver.py` 可复算已审核的线性 B31 梁时程，并分别输出位移、梁截面应力包络和绝对加速度 Exodus；`render_beam_fields.py` 生成对应三场 MP4。PR-RG-400gal-X 实例证据见 [`doc/prototypes/PR-RG-400gal-X-three-field-rerun.md`](doc/prototypes/PR-RG-400gal-X-three-field-rerun.md)。

## Abaqus CDP-compatible 本构开发

`tools/abaqus/abaqus_cdp.py` 已提供第一阶段 `.inp` 材料提取、四表校验、等效塑性应变
检查和可追溯 manifest 输出。它只使用 `.inp` 材料数据；ODB/CSV 保留为验证证据。
当前状态为 `prototype`。`AbaqusCDPStressUpdate` 已接入标准 MOOSE
`ComputeMultipleInelasticStress` 链，并通过七类单 HEX8 自洽诊断；生成的 `.i` 片段可
作为材料块接入该自有应用。由于 Abaqus 黄金材料点、自动稳定化、C3D8R 减缩积分/沙漏
控制及正式 C30 数据尚未完成验收，不能据此声称与 Abaqus 商业实现一比一等价。

使用方法和范围见 [`tools/abaqus/README.md`](tools/abaqus/README.md)。系统需求、设计、
任务和验收矩阵由 `damASR` 集成仓库维护。

### 计算节点隔离开发

路线 B 在计算节点 `/home/kevin/DamSafetyApp-cdp-dev` 的
`codex/abaqus-cdp-compatible` 分支开发和验证，正式部署目录
`/home/kevin/DamSafetyApp` 不直接修改。构建脚本会校验锁定的 BlackBear/MOOSE SHA，
并复用计算节点现有隔离编译环境：

```bash
./scripts/cdp/build-app.sh
./scripts/cdp/run-tests.sh
./scripts/cdp/run-unit-tests.sh
```

构建产物为 `DamSafetyApp-opt`。开发目录中的二进制不得直接替换正式求解器；只有已提交、
已推送并通过材料点和单单元门禁的 SHA 才能进入 C06 部署流程。

`CDPMaterialTable` 是路线 B 的 C++ 输入层：读取 Python 转换器生成的四张 CSV，
检查列名、有限数值、严格递增横坐标、正应力、损伤范围和等效塑性应变单调性，并提供
常值端点外推、分段线性插值以及节点左右导数。本构公式、局部积分、损伤/恢复、黏性、
事务式子步和算法切线由 `CDPConstitutiveModel`/`AbaqusCDPStressUpdate` 分层实现。
