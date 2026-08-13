# TASK-MOOSE-001 上游组合与环境基线

> 状态：`verified`
>
> 核查日期：2026-08-13

## 1. 候选源码组合

| 组件 | Remote | 固定提交 | 来源 |
|---|---|---|---|
| DamSafetyApp | `git@github.com:RedCreationTech/DamSafetyApp.git` | 首次提交后补充 | 项目远端 |
| BlackBear | `https://github.com/idaholab/blackbear.git` | `1c190fd3d2b5f06a3518923f550a0e0a90b015d4` | 官方 `devel` 快照，仅作候选 |
| MOOSE | `https://github.com/idaholab/moose.git` | `4bce02d91b56c7ed845a5747df4d24f415592504` | BlackBear `moose` gitlink |
| NEML | `https://github.com/Argonne-National-Laboratory/neml.git` | `a01a27b524a737b6746a840150f5acc2bace778e` | BlackBear `contrib/neml` gitlink |
| demo-process | `git@github.com:zhaoyul/demo-process.git` | `0cc241b8dbee6b7eb9594e9f7538fbd493d1f7a8` | 只读迁移来源 |

候选组合只有在 TASK-MOOSE-002 构建成功且 TASK-MOOSE-004 官方门禁通过后，才能由 `candidate` 改为 `approved`。

## 2. 许可证基线

- BlackBear：GNU Lesser General Public License 2.1 or later（`LGPL-2.1-or-later`）。
- MOOSE、NEML 及其二进制依赖：发布前按锁定提交和实际制品生成完整许可证/SBOM；本任务不替代 TASK-MOOSE-026 的供应链审查。
- DamSafetyApp：项目许可证待项目方确认；首次发布前必须完成。

## 3. 138 节点环境快照

- 主机：Linux `7.0.10-arch1-Watanare-T2-2-t2`，x86_64。
- 资源：8 个逻辑 CPU，约 31 GiB 内存，`/home` 约 223 GiB 可用。
- 系统工具：Git 2.54.0、GNU Make 4.4.1、GCC 16.1.1、Clang 22.1.5、CMake 4.3.3、Python 3.14.5。
- 用户级环境：`/home/kevin/miniforge3/envs/moose`。
- 候选构建工具链：Conda GCC 14.3.0、MPICH 5.0.0、`moose-petsc 3.24.6`、`moose-libmesh 2026.04.13_0185b8b`、`moose-wasp 2025.09.19_02960f1`。
- 非交互 Shell 默认 `PATH` 不含 Conda 环境，系统路径中找不到 `mpiexec`；Conda `moose` 环境内存在 `mpiexec`/`mpirun`。

上述 Conda 包能否与 2026-08-12 的候选 MOOSE 提交兼容，必须由 TASK-MOOSE-002 实际构建验证，不能仅凭包名判定。

## 4. 构建策略

- 所有上游源码置于 `/home/kevin/DamSafetyApp/.upstream/`，所有构建产物保留在隔离工作区并从项目 Git 排除。
- 不使用 `/home/kevin/gt/demo/mayor/rig/build/moose`，不修改其未跟踪 CDP 源码或既有产物。
- 首轮使用 `METHOD=opt`、最多 4 个并行编译任务；扩大并行度前先观察内存和负载。
- 通过 `/home/kevin/miniforge3/bin/conda run -n moose` 显式进入依赖环境，不依赖交互 Shell 初始化。
- 构建和测试必须保存三层 SHA、环境包清单、命令、退出码、日志与制品哈希。

## 5. 完成结论

- 三层候选 SHA 均已解析为完整提交，BlackBear 和 DamSafetyApp 远端可访问。
- 候选状态和构建环境已记录；没有把移动分支名作为运行版本。
- TASK-MOOSE-001 满足完成定义，状态记为 `verified`；版本组合尚未经过 G1 批准。
