# TASK-MOOSE-003 资源探测与运行时冒烟记录

> 状态：`verified`
>
> 完成日期：2026-08-13

## 1. 资源探测

- 逻辑 CPU：8。
- 总内存：32,659,696 KiB；探测时可用约 29,864,140 KiB。
- `/home`：464 GiB，总可用约 217 GiB。
- 隔离环境编译器：Conda GCC/G++ 14.4.0。
- MPI：MPICH 5.0.1，路径位于 `.build/env/bin`。

系统默认非交互 `PATH` 中没有 `mpiexec`，项目脚本统一通过 `conda run -p .build/env` 激活依赖。无需 `sudo` 或系统级安装。

## 2. 冒烟结果

使用官方 `test/tests/concrete_ASR_swelling/asr_confined.i` 执行：

| 模式 | 命令口径 | 结果 |
|---|---|---|
| 串行 | `blackbear-opt --check-input -i asr_confined.i` | 退出码 0，`Syntax OK` |
| MPI | `mpiexec -n 2 blackbear-opt --check-input -i asr_confined.i` | 退出码 0，两个 rank 均完成输入检查 |

`smoke-runtime.sh` 默认使用 2 个 rank，允许通过 `MPI_PROCS` 调整，但会拒绝非正整数和超过节点逻辑 CPU 数的请求。脚本中不存在固定 `-n 8`。

日志哈希：

- 环境探测：`9a05a74d0d50abe170b26ac891717d3cd90763cd060a6dfb8ea33710ffb51d2d`；
- 串行/MPI 冒烟：`4b509bee45d467f36f373bf317db6411b4f7c815b811d00338f24773afa7e83a`。

## 3. 已知告警

启动时 BlackBear 调用已弃用的 `MooseApp::addCapability()`，候选 MOOSE 输出 deprecation warning。该告警不影响本次输入检查，但应在 TASK-MOOSE-005 的 `DamSafetyApp` 注册入口中使用新 capability API，并在上游升级评审中持续跟踪。

## 4. 结论

串行和受控核数 MPI 路径均可用，资源探测和超配保护已脚本化。TASK-MOOSE-003 满足完成定义，状态记为 `verified`。
