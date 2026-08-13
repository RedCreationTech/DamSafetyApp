# TASK-MOOSE-002 隔离构建记录

> 状态：`verified`
>
> 完成日期：2026-08-13

## 1. 隔离边界

- 项目：`/home/kevin/DamSafetyApp`。
- 上游源码：`/home/kevin/DamSafetyApp/.upstream/blackbear`。
- 用户级依赖环境：`/home/kevin/DamSafetyApp/.build/env`。
- 构建日志：`/home/kevin/DamSafetyApp/.build/logs`。
- `.upstream/` 和 `.build/` 均由项目 `.gitignore` 排除。
- 未读取为构建依赖、未修改 `/home/kevin/gt/demo/mayor/rig/build/moose`。

## 2. 锁定依赖

| 包 | 版本/构建 |
|---|---|
| `moose-dev` | `2026.07.30`, `mpich` |
| `moose-build` | `2026.07.30`, `mpich` |
| `moose-petsc` | `3.25.2`, `mpich_0` |
| `moose-libmesh` | `2026.06.05_ab36c00`, `mpich_1` |
| `moose-libmesh-vtk` | `9.6.2`, `mpich_0` |
| `moose-wasp` | `2026.06.15_d8777a8`, `build_0` |
| `moose-tools` | `2026.06.16` |

环境由 MOOSE `4bce02d...` 的 `scripts/versioner.yaml` 和 `moose-dev` recipe 确定；没有复用版本较旧的 `/home/kevin/miniforge3/envs/moose`。

## 3. 可复现命令

```bash
cd /home/kevin/DamSafetyApp
./scripts/p0/fetch-upstream.sh
./scripts/p0/create-env.sh
JOBS=4 ./scripts/p0/build-blackbear.sh
```

脚本会校验 BlackBear、MOOSE 和 NEML 的完整 SHA，发现版本漂移时立即失败。

## 4. 构建结果

- `blackbear-opt`：`/home/kevin/DamSafetyApp/.upstream/blackbear/blackbear-opt`。
- 文件大小：57,640 bytes。
- SHA-256：`551d197b26c6161836c7515bfac6c6af3211d51b843195f01c4da16c6d399a70`。
- 重复增量构建前后可执行文件 SHA-256 一致。
- BlackBear、MOOSE、NEML 均无已跟踪源码改动。
- 空间占用快照：`.upstream` 约 3.8 GiB，隔离环境约 5.4 GiB。

日志哈希：

- `TASK-MOOSE-002-build-blackbear.log`：`2d5864761a9d0af54a24562b5412709cfe2609f36bd5f4f3fb00ddc0ce132bd0`；
- `TASK-MOOSE-002-environment-explicit.log`：`be91a45276fd69e1afe58a9e41492f7a9d6ac373e8f7f843bee8a6666e5525e0`；
- `TASK-MOOSE-002-repeat-build.log`：`ccb1b00d4401658b375bca2baf942c5a9edfb96efca48e243afc297fa45435b8`。

日志位于被 Git 排除的本机构建目录，以上哈希用于核对本次证据，不代表发布制品摘要。

## 5. 结论

锁定源码组合可在目标节点的独立用户级环境中构建 `blackbear-opt`，重复构建路径和版本检查已脚本化。TASK-MOOSE-002 满足完成定义，状态记为 `verified`。
