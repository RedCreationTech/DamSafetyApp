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
