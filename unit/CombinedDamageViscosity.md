# 组合损伤黏性候选的开发检查

`damage_viscosity_scheme` 默认为 `branch_legacy`；`combined_scalar` 在每个真实材料子步 h 上把骨架损伤按当前黏性有效应力权重 r_v 组合，再用原 μ 松弛独立 dv。塑性和原 DamageT/C 分支保持，dv 允许随恢复下降。`cdp_combined_damage` / `cdp_stiffness_factor` 始终对应实际应力；目标、旧公式重建和 r_v 单列观测。

新增 dv 为第17个状态，转移矩阵23×23；原0～15状态编号保持。材料属性通过 old/init/propagate/store 接入事务，失败重分区从原接受状态重新开始。旧检查点缺此属性，控制与试验均从0秒新启动，不切换模式续算旧检查点。

`CombinedDamageViscosity.*` 检查：实际 AD 种子含旧 dv/组合目标/一般恢复参数/r_v/应力乘积；固定路径塑性分支不变；μ0、h0、恢复下降；塑性拉压完整子步切线；部分成功后失败重分区及重试一致。AD 核对复用未改的骨架 Jacobian 作为输入种子，混合对角应力 r_v 由 AD 比值独立求导；不等于全坝组装 Jacobian 或重复谱处导数已修复。

`tools/dump_cdp_legacy_state.C` 同一源文件分别链接冻结44df4c4和候选库，比较原22输出/22输入块及拉压固定路径；stderr 保留原单轴非对称方向的诊断。重复主值和表拐点不保证中心差分唯一，不放宽原模型容差或新增谱优化。

节点锁定环境构建后，在源码根目录执行：

```sh
./unit/DamSafetyApp-unit-opt --gtest_filter='CombinedDamageViscosity.*:AbaqusCDPStateIntegrator.*:AbaqusCDPSubstepIntegrator.*:AbaqusCDPLocalIntegrator.*:AbaqusCDPEffectiveStrength.*:HourglassElasticScale.*'
```

以上均为开发检查，未执行 FE 时间步。不代替真实坝25s默认控制、候选单变量试验、七个冻结完整回归或50s精度验收。
