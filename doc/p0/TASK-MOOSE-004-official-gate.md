# TASK-MOOSE-004 BlackBear 官方门禁记录

> 状态：`verified`
>
> 完成日期：2026-08-13

## 1. 门禁范围

锁定组合：

- BlackBear：`1c190fd3d2b5f06a3518923f550a0e0a90b015d4`；
- MOOSE：`4bce02d91b56c7ed845a5747df4d24f415592504`；
- NEML：`a01a27b524a737b6746a840150f5acc2bace778e`；
- `blackbear-opt` SHA-256：`551d197b26c6161836c7515bfac6c6af3211d51b843195f01c4da16c6d399a70`。

门禁分为：

1. `asr_confined.i --check-input`；
2. `concrete_ASR_swelling` 与 `concrete_ASR_validation` 定向回归；
3. BlackBear 全量默认非 `HEAVY` 回归。

TestHarness 必须从 BlackBear 根目录启动，并显式设置 `MOOSE_DIR=.upstream/blackbear/moose`。这是项目运行脚本的环境适配，不是上游源码补丁。

## 2. 测试结果

| 集合 | 结果 | 用时 | 处置 |
|---|---|---:|---|
| ASR 定向回归 | 8 passed、1 skipped、0 failed | 6.5 s | 跳过项为官方 `HEAVY` 的完整约束模型 |
| BlackBear 全量默认回归 | 188 passed、4 skipped、0 failed | 143.8 s | 4 项均因官方 `HEAVY` 标签按默认策略跳过 |

全量默认回归覆盖 ASR、温湿传输、微裂损伤、损伤塑性、蠕变/收缩、钢筋粘结滑移、NEML 及 assessment 输入语法检查。它证明锁定上游组合在目标节点的默认回归基线可用，但不等同于 DamSafetyApp 工程模型 V&V。

## 3. 跳过项处置

以下测试未失败，而是被官方 TestHarness 标记为 `HEAVY` 并在默认门禁中跳过：

- `concrete_ASR_swelling.ASR_swelling/test_full`；
- `steel_creep_damage_oh.creepdamage/ad_steel_creep_damage`；
- `steel_creep_damage_oh.creepdamage/steel_creep_damage`；
- `steel_creep_damage_oh.vector/multicreeplaw`。

处置：P0 不将 `HEAVY` 测试设为阻断项；在 TASK-MOOSE-025 的单节点性能/长期任务验证中，按资源预算运行 ASR `test_full` 和与项目材料路线相关的重测试。正式发布门禁是否包含全部 `HEAVY` 测试由 TASK-MOOSE-020/026 决定。

## 4. 结果与日志哈希

日志：

- ASR 定向门禁：`c884721ee8f19e793919d370879bf37acde12078e35dd263e129952181c99b65`；
- 全量默认回归：`df84a0cfb5bd46e017bffc351d61fadedd1474b771137ae01bd7858cbfb41741`。

代表性 ASR 结果：

- `asr_confined_strip_out.e`：`e5d2206d22b1318fcc81e46f1dd3a4d5ad7f39f03150b65b9c13df9fa4ba8708`；
- `asr_validation_case1_out.csv`：`5447e7b6c358eb80025b004507298700b9ac86576b479548589bc6d3e254829d`；
- `asr_validation_case2_out.csv`：`2bd07e4d68e7b1fdeea609eeb0df0a38e140495c0bd03e6124ff6ebd1ed953b5`；
- `asr_validation_case3_out.csv`：`78841ed188bb2a471bd2c685f61371a180da9226c085867d9620fe22078e98e8`。

运行结果保留在被项目 Git 排除的上游测试工作区，未覆盖官方 `gold/`。BlackBear、MOOSE 和 NEML 均无已跟踪源码改动。

## 5. 已知限制

- 启动时存在 BlackBear `MooseApp::addCapability()` deprecation warning，不影响本次结果；
- `--version` 只显示 MOOSE 快照，尚不满足三层 SHA 输出要求，由 TASK-MOOSE-005 修复；
- assessment 中部分用例只执行 `SYNTAX PASS`，不能据此宣称完成工程验证；
- P0 未运行官方 `HEAVY` 测试，也未运行 DamSafetyApp 项目级 TestHarness。

## 6. G1 结论

候选组合已通过目标节点构建、ASR 定向回归和 BlackBear 全量默认回归。批准该组合进入 P1 应用开发，版本状态从 `candidate` 更新为 `approved-for-development`。

该批准只表示上游开发基线可用，不表示生产发布、工程安全结论或寿命预测模型已经 `accepted`。TASK-MOOSE-004 满足完成定义，状态记为 `verified`。
