# TASK-MOOSE-007-A 官方 ASR 受约束圆柱回归复现记录

> 状态：`verified`
>
> 完成日期：2026-08-17
>
> 对应 damASR 周报（2026-08-10—08-16）下周清单第 1 项：在锁定 BlackBear 上复现 `concrete_ASR_swelling/asr_confined.i` 的 `test_strip`

## 1. 复现范围与基线

锁定组合（脚本启动时逐层校验，不一致即退出）：

- DamSafetyApp：`3abd6bf3c03067455aa58161cac24e423d5c7aa4`（工作树含 2D 坝完整动力未提交改动，本次未触碰）；
- BlackBear：`1c190fd3d2b5f06a3518923f550a0e0a90b015d4`；
- MOOSE：`4bce02d91b56c7ed845a5747df4d24f415592504`；
- NEML：`a01a27b524a737b6746a840150f5acc2bace778e`。

复现对象：`test:concrete_ASR_swelling.ASR_swelling/test_strip`

- 输入：`.upstream/blackbear/test/tests/concrete_ASR_swelling/asr_confined.i`；
- 网格：`mesh_contact_strip.e`（官方仓库自带）；
- 判定：`Exodiff` 对照 `gold/asr_confined_strip_out.e`，自定义容差 `asr_confined.cmp`；
- 覆盖物理：轴对称混凝土圆柱、钢套约束与接触、ASR 反应度/体积应变/各向异性本征应变、应力与位移、ASR 微裂损伤。

## 2. 执行命令与结果

执行入口（可重复）：

```bash
cd /home/kevin/DamSafetyApp
./scripts/p1/run-asr-confined-strip.sh
```

脚本沿用 `scripts/p0/common.sh` 的隔离环境（conda 前缀 `.build/env`），从 BlackBear 根目录以 `MOOSE_DIR=.upstream/blackbear/moose` 运行 TestHarness，`-j 2`。

结果（2026-08-17 09:31，Asia/Shanghai）：

```text
test:concrete_ASR_swelling.ASR_swelling/test_strip ... OK
Ran 1 tests in 1.7 seconds.
1 passed, 0 skipped, 0 failed
```

Exodiff 在 `asr_confined.cmp` 自定义容差内通过，即本次求解结果与官方金标准数值一致。Exodus 为二进制格式且含运行时元数据，逐次运行文件哈希不同，等价性以 Exodiff 判定为准，不以文件字节相等为准。

## 3. 证据与哈希

| 对象 | 位置 | SHA-256 |
|---|---|---|
| 完整测试日志（含三层 SHA 与 `--version`） | `.build/logs/TASK-MOOSE-007-A-asr-confined-strip.log` | `369c0d30e33cbe3f4cc13f808a366cdd8cccf5dabe5708f6070c0f2efc8cdb53` |
| 本次求解输出 | `.upstream/blackbear/test/tests/concrete_ASR_swelling/asr_confined_strip_out.e` | `0da60f9581a78be0b4e7390b4bed54393da51f9b4ddc6cc4949122c9fb4c8c98` |
| 官方金标准（未改动） | `.upstream/blackbear/test/tests/concrete_ASR_swelling/gold/asr_confined_strip_out.e` | `88f58d554add7e8ee855fd314c0e427a00a74fc15d6ac1be5454a708b7d5d4f1` |

运行结果保留在被项目 Git 排除的上游测试工作区，未覆盖官方 `gold/`；BlackBear、MOOSE、NEML 均无已跟踪源码改动。

## 4. 已知事项与限制

- `MooseApp::addCapability()` deprecation warning 与 P0 门禁记录一致，属上游已知告警，不影响结果；
- 本次只复现 `test_strip`；`test_strip_kelvin`、`test_strip_isotropic`、`test_strip_in_tension` 变体和 `HEAVY` 的 `test_full` 已在 P0 定向门禁中通过（8 passed、1 skipped），本任务不重复运行；
- 本记录证明锁定基线上官方 ASR 回归可复现，是派生 `dam_thm_asr_baseline.i`（TASK-MOOSE-009/010）的物理基线之一，不构成工程模型验收。

## 5. 补充：官方 HEAVY 回归 test_full（2026-08-17）

应人工交互验证需求，补跑官方 `HEAVY` 用例 `test:concrete_ASR_swelling.ASR_swelling/test_full`（`asr_confined.i` + `mesh_contact.e` 完整接触网格，覆盖圆柱更大轴向区域；1025 节点、920 单元、2 个 element block、6 个时间步）：

```bash
cd /home/kevin/DamSafetyApp/.upstream/blackbear
MOOSE_DIR=moose ./run_tests -j 2 --heavy --re='concrete_ASR_swelling\.ASR_swelling/test_full$'
```

结果：**1 passed, 0 skipped, 0 failed**（5.8 s），Exodiff 通过官方金标准。

| 对象 | SHA-256 |
|---|---|
| 测试日志 `.build/logs/TASK-MOOSE-007-A-asr-confined-full.log` | `5821f80ca2864abadc1473bea495aca3faf58ead6d8173f11c717f121e7a082e` |
| 求解输出 `asr_confined_out.e` | `9d094cba29474d3a770dd5276dca85a020d3ed59ad8df52b1b0fa02ab3c0e78b` |
| 官方金标准 `gold/asr_confined_out.e` | `d7fa8a10ef2a0542554c3c38915d61582980f928687d449d096201038d661fc8` |

物理说明：均匀温度驱动下该算例空间场仍接近均匀（末步 `stress_yy` 单元间相对差异约 0.026%），其验证价值在于更大网格/接触语义的 Exodiff 一致性与时间演化，不适合作为空间梯度云图素材。P0 门禁中 `HEAVY` 延期项由此减少 1 项，剩余 3 项（钢蠕变损伤系列）仍按计划留待 TASK-MOOSE-025。
