# Abaqus 2D 坝体静力—动力复算记录

- 文档状态：`prototype`
- 记录日期：2026-08-14
- 原始算例：`/home/kevin/Abaqus/2d-dam/dam-2b.inp`
- 实施目录：`.build/cases/abaqus-2d-dam-p0/full-static-dynamic-v1/`
- 求解器：DamSafetyApp 锁定的 BlackBear/MOOSE

## 1. 上次未执行第二步的原因

原 Abaqus 输入包含三个分析步：

| 顺序 | Abaqus 步 | 作用 | 上次状态 |
|---|---|---|---|
| 1 | `Fre` | Lanczos 频率提取，50 阶 | 未迁移 |
| 2 | `Gra` | 重力与上游面静水压力静力平衡 | 已计算 |
| 3 | `EQ` | 50 秒基底加速度动力响应 | 只做了 0.02 秒、常量小加速度的链路冒烟 |

上次产物的目标是验证 P0 阶段的混合 `CPS4R/CPS3` 网格、二维边界、附加质量和 Exodus/MP4 产出链路，并非完整物理复算。旧转换报告虽记录了基本分析步与幅值数据，但缺少自动、无歧义生成两阶段输入所需的关键语义：

- 模型级初始 `*Boundary`；
- `*Dload, GRAV` 的方向和重力值；
- `*Dsload, HP` 的水头参数；
- `*Boundary` 的 `type=ACCELERATION` 与 `op=NEW`；
- `*Step` 参数以及 `*Amplitude, time=TOTAL TIME`。

因此旧流程没有把 Abaqus 的 `Gra → EQ` 状态传递实现为可复现的 MOOSE 两阶段计算。0.02 秒结果只能证明动力求解和附加质量对象能够运行，不能证明第二分析步已经完成。

## 2. 本次修复

转换器现将上述语义写入 JSON 报告；新增 `gen_dam_two_step_case.py` 生成：

1. `dam_2d_full_static.i`：施加重力与静水压力，得到静力平衡状态；
2. `dam_2d_full_dynamic.i`：从静力 Exodus 精确读取 `disp_x`、`disp_y` 和 `strain_zz` 初值，保留静力荷载，再施加完整地震加速度；
3. `dam_2d_eq_acceleration.csv`：根据 Abaqus `time=TOTAL TIME` 将原幅值时间平移到 MOOSE 动力步的全局时间轴；
4. X/Y 向 `NodalTranslationalInertia`：读取 `add.dat` 导出的方向性节点附加质量；
5. Rayleigh 阻尼：使用 Abaqus 材料中的 `alpha=1.95`、`beta=0.00113`；
6. Newmark 平均加速度积分：`beta=0.25`、`gamma=0.5`，动力步长 `0.01 s`。

模型采用小变形线弹性平面应力。此次修复没有声称 Abaqus CDP 或 Lanczos 模态步已经迁移。

## 3. 可复现基线

| 项目 | 版本或 SHA-256 |
|---|---|
| DamSafetyApp | `3abd6bf3c03067455aa58161cac24e423d5c7aa4` 加本记录所列未提交改动 |
| BlackBear | `1c190fd3d2b5f06a3518923f550a0e0a90b015d4` |
| MOOSE | `4bce02d91b56c7ed845a5747df4d24f415592504` |
| `dam-2b.inp` | `76af2cad2d6fdd48111f8c421c275badac7203dbf1af19ceaf510a146dc4b1f2` |
| `add.dat` | `5582b22b5a7837b525fa3fb0785032962c4a30562ef2b46daa4eeb7cdc2cb2e0` |
| 转换报告 | `e110f69a97b621ed4d7163da28b06e2cd8b7f6d84f82f99f083ba0f14ee5f00f` |
| 静力 `.i` | `ff30af208ad9b95a9adb7eeb88b173ac4b0c3439f3ed13e05bfb14f5b80c08ed` |
| 动力 `.i` | `e5553fd2b230013c2df11f1028585a04baf0f5e47fe3c0d3c7e6da0d4f30ff8c` |
| 加速度 CSV | `91c34fe2d8c3f6830270a9d05e3a541383b328ee3979b03aced1959a0f39f2ad` |

## 4. 验证记录

| ID | 检查 | 结果 |
|---|---|---|
| `TEST-DAM-2STEP-001` | 转换器与生成器单元测试 | 6 项全部通过 |
| `TEST-DAM-2STEP-002` | 静力、动力输入 `--check-input` | 通过 |
| `TEST-DAM-2STEP-003` | 静力平衡求解 | 通过，输出可读且数值有限 |
| `TEST-DAM-2STEP-004` | 静力末状态与动力首状态逐数组比较 | `disp_x`、`disp_y` 完全一致 |
| `TEST-DAM-2STEP-005` | 5000 步、50 秒动力求解 | 通过，退出码 0，501 个 Exodus 状态，所有字段为有限值 |

完整动力计算使用 4 个 MPI 进程，墙钟时间约 2354 秒。输出包含节点位移、速度、加速度和单元应力场。

## 5. 结果摘要

| 指标 | 全时程绝对峰值 |
|---|---:|
| `disp_x` | `0.268466057216437 m` |
| `disp_y` | `0.005865165025626482 m` |
| `accel_x` | `24.0146461919731 m/s²` |
| `accel_y` | `11.17260569230275 m/s²` |
| `vonmises_stress` | `3.8763713720903 MPa` |

原始地震幅值在 `48.11 s` 结束时仍为 `-0.001267427`，输入也未包含基线校正。当前 MOOSE 幅值函数在数据范围外保持末值，因而较大的末段水平位移漂移必须在与 Abaqus ODB 对比时单独审计；本次没有擅自修改原始地震输入。

## 6. 交付物

| 文件 | SHA-256 | 说明 |
|---|---|---|
| `dam_2d_full_static_dynamic.e` | `9044d072db62476b8395b5473ac3d60ba126149a94de5cda135f55f89d5af8dc` | 501 个真实时间状态，位移、速度、加速度与应力字段 |
| `dam_2d_full_static_dynamic.mp4` | `38f9a343377b8b710ce88af89e46fcfa975c7f0e4957cb15b38d4daee60a4784` | 200 个真实状态等间隔抽样，1600×900、25 fps、8 秒 |

服务器产物位于：

```text
/home/kevin/DamSafetyApp/.build/cases/abaqus-2d-dam-p0/full-static-dynamic-v1/
```

本地交付目录：

```text
/Users/a123/Desktop/proj/北建/0810/模型示例/2D坝体模型/DamSafetyApp-results/
```

## 7. 已知限制

- 当前复算为线弹性平面应力模型，未迁移 Abaqus Concrete Damaged Plasticity；
- 未执行 Abaqus `Fre` 步的 Lanczos 频率提取，不能宣称模态等价；
- 尚未取得 Abaqus ODB 的同位置、同时间数值进行误差量化；
- 原输入的地震幅值末值非零且未做基线校正，位移漂移需要作为输入审计项；
- 结果状态为 `prototype`，不是工程验收结论。
