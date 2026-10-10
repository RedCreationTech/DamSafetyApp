# CDP总应变率阻尼候选：三维公式诊断准备

2026-10-10 · TASK-DAM2B-CDP-RATE-048 · 本批源码已获用户提交推送授权，编译、节点验证及新Job尚待执行。

## 为什么先用小模型

目标仍是2D坝体。这里的一块HEX8只检查“混凝土开裂后，阻尼是否按未损伤刚度乘总应变率计算”，不与2D坝的11条曲线混用，也没有假定它是专家Abaqus算例。二维厚度处理尚未确定，因此明确拒绝2D输入，不能把它拿去直接跑H34或延长15/50秒。

系统审计及设计在damASR的 `docs/verification/evidence/2026-10-10-cdp-velocity-damping-semantics/` 和 `docs/designs/2026-10-10-CDP总应变率阻尼候选设计.md`。

## 代码范围与默认值

新增CDPVelocityDynamicPhysics/CDPVelocityRateMaterial/CDPVelocityStressDivergence及语法注册；未改现有CDP本构、H30弹性算法、质量、沙漏或正式Solver。

`enable_cdp_velocity_damping`默认false：委托原Dynamic Physics，不创建新阻尼材料或新变量。启用必须显式确认 `diagnostic_3d_model=true`，限3D Cartesian、增量SMALL、常量C0和常量非负阻尼系数。拒绝有限变形、Bbar、厚度变量、AD/Lagrangian、温度/本征应变、静态初始化和非零初始应变率。纯弹性控制须显式 `require_cdp_material=false`。

新阻尼不采用CDP损伤系数或Jacobian_mult；后者只用于原材料应力的导数。诊断输出 `cdp_total_velocity_rate` 是当前原始率，`cdp_velocity_damping_stress` 是HHT时间层加权后的阻尼张量，与材料S输出分开。初始时刻阻尼置零，允许原有非零初始加速度；本批实际初态为零。

## 输入与逐项验算

共同0～0.12s、5ms、1 HEX8（全积分）、E=3.04e10Pa、ν=.2、ρ=2500、HHT α=−.05/β=.275625/γ=.55、质量阻尼.7、刚度阻尼.03。CDP沿用已有单元测试表、黏性.0005及原参数；这是诊断模型参数，不宣称等同整坝材料标定。X加载、Y/Z有自由未知量，左/底/后分别消除刚体移动。24步只是预期，必须终态核实际步序。

| 输入 | 唯一用途 | 后处理门禁 |
|---|---|---|
| cdp-default-off.i / cdp-old-5ms.i | 新语法默认关闭与独立原Dynamic控制 | 规范化输入一致，实际U/V/A、S/E、损伤和步序一致 |
| cdp-new-zero.i / cdp-old-zero.i | 阻尼系数为零 | 新旧残差无阻尼差，实际全部物理输出一致 |
| elastic-new-5ms.i / elastic-old-5ms.i | 纯弹性控制，小幅加载 | 独立按接受U/V/A重建Newmark率和C0应力，核全位移交叉列 |
| cdp-new-5ms.i / cdp-old-5ms.i | 拉伸、卸载、压缩、卸载循环 | 必须实际出现塑性与拉压损伤；未激活则覆盖不足，不能算通过 |
| cdp-retry.i / cdp-retry-control.i | 收敛试步主动拒绝与相同接受序列控制 | 实际接受时刻/步长逐个配对，核全U/V/A、S/E、DamageT/C、κ，不能只看最终一点 |

所有8个节点U/V/A、平均S/E、CDP的DamageT/C/κ、候选的率/阻尼单独保存；完整Exodus保留积分点值，平均值不能代替每IP恒等式检查。共享include避免输入互相漂移。切线有限差分必须在可微分支检查；屈服/恢复切换点另记录，不能放宽阈值掩盖。

`tests`有9个运行、5个错误拒绝、2个全Jacobian检查。RunApp成功仅证明可执行；配对一致性、独立率闭合、损伤激活和重试恢复还需终态物理复核。这里没有对应的专家动态HEX8参考，因此即使内部门禁通过，也不能写成Abaqus物理精度已通过。

## 执行边界

本地无MOOSE环境，仅完成独立公式和输入准备检查。2026-10-10用户已明确批准本批提交推送与继续验证；随后依项目锁定环境隔离构建、记录SHA与实际libMesh依赖。`--check-input`等开发预检不代替正式Job。实际小模型通过LIMS Facade→C06登记，逐终态独立报告/审计；不直接经SSH跑有限元或代操作CAE。没有改C06/LIMS或重启服务。

本候选不能处理二维厚度；三维公式门禁通过后还要明确CDP平面应力下总厚度率/阻尼消元和全部导数，才进入真实坝体短窗、10.3s、15s及最终回归。
