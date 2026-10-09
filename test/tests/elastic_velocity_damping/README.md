# 默认关闭的纯弹性速度阻尼候选

状态：源码WIP、小模型输入已准备，未在MOOSE环境编译、验算或提交LIMS Job。设计ID：DES-DAM2B-VELOCITY-DAMPING-039。

## 修改内容与边界

`ElasticVelocityDynamicPhysics` 的 `enable_velocity_rate_damping` 默认为false；关闭时使用原DynamicSolidMechanicsPhysics的实体核，没有新增速率变量／材料／方程。开启时需要显式声明严格恒定线弹性、SMALL incremental、WEAK_PLANE_STRESS、Newmark参数、弹性常数，以及与原strain_zz相同插值和块的非线性`thickness_rate`变量。

只替换实体刚度比例阻尼的变化率：由当前位移和上一接受V/A重建Newmark速度；保留原静态应力、质量比例阻尼、地面运动和输入阻尼系数。与B不同；本目录诊断模型没有沙漏核，不能据此判断H28沙漏／全坝效果。输出`elastic_velocity_damping_stress`与`elastic_velocity_rate_stress`独立于原S输出。

禁止CDP、有限应变、Bbar、温度／本征应变、静态初始化与非零初始应变率；初始非零加速度保留。已知CDP属性和实际弹性／切线张量会校验。该运行校验不替代输入审计，不能排除所有未知自定义材料；只通过Physics入口应用于已审计纯弹性组合。

## 小模型与准备的门禁

矩形2×1，计划节点0=(0,0)、1=(2,0)、2=(0,1)、3=(2,1)，底部X/Y固定，上部X/Y自由；实际Job必须核实节点编号。E=100、ν=.2、ρ=2.5，α=−.05、β=.275625、γ=.55、η=.7、ζ=.03。它是为了测算子／导数的诊断模型，既不是全坝，也不是新Abaqus试验。

- `upstream.i`与`default-off.i`：同一弱平面应力模型，原Physics与新Physics默认关逐步U/V/A/S比对；RunApp仅证明运行，不证明解相同。
- `velocity-10ms.i`／`velocity-5ms.i`：显式开启，0～.2s，完整4点积分用于隔离全Jacobian／弱约束问题。TIMESTEP_END记录每节点U/V/A/z/r；两独立速率应力诊断为单元量。
- `velocity-retry.i`：变步并故意拒绝.04s试算，对独立参考与未拒绝的相同接受序列比较，不只检查拒绝日志。
- `velocity-5ms.i`启用checkpoint，节点TestHarness的恢复检查及正式Job恢复证据仍待执行。
- `upstream-reduced-q1.i`／`velocity-reduced-q1.i`：保持原FIRST厚度空间，改成H28的一点积分；这是数值秩风险诊断，未加入任何惩罚、沙漏或厚度边界来掩盖奇异性，不作为预设成功的RunApp测试。
- TestHarness包含完整Jacobian及未确认模型、有限应变、厚度空间不一致、弹性常数不匹配的负对照。数学helper单元测试覆盖变步重试、初态、厚度两列导数及零步长拒绝。

## 一点积分厚度空间的风险

独立实际Q1矩阵：4点积分Kzz秩4，一点积分Kzz秩1（4个厚度自由度中3个不可观测）。原静态厚度也有此风险；新增同空间代数速率会再增加无约束方向。不能用额外刚度、换插值、清零历史或线性求解器成功退出宣称问题解决。

必须先核查节点的矩阵／解唯一性与可观測量、拒绝／重启、求解成本；若目标减缩积分空间不稳，应在同一弱空间的可观测子空间上做无惩罚消元／约束，再单独更新设计。当前候选没有跨过这一门禁，不能直接复制为H28整坝验证。

## 验证与归属

本地没有MOOSE工具链，不在本地编译。用户已明确授权本批源码提交推送；在本地所属worktree提交推送，再用节点锁定MOOSE/BlackBear/libMesh和独立release构建；不改共享MOOSE和生产默认版本。静态文本／独立Python矩阵通过不等于C++或MOOSE测试通过。

开发单元／TestHarness按仓库规则在节点做；正式有限元算例只从LIMS Facade提交并取得job_id，每Job保留独立报告。小模型报告使用ground-motion-quad4图和独立矩阵参考，不能出现C3D8图或把独立曲线标为Abaqus。七冻结回归与性能／重复性在候选推广前必须通过。
