# DS04G 2597 材料入口被动观测候选

基于 a9 源码，尚未在锁定 MOOSE 环境编译或验证，**不准直接用于生产 Job**。默认没有启用观测；只有显式配置 `CDPObservationProblem` 和 `observation_enabled=true` 才捕获。观测与原始积分器共享求解路径，但是否路径中性必须靠 OFF/ON 的测试实证，不能仅看代码断言。

从已哈希的 Q0 主输入生成专用候选（输出路径须不存在，输入附件要原样保留）：

```sh
python3 tools/diagnostics/prepare_ds04g_2597_observer.py /path/to/ds04g_sp1_q0_eq2p56.i /path/to/new-observer.i
```

工具拒绝 Q0 SHA 不符、已有 `Problem`、或时间步结构不符的输入；只插入一个 `Problem` 并将 `TimeSequenceStepper` 换为其被动事件子类 `CDPObservationTimeStepper`，不修改原 time_sequence。输入必须放在与附件兼容的作业快照中，不能只上传单个 `.i`。

按专家 2597 四节点坐标唯一映射，仅 QP0，总时间 3.018–3.202（EQ 2.018–2.202）秒，不强制新时间点。时间步器只转发原有 `step`/`acceptStep`/`rejectStep`，并记录事件；`accepted_end` 是调用基类 acceptStep 之前的事件标签，不单凭它证明已提交状态。逐 rank CSV 包含各试算的 evaluation、call、partition/N/substep 以及字段；`residual_solution_vector` 的向量身份也不足以证明该调用最终被接受。必须将收敛事件、accepted 状态、解向量和记录关联；不唯一或缺失时停止解释。

入口捕获 old elastic、raw increment、old/new total、old16、dt；成功捕获应力、新状态、分区统计与 6×6 列主序切线。局部积分器直接记录计算出的原生 trial stress、ft/fc、stress scale、trial_yield 与实际弹性判断阈值，区分子步及重试；这是 MOOSE 原生 trial，不是 Abaqus 材料驱动。`local_solve` 为各 rank 上真实局部积分调用的递增序号；应与 rank、evaluation/call/partition/substep 联合使用。分区成功或失败单独标记，记录顺序无法关联接受态时仍属证据不足。QP 为张量应变 xx/yy/zz/xy/yz/xz，xy 为张量剪切；应力 Pa、时间 s。

没有安装 SNES monitor，也不为残差求额外 MPI 全局范数；初始化的四节点唯一匹配使用一次按标签、所有 rank 共同执行的 `comm().sum`，要求配置在所有 rank 一致。记录 `SNES` 现有迭代号与求解终态，不表示有逐次 accepted iterate monitor。记录文件和各 rank 状态文件须全部存在，检查截断、I/O、history_incomplete 和唯一单元匹配。没有完整状态或关联不得做重放。

目标机流程：先由协调者审查并将本独立分支提交推送，锁定环境中构建隔离候选；运行 C++ 单元、既有 SP1 与受控 ON/OFF 检验，含 N1/N>1、失败重试、MPI、路径与输入完整性。比较所有物理字段、实际接受步、残差、应力、状态和切线，完成原始材料输入重放；再由协调者判定是否提交 C06 观测 Job。不得将本地静态测试当成编译通过或路径中性证明。
