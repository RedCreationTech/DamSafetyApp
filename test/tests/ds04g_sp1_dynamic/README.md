# DS04G 测试专用动态装配

三真实几何×四区间，12组规定历史，0接受步。当前U来自SP1原生快照；V/A缩放1e-3，预应力缩放0.05且S33投影到弹性平面应变。old/older、κ=0及F为制造历史，不是真实CDP状态；原生产材料与所有内核公式不改。

`DS04GSP1DynamicAuditProblem`仅允许serial/1thread/1QUAD4/8DOF，实际切换内核活动缓存逐项装配R/J与三个尺度FD。拒绝检查先装配试步，再调用同FEProblem.restoreSolutions API，校验接受old/older和上次接受区间未变、R/J重放相同；不冒称自动Executioner拒绝或整坝真实CDP异常恢复。

结果CSV含base/FD/R/J/活动组/材料历史及每H状态；分析按锁定独立Q4质量、当前重力、HG和实体Taylor假设核验。动态分项不匹配需如实报告，不改变生产残差救门禁。本夹具不代表整坝响应改善，不触发七冻结回归。

## 可选两接受区间历史

`e*_two_interval.i`只在原重力补项输入上开启`sp1_two_interval_history=true`，默认false。实体/HG共用上一接受dt，当前材料系数不变，对应有效应力的几何项同步。`dynamic_interval_weights.csv`补充两核在干净/拒绝/重放时的h/hp与权重；原诊断文件列不变。

`accepted_interval_reject.i`为独立1D/1单元/无材料测试，用真实Transient＋ConstantDT的测试派生Stepper强制一次失败，观察hp经减步及恢复后保持；不代表CDP自然失败事务或坝体计算。总时间.007秒、最多8接受步；不得将该测试Stepper用于生产输入。
