# DS04G 测试专用动态装配

三真实几何×四区间，12组规定历史，0接受步。当前U来自SP1原生快照；V/A缩放1e-3，预应力缩放0.05且S33投影到弹性平面应变。old/older、κ=0及F为制造历史，不是真实CDP状态；原生产材料与所有内核公式不改。

`DS04GSP1DynamicAuditProblem`仅允许serial/1thread/1QUAD4/8DOF，实际切换内核活动缓存逐项装配R/J与三个尺度FD。拒绝检查先装配试步，再调用同FEProblem.restoreSolutions API，校验接受old/older和上次接受区间未变、R/J重放相同；不冒称自动Executioner拒绝或整坝真实CDP异常恢复。

结果CSV含base/FD/R/J/活动组/材料历史及每H状态；分析按锁定独立Q4质量、当前重力、HG和实体Taylor假设核验。动态分项不匹配需如实报告，不改变生产残差救门禁。本夹具不代表整坝响应改善，不触发七冻结回归。
