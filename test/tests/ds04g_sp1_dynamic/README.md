# DS04G 测试专用动态装配

三真实几何×四区间，12组规定历史，0接受步。当前U来自SP1原生快照；V/A缩放1e-3，预应力缩放0.05且S33投影到弹性平面应变。old/older、κ=0及F为制造历史，不是真实CDP状态；原生产材料与所有内核公式不改。

`DS04GSP1DynamicAuditProblem`仅允许serial/1thread/1QUAD4/8DOF，实际切换内核活动缓存逐项装配R/J与三个尺度FD。拒绝检查先装配试步，再调用同FEProblem.restoreSolutions API，校验接受old/older和上次接受区间未变、R/J重放相同；不冒称自动Executioner拒绝或整坝真实CDP异常恢复。

结果CSV含base/FD/R/J/活动组/材料历史及每H状态；分析按锁定独立Q4质量、当前重力、HG和实体Taylor假设核验。动态分项不匹配需如实报告，不改变生产残差救门禁。本夹具不代表整坝响应改善，不触发七冻结回归。

## 可选两接受区间历史

`e*_two_interval.i`只在原重力补项输入上开启`sp1_two_interval_history=true`，默认false。实体/HG共用上一接受dt，当前材料系数不变，对应有效应力的几何项同步。`dynamic_interval_weights.csv`补充两核在干净/拒绝/重放时的h/hp与权重；原诊断文件列不变。

`accepted_interval_reject.i`为独立1D/1单元/无材料测试，用真实Transient＋ConstantDT的测试派生Stepper强制一次失败，观察hp经减步及恢复后保持；不代表CDP自然失败事务或坝体计算。总时间.007秒、最多8接受步；不得将该测试Stepper用于生产输入。

`accepted_interval_observer.i`复用同一拒步注入，改为无材料的单QUAD4、`CDPObservationProblem`和`CDPObservationConstantDT`调用链，只检查stepper retry/accept及独立TIMESTEP_END execute/output身份，不宣称覆盖材料或swapBack。`observation_enabled=false`与`trace_state_swap=false`均为默认；标准`FEProblem + ConstantDT`不经过observer。诊断输入须同时使用后两个observer类型，OFF覆盖直接调用基类。

`cdp_accepted_output_observer.i`复用`cdp_material_table`数据与单单元CDP加载方式，在真实单QUAD4中小应变求解并开启`trace_state_swap=true`。同一测试Stepper注入一次拒步，终态以真正CSV解析检查材料入口/出口、local trial、实际partition/substep、swapBack候选与后快照、accept/output；候选键必须回连到`material_success`，但仍不据此指定唯一accepted调用。状态文件须无截断、I/O、历史或swap不完整标志。本用例未可靠注入真实材料异常，`material_exception`仅有静态源码检查；完整OFF/OFF与OFF/ON数值对照也未在此声称覆盖，留待后续C06作业。

observer所有CSV固定为`rank,record,time,dt,dt_old,step,attempt,evaluation,execution,output,call,global_iteration,context,exec_flag,element,source_element,qp,partition,N,substep,local_solve,local_iteration,event,field,values,note`。`evaluation`只在真实residual/Jacobian内非零并由RAII恢复；`execution`和`output`分别计数；`call`按rank的真实材料调用单调递增；`attempt`含拒步重试。材料`old16/new16`顺序为backbone塑性6、κt/κc、黏性塑性6、损伤t/c，rotation为row-major九分量；出口末项是最后局部子步backbone残差而非全局SNES残差。异常仅写`material_exception`后沿原路径抛出，不重复积分。详细记录仅目标单元QP0和时间窗；`observation_byte_limit`是每rank全部observer CSV的合计上限，包含各文件表头和attempt-events，完整行生成后若无法容纳则不写并立即标记truncated，独立status JSON不计入。候选队列只表示已完成调用，不单独判定唯一accepted调用。

## 可选HG材料β排除

`sp1_exclude_hourglass_stiffness_damping`默认false，仅`sp1_center_hg`允许开启。开启时只让两个独立HG核使用ζ=0；保留HG弹性Kh和HHT α、实体β、质量α、CDP与全部其他控制。`e*_exclude_hg_beta.i`与`e*_two_interval.i`仅此一行差异。三个几何×四制造历史用于实际R/J、old保持与拒绝重放；不是整坝自然失败、响应改善或完整20秒验收。

`dynamic_interval_weights.csv`新增`zeta_for_weight_check`：HG取实际核参数，实体沿夹具冻结的材料ζ=.003604作独立权重核对；并非新增生产观测。所有原字段及关闭模式原值保持。完整物理门禁由集成仓独立参考比较，不以进程退出成功代替。

## SP1实体C0速度率候选（2026-10-29）

`e*_entity_rate.i`复用既有`enable_finite_newmark_damping=true`，要求SP1相容切线；仅实体散度ζ=0，再由参考构形AD速度率核加βR C0:D。独立HG ζ=.003604、Kh/HHT、质量/重力和材料四表保持。两个新核必须列入TimePeriod静力禁用清单，否则输入被拒绝。默认false及已有LEGACY ON接线保持。

夹具仅在检测到两个速度率核时增加`entity_rate`和`entity_total`分组及总力；原OFF分组/列/材料记录保持。`dynamic_interval_weights.csv`现在从实体实际参数读ζ（ON=0），HG仍取实际参数。拒绝/重放比较包括新核R/J；不是整坝动力计算。

P01/P06真实损伤塑性材料输入另在既有单元测试中检查C0速度率与CDP历史分离，不能把这两点材料核验冒充12组制造历史已成为真实损伤轨迹。有限速度率helper及生产CDP源代码不改。
