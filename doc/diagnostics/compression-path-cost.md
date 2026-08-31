# 压缩 CDP 路径与成本诊断

关联集成仓 `TASK-CDPC-UC-D04-TAN / D06-COST`。诊断默认关闭，不改变材料、容差、时间推进或返回映射。正式计算必须通过 LIMS 提交；不能直接以 SSH 启动求解器。

在 `AbaqusCDPStressUpdate` 材料块内启用：

```text
enable_path_diagnostics = true
diagnostic_trace_elements = '54 550 954'
diagnostic_time_begin = 0.015
diagnostic_time_end = 0.05
diagnostic_max_trace_calls = 2
diagnostic_max_tangent_checks = 2
diagnostic_max_failure_samples = 4
```

三个数量上限均按材料对象 / MPI rank / thread 计数，而非整个 Job 共享。元素为 MOOSE 零基索引。轨迹和切线样本是求解过程中真实材料调用，可能是全局试探状态，不自动等于接受状态。应先用更小的单元任务核查数值不变性，再对原结构采样。

输出：

- `cdp_cost_<material>_rankN_threadN.csv`：全运行累计调用、异常失败调用、包含子调用的耗时及失败耗时，包含此前失败分区。每 1000 次材料调用、材料异常及对象销毁时刷新；C06 将 CSV 收入结果目录。非正常强杀可能丢失最后一个刷新周期，不能把残缺记录视为完整成本。
- `logs/cdp_trace_<material>_rankN_threadN.jsonl`：`material_input`、`substep`、`material_failure`、`local_failure_jacobian`、`material_tangent`。随 LIMS Job 日志归档。
- 原有材料属性保持原语义：它们仍是一次最终成功分区的属性，不用它们替代新的累计账。

成本层级：`material → partition → state_linearized → local_linearized → local → newton_step`，其中 residual、AD/FD Jacobian、普通 / AD 主应力计算、factor、backsolve 为细分项；`state_chain` 为跨子步历史导数传播。时间为**包含子调用的时间**，不得直接相加。分 rank 统计后结合 Job 完整墙钟判断负载，不把 rank 时间和作为墙钟。`diagnostic_recompute` 单独记录并从所有外层计时中扣除；轨迹 I/O 单列，但输出、包装、MPI 和其他 MOOSE 成本不由此完整覆盖。

失败 Jacobian 记录实际 9 维未知量、旧骨架状态、应变、缩放、残差、AD 矩阵及 4 个扰动尺度的中心 / 前向差分。主应力过零、重复主值、表节点和非负约束边界可能导致差分跨分支；不能只看一个矩阵误差就判定 AD 错误。诊断异常单独记录，不替换原本的求解异常。

完整切线按固定旧状态、物理张量剪切分量，对 6 个应变分量做 3 个尺度的正负扰动，记录前向 / 后向 / 中心导数、分区数和末分支。仅在分区一致且路径光滑时应用相对误差门禁；末分支一致并不证明所有子步分支一致。每次扰动使用状态副本，并在最后重复原调用核对应力和分区。所有这些诊断重算不计入生产调用计数。

验证要求：既有材料回归 + 失败成本 / 状态不变性回归；LIMS 中关闭 / 开启开关配对；原网格短程真实失败样本。未经这些验证，不得宣称新诊断 release 已完成工程验收或性能优化。
