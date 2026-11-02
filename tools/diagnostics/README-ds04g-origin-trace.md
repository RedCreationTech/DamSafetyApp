# DS04G 目标材料调用来源：候选说明

本候选在状态推进观测器上增加默认关闭的 `trace_call_origin`。必须同时启用 `observation_enabled`，且选中的单元/QP 位于观测窗口，才对材料入口记录 native stack；不增加有限元残差、Jacobian、材料求解或 MPI 全局归约。记录按 rank/evaluation/call 分开，每次包含模块路径、**未调整的相对返回地址**及动态符号；调用点解释须按模块 SHA、载入基址和目标平台 ABI 在线下验证，不假装 `dladdr` 的动态符号能解析全部栈。

`trace_calls_per_evaluation` 默认 8，上界 1000；`trace_frame_limit` 默认 32，上界 64。真实调用超过上限标 `origin_trace_capped`，多于最大帧标 `origin_trace_truncated`，系统无采样接口、捕获/解析/序列化异常标 `origin_trace_failed`，模块基址不可解析计 `origin_trace_unresolved`，动态符号缺失计 `origin_trace_symbol_missing`。任一标志或帧/调用缺失、I/O/总字节截断、MPI rank 记录丢失时均不能据此宣布来源或接受调用唯一；可通过保留的模块+相对地址和对应未剥离符号表离线补充分析，但仍须证明解析完整，不能以近似名字补值。非 Linux 标失败，没有在运行时引入替代求解。

该记录只解释对应 evaluation 的调用上下文。哪一次产生最终接受材料状态仍必须联合完整状态推进、求解返回的向量候选、所有 rank 及碰撞/多候选核查；不同调用栈并不自动证明某个调用被提交。求解返回后的 SNES 不重新创建；堆栈采样可能触发动态加载、分配内存和额外 I/O，应测量成本，不承诺零耗时或逐位数值不变。

当前仅为候选：本地 Python 静态测试与注入式 C++ 故障测试不代替锁定节点编译及 C06 实测。部署前依次验证库链接、默认关闭路径、目标过滤、栈/调用数上限、未解析符号与输出失败、MPI 与同二进制 OFF/ON，逐轮保留快照和独立报告；无可辨符号、相同栈或无法唯一关联时保留 `inconclusive`，不得扩展到目标 EQ2秒观测轮。
