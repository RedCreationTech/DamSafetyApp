# DS04G 仅观测候选

仅记录已有计算量，不增加残差、Jacobian、材料积分或扰动计算。默认 `observation_enabled=false`，不建立全局监听器或文件。支持单个Problem、可配置二维QUAD4目标列表、完整6分量平面应变。本轮验证口径仍为4 MPI×1线程；历史默认目标1881/767保留。原数学公式、材料表、HHT/Newmark和步长策略未更改。

`prepare_ds04g_observation.py`以原OFF20输入哈希、专家节点/连接CSV冻结累计2.720秒的控制/禁用/启用三份输入。节点坐标来自ODB单精度导出，源INP舍入最大差小于5e-8m，运行时采用1e-7m唯一物理连接匹配。标签不是运行时ID；每个目标必须全MPI恰好命中一次。

新增对象：`CDPObservationProblem`、`CDPObservationTimeStepper`、`CDPObservationMaterial`。详细窗口按累计时间筛选；默认2.698～2.720秒。材料端回调仅复制原计算量，Old/Older均从已分配存储读取，不通过新getter创建历史。TimeStepper各方法恰好委托一次原基类。仅在候选隔离目录验证，不改变正式registry或默认Solver。

## 输出与身份

每类输出为`<prefix>-<stream>-rankN.csv`，避免MPI共享写文件。公共列包括累计/EQ时间、dt/dt_old、step、独立attempt、evaluation、global_iteration、kind、源/运行单元、QP、material call、partition/N/substep、local solve/iteration。`values`是带引号的空格分隔数值数组，CSV标准双引号转义；17位数字可往返恢复double。Rank取文件名，Job身份由LIMS制品manifest关联，禁止为开发检查冒造job_id。

- `probe-map`：真实节点连接、装配参考/物理QP与JxW×coord；四点权重不假设相等。
- `qp-updates`：真实材料输入、骨架/黏性old、应变、F/转动、应力、最终DT/DC/κ。
- `local-newton`：真实9维unknown、归一化残差、尺度、线搜索候选/拒绝、分支和异常。
- `local-jacobian`：实际已计算9×9，**行优先**，区分AD、FD回退和收敛根线性化；弹性分支不填假矩阵。
- `substep-events`：成功和放弃分区、工作历史、N与dt/N、黏性松弛、原表响应/导数；表段/表尾可由同次κ与冻结原表离线标注。
- `algorithmic-tangent`：完整链式6×6，**列优先**；真正Jacobian装配时另记81项RankFour。两者按同一evaluation/point/call关联；call=0表示未关联，不宣称刷新有效。
- `attempt-events`：原尝试/拒绝前后/接受事件、全局SNES/KSP原因及范数；窗口内残差回调记录目标节点实际输入U。`residual_solution_vector`/`residual_trial_vector`只表示PETSc向量身份，不能单独证明线搜索接受。SNES monitor的接受是Newton迭代接受，不是时间步接受。
- `history-events`：已分配存储的state0/1/2，state0标`stored_current_unverified`，不能把拒绝后残留current直接当已接受状态；节点记录node_id/dof/current/old/older，包括U和已分配V/A。
- `status-rankN.json`：写出耗时/字节/分类计数与缺项标记。任一`capture_truncated`、`io_failed`、`history_incomplete`使完整性审计失败。默认每rank65,536,000字节软上限（四rank总约250MiB）；达到上限停止详细转储并继续原求解/轻量事件。文件I/O错误仅记失败，不篡改材料求解结果。

6分量统一`xx,yy,zz,xy,yz,xz`，剪分量为tensor。9维unknown依次为有效应力/应力尺度、塑性乘子/应变尺度、κt/c增量/应变尺度。9维残差不是Pa，6×6导数不是全局Jacobian；矩阵不一定对称，不能对称化再判断。

## 开发验证与正式任务的边界

隔离构建后运行材料/观测单元测试、三份输入`--check-input`和4MPI零时间步初始化。初始化仅用`Executioner/num_steps=0`；为检查接口可将观测窗口暂设为0，但不得把该测试覆盖范围当成EQ窗口证据。材料失败用已有C++单测构造，正式输入不强造失败。

`audit_ds04g_observation.py <目录> --development --prefix trace`检查零步观测文件；截断试验应返回不完整。正式任务取得真实Job后使用`--job-id <真实ID>`，输出`diagnostic-manifest.json`。审计只检查结构和关联，不能据此认定切线正确或原损伤路径不变。

本次开发验证不提交正式地震Job。EQ1.700～1.720秒实际失败恢复、原基线路径保持、动态观测开销及工程门限仍需后续正式任务与专家确认。全局方向导数/材料扰动重放不属于本候选推进路径。

## 0.5ms 三单元准备（2026-10-21）

新工具分支基于8a，仅修改Python工具/测试/本文；运行二进制仍锁8a，不能将工具提交SHA填成求解器SHA。优先以原0.5ms完整输入及其SHA生成新配置，不重建时间序列：

```bash
python3 tools/diagnostics/prepare_ds04g_observation.py parent.i mesh-nodes.csv mesh-elements.csv prepared \
  --source-sha b1b7575cff284cc2ffae268bc1c891efd9267d5f4dac8574783897cd41f38ff8 \
  --preserve-source --labels 833 1312 959 \
  --observation-start 2.680 --observation-end 2.720 --byte-limit 4294967296
```

保留控制输入逐字节相同；另生成仅在Problem观测块不同的关闭/开启配置。网格、四表和波直接复制。`--preserve-source`不会改时长、步长、输出同步、OFF或材料；source SHA不符即拒绝。生成后仍须核六附件与原manifest，执行节点输入检查和零步验证，不是直接提交命令。

未来正式原文结构审计：

```bash
python3 tools/diagnostics/audit_ds04g_observation.py results --prefix ds04g-observation \
  --job-id ACTUAL_JOB_ID --labels 833 1312 959 --allow-no-reject \
  --window 2.680 2.720 --nominal-dt 0.0005
```

参数中的ACTUAL_JOB_ID必须来自真实任务。该审计要求12QP唯一owner、实际参考位置、矩阵同调用关联、逐rank记录连续、状态计数及81个目标覆盖。默认旧8QP仍要求出现拒绝；新细步允许无拒绝，但`rollback_observed=false`，不能据此说回退已验证。`--development`只检查零步结构，不声称动态窗口或接受链已核。完整材料接受链、黏性恒等式和分区归因需终态原文再做；旧rank1/两单元分析器不可直接用于新12QP。
