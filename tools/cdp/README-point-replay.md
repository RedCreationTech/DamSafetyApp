# CDP失败材料点采集与精确缓存候选

TASK-AQ-PERF-020/022/027。当前为本地源码准备，未编译、未节点验证、未证明提速。
不修改正在运行的4304a9f完整渡槽Job；不得为了本工具在运行中的121构建或追加计算。

## 默认与三种模式

`AbaqusCDPStressUpdate`新增参数默认：

```text
enable_exact_result_cache = false
observe_exact_result_cache = false
exact_result_cache_capacity = 65536
exact_cache_diagnostics_prefix = ''
failure_point_capture_prefix = ''
failure_point_capture_limit = 1
failure_point_capture_start = 0
failure_point_capture_end = <maximum finite double>
```

- off：两个缓存开关均false，执行原`integrateLinearized`，不分配缓存。
- observe：仅observe为true，计潜在命中但每次仍执行原积分。
- reuse：仅enable为true，精确命中才复用完整成功结果；两开关同时true拒绝输入。

缓存每材料实例/线程独立，有界、每elem/side/QP一个成功结果；volume使用无side标记，MOOSE的volume/face/neighbor材料实例独立持有缓存。键为旧/新总应变12值、dt、16旧状态，共29个double原始位；参数/表在该材料实例内固定。不近似比较、不跨点共享、不缓存失败。timestepSetup、residualSetup、meshChanged和异常路径清空；旧状态仍每次读取。预期作用域为最近残差求值到Jacobian，框架生命周期/重初始化和拒绝恢复须后续真实检查，不能只凭这段描述标验证通过。

`exact_cache_diagnostics_prefix`启用时写每材料/rank/thread的CSV，生命周期结束给出requests、potential_hits、reused_hits、实际积分次数、容量跳过、失效、峰值entry数/entry payload字节和两类累计微秒。微秒是每rank/thread材料调用累计，不是全Job墙钟；payload不含map节点/分配器。未正常结束只有表头时不伪造统计。

缓存返回值内的局部迭代/分区/分解计数保留原结果语义，不能累加后当作命中调用实际执行工作；实际工作以独立CSV的integrations/reused_hits及计时核查。诊断CSV在材料构造时创建表头，因此后续`--check-input`必须使用独立新目录，不能把预检留下的CSV带入正式Job运行目录。正式提交前须核诊断路径无旧文件。

## 最小失败包

只在显式非空`failure_point_capture_prefix`及时间窗内工作。每rank/thread/材料最多指定次数，I/O失败仅打印诊断，不替换原MooseException，不更新材料状态。包名为：

```text
<prefix>.<material>.rank<R>.thread<T>.record<N>.txt
```

`AQ_CDP_POINT_V1`包含elem/side/QP/位置/rank/thread、真实step/time/dt，原局部/状态/分区参数、29个积分输入double、四张原表完整数值和原异常。浮点值以IEEE754 binary64的16位十六进制保存。包仅记录该rank/thread遇到的首个合格失败，不声称为所有rank的全局时间顺序首个。初步窗口建议总1.19—1.22s，另可为首次分叉收窄；先用候选关闭缓存的控制采集，不与缓存复用混在同一机制轮。新候选采集不能声称恢复旧4304轮未保存的原失败点。

C++读取器对ID、迭代数、分区数和表行数使用十进制无符号范围检查，负数、符号和溢出不接受；不得用无符号回绕后的ID定位失败点。新增的此项C++检查仍待锁定环境编译执行。

后续必须确认C06清单/下载接口实际发布这些TXT、缓存CSV，或在同一Job归档中提供下载引用；缺制品发布不得先跑完再冒称已采集。当前完整Job没有这些新输出，不能从旧舍入日志补造包。

## 本地打包（没有求解）

得到真实包及对应LIMS Job/节点binary身份后：

```sh
python3 tools/cdp/prepare_point_replay.py \
  --packet "/absolute/path/capture.txt" \
  --output "/absolute/new/replay-package" \
  --source-commit <full-capture-source-SHA> \
  --binary-sha256 <actual-capture-binary-SHA256> \
  --job-id <capture-LIMS-Job-ID>
```

打包器不会验证节点身份真实性；manifest明确其为声明值，执行前须对Job manifest/节点SHA核实。只校验包结构/范围、全状态、四表数值的位保持及文件SHA，不执行CDP。本地Python测试只用标记为synthetic的传输样本，不代表真实失败点。

## 后续隔离重放（本轮不执行）

待当前Job终态、用户通知并完成归档，且授权节点测试后，用锁定环境编译同候选的unit二进制，在应用源码根目录执行：

```sh
CDP_REPLAY_PACKET="/absolute/replay-package/capture.txt" \
CDP_REPLAY_OUTPUT="/absolute/new/replay-output" \
  ./unit/DamSafetyApp-unit-opt --gtest_filter=AbaqusCDPPointReplay.ExternalCapturedFailure
```

这属于材料点C++诊断，不是正式FE Job。任何结构计算仍只能通过LIMS提交。
重放器拒绝已存在的输出目录，先无观察、再附加有界trace运行**原**local/state/substep积分器，核原异常逐字一致；否则判失败，不默认“差不多复现”。输出local-trace.csv含每次局部调用、原9未知量/残差、AD/FD候选与接受/非负筛选，trace上限100000；丢弃计数必须为零，达到上限即诊断不完整。sink故障不会改变积分路径，但会使诊断门禁失败。原误差串和重放误差串另存replay-receipt.txt。

原C++单位测试与新测试本轮均未执行；本地没有MOOSE工具链。后续顺序和full0.1/full0.35门禁由damASR集成验证计划约束。七冻结CDP回归继续暂缓，未运行不标通过。
