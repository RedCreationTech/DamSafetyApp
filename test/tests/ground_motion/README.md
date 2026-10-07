# 规定运动诊断

本目录只验证运动、惯性、质量比例阻尼和一次试步拒绝恢复，不代表坝体或CDP验收。

`IntegratedAcceleration`读取原始无表头两列CSV（time,acceleration，SI单位），A线性分段，V/U解析积分，不记录接受步历史；不允许域外查询。主变量q_x相对基底，绝对量显式加同一函数。初始相对加速度内域+0.004、固定底部0，对应原绝对内域0、底部−0.004。质量校正需同时传ground_acceleration/ground_velocity以恢复不规则QUAD4的非零行和载荷。

两个输入为小型自由求解弹性单元，所有Y固定，底部q_x固定。保留HHT和两类Rayleigh阻尼；retry输入主动拒绝一次已收敛试步，接受记录打印17位。`GroundMotionValidationStepper`必须显式validation_only=true，不能用于正式地震计算。

正式FE执行按集成仓库约束经LIMS Facade→C06。独立矩阵解、输入身份和逐Job报告保存在damASR，不在本目录制造金标准或修改冻结CDP基线。
