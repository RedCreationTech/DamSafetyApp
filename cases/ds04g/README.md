# DS04G 历史输入快照

本目录仅归档原 S4 与 LC05 OFF20 的文本输入，不代表当前分支的默认算例。
数值基线锁定应用源码 `897671fa429da508a19eb02addd3f09d5b77e26f`，
显式 Solver ID 为 `dam-safety-app-ds04g-newmark-897671f-mpi4`。
本目录提交所在分支含后续观测/测试代码，不能据目录位置替代该基线源码。

完整可提交输入包以 damASR 仓库 `examples/ds04g-lc05-off20-v1/`
及 LIMS `api/damage_models/dam-cdp-2d-gravity-v0/` 为准。
本仓库按既有规则不提交 Exodus 网格；运行前从受控模型包取得
`DamSeismic_04g_si.e`，SHA-256：
`e6200cc785062f8b98f5fffaade5602db04fa08fff45420275484b07c1936ebe`。

- `lc05-off20-baseline-20261020/`：自重1秒＋EQ20秒，主输入 SHA-256
  `9813ae91f1433b7798fec4233be0efdb23739009e544405247f04cac9a9da804`。
  除网格外已含地震波和本构四表。`candidate-identity.json` 为当时提交身份快照，
  其中IP、构建路径及二进制哈希属于历史节点，不能用作新节点配置。
- `s4-eq1p5-20261020/`：历史1.5秒输入与地震波；本构四表与OFF20一致，
  准备实际任务时须从受控包复制四表和网格到同一输入目录。

不得把此目录作为可以绕过 LIMS 直接运行求解器的入口。
本次归档不改输入数值，不重新运行计算，不晋级任何候选。
