# MOOSE 4bce02d 弱平面应力 × B-bar Jacobian 修正(Q2-009 失稳根因修复)

## 背景

dam-2b 纯 EQ 0~5 s Job(D2EQ-Q4-BBAR-WEAKPS-5S-BASE-009)在 t≈4.85 s 发生
`DIVERGED_LINE_SEARCH` 失稳。定位链(详见 damASR
`docs/verification/evidence/2026-09-25-dam2b-ls010-expert-package/我方工作项-Q2源码修复与回归计划.md`):

1. PJFNK 对照证明**组装 Jacobian 与残差不一致**(行为学);
2. PETSc `-snes_test_jacobian`:无缩放 3.45%、有 automatic_scaling 20.5%(定量);
3. 纯弹性 / 去 Rayleigh / 静态化 / 零密度对照逐一排除本构、阻尼、HHT、惯性项;
4. 源码通读定位:`Compute2DSmallStrain` 在 `volumetric_locking_correction=true` 时把
   应变迹替换为单元体平均值,而 **strain_zz 耦合的 Jacobian 入口**
   (`StressDivergenceTensors` 的 out_of_plane off-diagonal、`WeakPlaneStress` 的
   diagonal 与 displacement off-diagonal)未携带对应的 `(avg - local)/3` 修正项。

## 本补丁内容

`apply_moose_zz_bbar_fix.py <moose-root>` 对 pristine MOOSE `4bce02d` 打四类修正:

- `WeakPlaneStress`:新增 `volumetric_locking_correction` 参数(由
  SolidMechanics Physics 动作自动前递)、体均值形函数辅助函数,并在 diagonal 与
  displacement off-diagonal Jacobian 中补 B-bar 修正项;
- `StressDivergenceTensors`:新增体均值形函数值辅助函数,并在 out_of_plane_strain
  的 off-diagonal Jacobian 中补 trial 侧与 test 修正项侧的 B-bar 项。

幂等;锚点针对 pristine `4bce02d` 校验。

## 构建方式(计算节点)

```bash
rsync -a --exclude=.git <pristine-moose>/ <patched-moose>/
python3 patches/moose-4bce02d-zz-bbar/apply_moose_zz_bbar_fix.py <patched-moose>
make -j4 MOOSE_DIR=<patched-moose> LIBMESH_DIR=<env> WASP_DIR=<env>
```

## 验证链(门禁顺序)

1. 探针 `-snes_test_jacobian`:ratio 0.0345 → 期望 ~1e-8;
2. 2D 坝 0~5 s:0~2.5 s 与冻结基线身份一致且 4.85 s 失稳消除(2D 有效性报告,
   交项目内部确认);
3. 确认后跑 7 个冻结 CDP 标准算例回归(同一候选二进制)。

未验证前不得将此二进制用于任何对照轮或部署为正式 release。
