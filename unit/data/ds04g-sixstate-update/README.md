# DS04G six actual inputs: test-only fixed-N complete update

来自job_20261021_203401_5ygerd，P01/P02=833/Q1 EQ1.7095/1.710；P03/P04=1312/Q3 1.693/1.6935；P05/P06=959/Q0 1.695/1.6955。

每点独立从真实old开始，原17位dt/输入和16维old/出口状态保持。文件依次为id、old strain6、target strain6、old state16、dt、expected core stress6、expected state16、recorded C6 column-major36、backboneDT/DC、branch。状态顺序见source-identity.json，tensor剪不作工程剪转换。

仅测试公开生产接口；不改生产本构、μ或选择器。先P01的G0计量，再其余G0；全部六点G0标记和结果经外部审查后才允许对照。G0是原adaptive linearized/primal及fixed1；对照从同old做fixed2再fixed1，F1结果复用G0的17位文件，所以名义总36次State更新，不串接六点的新轨迹。

未来测试命令需DS04G_SIXSTATE_OUTPUT指向新空目录。比较测试同时要求该目录六个G0通过文件和显式DS04G_SIXSTATE_GATE=G0_ALL_VERIFIED（由审查流程提供）；不绕过G0。仅材料点诊断，不是LIMS Job或整坝验收。
