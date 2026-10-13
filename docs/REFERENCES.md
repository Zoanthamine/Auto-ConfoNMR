# 来源、引用与核对范围

本流程的方法依据、工作簿来源和软件文档如下。计算条件与原文的差异在各节说明。

## DP4+

Grimblat, N.; Zanardi, M. M.; Sarotti, A. M. **Beyond DP4: an Improved Probability for the Stereochemical Assignment of Isomeric Compounds using Quantum Chemical Calculations of NMR Shifts**. J. Org. Chem. 2015, 80, 12526–12534. [DOI:10.1021/acs.joc.5b02396](https://doi.org/10.1021/acs.joc.5b02396)。应同时查看Supporting Information的参数表及方法选择要求。

核对位置：正文Experimental Section，印刷12533页/PDF第8页；SI S-2。Gaussian09、气相B3LYP/6-31G*几何；NMR考察B3LYP/mPW1PW91、多基组、气相/PCM氯仿；298K对应NMR单点电子能加权。本版气相色散几何和ORCA/Shermo权重属于自定义差异。

## STS

Li, J.; Liu, J.-K.; Wang, W.-X. **GIAO 13C NMR Calculation with Sorted Training Sets Improves Accuracy and Reliability for Structural Assignation**. J. Org. Chem. 2020, 85, 11350–11358. [DOI:10.1021/acs.joc.0c01451](https://doi.org/10.1021/acs.joc.0c01451)。

核对位置：几何/模型正文11352页/PDF第3页，Experimental Methods 11356页/PDF第7页。IEFPCM、B3LYP-D3(BJ)/TZVP几何/Gibbs权重、wB97X-D/6-31G*屏蔽，Gaussian09 D01或Gaussian16加g09defaults。STS处理13C化学环境分类，不是位移排序匹配。

## 上游实现

[chemicalchange/EasyCalcNMR_for_OrgSynLab](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab)，核对固定提交 [095dbc16db04b40371edfaa61a03aa581e5c7d93](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/tree/095dbc16db04b40371edfaa61a03aa581e5c7d93)。参考其DP4+/STS表格与对接方式。本版calibration.json保留提取来源与工作簿哈希，仓库及运行包中的接口表位于confonmr-resources/upstream。上游也有自己的加速/权重近似；不能把其CENSO出现的SMD误判为NMR屏蔽步骤遗漏PCM/G09。

该固定版本的NMR屏蔽采用PCM/G09兼容设置；CENSO阶段的SMD设置属于另一计算步骤。

## 相关讨论

[计算化学公社帖子41199，回复pid270367](http://bbs.keinsci.com/forum.php?mod=viewthread&tid=41199&page=1&extra=#pid270367)。方法参数与计算条件以本页列出的原始论文、软件文档和实际输入为准。

## 软件文档与引用

软件引用需使用实际安装版本的官方citation输出/手册，不把本版脚本作者当作这些程序的作者。

|组件|官方资料|
|---|---|
|Gaussian|[官方默认值与版本说明](https://gaussian.com/relnotes/)、[NMR](https://gaussian.com/nmr/)、[SCRF](https://gaussian.com/scrf/)|
|ORCA|[FACCTs ORCA手册](https://www.faccts.de/docs/orca/6.1/manual/)|
|CREST|[CREST官方文档](https://crest-lab.github.io/crest-docs/)|
|xTB|[xTB官方文档](https://xtb-docs.readthedocs.io/en/latest/)|
|GROMACS|[官方手册](https://manual.gromacs.org/)|
|Open Babel|[官方文档](https://openbabel.org/docs/)|
|VMD|[官方用户指南](https://www.ks.uiuc.edu/Research/vmd/current/ug/)|
|Molclus/isostat|[开发者页面](http://www.keinsci.com/research/molclus.html)|
|Sobtop|[开发者页面](http://www.keinsci.com/research/sobtop.html)|
|Shermo|[开发者页面](http://www.keinsci.com/research/shermo.html)|
|Multiwfn|[开发者页面](http://sobereva.com/multiwfn/)|

引用计算研究时还应明确报告实际采样、优化、单点、布居、NMR和统计条件，包括与训练协议的差异。仅写“使用ConfoNMR”不足以描述方法。

## J-DP4 与两种J标度参考

- J-DP4 原始方法：[Org. Lett. 2019，DOI 10.1021/acs.orglett.9b01193](https://doi.org/10.1021/acs.orglett.9b01193)；[作者操作指南](https://sarotti-nmr.weebly.com/uploads/2/4/2/9/24298506/instructive_j-dp4.pdf)。指南第1/3/5页说明默认B3LYP/6-31G**//MMFF、FC与GIAO、平均及J线性标度。本版DFT几何和自由能权重不等同该校准条件。
- Bally, T.; Rablen, P. R. *Quantum-Chemical Simulation of ¹H NMR Spectra. 2. Comparison of DFT-Based Procedures for Computing Proton–Proton Coupling Constants in Organic Molecules*. J. Org. Chem. 2011, 76, 4818–4830. [DOI 10.1021/jo200513q](https://doi.org/10.1021/jo200513q)。0.9155对应H-only去收缩/核心增广FC协议；普通6-31G(d,p)的参考外推未具备该原协议的误差保证。不是²J专用校准。
- ECD展宽：[Nugroho与Morita，DOI 10.1007/s11418-013-0768-x](https://doi.org/10.1007/s11418-013-0768-x)。FWHM与1/e半宽明确换算；不将谱宽当成高斯标准差。
