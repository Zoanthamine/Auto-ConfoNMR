> 当前版本 V2.1，整合 V1.2 接口与两张 J 参考页。操作入口见 [README](../README.md)。

# V2 模块、输入文件与计算依据

## 1. 模块表

| 基础/附加 | 选项 | 实际计算 | 实际输出 |
|---|---|---|---|
| 基础0 | 默认NMR | 现有采样、DFT优化/频率、ORCA、Shermo，01/02 DP4+、03 STS | V1全部输出，另加构象坐标/布居表 |
| 基础1快速 | `conformers` + `gfn2_energy` | 保留RESP2/MD/GFN1/GFN2采样；跳过DFT和谱学性质 | XYZ与GFN2电子能估算布居 |
| 基础1精修 | `conformers` + `free_energy` | 采样、DFT优化/正频率、ORCA、Shermo；跳过三套NMR | XYZ与自由能布居 |
| 附加2 | `coupling = true` | 独立气相B3LYP/6-31G(d,p) FC-only和GIAO，G09Defaults | 每构象J、H–H矩阵、独立屏蔽、C/H/J工作簿 |
| 附加3 | `ecd = true` | 旧三套TDDFT，默认30态/SMD甲醇 | 逐态表、逐构象及加权曲线、PNG/SVG |

0/1二选一；J、ECD可同时选择，也可独立附加到精修基础1。快速档不与J/ECD混用。未选择的模块不会执行。附加模块复用本项目已完成优化和布居的构象，仍受DFT预算和累计布居限制。

## 2. 只改 project.ini

`[模块]`属于计算身份：基础模式、布居模型、J/ECD开关、ECD方法/态数/溶剂。提交后修改这些项必须另建项目；不能对活动项目打开新模块。

`[模块后处理]`不重跑量化：导出数量、J探索比较开关、实验表路径、ECD展宽/规范/显示平移/能量缩放/波长范围。计算结束后修改，再运行03重新导出；04刷新集中查看文件。先保留旧结果版本。

`top_n=6`表示导出权重最高的6个，**不是DFT预算**。增加`top_n`不会补算被预算截掉的构象。`max_conformers=6`才是送DFT的预算，0为全部。

## 3. 构象比例如何读

`最终结果/附加模块/conformers/候选/populations.csv`包含全集布居、是否导出、子集内比例及累计覆盖；`summary.json`记录温度、能量来源、可用构象数和前N覆盖。坐标分别为`cXXXX.xyz`及多构象`top_conformers.xyz`。

快速档用GFN2电子能做Boltzmann估算，未计熵和热校正；精修档用ORCA单点与Shermo热校正的G。两档都只代表当前采样/计算集合。例：前6个累计占80%，子集混合表可归一到100%，但全集覆盖仍是80%。

PDB输入走现有采样。高级XYZ输入可跳过采样；快速档必须在`job.json`显式声明`input_energy_model: "gfn2"`且XYZ注释为Hartree能量。不可用任意软件能量冒充GFN2。

## 4. J与dJ-DP4的实验输入

没有实验表也能提交J计算并得到矩阵。之后准备以下四份表；从`.template.csv`复制成同名正式CSV并填写，UTF-8编码。

### experiment.csv 和 assignments.csv

J分支使用 `[对比] experiment/assignments` 显式指定的C/H表，以及 `[模块后处理]` 指定的J表；同一甲基3H一组。这里**使用另算的B3LYP/6-31G(d,p)气相屏蔽**。V1.2整合后的多样品扫描仅用于两组DP4+与STS，不自动将J/ECD匹配到每个样品。多个实验点的J数据不同，需分别保存并修改显式文件路径，保留各轮模块快照；无真实实验J时仅输出计算矩阵和工作簿。

### experiment_j.csv

| 列 | 填写内容 | 示例（演示，不是真实实验） |
|---|---|---|
| signal | 唯一耦合信号名称 | J_Ha_Hb |
| j_hz | 实验J，Hz | 7.5 |
| mode | `signed`保留符号；`magnitude`比较绝对值 | magnitude |
| bond_count | 本版dJ只允许3键的3JHH | 3 |

### assignments_j.csv

| 列 | 填写内容 |
|---|---|
| candidate | 与项目候选ID相同 |
| signal | 与实验J信号相同 |
| atoms_a | 第一个氢原子或等价氢组；用分号连接，例如14;15;16 |
| atoms_b | 第二个氢原子或等价氢组，例如20 |

**编号以本次`sampling/atom_map.json`中的calculation_atom为准。** 每候选、每信号都填写。软件核查H元素、原子组不重叠、键连距离为3、同一耦合对不重复。先平均带符号J，再按`magnitude`取绝对值；不会先取绝对值再平均。采样的`molecule.mol2`提供键连证据；XYZ跳过采样时，本开发版只导出J矩阵，不支持自动dJ比较。

`coupling_allow_exploratory=false`默认导出对比输入表；审核后改true再运行03，会输出H、C、H+C、J、all概率和逐项误差。每个用于回归的C/H核种至少3个不同信号。概率为候选集合内的相对概率，不是绝对正确率。

### 原文与本版差异

原J-DP4文献：Grimblat等，*Org. Lett.* 2019, 21, 4003–4007，[DOI](https://doi.org/10.1021/acs.orglett.9b01193)。默认统计表针对B3LYP/6-31G**//MMFF，耦合只取FC。上游[EasyCalcNMR dJ_DP4](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/tree/095dbc16db04b40371edfaa61a03aa581e5c7d93/dJ_DP4)使用独立FC与屏蔽单点，含G09兼容关键词。本版按同一FC/屏蔽计算等级新建分支，继续用用户DFT几何及ORCA/Shermo布居，所以必须保留自定义条件标记。

数值实现核对作者工作簿：J斜率0.9509、截距−0.1405、σ=0.992、ν=3.06；C/H统计参数保存在`coupling.PARAMETERS`。传统Excel `TDIST`将自由度截为整数，本版明确复现该行为并在元数据记录。采用对数似然避免连乘下溢。dJ-DP4独立于DP4+；本版只提供dJ分支，iJ约束构象搜索和iJ/dJ组合尚未实现。

`dJ_DP4_tables.xlsx`是本项目新建工作簿，含Main、各候选原始J矩阵、各候选“经验标度”页、条件和来源；启用探索比较后另含概率页。标度页按 `(J_raw + 0.1405)/0.9509` 对整张H–H矩阵非对角元输出一次标度后的带符号数值，不自动筛选键距。页首明确注明该经验标度仅适用于³JHH，使用者应按原子归属自行挑选；²J及长程J的本页代入值未经过对应校准验证。对角线为空。等价氢组先平均原始值，再按实验符号约定做比较；不能将标度页数据作为原始J再次输入统计程序。新增页不会改动原始CSV、构象权重或概率算法，仍标明自定义几何/权重的校准差异。

`dJ_DP4_input.csv`的Type/Exp/候选列可用于向作者原表粘贴，Signal为额外追溯列。dJ工作簿由本脚本生成；随附DP4+/STS接口表的来源和许可见辅助资源目录。

## 5. ECD的输入、方法与图片

ECD默认保留旧流程参数：

| ECD编号 | TDDFT方法 | 溶剂与激发态 |
|---|---|---|
| 01 | PBE1PBE（PBE0）/TZVP | SMD甲醇，30态 |
| 02 | CAM-B3LYP/6-311G* | 同上 |
| 03 | B3LYP/6-311G* | 同上 |

三套均保留旧`Int=Fine`和`EM=GD3BJ`，不继承NMR的PCM或G09Defaults。只支持闭壳层单重态TDDFT；其他情况需扩展并重新验收。共用气相几何和气相自由能仍有溶剂条件差异；这不是甲醇中完整构象自由能研究。

默认读取速度规范，兼容旧Multiwfn选项；原始表同时保留length和velocity值。可在模块后处理选`length`；不要为了匹配结论反复切换。展宽采用高斯函数，参数`ecd_fwhm_ev=0.5`为**半高全宽FWHM**，不能与文献sigma或1/e半宽直接等同。公式与换算见[Nugroho与Morita，2014，式1](https://doi.org/10.1007/s11418-013-0768-x)；FWHM=2√ln2×1/e半宽。

实验可选：`experiment_ecd.csv`两列`wavelength_nm,delta_epsilon`，后者单位L mol⁻¹ cm⁻¹。若仪器只有mdeg，先按实测浓度与光程换算，不能直接把mdeg改列名当Δε。未提供实验曲线仍会生成计算图。

每候选/方法输出`transitions.csv`、`curves.csv`、`ecd_raw.png/.svg`；额外有同方法四候选叠加图。CSV列包含原始加权曲线、显式显示修正后曲线和各构象曲线。曲线先按每个构象的原始跃迁展宽再混合，不把不同构象的第n个激发态强行对齐平均。

默认能量缩放1.0、波长平移0nm。手动修正时保留原始图和调整后图，并在metadata中记录；可选归一化图使用任意单位。软件不自动调偏移/宽度，不给出伪ECD候选概率。30态可能不足以覆盖短波段，需检查跃迁覆盖和态数收敛。

## 6. 集中查看、迁移与复现

01/02/03/04入口不增加。新建项目的真实Gaussian中间文件位于`过程文件/候选/conformers/cXXXX/coupling`、`coupling_shielding`及`ecd01/02/03`；文件名含软件、阶段、候选和构象，例如`Gau_ecd02_001RR_C0001.log`。ORCA、checkpoint、采样、Shermo也直接保存在项目内；所有模块结果位于`最终结果/附加模块`。`运行日志/run_*.log`保留阶段开始、结束、路径与返回值，量化正文保留独立软件日志；搜索与聚类输出仍写入总日志。

原集群资源配置与模块化环境按MIGRATION填写。ECD额外需要matplotlib，安装`requirements-ecd.txt`；未启用ECD时不导入绘图库。J使用已配置Gaussian，不需要增加ORCA J模块或CENSO。辅助工具随资源包提供；Gaussian等完整安装按各自许可另行配置。

提交后项目固定workflow及哈希；修改计算参数或脚本应新建项目。失败保存日志和输入，不删除receipt/锁文件绕过保护。新V2项目不读取或覆盖旧生产作业。

## Bally–Rablen H–H参考页

工作簿同时保留原始矩阵、dJ-DP4经验标度页和Bally参考页。Bally页由原始带符号FC矩阵独立生成：`J_reference = 0.9155 * J_raw`，没有截距；不是先做dJ标度再乘因子。整页不筛选键距。该因子来自Bally–Rablen（JOC 2011，DOI https://doi.org/10.1021/jo200513q）推荐的H-only去收缩/核心增广`B3LYP/6-31G(d,p)u+1s[H]`协议及其probe-set/Figure 4说明。它属于总体H–H校准，包含偕位耦合案例，不是单独拟合的²J系数。

当前任务实际算的是普通6-31G(d,p)，并未加紧s函数；因此本页明确标记跨基组外推、未经当前体系验证。新增页不更改Gaussian模板、不重算量化、不修改构象权重、原始CSV或dJ概率算法。不得把本页标成实际用+1s计算，或用文献原协议的误差指标代表本页精度。

## dJ 的 C/H 文件引用规则（V2.1复核）

dJ读取 `[对比] experiment` 和 `assignments` 指定的一组位移表，默认 `experiment.csv`、`assignments.csv`；另读取模块后处理指定的 J 表。不自动循环采用全部 `experiment_样品名.csv`。如果只准备了多样品表，须明确指定其中一组供J使用；J/ECD与NMR批次独立。相对dJ概率至少两个候选，单候选仍可输出原始J矩阵和标度参考页。

长候选名工作表使用确定性摘要后缀，完整身份保留在来源元数据中；短名称保持可读。两张标度页都从原始矩阵独立生成，不能把其中一页作为另一页输入。
