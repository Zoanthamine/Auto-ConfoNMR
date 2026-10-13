> 本页保留默认NMR基础说明。V2.1模块选择见 [QUICKSTART_V2.md](QUICKSTART_V2.md)，V1.2整合后的多样品归属、检查及批次收集见 [MULTISAMPLE.md](MULTISAMPLE.md)。02已使用实时阶段回执，03支持 --check-only；旧上下样品字段只作兼容。

# ConfoNMR 详细使用说明

本手册按实际操作顺序说明需要什么文件、输入什么、如何验收。安装先读 [迁移配置](MIGRATION.md)，日常简版见 [快速使用](QUICKSTART.md)。本版方法编号为 01/02 DP4+、03 STS；旧版项目有不同编号，不能只改旧结果标签。

## 一、开始前准备清单

|材料|要求|由谁确认|
|---|---|---|
|每个候选的三维 PDB|单个分子、完整H、正确元素列和CONECT；文件名英文|研究者|
|候选集合|相同待研究物种，电荷/质子化/互变异构状态明确|研究者|
|总电荷、多重度|当前简易入口所有候选共用|研究者|
|DFT预算与聚类档位|预算6是现有默认，0=全部；根据资源与分子柔性决定|研究者|
|实验NMR表|信号标识、核种、ppm；记录真实溶剂、温度、文献来源|研究者|
|显式原子归属|每个候选每个信号对应的计算原子号；等价组与化学分类明确|研究者|
|集群配置|config/cluster.json、外部资源、软件环境、队列|安装者/管理员|

默认 NMR 使用氯仿/CDCl3 参数。实验不是该溶剂时，不能仅修改显示标签后沿用此表。现有简易入口不提供任意溶剂统计模型。

## 二、建项目与选项

在模板目录运行 `bash 00_new_project.sh`，按快速说明填写名称、输入目录、yes、电荷、多重度、预算、聚类选项、NMR选项。建项目不提交；它复制输入与独立程序、资源、配置和SHA256清单。

|NMR选项|实际操作|适用|
|---|---|---|
|1 / legacy_pcm_g09|SCRF=(IEFPCM,Solvent=Chloroform)，G09Defaults，GIAO|沿用恢复后的旧屏蔽方案|
|2 / pcm_g16|同一IEFPCM，G16默认值|单独研究版本默认值差异|
|3 / smd_g16|SMD氯仿，G16默认值|保留历史探索条件，不能冒充PCM训练条件|

三种选项只改变 NMR 屏蔽步骤。共同气相几何、单点和 Shermo 权重保持不变。STS 原文明确允许 Gaussian09 或 Gaussian16 加 g09defaults；DP4+ 原文用 Gaussian09，因此 G16 的兼容设置是本项目复现旧默认值的选择。它不能保证所有版本逐位相同。

|聚类|阶段|Edis/kcal mol-1|Gdis/Å|显示温度/K|Nout|Eout/kcal mol-1|
|---|---|---:|---:|---:|---:|---:|
|coarse|GFN1|0.5|0.5|300|100|100|
|coarse|GFN2、DFT|0.3|0.2|300|100|100|
|fine|GFN1|0.2|0.2|350|100|100|
|fine|GFN2、DFT|0.2|0.1|350|100|100|

isostat 温度只影响其布居显示；最终权重由 Shermo 的温度决定。细分类可能留下更多相似构象，并非自动提高精度；预算仍为6时，更多构象可能没有进入DFT。

## 三、修改 project.ini

用文本编辑器保存 UTF-8，保留三个节名 `[方案]`、`[计算]`、`[对比]` 和所有字段，不改字段拼写。

|通常需要检查|含义与限制|
|---|---|
|charge / multiplicity|总电荷/自旋多重度|
|max_conformers|低GFN2能构象送DFT的数量上限，不是MD帧数|
|population_coverage|完成DFT并通过验收的集合内累计布居阈值；默认0.99|
|cores / md_threads|量化核数/MD线程；PBS申请与cores同步|
|memory_gb / orca_maxcore_mb|Gaussian总GiB / ORCA每进程MB；后者×核数需留余量|
|output|固定为过程文件；存储位置通过建立项目时的项目根目录选择，运行后禁止修改|
|temperature / pressure|Shermo K/atm，不是MD退火温度|
|allow_exploratory|默认true，输出含自定义条件差异的探索比较；false时保留输入导出并提示待审核，不会自动修正文献条件|
|experiment / assignments|单点实验与共用归属表|
|upper_experiment / lower_experiment|两个实验点对应CSV|
|upper_assignments / lower_assignments|留空使用共用表；不同归属时填写专用文件名|

通常不要编辑 workflow/ 或 job.json。更高级的 MD种子、Shermo缩放、模板路线等属于维护者变更：在未运行的新模板上修改、测试、记录，再新建项目。`计算方案.md` 和 `参数预览` 是实际配置的摘要，首次提交前务必阅读。

## 四、提交、监测、日志

项目内运行 `bash 01_submit.sh`。调度器不可用、记录作业仍活动、配置变化、输入或资源不完整时会停止。提交日志、作业号、提交历史保留在本项目。候选和构象串行推进，Gaussian/ORCA/CREST各自可使用配置核数；不代表四候选同时各占一节点。

运行 `bash 02_status.sh` 查看状态。`运行日志/run_日期_作业号.log` 保留进度、搜索/聚类及工具输出；Gaussian（含RESP）与ORCA详细输出只写独立文件。总日志保留它们的起止、返回值、工作目录和日志路径，便于定位失败。

|查什么|去哪里|
|---|---|
|完整运行过程/PBS结束码|运行日志/run_*.log及调度器输出|
|MD是否结束、LINCS/NaN|过程文件中的采样目录/md*.log、md*_run.log、md*_check.log|
|Sobtop补参警告|sobtop.log、topology_quality.json|
|GFN1/GFN2、聚类结果|采样目录/gfn1、gfn2及isostat.log|
|实际Gaussian输入/输出/失败优化|过程文件/候选/conformers/c0001/opt1、opt2、opt3、freq、nmr01等|
|ORCA输入/单点输出|过程文件/候选/conformers/c0001/sp/|
|最终构型|过程文件中的构型汇总、Gaussian输出末尾几何、原始manifest|
|每个构象比例|Shermo权重汇总CSV与原始conformers.json中的weight/nmr_weight|
|原子号对应关系|采样目录/atom_map.json、input.xyz、原始PDB|
|不再选用的构象|conformers.json的exclusions及DFT去重记录|

新项目过程文件目录直接保存实际记录。02/04只更新索引、构型和汇总，不建立外部链接；文件名规则与目录示意见OUTPUTS.md。

## 五、实验表与归属表：逐列填写

实验表例子（仅示意，不能当作真实结构数据）：

```csv
signal,nucleus,ppm
C1,C,15.2
C2,C,40.3
C3,C,130.4
H1,H,1.2
H2,H,2.4
H3,H,7.1
```

归属表使用与输入文件名归一化后一致的候选ID。例如 candidate_A.pdb 对应 candidate_A：

```csv
candidate,signal,atoms,hybridization,carbon_type,heavy_atom
candidate_A,C1,1,sp3,sp3-CH3,None
candidate_A,C2,2,sp3,sp3-CH2/CH/C,None
candidate_A,C3,3,sp2,sp2-CH,None
candidate_A,H1,4;5;6,sp3,,
candidate_A,H2,7,sp3,,
candidate_A,H3,8,sp2,,
```

每个其他候选也必须提供全部信号；不能漏填一个候选或按预测位移给它重新排序配对。分号分隔等价原子；CSV列仍由逗号分隔。每种参与拟合的核至少三条非退化信号。氢的 sp2/sp3 按其所连原子的化学环境核对，不能从截图的空白/X机械推断。

DP4+ 为C/H分别填写 `hybridization`。STS为每个C填写 `carbon_type`，**精确合法值以 calibration.json 的 sts.CDCl3.scale 键为准**，可用以下只读查看命令（高级用户）：

```bash
python3 -c 'import json; print("\n".join(json.load(open("workflow/calibration.json"))["sts"]["CDCl3"]["scale"]))'
```

类别含sp3的甲基/其他碳、sp2-CH/其他、羰基等十类，具体拼写不自行翻译。`heavy_atom` 没有对应修正填 `None`；有重原子时必须在参数表支持的类别组合内，不能随意写元素标签让程序猜。当前Excel导出要求同一信号在不同候选中的分类兼容，否则会停止；不能为绕过此检查篡改真实分类。

赋值必须从 atom_map.json 核对 `input_serial` 与 `calculation_atom`。PDB序号与列表行号并不总一致。同一甲基的3H先按原子屏蔽平均，再做构象加权，以一个信号参与概率。区分“构象混合比例”和“候选统计概率”；后者不是实验样品的异构体浓度。

## 六、实算完成后的验收

03内部会复核manifest、屏蔽文件哈希、原子集合、几何和布居，但不能代替化学审核。正式解释结果前检查：

1. 所有候选state.json为complete；PBS返回0；无残余活动同批作业。
2. 实际Gaussian优化/频率/NMR、ORCA正常结束；优化失败分支保留且有明确后续成功分支。
3. 频率模式数完整，有限且严格正值；被拒绝构象的能量与原因可追溯，不能悄悄删除低能异常结构。
4. Shermo设置、ORCA能量和频率文件按构象ID对应；weight和nmr_weight各自归一；覆盖只针对计算集合。
5. MD有完整退火周期，日志无LINCS/NaN；核对Sobtop缺失项、键连、原始与优化后立体构型。代码的体积符号检查是有限启发式，不能证明完整化学同一性。
6. 程序、资源、输入、方法布局和实际关键词与本项目记录一致；实验归属、等价组与分类独立审核。

## 七、比较和导出

运行 `bash 03_compare.sh`，然后 `bash 04_collect.sh`。若存在single、upper、lower三个实验文件，三个都会分别比较。每次结果有时间戳，不覆盖之前的点版本。

查看 comparison_results.json 的统计概率与误差、weighted_shielding.csv/per_conformer_shielding.csv 的原始贡献、实验映射和参数差异记录，以及 dp4plus_method01/02、sts_method03 工作簿。文件名可能带 exploratory。Excel模板保存公式，由Excel打开重算；脚本JSON不依赖Excel公式缓存。最高概率仍需与MAE/RMSE、异常信号、候选集合和外部结构证据一起判断。

## 八、失败、续跑和备份

软件非正常结束时保留原始输入/输出、receipt、state和日志。先用02和调度器确认已停止，再定位问题。不改原因盲目重提无意义。完整、哈希一致的量化阶段可复用；改了输入输出会被拒绝。**不完整采样阶段当前不支持自动原地续接**：出现“Incomplete sampling”时保留目录，用新项目/新run并记录原因，不删除中间文件绕过保护。

修改科学条件或候选结构必须建新项目；仅改实验数据/显式归属可以复用屏蔽结果再运行03。新项目备份整个项目目录（含过程文件）；04更新汇总，不代替备份。历史旧项目仍需另外保存其完整output目录及链接目标。

## 附表：STS合法类别和重原子修正标签

精确名称来自本版固定参数表。归类化学依据仍应核对原文/原始表；不按位移大小猜类别。

|carbon_type（原样填写）|
|---|
|sp3-CH3|
|sp3-CH2/CH/C|
|sp2-CH2|
|sp2-CH|
|sp2-C|
|C=O|
|X-C=O|
|COOH|
|sp-CH|
|sp-C|

|heavy_atom标签|支持的carbon_type|
|---|---|
|C-S|sp3-CH3, sp3-CH2/CH/C, sp2-CH2, sp2-CH, sp2-C, X-C=O|
|C-Cl|sp3-CH3, sp3-CH2/CH/C, sp2-CH2, sp2-CH, sp2-C, X-C=O|
|C-Br|sp3-CH3, sp3-CH2/CH/C, sp2-CH2, sp2-CH, sp2-C, X-C=O|
|C-OSO2R|sp3-CH3, sp3-CH2/CH/C|

无对应修正填None。超出支持组合时停止并另行研究适用模型，不用错误类别绕过检查。
