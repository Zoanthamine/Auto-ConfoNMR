# Auto ConfoNMR

自动化构象加权核磁计算与候选结构评估 · **V1.1**

Auto ConfoNMR 从三维候选结构出发，串联构象搜索、几何优化、频率分析、单点能和 NMR 屏蔽计算，导出构象加权结果，并根据实验归属进行 DP4+/STS 比较。适用于 Linux 集群，目前提供 PBS 单节点提交方式。

日常操作使用一个建项目脚本和四个项目内脚本。计算条件写在带中文注释的 `project.ini` 中；软件路径、module 和队列设置集中在 `config/cluster.json` 中。

本仓库提供脚本、输入模板和统计参数。Gaussian、ORCA 等外部软件及部分辅助资源需要另行配置。首次安装请先完成后文的[迁移配置](#迁移配置)。

## 日常使用

### 1. 准备结构并建立项目

每个候选准备一个完整的三维 PDB，包含氢原子、元素信息及正确的 `CONECT` 键连。检查质子化状态、电荷和立体构型，一个文件只放一个模型。文件名使用英文、数字或下划线，例如 `candidate_A.pdb`。将本次所有候选放进同一输入目录，避免混入其他结构。

上传到集群后，在程序目录运行：

```bash
cd /path/to/Auto_ConfoNMR
bash 00_new_project.sh
```

按提示填写：

|提示|输入或选择|说明|
|---|---|---|
|新项目名称|例如 `NMR_demo`|英文、数字、下划线或短横线；目录不能已经存在|
|PDB 所在目录|例如 `/home/user/inputs/demo`|填写集群上的路径|
|确认候选|核对文件列表后输入 `yes`|列出的每个 PDB 都会参与计算|
|总电荷|中性分子填 `0`|填写分子总电荷|
|自旋多重度|闭壳层通常为 `1`|按实际电子态确认|
|DFT 最多构象数|默认 `6`；`0` 为全部|按 GFN2 能量选取，是计算预算，不代表构象搜索已收敛|
|聚类档位|`1` 粗分类；`2` 细分类|默认 `1`；具体阈值见后文|
|NMR 方案|`1` IEFPCM/G09 兼容；`2` IEFPCM/G16；`3` SMD/G16|默认 `1`；方案 2、3 用于独立对照项目|

程序在配置的 `projects_root` 下创建独立项目，显示其实际路径。此时只准备文件，尚未提交计算。

### 2. 核对参数并提交

进入新项目，用文本编辑器打开 `project.ini`，核对电荷、多重度、构象预算、核数和内存。文件保存为 UTF-8，Shell 脚本使用 LF 换行。

|文件|提交前检查什么|
|---|---|
|`inputs/`|候选结构齐全，内容正确|
|`project.ini`|计算条件、资源用量、后处理文件名|
|`计算方案.md`|方法、聚类阈值及与校准条件的差异|
|`参数预览/`|三组 NMR 的实际关键词；其中坐标占位符仅供预览|

```bash
cd /path/to/projects/NMR_demo
bash 01_submit.sh
```

看到 `Submitted 作业号` 后，检查 `提交记录.json`。提交脚本会检查配置、软件、输入和活动作业；失败时先查看 `运行日志/submit_*.log`。

默认申请 16 核，Gaussian 总内存 60 GB，ORCA 每进程内存 3000 MB，MD 使用 4 线程。请按节点限制调整；ORCA 的每进程内存乘核数后还应留出系统余量。运行后计算条件与输入会被记录，改变候选结构或科学条件应另建项目。

### 3. 查看进度

```bash
bash 02_status.sh
```

02 直接读取已完成的阶段回执及软件启动/结束记录，不必等待候选全部算完。它显示 GFN2 构象数、DFT 输入数、优化/正频率/ORCA 完成数，以及 NMR01/02/03 分别完成的任务数、当前构象、阶段和独立日志更新时间。

NMR 的“任务数”与“构象数”分别显示：一个构象跑完三套方法是 3 个任务、1 个完整 NMR 构象。覆盖率和性质计算目标总数在最终汇总前显示“待最终汇总”，不把空白理解成未计算。频率剔除不会计入正频率数，多次优化尝试不会重复算成多个构象。

02 刷新 `计算总览.html`、`实时进度.html` 和文件索引，不提交计算。HTML 是最近一次读取的快照，再运行 02 后刷新页面；判断是否还在运行应同时看队列。总日志 `运行日志/run_*.log` 显示阶段进度、搜索和聚类输出；Gaussian/ORCA 的完整输出保存在 `过程文件/` 下的独立日志中。

### 4. 填写实验数据与原子归属

可以在计算期间准备实验表。复制项目中的 `experiment.template.csv` 为 `experiment.csv`，每个实验信号一行：

```csv
signal,nucleus,ppm
C1,C,15.2
C2,C,40.3
C3,C,130.4
H1,H,1.2
H2,H,2.4
H3,H,7.1
```

复制 `assignments.template.csv` 为 `assignments.csv`，为**每个候选的全部信号**填写归属。下面仅示意字段格式，不是真实分子的归属：

```csv
candidate,signal,atoms,hybridization,carbon_type,heavy_atom
candidate_A,C1,1,sp3,sp3-CH3,None
candidate_A,C2,2,sp3,sp3-CH2/CH/C,None
candidate_A,C3,3,sp2,sp2-CH,None
candidate_A,H1,4;5;6,sp3,,
candidate_A,H2,7,sp3,,
candidate_A,H3,8,sp2,,
```

`candidate` 必须与项目中的候选 ID 一致。`atoms` 使用计算原子号，应通过 `过程文件/<候选>/sampling/atom_map.json` 的 `input_serial` 和 `calculation_atom` 核对到原 PDB。PDB 序号不一定等于计算原子号，也不一定等于文件行号。

等价原子用分号分隔，例如甲基三个氢写为 `4;5;6`，作为一个实验信号参与比较。不能复制成三条独立证据，也不能按计算位移大小排序配对。每种参与拟合的核至少需要三条非退化信号。

DP4+ 使用 `hybridization` 的 `sp2`/`sp3` 分类；氢按其所连原子的环境判断。STS 使用碳的 `carbon_type`，合法名称为：

|类别|`carbon_type` 原样填写|
|---|---|
|饱和碳|`sp3-CH3`、`sp3-CH2/CH/C`|
|不饱和碳|`sp2-CH2`、`sp2-CH`、`sp2-C`|
|羰基等|`C=O`、`X-C=O`、`COOH`|
|炔碳|`sp-CH`、`sp-C`|

无对应重原子修正时 `heavy_atom` 填 `None`。`C-S`、`C-Cl`、`C-Br` 支持 `sp3-CH3`、`sp3-CH2/CH/C`、`sp2-CH2`、`sp2-CH`、`sp2-C`、`X-C=O`；`C-OSO2R` 支持前两类饱和碳。超出参数表的组合需要另行评估。工作簿导出要求同一信号在不同候选中的分类兼容，不能为通过检查而改变真实分类。

比较两个实验样品时，分别保存为 `experiment_upper.csv` 和 `experiment_lower.csv`；若归属不同，在 `project.ini` 中填写 `upper_assignments`、`lower_assignments`，留空则共用 `assignments.csv`。程序会处理所有存在且在配置中指定的实验表，因此同一数据不要同时保留为 `experiment.csv`。

### 5. 验收、比较与保存

比较前检查：

1. 全部候选的 `state.json` 为 `complete`，PBS 返回码为 0，本批次无活动作业。
2. Gaussian 和 ORCA 正常结束，频率模式完整且严格为正；失败或被拒绝的构象有可追溯原因。
3. 频率、单点能和 Shermo 文件按构象 ID 对应，构象权重与 NMR 权重各自归一化。
4. MD 无未解释的 LINCS/NaN，核对 Sobtop 缺失参数、键连及搜索前后的立体构型。
5. 实际计算关键词、输入与程序哈希一致，实验归属和分类已独立审核。

随后运行：

```bash
bash 03_compare.sh
bash 04_collect.sh
```

03 读取现有屏蔽结果，不重跑量化计算；每个实验表生成独立的 `results/<点名>_<时间戳>/`，保留之前的结果。至少需要两个候选才能进行 DP4+ 相对比较。单候选仍可完成计算、加权导出和文件整理。

查看 `comparison_results.json` 中的概率和误差，以及工作簿、原子映射和条件差异记录。Excel 工作簿需打开并重新计算公式；JSON 中的概率由脚本计算，不依赖 Excel 的公式缓存。最高概率应结合 MAE/RMSE、异常信号、候选集合和其他结构证据解释，不能当作实验样品的异构体浓度。

04 更新构型、构象分布、索引和结果汇总。备份整个项目目录，包括 `过程文件/`；04 不代替备份。

## 文件组织

```text
NMR_demo/
├── 01_submit.sh … 04_collect.sh
├── project.ini                 # 日常编辑的配置
├── job.json                    # 程序生成的实际配置
├── inputs/                     # 原始候选结构
├── workflow/                   # 本项目固定的程序与资源
├── 参数预览/、计算方案.md
├── 运行日志/                   # 提交、总进度和比较日志
├── 过程文件/
│   └── candidate_A/
│       ├── sampling/           # RESP、拓扑、MD、GFN 搜索与原子映射
│       ├── conformers/c0001/   # 优化、频率、单点与三组 NMR
│       └── shermo/             # 热校正及构象权重
├── 最终结果/                   # 汇总与最新结果位置
├── results/                    # 每次实验比较的独立结果
├── 构象分布与混合比例.csv
├── 文件索引.csv、计算总览.html
└── 提交记录.json、程序版本SHA256.json
```

过程文件直接写入项目目录。Gaussian 和 ORCA 文件名包含软件、阶段、候选和构象，例如 `Gau_opt_candidate_A_C0001_A01.log`、`Gau_nmr02_candidate_A_C0001.log`、`Orca_sp_candidate_A_C0001.gbw`。优化回退步骤用 `A01/A02/A03` 区分，检查点与输入使用相同前缀。

`weighted_shielding.csv` 保存加权屏蔽，`per_conformer_shielding.csv` 保存各构象贡献；构象权重与实际采用的 NMR 混合比例可在构象分布表中查看。截取部分构象后，混合比例会在保留集合内重新归一化。

## 计算方法与自动化过程

### 从结构到加权核磁

流程依次完成：

1. **结构准备与电荷**：Open Babel 转换结构；Gaussian 在 B3LYP/6-31G 级别计算气相与水相静电势，Multiwfn 拟合电荷，按 `delta=0.5` 混合。这里采用两相电荷混合思路，电子结构条件不等同于发表的 RESP2 标准方案。
2. **动力学搜索**：Sobtop 生成拓扑，GROMACS 执行三组固定种子的退火 MD，每组 200 ps，VMD 提取轨迹。当前实现按 md03、md02、md01 的顺序选择最后一组成功轨迹，不合并三组轨迹。
3. **构象精修与筛选**：CREST/xTB 依次进行 GFN1、GFN2 精修，isostat 按能量和 RMSD 聚类，再按 GFN2 能量及 DFT 预算选取构象。这不是完整的 CREST iMTD-GC 搜索。
4. **量化计算**：Gaussian 做几何优化及频率；ORCA 计算单点电子能；Shermo 组合频率热校正与单点能得到自由能及构象分布。优化主方案失败时有两个回退输入，保留每次尝试的输出。
5. **NMR 与比较**：按累计权重选择构象，进行三组 GIAO 屏蔽计算；先平均等价原子，再按构象权重平均，最后用显式实验归属进行 DP4+/STS 评估。

构象权重满足 Boltzmann 关系 `p_i ∝ exp[−(G_i−G_min)/(RT)]`，在已计算的有效集合内归一化。默认累计权重阈值为 0.99；它描述这个集合内的覆盖，不能说明已覆盖整个构象空间的 99%。构象较多或柔性较强时，需要检查预算和搜索收敛。

|步骤|默认条件|
|---|---|
|几何优化、频率|气相 B3LYP-D3BJ/6-31G(d)|
|单点能|气相 ORCA wB97M-V/def2-TZVP，RIJCOSX、def2/J|
|热力学|Shermo，298.15 K、1 atm，低频熵插值阈值 100 cm⁻¹，缩放因子 1|
|01 屏蔽 / 比较|B3LYP/6-31G(d) / DP4+|
|02 屏蔽 / 比较|mPW1PW91/6-311+G(d,p) / DP4+|
|03 屏蔽 / 比较|wB97XD/6-31G(d) / STS，仅 ¹³C|
|NMR 溶剂与默认值|IEFPCM 氯仿，GIAO，G09Defaults|

|聚类档位|GFN1 Edis / Gdis|GFN2、DFT Edis / Gdis|isostat 温度|
|---|---|---|---|
|coarse 粗分类|0.5 / 0.5|0.3 / 0.2|300 K|
|fine 细分类|0.2 / 0.2|0.2 / 0.1|350 K|

Edis 单位为 kcal/mol，Gdis 为 Å。两档的 Nout、Eout 均为 100；isostat 温度仅用于其分布显示，最终权重使用 Shermo 温度。细分类不会自动提高 DFT 预算。

### DP4+/STS 的适用条件

DP4+ 按核种及杂化类别使用固定误差分布，综合缩放与未缩放误差计算候选相对概率。STS 按碳环境分组使用训练参数，并对支持的重原子环境作修正。参数固定在 `calibration.json` 中，保留上游工作簿来源与哈希；不能根据期望结论重新调整参数。

DP4+ 文献的几何与能量加权方案、STS 文献的溶剂中 B3LYP-D3BJ/TZVP 几何及 Gibbs 能，与本流程共同采用的气相几何和 ORCA/Shermo 权重不同。因此本版比较保留 `exploratory` 条件差异标记，不能视为完整复现原文校准体系。

STS 文献使用 IEFPCM，Gaussian 09 或 Gaussian 16 的 G09Defaults 兼容设置。PCM 与 SMD、Gaussian 09 与 Gaussian 16 默认值均不可直接互换。默认 NMR 方案采用 IEFPCM/G09Defaults；另两种设置用于对照，不表示相同校准参数已经适用于这些设置。`G09Defaults` 也不保证跨版本结果逐位相同。

Sobtop 若缺少扭转参数，可能将势垒置零；缺少 improper 时可能采用替代项。实际验证体系出现过这两类提示，虽未在已检查轨迹和构象中发现键连或可区分四面体构型变化，仍会限制采样可信度。每个新体系都应检查其拓扑日志和结构，不能只看任务是否正常退出。

### 程序如何记录与续跑

建项目时复制程序、模板和资源，生成固定快照与 SHA256 清单。提交后，各阶段通过状态文件和回执记录输入、输出及执行结果；符合原始身份和哈希的已完成阶段可以复用。活动作业检查与输出锁用于防止同项目重复运行。

失败后保留输入、输出和提交记录，用 02 及调度器确认作业停止，再根据独立软件日志诊断。修复运行环境后可对同一计算条件续跑；输入或科学条件变化应新建项目。不完整采样阶段目前不支持自动原地续接，出现 `Incomplete sampling` 时应保留现场，在新项目中重新开始采样。遗留锁目录也不能在确认作业仍活动时删除。

|程序|职责|
|---|---|
|`project_console.py`|交互建项目、INI 同步、提交、查询和分点比较|
|`cluster_config.py`|环境配置、辅助资源检查、Bash 环境生成与 PBS 参数|
|`sampling.py`|电荷、拓扑、动力学与 GFN 构象筛选|
|`nmrflow.py`|量化调用、解析、权重、状态、回执与锁|
|`protocol_options.py`|方法分支、NMR 方案、聚类阈值和条件差异|
|`postprocess.py`|显式映射、加权、统计比较与工作簿导出|
|`project_views.py`|文件索引、HTML 总览与构象分布汇总|
|`init_project.py`、`submit.sh`、`run_nmr.pbs`|初始化与调度内部入口，由日常脚本调用|
|`extract_calibration.py`|维护者从上游表提取参数，非日常操作|
|`test_*.py`|解析、统计、身份、调度保护、文件输出与配置测试|

输入模板位于 `resources/`：`03_tem1/2/3.gjf` 为优化主方案和两个回退，`04_tem.gjf` 为频率，`05_tem.inp` 为 ORCA 单点，`06_tem01/02/03_nmr.gjf` 对应三组屏蔽，`md01/02/03.mdp` 为动力学设置。修改这些模板需要方法审核和新项目验证。

## 迁移配置

### 1. 准备外部环境

下载源码：

```bash
git clone https://github.com/Zoanthamine/Auto-ConfoNMR.git Auto_ConfoNMR
cd Auto_ConfoNMR
git checkout V1.1
```

源码不包含化学软件安装。目标节点需要访问以下环境：

|软件或资源|用途|配置位置|
|---|---|---|
|Linux、Bash、PBS|运行与排队；使用 qsub/qstat|本版仅提供 PBS，Slurm/LSF 需另行适配提交器|
|Python 及 numpy、openpyxl|执行脚本、读写数值和工作簿|`python`；依赖见 `requirements.txt`|
|Gaussian、同安装的 formchk|静电势、优化、频率、NMR|`software.g16/formchk`、`gaussian_profile`|
|ORCA、发行包要求的 OpenMPI|单点能|`software.orca/openmpi`，MPI 填安装根目录|
|CREST、配套 xTB 与数据/库|GFN 精修|`software.crest`，xTB 的 bin 目录加入 PATH|
|GROMACS|动力学|`software.gmx` 或对应 module|
|VMD 完整安装|轨迹提取|`software.vmd`，填写实际执行文件|
|Open Babel|结构格式转换|`software.obabel`|
|Multiwfn 无图形版本及配置|电荷拟合|`software.multiwfn`|
|Shermo|热校正与构象分布|`software.shermo`|
|Molclus 辅助程序、Sobtop、力场数据|聚类、拓扑及格式转换|`resource_bundle`|
|上游 DP4+/STS 工作簿|可追溯工作簿导出|资源包内 `upstream/`|

软件须合法取得，并与自动输入参数及菜单兼容。VMD 不能只复制一个二进制，配置示例使用 `lib/vmd/vmd_LINUXAMD64`；Multiwfn 菜单变化可能使自动输入失效。ORCA 的 OpenMPI 应匹配其发行要求；GROMACS 的 FFTW、编译器运行库等依赖按实际构建配置，不能任意拼接版本。

已验证环境包含 Gaussian 16、ORCA 6.1.0/OpenMPI 4.1.8、CREST 2.12/xTB 6.6.1、GROMACS 2020.1、VMD 1.9.3、Shermo 2.4.1、Multiwfn 3.8 开发版。更换版本后应在目标计算节点做小规模验证。

### 2. 填写集中配置

```bash
cd /path/to/Auto_ConfoNMR
cp config/cluster.example.json config/cluster.json
```

将所有 `CHANGE_ME` 替换为实际值，JSON 不支持注释。

|字段|填写内容|
|---|---|
|`python`|安装了依赖的 Python 可执行文件绝对路径|
|`projects_root`|新项目父目录；计算节点必须可读写|
|`scratch_root`|暂存根目录，配置校验仍要求填写；当前项目过程文件直接放在项目内|
|`resource_bundle`|下一节辅助资源包的绝对路径|
|`software`|各程序可执行文件；`openmpi` 为安装根目录|
|`modules`|按加载顺序填写本站 module 名称，无需 module 则填 `[]`|
|`gaussian_profile`|Gaussian 的 `g16.profile` 路径；若其他初始化已完整提供环境可留空并验证|
|`path_prepend`|额外 bin 目录，例如 xTB；填路径数组|
|`ld_library_path_prepend`|实际需要的运行库目录；无额外库填 `[]`|
|`scheduler.queue/walltime`|本站队列与时间限制|
|`scheduler.node_resource`|PBS 资源格式，须包含 `{cores}`；可使用 `{memory_gb}`|
|`scheduler.extra_qsub_args`|本站额外 qsub 参数，按独立参数填数组|

站点路径和队列通常只改这一份文件。项目建立后会冻结配置和环境脚本；修改模板配置不会自动改变已建项目。科研参数仍在 `project.ini` 和计算模板中，不应为了迁移路径而更换方法。

module 的作用是加载集群已经安装的软件环境，不负责安装。已验证站点使用过以下名称，迁移时应换成本站提供的模块：

```bash
module load apps/anaconda3/2021.05
module load apps/gromacs/2020_1
module load 7/mathlib/fftw3-double_3.3.8
module load 7/mathlib/fftw3-double-openmpi_3.3.8
```

其中 Anaconda 模块提供基础环境，脚本实际使用 `python` 指定的独立环境；GROMACS 模块提供 gmx，其依赖可能自动加载。两套 FFTW 是该站点的配置，不要求在其他集群全部照搬，也不能替代 ORCA 所需 MPI。

没有 module 时，先独立安装软件及依赖，设置 `modules: []`，填写绝对路径，将必要的 bin 和库目录加入两个路径数组。环境生成顺序为 module、PATH/库路径、Gaussian profile、Multiwfn 路径。若有 module 但批处理内不可用，需要在本站 Shell 初始化中配置 Environment Modules/Lmod。

### 3. 准备辅助资源包

按以下结构放置匹配的完整资源，保留二进制执行权限：

```text
confonmr-resources/
├── isostat、sobtop、sobtop.ini、atomtype、mol2_xyz、xyz2QC
├── ATOMTYPE_AMBER.DEF、ATOMTYPE_GFF.DEF、CORR_NAME_TYPE.DAT
├── LJ_param.dat、assign_AT.dat、bondcrit.dat、bonded_param.dat
└── upstream/
    ├── DP4plus.xlsx
    └── STS.xlsx
```

每个名称均为独立文件。辅助程序和参数文件应来自兼容版本，不能只准备可执行文件。上游工作簿须与 `calibration.json` 记录的来源匹配；修改训练参数或工作簿应单独审核。本公共包不携带这些第三方二进制、力场资源和工作簿。

### 4. 检查并试运行

使用配置中填写的 Python 安装依赖并检查配置，下面的路径须替换：

```bash
/path/to/env/bin/python -m pip install -r requirements.txt
/path/to/env/bin/python cluster_config.py check
/path/to/env/bin/python -m unittest discover -v
```

配置检查验证字段、路径及资源完整性，不等于所有软件都能在节点启动。随后按日常步骤建立一个小规模项目，限制 DFT 构象预算，检查真实软件调用、频率、权重、命名和输出。至少两个候选及审核后的实验归属才能验证相对结构比较。环境验证成功后再开展完整采样。

没有安装上游工作簿时，部分工作簿测试会跳过；Windows 上 POSIX 锁测试也会跳过。Linux/PBS 是正式运行环境。项目路径建议避免空格，不要在登录节点直接启动量化主流程。

迁移新计算时复制程序目录并重新配置即可。迁移已经计算一半的项目涉及固定快照、绝对路径、状态和软件版本，本版不提供自动跨集群断点迁移；应完整备份后另行评估。

## 验证范围

V1.1 已完成单候选、两个 DFT 构象的真实集群验证：动力学搜索、GFN 精修、优化、正频率、ORCA 单点、Shermo 分布、三组 NMR、加权导出与文件整理均完成，PBS 返回码为 0。两个构象是该次验证的预算，不是通用推荐采样数量。

DP4+/STS 的解析和统计有回归测试，但上述单候选验证不构成多候选真实数据比较的端到端验收。四候选完整实算验收仍在进行，软件运行验证与科学适用性验证分别记录。

## 参考与许可

- DP4+：Grimblat, N.; Zanardi, M. M.; Sarotti, A. M. *Beyond DP4: an Improved Probability for the Stereochemical Assignment of Isomeric Compounds using Quantum Chemical Calculations of NMR Shifts*. J. Org. Chem. **2015**, 80, 12526–12534. [DOI](https://doi.org/10.1021/acs.joc.5b02396)。
- STS：Li, J.; Liu, J.-K.; Wang, W.-X. *GIAO 13C NMR Calculation with Sorted Training Sets Improves Accuracy and Reliability for Structural Assignation*. J. Org. Chem. **2020**, 85, 11350–11358. [DOI](https://doi.org/10.1021/acs.joc.0c01451)。
- 参数表与工作簿接口参考 [EasyCalcNMR_for_OrgSynLab](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab)，固定来源提交 `095dbc16db04b40371edfaa61a03aa581e5c7d93`；参数来源与哈希保存在 `calibration.json`。
- 流程相关讨论：[计算化学公社帖子 41199](http://bbs.keinsci.com/forum.php?mod=viewthread&tid=41199&page=1&extra=#pid270367)。
- 软件说明：[Gaussian](https://gaussian.com/)、[ORCA](https://www.faccts.de/docs/orca/6.1/manual/)、[CREST](https://crest-lab.github.io/crest-docs/)、[xTB](https://xtb-docs.readthedocs.io/)、[GROMACS](https://manual.gromacs.org/)、[VMD](https://www.ks.uiuc.edu/Research/vmd/)、[Open Babel](https://openbabel.org/docs/)、[Molclus](http://www.keinsci.com/research/molclus.html)、[Sobtop](http://sobereva.com/soft/Sobtop/)、[Shermo](http://sobereva.com/soft/shermo/)、[Multiwfn](http://sobereva.com/multiwfn/)。使用时应同时引用所用方法与软件的原始文献。

本版本暂未附带开源许可证；第三方软件、资源及工作簿分别受其自身许可约束。详见 [LICENSE_NOTICE.md](LICENSE_NOTICE.md)。公开仓库不应包含私有集群配置、凭据、软件二进制和研究数据。
