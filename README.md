# Auto ConfoNMR

自动化构象加权核磁计算与候选结构评估 · **V2.1**

Auto ConfoNMR 从三维候选结构出发，完成构象搜索、几何优化、频率分析、单点能和构象加权 NMR 计算，并根据实验归属进行 DP4+/STS 比较。V2.1 保留 V1.2 的多样品比较和四个日常脚本，增加单独构象搜索、H–H 耦合常数及 ECD 模块。适用于 Linux 集群，目前提供 PBS 单节点提交方式。

日常使用下载 [V2.1 运行包](https://github.com/Zoanthamine/Auto-ConfoNMR/releases/download/V2.1/Auto_ConfoNMR_V2.1_runtime.zip)，解压后按后文配置。也可以克隆本仓库。两种方式都包含辅助程序、配套参数和 DP4+/STS 工作簿；Gaussian、ORCA 等完整软件环境需另行准备。

|下载文件|用途|内容|
|---|---|---|
|`Auto_ConfoNMR_V2.1_runtime.zip`|日常计算|运行脚本、计算模板、说明、辅助程序和工作簿|
|`Auto_ConfoNMR_V2.1_source.zip`|查看或维护源码|运行包内容，以及 `tests/` 和参数提取工具|
|`Auto_ConfoNMR_V2.1_linux_resources.tar.gz`|单独配置辅助资源|辅助程序、参数、工作簿、许可及对应源码；前两个包已包含，无需重复下载|
|`SHA256SUMS.txt`|检查下载文件|三个压缩包的 SHA256|

软件路径、module 和队列写在 `config/cluster.json`，每次计算的条件写在带中文注释的 `project.ini`。辅助资源沿用 V1.2 的版本，目录与来源见 [confonmr-resources/](confonmr-resources/)。

## 选择计算内容

|建项选项|计算内容|主要结果|
|---|---|---|
|基础 `0`：NMR，默认|搜索、DFT 精修、自由能布居、三套屏蔽|加权屏蔽、构象坐标与比例；填写实验表后生成 DP4+/STS 比较|
|基础 `1`，布居 `1`：快速构象|完整搜索、GFN2 电子能布居|默认前 6 个构象坐标、比例及覆盖|
|基础 `1`，布居 `2`：严谨构象|完整搜索、DFT 优化/频率、ORCA/Shermo 自由能布居|默认前 6 个构象坐标、自由能比例及覆盖|
|附加 J|对所选构象另算 FC 耦合和 GIAO 屏蔽|H–H 矩阵、逐构象 J、权重和 dJ 工作簿|
|附加 ECD|对所选构象运行 TDDFT|跃迁、加权曲线、PNG/SVG；可叠加实验曲线|

基础 0/1 二选一。J 和 ECD 可以加到 NMR 或严谨构象模式，快速档不支持。严谨构象加 J 时，只运行 J 所需的一套屏蔽，不额外运行三套 NMR。

快速与严谨档的退火搜索时长相同；快速档省去后续 DFT、频率、ORCA 和 Shermo。导出前几个构象由 `top_n` 控制，送入 DFT 的数量由 `max_conformers` 控制，两个设置分别生效。

## 日常使用

### 1. 准备结构并建立项目

每个候选准备一个完整三维 PDB，包含氢、元素列和正确的 `CONECT` 键连，一个文件只放一个模型。检查质子化状态、电荷及立体构型。文件名使用英文、数字、下划线或短横线，避免只有大小写不同的名称。本次所有候选放进同一输入目录。

在程序目录运行：

```bash
cd /path/to/Auto_ConfoNMR
bash 00_new_project.sh
```

|提示|输入或选择|
|---|---|
|新项目名称|例如 `NMR_demo`，目录不能已经存在|
|PDB 所在目录|集群上的完整路径，例如 `/home/user/inputs/demo`|
|确认候选|核对列表后输入 `yes`|
|总电荷 / 多重度|按实际电子态填写，中性闭壳层通常为 `0` / `1`；所有候选共用|
|基础模块|`0` NMR；`1` 单独构象|
|构象布居，基础 1 时询问|`1` 快速；`2` 严谨，默认 2|
|DFT 最多构象数|默认 `6`，`0` 为全部；快速档不使用|
|导出构象数量|默认 `6`，可填其他正整数|
|增加 J / ECD，精修模式时询问|分别回答 `y` 或 `n`，默认 `n`|
|聚类档位|`1` 粗分类；`2` 细分类，默认 1|
|NMR 方案|`1` IEFPCM/G09 兼容；`2` IEFPCM/G16；`3` SMD/G16，默认 1|

最后一项只影响 NMR 屏蔽，单独构象模式保留该提示但不运行 NMR。J/ECD 使用各自的输入模板。程序在 `projects_root` 下创建独立项目并显示路径，此时尚未提交计算。

### 2. 核对参数并提交

进入新项目，编辑 `project.ini`，核对计算内容、电荷、多重度、预算、核数及内存。保存为 UTF-8，Shell 脚本使用 LF 换行。

|配置区域|可修改的内容|
|---|---|
|`[方案]`|粗/细聚类、NMR 溶剂与版本方案|
|`[计算]`|电荷、多重度、DFT 预算、性质覆盖、CPU/内存、Shermo 温度/压力|
|`[模块]`|基础模式、J/ECD 开关、ECD 方法、态数及溶剂|
|`[对比]`|位移实验表、原子归属、NMR 探索比较开关|
|`[模块后处理]`|导出数量、J 实验表及探索开关、ECD 展宽与显示设置|

```bash
cd /path/to/projects/NMR_demo
bash 01_submit.sh
```

看到 `Submitted 作业号` 后检查 `提交记录.json`。01 先查询活动作业，再同步配置、预检和提交。失败时查看 `运行日志/submit_*.log`，已有作业号时先核实队列再操作。

默认 16 核、Gaussian 总内存 60 GB、ORCA 每进程 3000 MB、MD 4 线程。按节点限制调整，ORCA 每进程内存乘核数后需留出系统余量。`计算方案.md`、`模块选择.json` 和 `参数预览/` 在建项及提交同步时生成，实际计算输入保存在过程目录。

项目提交后保留计算条件和输入。改变结构、计算模块或科学条件时另建项目；实验数据和导出数量可在计算结束后修改。

### 3. 查看进度与中间文件

```bash
bash 02_status.sh
```

02 显示 GFN2 构象数、DFT 输入数、优化/正频率/ORCA 完成数，以及所选 NMR、J、ECD 任务的进度、当前构象和独立日志更新时间。正在计算的单点尚不计为完成；覆盖率和目标总数在汇总前显示“待最终汇总”。

`计算总览.html`、`实时进度.html` 是最近一次读取的快照，再运行 02 后刷新页面。总日志 `运行日志/run_*.log` 显示阶段、搜索、聚类和错误位置，Gaussian/ORCA 的完整输出在 `过程文件/` 下各自的日志中。判断作业是否仍在运行时，同时查看队列。

### 4. 填写位移实验表与归属

仅导出构象、计算 J 矩阵或画计算 ECD 时，不需要填写 NMR 实验表。

需要 DP4+/STS 时，每个实验样品各准备 `experiment_A.csv` 和 `assignments_A.csv`，与 `03_compare.sh` 放在同一目录。可以从项目生成的 `.template.csv` 复制，以下数值仅作格式示例：

```csv
signal,nucleus,ppm
C1,C,15.2
C2,C,40.3
C3,C,130.4
H1,H,1.2
H2,H,2.4
H3,H,7.1
```

归属表为每个候选填写完整的一组信号。下面展示一组，另一个候选按其自己的计算原子号填写：

```csv
candidate,signal,atoms,hybridization,carbon_type,heavy_atom
candidate_A,C1,1,sp3,sp3-CH3,None
candidate_A,C2,2,sp3,sp3-CH2/CH/C,None
candidate_A,C3,3,sp2,sp2-CH,None
candidate_A,H1,4;5;6,sp3,,
candidate_A,H2,7,sp3,,
candidate_A,H3,8,sp2,,
```

`candidate` 与项目 ID 一致，`signal` 与实验表一致。`atoms` 使用从 1 开始的计算原子号，通过 `原子编号对照.csv` 或 `sampling/atom_map.json` 对照到原 PDB。等价原子用分号分隔，非等价氢分别填写，不按位移大小排序配对。参与拟合的每种核至少三条信号，位移不能全部相同。

DP4+ 的 `hybridization` 填 `sp2` 或 `sp3`，氢按所连原子环境判断。STS 的碳分类见下表；无重原子修正时 `heavy_atom` 填 `None`。

|类别|`carbon_type`|
|---|---|
|饱和碳|`sp3-CH3`、`sp3-CH2/CH/C`|
|不饱和碳|`sp2-CH2`、`sp2-CH`、`sp2-C`|
|羰基等|`C=O`、`X-C=O`、`COOH`|
|炔碳|`sp-CH`、`sp-C`|

支持的重原子组合及完整填表规则见 [详细 NMR 手册](docs/USER_GUIDE.md)。分类应依据真实结构，表头和类别名称需原样填写。

多个样品继续增加 `experiment_B.csv`、`experiment_product3.csv` 等，程序逐个比较。优先采用同名 `assignments_<样品>.csv`，缺少时采用共用 `assignments.csv`。单样品也支持 `experiment.csv`。同一份数据只保留一个正式文件名，`.template.csv` 不参与比较。

```bash
bash 03_compare.sh --check-only
bash 03_compare.sh
bash 04_collect.sh
```

`--check-only` 检查位移表的格式、元素、原子号和分类，保存检查记录，不计算概率或改写计算配置。需要项目已经产生原子映射或构象记录。正式 03 等全部候选完成、队列无活动作业后运行，逐样品保存输入快照、结果和日志；04 收集所有已完成样品及历史批次。

每样品查看 `comparison_results.json`、三份工作簿、`信号误差.csv` 和 `数据与结果复核.md`。Excel 工作簿打开后重算公式；信号误差是当前数据拟合后的残差。DP4+ 至少两个候选，概率反映所提供候选间的比较。

多样品文件匹配和结果目录见 [MULTISAMPLE.md](docs/MULTISAMPLE.md)。

### 5. 使用 J 模块的结果

启用 J 后，脚本另算气相 B3LYP/6-31G(d,p) 的 FC-only 耦合和 GIAO 屏蔽，带 G09Defaults。每个候选输出原始 H–H 矩阵、逐构象 J、权重及 `dJ_DP4_tables.xlsx`。没有实验数据时也可查看这些结果。

|工作表|内容|
|---|---|
|`J_<候选>`|加权原始 FC J，Hz，保留正负号|
|`经验标度_<候选>`|整页非对角元 `(J_raw+0.1405)/0.9509`；该校准用于 ³JHH，²J和长程项需自行区分|
|`Bally参考_<候选>`|原始值独立乘 `0.9155`；系数来源为 H-only 去收缩/紧 s 基组协议，普通基组上仅供参考|

两张参考页分别从原始矩阵生成。当前计算使用普通 6-31G(d,p)，没有加入 `+1s`。工作簿页首写明公式和适用条件，原始 CSV 保持原始数值。

需要 dJ-DP4 探索比较时，准备一组 C/H 位移表和两张 J 表：

|文件，默认名|列名|
|---|---|
|`experiment.csv`|`signal,nucleus,ppm`|
|`assignments.csv`|与 NMR 归属表相同|
|`experiment_j.csv`|`signal,j_hz,mode,bond_count`|
|`assignments_j.csv`|`candidate,signal,atoms_a,atoms_b`|

J 的 `mode` 填 `signed` 或 `magnitude`，`bond_count` 填 `3`，两端原子组用分号分隔。程序根据采样 MOL2 核对三键路径，先平均带符号的原始 J，再按实验符号约定比较。直接 XYZ 输入目前支持矩阵导出，不支持 dJ 的键连审核。

J 使用 `[对比] experiment/assignments` 指定的一组 C/H 表，默认是 `experiment.csv/assignments.csv`；使用多样品文件时，在这里明确指定其中一组。J/ECD 不随 NMR 样品自动循环。

默认 `coupling_allow_exploratory=false`，只输出计算表。确认几何和布居与原文 MMFF 校准的差异后，可改为 `true` 再运行 03。相对 dJ 概率至少需要两个候选，工作簿手动改值不会更新脚本统计。

### 6. 使用 ECD 模块的结果

默认三套 TDDFT：01 PBE0/TZVP、02 CAM-B3LYP/6-311G*、03 B3LYP/6-311G*，SMD 甲醇、30 个激发态。ECD 溶剂和态数在 `[模块]` 中修改，与 NMR 设置分开。

每候选、每方法输出 `transitions.csv`、`curves.csv`、`ecd_raw.png/.svg`，另外生成同方法的候选叠加图。默认速度规范、FWHM 0.5 eV，能量缩放 1、波长平移 0；修改显示设置时保留原始曲线。30 态是否覆盖目标波段需检查跃迁表。

实验曲线可选，`experiment_ecd.csv` 两列为 `wavelength_nm,delta_epsilon`，单位分别是 nm 和 L mol⁻¹ cm⁻¹。仪器若给出 mdeg，需用实际浓度和光程换算后填写。没有实验表也会生成计算图。

J/ECD 的完整列名、原理和设置见 [MODULES.md](docs/MODULES.md)。

### 7. 保存和下载

```bash
bash 04_collect.sh
```

|位置|内容|
|---|---|
|`过程文件/<候选>/sampling/`|电荷、拓扑、MD 轨迹、GFN、原子映射和质量警告|
|`过程文件/<候选>/conformers/cXXXX/`|每个构象的优化、频率、ORCA和所选性质输入输出|
|`过程文件/<候选>/shermo/`|自由能及布居记录|
|`最终结果/附加模块/conformers/`|前 N 坐标、全集比例和导出子集比例|
|`最终结果/附加模块/coupling/`|J 矩阵、权重、工作簿和比较记录|
|`最终结果/附加模块/ecd/`|跃迁、曲线及图片|
|`results/comparison_<时间>/`|各样品比较结果、输入快照和批次记录|
|`最终结果/实验比较/`|04 集中的已完成样品及历史结果|

实际过程文件直接保存在项目内。文件名包含软件、阶段、候选和构象，例如 `Gau_opt_001RR_C0001_A01.log`、`Orca_sp_001RR_C0001.log`。Gaussian 检查点、ORCA 波函数及大轨迹也保留在过程目录。长期保存建议备份整个项目，单独下载最终结果会遗漏这些文件。

修改 `top_n` 后运行 03、04 可重新导出坐标和比例，无需重跑量化。计算或后处理失败时保留日志、回执和提交历史，先确认作业停止再诊断；不完整采样阶段需保留现场并在新项目重新开始。常见问题见 [简明 SOP](docs/QUICKSTART_V2.md)。

## 计算流程与参数

PDB → 两相 RESP 混合电荷 → Sobtop 拓扑 → 三套退火 MD → VMD 提取 → GFN1/GFN2 优化聚类。快速档到这里计算电子能布居；精修档继续 DFT 优化、去重、频率、ORCA 和 Shermo，再运行所选性质。

三套 MD 各 200 ps，最高温度分别为 1000、1500、2000 K，并完整降温；从最后成功完整轨迹按 md03→md02→md01 的顺序选择。GFN 使用 CREST `-mdopt` 精修已有构象。电荷采用旧流程的两相 RESP 混合，气相/水 B3LYP/6-31G MK，delta=0.5。

|步骤|默认条件|
|---|---|
|优化 / 频率|气相 B3LYP-D3BJ/6-31G(d)|
|单点能|气相 ORCA wB97M-V/def2-TZVP，RIJCOSX、def2/J|
|热力学|Shermo，298.15 K、1 atm，低频熵处理阈值 100 cm⁻¹，缩放因子 1|
|NMR 01 / 比较|B3LYP/6-31G(d) / DP4+|
|NMR 02 / 比较|mPW1PW91/6-311+G(d,p) / DP4+|
|NMR 03 / 比较|wB97XD/6-31G(d) / STS，仅 ¹³C|
|NMR 溶剂 / 默认值|IEFPCM 氯仿，GIAO，G09Defaults|

|聚类档位|GFN1 Edis / Gdis|GFN2、DFT Edis / Gdis|isostat 温度|
|---|---|---|---|
|coarse 粗分类|0.5 / 0.5|0.3 / 0.2|300 K|
|fine 细分类|0.2 / 0.2|0.2 / 0.1|350 K|

Edis 单位 kcal/mol，Gdis 单位 Å，Nout/Eout 均为 100。粗分类通常保留较紧凑的代表集合；细分类保留更多相近构象的差别，两档搜索时长相同。细分类不会自动增加 DFT 预算，isostat 温度只用于其布居显示。

自由能布居按 `p_i ∝ exp[−(G_i−G_min)/(RT)]` 在合格的已计算集合中归一，快速档将 G 换为 GFN2 电子能 E。默认累计 0.99 从精修集合选择性质构象，选中部分再归一后混合；前 N 坐标另按布居排序导出。这里的覆盖率只描述该集合，不是整个构象空间的覆盖。

DP4+、STS 和 dJ 原文采用的几何及能量加权条件与本流程有差异，比较结果保留 `exploratory` 标记。PCM/SMD、G09/G16 默认值需要分别处理；G09Defaults 也不保证跨版本逐位相同。方法说明见 [THEORY.md](docs/THEORY.md) 和 [V2_INTERNALS.md](docs/V2_INTERNALS.md)。

Sobtop 缺少扭转参数时可能将该项势垒置零，缺少 improper 时可能使用替代参数。这会影响早期采样的势能面，后续量化优化也未必能补回漏掉的构象。查看 `sobtop.log` 和 `topology_quality.json`，并检查 MD 稳定性、键连及立体构型。已验证体系出现过 10 类扭转缺项和 2 类 improper 替代提示。

## 迁移配置

### 1. 下载并准备环境

运行包和克隆仓库任选一种。克隆固定版本：

```bash
git clone https://github.com/Zoanthamine/Auto-ConfoNMR.git Auto_ConfoNMR
cd Auto_ConfoNMR
git checkout V2.1
cp config/cluster.example.json config/cluster.json
```

随附辅助程序适用于 Linux x86-64。Gaussian/formchk、ORCA/匹配 OpenMPI、GROMACS、VMD、Open Babel、Python 和 PBS 仍需自行配置。Python 安装 `requirements.txt`，启用 ECD 再安装 `requirements-ecd.txt`。计算节点需能访问程序和项目目录。

### 2. 填写集中配置

用编辑器替换 `config/cluster.json` 中所有 `CHANGE_ME`，JSON 不支持注释。

|字段|填写内容|
|---|---|
|`python`|安装了依赖的 Python 可执行文件绝对路径|
|`projects_root`|新项目父目录，在程序目录外且节点可读写|
|`scratch_root`|保留有效暂存根路径，当前过程文件直接写入项目|
|`resource_bundle`|`/实际程序目录/confonmr-resources`|
|`software`|各程序实际执行文件；`openmpi` 填安装根目录|
|`modules`|本站需要加载的 module，独立安装时填 `[]`|
|`gaussian_profile`|本机 `g16.profile`，或已由其他初始化完整设置环境时留空|
|`path_prepend`|额外 bin 目录，包含资源包的 `bin`|
|`ld_library_path_prepend`|实际需要的库目录，无额外库填 `[]`|
|`scheduler.queue/walltime`|本人队列权限和本站时限|
|`scheduler.node_resource`|单节点资源格式，含 `{cores}`，可含 `{memory_gb}`|
|`scheduler.extra_qsub_args`|账户等额外 qsub 参数，逐项填写，无则 `[]`|

使用随附资源时，CREST 填 `confonmr-resources/bin/crest`，Shermo 填 `bin/Shermo`，Multiwfn 填 `multiwfn/Multiwfn_noGUI`，均使用完整绝对路径。建项时自动填写 Sobtop 的 Multiwfn/Open Babel 路径，环境脚本自动设置随附 xTB 参数目录。

module 用于加载集群已有软件，不负责安装。将本站模块写入数组；无 module 时填写独立安装路径及必要库目录。ORCA/OpenMPI 需按发行要求配套。PBSPro 可使用 `select=1:ncpus={cores}:mem={memory_gb}gb`，Torque 可使用 `nodes=1:ppn={cores}`；Slurm/LSF 需要另写提交适配。

#### 不熟悉集群配置时

如果对 Linux、module 或 PBS 不太熟悉，可以请能访问文件和终端的 AI 助手协助配置。让它结合本 README、`config/cluster.example.json` 和集群实际环境，查找已安装的软件、填写路径和队列，并解释配置检查中的报错。可以把下面这段话交给它：

> 请帮我在这台 Linux/PBS 集群配置 Auto ConfoNMR。先阅读 README 和配置示例，检查已有软件、module、Python 环境及可用队列，再填写 config/cluster.json。不确定的路径和版本请先核实，计算方法、溶剂、G09Defaults 和统计参数保持原设置。完成后检查配置，并协助建立一个小规模项目试运行。

集群登录在本机完成认证，不需要在聊天中提供密码或私钥。配置完成后，按下一节检查资源并试运行。

### 3. 检查资源并试运行

```bash
cd /path/to/Auto_ConfoNMR/confonmr-resources
sha256sum -c SHA256SUMS.txt
chmod +x bin/crest bin/xtb bin/Shermo multiwfn/Multiwfn_noGUI \
  isostat sobtop atomtype mol2_xyz xyz2QC optional/molclus/molclus
cd ..
/path/to/env/bin/python -m pip install -r requirements.txt
# 使用 ECD 时：
/path/to/env/bin/python -m pip install -r requirements-ecd.txt
/path/to/env/bin/python cluster_config.py check
```

源码包或仓库可运行 `python -m unittest discover -s tests -v`，运行包不含开发测试。随后用小规模项目检查目标节点的真实软件调用、文件输出和权重，再开展正式计算。完整资源版本、module 示例、菜单兼容及跨集群注意事项见 [MIGRATION.md](docs/MIGRATION.md)。

模板更新只影响后续新项目，已有项目保留自己的 `workflow/`。计算中的项目涉及绝对路径和状态文件，跨集群续跑需单独评估。

## 与 V1.2 的区别和验证范围

|内容|V1.2|V2.1|
|---|---|---|
|三套 NMR、多样品比较、01—04 操作|已有|保留|
|实际过程目录、独立软件日志、实时进度|已有|保留，进度按所选模块显示|
|快速/严谨单独构象|未提供|新增，默认导出前 6 个，可改数量|
|J / ECD|未提供|建项时可选，J 含原始及两种参考标度页|
|辅助程序、力场与工作簿|随包提供|沿用同一资源版本|
|路径配置|集中配置|保留，补齐 Sobtop 与 xTB 资源适配|
|提交和数据检查|已有保护|加强排队配置、防重提和只读检查|

程序回归、四种建项预检、Bash 语法及 ShellCheck 已通过。真实有机分子完成过快速搜索、严谨构象和耦合全流程；ECD 的真实节点测试采用 H2O2，完成三套 TDDFT、曲线导出及绘图。具体代码版本与测试记录见 [VALIDATION_V2.md](docs/VALIDATION_V2.md)。

有机多构象 ECD、真实实验 J 的预测误差及 dJ 结构鉴别准确率尚未验证。默认 DFT 预算和已有力场缺项也需要结合具体分子评估。

## 参考与许可

- [DP4+](https://doi.org/10.1021/acs.joc.5b02396)、[STS](https://doi.org/10.1021/acs.joc.0c01451)、[J-DP4](https://doi.org/10.1021/acs.orglett.9b01193)及其[作者操作指南](https://sarotti-nmr.weebly.com/uploads/2/4/2/9/24298506/instructive_j-dp4.pdf)。
- [Bally–Rablen H–H 耦合标度](https://doi.org/10.1021/jo200513q)、[ECD 展宽方法](https://doi.org/10.1007/s11418-013-0768-x)。
- 工作簿接口与统计参数参考 [EasyCalcNMR](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/tree/095dbc16db04b40371edfaa61a03aa581e5c7d93)，来源和哈希保存在校准文件中。
- 相关讨论：[计算化学公社帖子 41199](http://bbs.keinsci.com/forum.php?mod=viewthread&tid=41199&page=1&extra=#pid270367)。完整文献及软件引用见 [REFERENCES.md](docs/REFERENCES.md)。

Auto ConfoNMR 暂未附带开源许可证。随附辅助程序、参数和工作簿保留原作者许可及引用要求，详见 [资源说明](confonmr-resources/README.md)和 [LICENSE_NOTICE.md](LICENSE_NOTICE.md)。Gaussian、ORCA 等完整安装需另行准备；私有集群配置、凭据和实际研究数据不随包发布。
