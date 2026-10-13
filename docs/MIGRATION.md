# 在自己的集群配置 ConfoNMR

## 先说清楚：复制项目以后，还缺什么？

ConfoNMR通过脚本调用多个计算软件。V2.1运行包、源码包和仓库均附带confonmr-resources辅助资源；Gaussian、ORCA/OpenMPI、GROMACS、VMD、Open Babel、Python和PBS仍需在目标集群准备。

迁移有两种情况：

- **换一台集群，准备开展新计算**：复制 ConfoNMR 模板，准备下面的外部软件和资源，修改 `config/cluster.json`，再用 `00_new_project.sh` 建新项目。这是本说明的主要用法。
- **想接着跑一个已经算到一半的项目**：不能只改几个路径就接着提交。旧项目的固定程序、提交脚本、状态文件和构象清单里有原集群绝对路径，还可能有软件版本差异。本版没有“跨集群断点迁移”按钮。先保留完整旧项目，通常在新环境另建项目；不要覆盖活动计算。

## 1. 哪些要自己准备，哪些已经在项目里？

|要准备的东西|做什么用|迁移时怎么处理|
|---|---|---|
|Linux、Bash、PBS|运行脚本、排队提交|目标集群须提供PBS和qsub/qstat；Slurm/LSF需要另写提交适配，改一个type字段不能解决|
|Python环境|运行ConfoNMR、读CSV/Excel、数值处理|自己建环境，安装requirements.txt中的numpy和openpyxl；这两个包不会替你安装下方化学软件|
|Gaussian及formchk|RESP静电势、DFT优化/频率、三套NMR屏蔽|使用目标集群合法可用的安装，formchk来自同一套Gaussian；还要初始化Gaussian运行环境|
|ORCA及匹配的OpenMPI|高精度单点电子能|分别准备ORCA和其发行包要求的MPI；不能把任意MPI版本拼上去|
|CREST及配套xTB环境|GFN1/GFN2构象精修|安装兼容本脚本参数的版本，并准备所需程序/数据/运行库；本版没有单独的software.xtb字段，需要把配套bin目录加入PATH|
|GROMACS|退火动力学搜索|可使用集群module，也可自己安装；本脚本使用CPU单进程/多线程命令，不依赖GPU|
|VMD|从轨迹提取构象|保留完整VMD安装，不能只拷一个二进制；配置中填写真实VMD执行文件，见下文|
|Open Babel|PDB/MOL2/XYZ转换|可由独立安装或自己的Conda环境提供obabel|
|Multiwfn|RESP电荷拟合|准备无图形程序及配套设置/库；菜单版本变化可能导致自动输入失效|
|Shermo|热校正及构象布居|准备程序；settings.ini由项目按固定参数生成|
|Molclus辅助程序、Sobtop和力场数据|聚类、拓扑、格式转换|准备匹配的完整资源包，文件清单见第4节；不是装一个Python包就有了|
|DP4+/STS上游Excel工作簿|导出可追溯对比工作簿|随附资源包已含固定版本，保留原来源及许可|
|编译器运行库、FFTW等依赖|让以上二进制能在节点启动|按所用程序构建版本准备，可能由module提供，也可能随安装提供；出现缺库时检查实际缺哪一项|

源码包括计算输入模板、协议选项、`calibration.json`、后处理代码和测试；运行包省去开发测试和提取工具。迁移环境时不需要为了“适应路径”去改泛函、基组、溶剂或统计参数。

## 2. 原集群到底加载了哪些 module？

以下是已核对的原集群私有配置，不是要求所有集群都有这些名称：

```bash
module load apps/anaconda3/2021.05
module load apps/gromacs/2020_1
module load 7/mathlib/fftw3-double_3.3.8
module load 7/mathlib/fftw3-double-openmpi_3.3.8
```

|原配置中的module|作用与迁移处理|
|---|---|
|apps/anaconda3/2021.05|加载原站点Anaconda环境；实际运行脚本的Python另外指向自己的tjenv环境。换集群后可改用自己的Python环境，不要求保留这个module名称|
|apps/gromacs/2020_1|提供原流程使用的gmx命令。新集群有兼容GROMACS就改填当地名称；没有则安装自己的GROMACS并填写可执行路径|
|7/mathlib/fftw3-double_3.3.8|原环境加载的双精度FFTW库。新环境是否需要单独加载，取决于你的GROMACS/其他二进制如何构建|
|7/mathlib/fftw3-double-openmpi_3.3.8|原环境加载的FFTW/OpenMPI配套库。它不等于ORCA配置里的OpenMPI安装；不要用它随意替换ORCA需要的MPI|

原配置还把自己的Python环境bin目录、xTB配套bin目录加入PATH，并额外加入Intel编译器运行库目录。这些都属于原站点环境，**新集群不要照抄原用户目录和库目录**。Gaussian另行source自己的g16.profile，ORCA单点另外使用`software.openmpi`指定的MPI环境。

module可以理解为“把集群已经装好的软件接到当前终端”。`module load`本身不下载安装软件。

### 新集群也有 module

先向管理员确认名称和版本，或用`module avail`查看。将你实际需要加载的名称按顺序写入`config/cluster.json`的`modules`数组。GROMACS module若自动加载匹配的FFTW，就不必机械照抄原环境的两套FFTW名称；最终以节点实测为准。

```json
"modules": ["gromacs/your_site_version"]
```

这是一个字段片段，不是完整JSON配置；版本名需要替换。原流程使用GROMACS 2020_1，并不保证任意新版本的行为和输出都一样。

### 新集群没有 module

**可以使用ConfoNMR，但要先把需要的软件和依赖装好。** 管理员安装或用户目录安装均可，前提是计算节点能访问。把`modules`写成`[]`，在`software`中填写程序的绝对路径；有需要的bin目录加入`path_prepend`，动态库目录加入`ld_library_path_prepend`。例如：

```json
"modules": [],
"path_prepend": ["/home/yourname/envs/confonmr/bin", "/home/yourname/apps/xtb/bin"],
"ld_library_path_prepend": ["/home/yourname/apps/required-runtime/lib"]
```

只把真实存在且确实需要的目录写进去；没有额外库就填`[]`。FFTW和编译器运行库是否要另装，取决于软件发行包及其动态链接依赖，不能只看有没有module。ORCA的MPI仍单独填写安装根目录。

如果登录终端能用module、批作业却提示`module: command not found`，说明批作业没有初始化Modules/Lmod。需要按目标集群规定在非交互shell中初始化它，或使用不依赖module的独立安装。本版配置数组只支持module名称，不支持把`source ...`命令塞进去执行；请管理员提供合适的批环境初始化方式，或由维护者在`cluster_config.py`的`runtime_text()`中明确适配并测试。

## 3. 到底改哪个文件？通常只改这一份

在ConfoNMR模板目录复制：

```bash
cp config/cluster.example.json config/cluster.json
```

用编辑器把`config/cluster.json`中的`CHANGE_ME`替换成自己的值。JSON用双引号，不支持中文注释，不要留尾逗号；不使用的数组填`[]`，字段不要删。`cluster.json`是私有机器配置，不要上传公共GitHub。

|字段|你要填什么|容易弄错的地方|
|---|---|---|
|python|自己的Python可执行文件，如/home/yourname/envs/confonmr/bin/python|不是环境目录；00入口还需要能调用python3读取JSON，或预先设置CONFONMR_PYTHON|
|projects_root|新计算项目存放根目录|应在源码目录外；必须所有计算节点可读写且空间充足|
|scratch_root|保留一个合法绝对路径，供历史配置兼容|V1.1新项目实际输出在项目/过程文件，不再靠此字段决定输出位置|
|resource_bundle|第4节资源包目录|不是ConfoNMR/resources占位模板目录|
|software.g16|真实的.../g16/g16|不能填shell alias；程序和安装目录要区分|
|software.formchk|同一Gaussian安装中的formchk|不能拿不匹配的另一安装代替|
|software.orca|ORCA主程序绝对路径|不是整个安装目录；完整ORCA安装要保留|
|software.openmpi|与ORCA匹配的MPI安装根目录|填.../openmpi，不是.../openmpi/bin/mpirun；其下须有bin和lib|
|software.crest|CREST可执行文件|xTB配套路径另外放入PATH|
|software.shermo|Shermo可执行文件|核对大小写和执行权限|
|software.multiwfn|Multiwfn无图形执行文件|程序自动用父目录设置Multiwfnpath；保留配套文件|
|software.obabel|Open Babel执行文件|例如自己环境的bin/obabel|
|software.vmd|真实VMD执行文件，如.../lib/vmd/vmd_LINUXAMD64|当前代码用父目录设置VMDDIR；不要直接填写跨目录启动包装脚本|
|software.gmx|gmx命令或绝对路径|无module时建议填绝对路径；当前调用不是gmx_mpi多进程方案|
|modules|目标集群真实module名称列表|没有module填[]；别照抄原站点名称|
|gaussian_profile|自己的.../g16/bsd/g16.profile|程序从g16路径推导g16root；若module已完整初始化Gaussian，可填空字符串并实测|
|path_prepend|运行所需额外bin目录列表|至少核对自己的Python/Open Babel与xTB环境是否能找到|
|ld_library_path_prepend|确实需要的运行库目录列表|不是程序目录；不要无差别加入多套MPI以免冲突|
|scheduler.type|pbs|Slurm/LSF尚未适配|
|scheduler.queue|自己有权限的队列名|向目标集群确认，不要求沿用原队列|
|scheduler.walltime|队列允许的运行时限|例如72:00:00，按实际限制填写|
|scheduler.node_resource|申请一台节点的写法|Torque可用nodes=1:ppn={cores}；PBSPro可用select=1:ncpus={cores}:mem={memory_gb}gb|
|scheduler.extra_qsub_args|其他qsub参数列表|例如["-A","your_account"]，每个参数独立；没有填[]|

**不需要把自己的路径散改到每个.py/.sh文件里。** 00建项目时自动把这份配置复制进新项目workflow，生成`workflow/cluster_runtime.sh`，并生成项目的01/02/03/04入口。脚本会先加载module，再补充PATH/动态库目录，再初始化Gaussian与Multiwfn环境。ORCA运行时进一步设置自己的MPI路径。

|你想修改什么|应该改哪里|
|---|---|
|软件位置、module、队列、环境库|模板config/cluster.json；然后建新项目|
|本次电荷、多重度、核数、内存、DFT预算、已有科学方案选项|新项目project.ini，首次提交前核对|
|实验数据和原子归属|项目experiment*.csv与assignments*.csv；计算后也可以修改并重新比较|
|Slurm/LSF提交机制|需要开发适配cluster_config.py、project_console.py及submit.sh/run_nmr.pbs，不能只改路径|
|软件新版改了交互菜单或输出格式|由维护者适配sampling.py的resp2/sample或nmrflow.py相关解析，并做小分子实测|
|VMD安装布局不同|优先填写匹配的真实二进制路径；确实不同才适配sampling.py中VMDDIR设置|

正在运行或已有计算的项目继续使用它自己的固定workflow。修改模板的cluster.json不会回头更新旧项目；不要覆盖旧workflow/cluster_runtime.sh来“热更新”。

## 4. 外部辅助资源放在哪里？

随附资源位于程序目录的confonmr-resources/，包含下列辅助程序、力场数据和工作簿。也可按原许可使用自己准备的兼容资源目录：

```text
/home/yourname/confonmr-resources/
  isostat
  sobtop
  sobtop.ini
  atomtype
  mol2_xyz
  xyz2QC
  ATOMTYPE_AMBER.DEF
  ATOMTYPE_GFF.DEF
  CORR_NAME_TYPE.DAT
  LJ_param.dat
  assign_AT.dat
  bondcrit.dat
  bonded_param.dat
  upstream/
    DP4plus.xlsx
    STS.xlsx
```

把此目录填写到`resource_bundle`。Sobtop及参数要使用配套版本；不要把几个发行版里的同名文件混起来。isostat、sobtop等二进制要能在目标Linux运行，并有执行权限。

上游工作簿对应[已核对的EasyCalcNMR提交](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/tree/095dbc16db04b40371edfaa61a03aa581e5c7d93)。本版预期DP4plus.xlsx的Main页、STS.xlsx的CDCl3页及原布局；不能把网上其他表直接改名来替代。DP4plus.xlsx是从上游custom-dp4plus.xlsx透明解密后保存的本地模板，文件哈希会不同；获取原表后按`extract_calibration.py`及来源记录核对，不要仅改文件名。重新提取参数的维护工具额外需要msoffcrypto-tool，日常运行不用它。

00会把资源复制进每个新项目的workflow/resources和workflow/upstream，并记录哈希。这让项目固定所用资源，但不等于Gaussian等外部安装也变成了可移植包。

## 5. 安装Python依赖，再检查环境

以下是示例，把路径换成自己的。公共安装按Python 3.8+准备；目标版本仍须运行回归和小分子验证。

```bash
cd /home/yourname/apps/ConfoNMR
python3 -m venv /home/yourname/envs/confonmr
/home/yourname/envs/confonmr/bin/python -m pip install -r requirements.txt
/home/yourname/envs/confonmr/bin/python cluster_config.py check
```

也可以使用已有Conda环境，把其Python绝对路径写入配置。`cluster_config.py check`只检查配置格式、Python路径和资源包是否齐全，不会检查所有软件真正能算。

新环境中还需要Bash 4+以及qsub、qstat、flock、readlink、tee等系统命令。原环境已使用ORCA 6.1.0/OpenMPI 4.1.8、CREST 2.12（位于xTB 6.6.1配套目录）、Multiwfn 3.8 dev、Shermo 2.4.1；更换版本不能只看可执行文件存在，还要检查参数、交互菜单和输出解析。

若要跑完整源码测试，把合法取得的两个工作簿复制到模板私有upstream目录，再执行：

```bash
mkdir -p upstream
cp /home/yourname/confonmr-resources/upstream/*.xlsx upstream/
/home/yourname/envs/confonmr/bin/python -m unittest discover -p 'test_*.py' -v
```

然后用00建立一个已知小分子的验证项目。在这个**新项目目录**里做不提交计算的预检：

```bash
source workflow/cluster_runtime.sh
/home/yourname/envs/confonmr/bin/python workflow/nmrflow.py check job.json
command -v xtb
```

必须使用新项目的workflow/nmrflow.py：公共模板本身不含已安装的第三方资源，不能拿公共目录直接替代固定项目预检。若module来自目标集群，先确认此source过程在批作业环境也能成功。

常见报错可以这样理解：

|报错|通常需要检查|
|---|---|
|command not found / Software not executable|路径填错、module没加载、执行权限缺失，或安装目录没挂到计算节点|
|error while loading shared libraries|缺程序依赖库；按报错查具体库，用ldd检查可信的本地程序，安装匹配库或补真实库目录|
|module: command not found|批shell没有初始化Modules/Lmod；不是简单把module名称改一下就行|
|Gaussian初始化失败 / GAUSS_EXEDIR缺失|g16路径、g16root推导、profile位置或站点Gaussian环境不正确|
|ORCA串行能启动、并行失败|MPI版本、bin/lib布局、ORCA配套文件及节点权限；不要随便换一套MPI|
|Multiwfn没有生成chg|输入文件、工作目录、菜单版本或运行库不匹配，检查独立日志|
|qsub队列或资源参数错误|实际PBS版本、队列权限、node_resource和账号参数|

最后通过项目01提交一次小分子完整实算，在计算节点核对MD、Gaussian、ORCA、频率、Shermo、原子映射及DP4+/STS导出。登录节点预检和模拟测试通过，不能代替这次实算。不要在登录节点直接运行nmrflow run。

## 6. 迁移时带走哪些文件？

- **发布或给别人安装**：源码、README、requirements、示例配置、计算模板、统计参数及测试。把config/cluster.json、.site-resources、真实工作簿、研究输入和结果留在私有位置，外部资源另外按许可处理。
- **保存自己的研究记录**：备份整个研究项目，包括inputs、project.ini/job.json、workflow、过程文件、运行日志、results、哈希清单和提交记录。旧布局还要保存外部scratch及链接目标。
- **换环境后开始新任务**：不要沿用旧集群绝对路径；准备新环境、改模板配置、重新00建项目。共享目录需容量充足；若在临时scratch，要了解自动清理周期。

这套脚本自动化的是计算和文件衔接。目标软件、许可证、节点环境、原子归属和结果审核仍需要使用者自己确认。

## V2.1随附资源配置

resource_bundle填写程序目录下confonmr-resources的绝对路径；software.crest为bin/crest，software.shermo为bin/Shermo，software.multiwfn为multiwfn/Multiwfn_noGUI，均补全绝对路径。将资源bin加入path_prepend。

建项时只复制运行所需的辅助工具、参数和工作簿到固定workflow，不复制整份辅助源码档。Sobtop的Multiwfn/Open Babel命令从集中配置写入项目副本，原资源不改；环境脚本在share/xtb存在时设置XTBPATH。来源、版本、许可及校验表见程序目录confonmr-resources/README.md。
