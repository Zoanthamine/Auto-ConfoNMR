> 本页保留默认NMR基础说明。V2.1模块选择见 [QUICKSTART_V2.md](QUICKSTART_V2.md)，V1.2整合后的多样品归属、检查及批次收集见 [MULTISAMPLE.md](MULTISAMPLE.md)。02已使用实时阶段回执，03支持 --check-only；旧上下样品字段只作兼容。

# 脚本、配置与输入输出参考

## 面向日常用户

|脚本|运行位置|输入|输出/作用|
|---|---|---|---|
|00_new_project.sh|模板根目录|交互填写项目、PDB目录、电荷/多重度、预算、两个方案选项|独立项目，尚未提交|
|01_submit.sh|项目根目录|project.ini/输入/集中环境|检查后提交一次，作业号和提交日志|
|02_status.sh|项目根目录|job.json及现有输出|状态、总览、实际过程文件索引；不提交|
|03_compare.sh|项目根目录|实验CSV/显式归属及完成屏蔽结果|各点独立时间戳results目录|
|04_collect.sh|项目根目录|已产生的中间文件|真实文件索引、构型/权重和结果汇总|

旧00_新建项目.sh仅保留为入口别名。新项目统一使用英文短脚本名，避免文档、终端和文件系统中的命名不一致。

## 内部程序

|文件|职责|维护时注意|
|---|---|---|
|project_console.py|交互、建项目、INI同步、活动作业检查、四入口、分点比较|不要绕过queue_guard和identity|
|cluster_config.py|严格读取site配置、安装外部资源、生成quoted Bash环境、PBS参数|只管环境，不改变科学模板|
|init_project.py|显式输入列表复制、候选ID和CSV空模板|内部使用；直接调用不会生成完整便捷项目|
|run_nmr.pbs|PBS固定路径、全程tee、环境、锁、run与退出整理|模板内路径占位符由00替换|
|submit.sh|软件预检、qsub参数、提交记录|内部使用，不直接提交模板示例配置|
|sampling.py|两相RESP、Sobtop、MD/VMD、GFN优化和采样过滤|菜单/软件版本改变须小分子验证|
|nmrflow.py|配置校验、结构/输出解析、DFT/频率/ORCA/Shermo/NMR、状态与锁|构象身份不可按文件排序替代|
|protocol_options.py|三NMR方案、粗细阈值、方法编号、预览/差异记录|旧编号兼容不能删除后假装历史数据改成新方法|
|postprocess.py|映射/加权导出、DP4+/STS统计、工作簿|两类拟合方向不同，参数不可按实验重拟合|
|project_views.py|新项目实际文件索引、HTML、构型/权重表；兼容旧项目链接/复制|不是完整scratch备份工具|
|extract_calibration.py|维护者从原始上游表提取校准参数|不属于日常运行；更换表需版本与数值审计|
|test_*.py|解析/统计/身份/队列/视图/迁移回归|模拟调用不是端到端量化验收|

## 配置分工

|文件|谁改|何时改|
|---|---|---|
|config/cluster.json|安装者|安装/迁移；新项目读取并冻结|
|config/cluster.example.json|维护者|公共占位格式，不可直接当真实环境|
|project.ini|研究用户|首次提交前；实验表引用可在计算后改|
|job.json|程序生成，维护者高级审核|不要与INI同时手动改出矛盾|
|resources/*.gjf、05_tem.inp、md*.mdp|维护者|在新模板中经审核测试后变更|
|resources/shermo_settings.ini与job中的shermo|维护者/研究者|全部设置以实际生成的Shermo settings.ini和manifest核对|
|calibration.json|维护者|固定方法训练参数，不因想要某个候选胜出而调整|
|assignments.csv / experiment*.csv|研究者|量化计算不依赖它们，可后填；比较前审核|

## 关键科学与资源默认值

|项|当前默认|
|---|---|
|候选循环|串行；同候选内量化阶段串行|
|cores / md_threads|16 / 4，需按站点资源修改|
|Gaussian memory_gb / ORCA maxcore|60GiB / 每进程3000MB|
|DFT预算 / 覆盖|6 / 已计算集合内0.99|
|MD seeds|1235761163、-700411449、963984578|
|grompp_maxwarn|0|
|RESP混合|B3LYP/6-31G、气相/水、delta0.5|
|DFT优化模板1|B3LYP/6-31G* em=GD3BJ opt=maxcyc=50 int=fine|
|优化模板2|opt(maxstep=5,notrust,maxcyc=40,gdiis) scf=xqc，其余同级别|
|优化模板3|opt(maxstep=2,notrust,maxcyc=50,gdiis) scf=xqc，其余同级别|
|频率|同DFT气相级别，严格正频率|
|ORCA|wB97M-V/def2-TZVP，RIJCOSX def2/J，气相|
|Shermo|298.15K、1atm、ilowfreq2、100cm-1插值阈值、缩放1|
|NMR|01 B3LYP/6-31G(d)；02 mPW1PW91/6-311+G(d,p)；03 wB97XD/6-31G(d)|
|NMR默认环境|GIAO、IEFPCM氯仿、G09Defaults|

这些是现有实现的记录，不是通用最优建议。改核心科学设置需新项目，重新验收参数适用性。

## 高级命令行

量化输入模板的文件名与用途：`03_tem1.gjf/03_tem2.gjf/03_tem3.gjf`对应优化主方案和两个回退；`04_tem.gjf`是频率；`05_tem.inp`是ORCA单点；`06_tem01_nmr.gjf/06_tem02_nmr.gjf/06_tem03_nmr.gjf`对应固定01/02/03屏蔽分支。旧文件名保留是为了输入/调用可追溯，日常用户无需逐个编辑。

非交互准备只适合熟悉参数的用户，例如：

```bash
python3 project_console.py prepare /home/user/projects/demo --inputs /home/user/inputs/A.pdb /home/user/inputs/B.pdb --charge 0 --multiplicity 1 --cap 6 --nmr-profile legacy_pcm_g09 --clustering-profile coarse
```

先保证python及集中配置可用。生产运行不从登录节点直接执行nmrflow run；使用项目01交给PBS。nmrflow check不计算；nmrflow post使用现有量化数据，但日常应使用03，以保留分点版本和日志。
