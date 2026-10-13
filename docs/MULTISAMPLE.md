# V2 多样品核磁对比（整合 V1.2）

此功能用于基础模块 `nmr`，沿用显式原子归属，不按位移大小排序配对。单独构象模块不需要实验表。一个实验样品可以对比多个候选；DP4+至少需要两个候选。

1. 每个实验样品准备 `experiment_样品名.csv`，例如 `experiment_upper.csv`、`experiment_lower.csv`。样品名使用英文、数字、下划线、短横线。表头为 `signal,nucleus,ppm`，核种只写 C 或 H，同一信号在各候选中含义相同。单表 `experiment.csv` 仍兼容。不要把同一数据重复保存成多个实验表。
2. 从 `assignments.template.csv` 准备 `assignments_样品名.csv`；若所有样品使用相同归属，可以只准备 `assignments.csv`。表头为 `candidate,signal,atoms,hybridization,carbon_type,heavy_atom`。每个候选必须恰好覆盖所有实验信号，计算原子号从1开始，等价组用分号分隔。原子号对照项目的 `原子编号对照.csv`、`sampling/atom_map.json` 与三维结构，不能直接套用绘图软件编号。
3. 按实验依据填写 sp2/sp3、STS碳分类和重原子字段；示例见根目录的两份 example CSV。非等价氢分别归属，甲基等价3H可合为一个信号。缺失峰不编造，未纳入信号的原因应记录在实验笔记中。不能为通过检查改变真实分类或归属。
4. 放在项目根目录后运行 `bash 03_compare.sh --check-only`。它检查表头、数字、重复/遗漏、元素、原子号及分类；不进行统计比较。自动检查不能证明化学归属正确。发现问题先修正真实输入，不要改统计参数追求预期结构。
5. 计算完成后运行 `bash 03_compare.sh`。脚本先检查全部表，再逐样品输出两组DP4+与一组STS，保存独立输入快照、哈希和批次清单。某样品失败时保留已经完成的结果和故障记录，禁止覆盖旧批次。打开每样品的 `数据与结果复核.md`、`信号误差.csv` 和输入表复核。拟合后的残差不是外部预测误差。
6. 运行 `bash 04_collect.sh`，在 `最终结果/实验比较/` 查看所有已完成样品及历史批次；`各样品最新结果.csv` 给出最新结果位置。工作簿打开时重算公式，并与JSON核对。

`project.ini` 的 `[对比]` 可修改 `experiment_pattern` 和公共 `assignments`。旧 upper/lower 配置仍可读取。CSV保存UTF-8，数字不带单位和首尾空格，每种参与拟合的核至少三条且位移不全部相同。

J与ECD使用不同格式，分别由 `[模块后处理]` 的 `j_experiment`、`j_assignments`、`ecd_experiment` 指定。它们及 `*.template.csv` 不参与核磁 `experiment_*.csv` 样品扫描，避免混用。附加模块仍按各自显式指定的实验表独立导出，不自动把J/ECD复制匹配到每个核磁样品。J默认输出原始带符号FC矩阵和工作簿，无实验数据不输出比较概率。

基础为 `conformers` 或没有核磁实验表时，03普通运行只刷新已有计算的原始导出和所选模块；`--check-only` 提示无核磁实验表且不刷新模块。01才提交计算，02查看，03导出/比较，04收集。运行中不应执行03。
