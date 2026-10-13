# 计算原理与适用边界

本文件把**当前实现**与**原文协议**分开说明。程序可自动执行多个软件，并不意味着采样完整、归属自动正确或概率具有绝对置信度。参考方法与版本来源见 [参考资料](REFERENCES.md)。

## 1. 从候选结构到构象集合

同一个候选结构存在不同构象，溶液实验通常反映快速交换构象的平均核磁响应。不同候选结构则是被比较的结构假说。构象权重属于一个候选内部；DP4+/STS概率在候选集合之间归一，含义不同。

PDB经Open Babel转为MOL2/XYZ。当前电荷步骤在原始几何上分别进行气相和水连续介质 B3LYP/6-31G MK静电势计算，再用Multiwfn进行两次RESP拟合。按 q=(1-delta)q_gas+delta*q_solvent 混合，delta=0.5。代码标记为 legacy-level two-phase RESP mixing，**并非声称完整复现发表的RESP2电子结构协议**。这里的water用于采样电荷，不是NMR屏蔽的chloroform，不能把两种溶剂混为一谈。

Sobtop根据这些电荷和力场参数生成GROMACS拓扑。旧体系曾出现10类缺失扭转势垒置零、2类improper替代项；因此MD高温采样存在势能面近似局限。没有LINCS/NaN只能说明数值运行稳定，不能证明力场准确或没有错过构象。topology_quality.json保留警告，研究者需检查键连和可区分四面体构型变化。

三套MD种子固定，步长0.001ps、200000步，各200ps。退火时间点0/30/80/200ps，温度分别为0→1000→1000→0、0→1500→1500→0、0→2000→2000→0K。这是高温构象搜索模板，不能当作室温平衡动力学或实验温度。

目前执行三条轨迹，但只从**最后一条成功完整轨迹**提取构象，优先md03→md02→md01，没有合并三条轨迹。VMD提取后用CREST的-mdopt分别做GFN1、GFN2优化，再用isostat去重；不是额外执行一次完整CREST iMTD-GC独立搜索。过滤包括原子顺序和可识别立体中心的有向体积符号，属于有限保护，不是完整CIP或图同构判定。

Edis是isostat能量去重阈值（kcal/mol），Gdis是几何阈值（Å）。输入XYZ注释中的能量使用Hartree，二者单位不同。coarse/fine阈值表见详细手册。聚类温度用于工具显示，不能代替Shermo最终混合温度。

## 2. DFT精修与热力学

GFN2集合按能量选最多max_conformers个构象，保留稳定ID。Gaussian气相B3LYP-D3(BJ)/6-31G(d)优化，失败时按三个固定模板逐级重试；成功后DFT聚类、频率计算。最低点要求模式数完整且所有频率严格大于零；有虚频的构象被记录并排除。

ORCA在该几何计算气相wB97M-V/def2-TZVP单点（RIJCOSX、def2/J）。Shermo从Gaussian频率得到热力学校正，以该ORCA电子能替代原低级别电子能：

```text
G_i = E_ORCA,i + G_thermal_correction,i
f_i = exp[-(G_i-G_min)/(R*T)]
w_i = f_i / sum_j(f_j)
```

R与能量单位要一致；实现中Hartree差乘627.509474转换kcal/mol，R=0.00198720425864083 kcal mol-1 K-1。默认T=298.15K、P=1atm、ilowfreq=2、intpvib=100cm-1，各缩放因子1。准RRHO用于低频熵处理，但仍存在柔性模式/自由转子近似。

按G由低到高累计选取至population_coverage=0.99；入选NMR的权重重新归一为 w_i/sum_selected(w)。这是**实际完成并验收的DFT集合内的99%**，不保证真实构象空间99%。预算6、聚类和采样方法都可能改变结果。

## 3. GIAO屏蔽与平均

每个入选构象在同一已优化几何做三种NMR单点。默认三者均NMR=GIAO、IEFPCM氯仿、G09Defaults，泛函/基组见README。原子a的加权屏蔽为 sigma_bar(a)=sum_i(w_i_NMR*sigma_i(a))。等价组内先平均原子，再构象加权；不把相同甲基实验信号当三条独立证据。

屏蔽常数sigma不是实验化学位移delta。DP4+使用对应参数中的TMS参考做未标度转换 delta_u=sigma_TMS-sigma_bar；内部拟合再给出标度位移。STS则按碳化学环境分段由屏蔽换算，带支持范围内重原子修正。两种统计模型的拟合方向不一样，不能互换。

## 4. DP4+ 当前算法

C和H分开处理。将计算未标度位移对实验位移拟合：delta_u=m*delta_exp+b，再反变换 delta_s=(delta_u-b)/m。未标度误差为delta_u-delta_exp，标度误差为delta_s-delta_exp。

按相应泛函/PCM/基组参数，对未标度C/H再分sp2/sp3，使用各自mu、sigma和自由度；标度误差使用对应核种参数。以t分布尾概率形成逐信号似然，在log空间相加避免数值下溢，分别组合H、C及all，最后在候选集合间归一。实现保留原表Excel TDIST对自由度取整的约定。参数是固定训练模型，不随本次实验结论调整。

DP4+原文使用Gaussian09、气相B3LYP/6-31G(d)几何，并考察气相/PCM和多套屏蔽级别；布居来自对应NMR级别电子能，在298K加权。当前使用带D3BJ气相几何和ORCA+Shermo Gibbs权重，因此即使屏蔽级别与参数键匹配，仍有协议差异。

## 5. STS 当前算法

STS意为Sorted Training Sets，是按**碳化学环境划分训练参数**，不是将实验位移大小排序后匹配。当前仅13C、CDCl3参数，每个C按carbon_type选择分段系数；有支持的重原子组合时先修正屏蔽。得到预测位移后，将实验位移对预测位移做线性拟合，误差为 m*delta_pred+b-delta_exp。

每个信号用该类sigma/自由度求双侧t尾概率；Pmean为逐信号概率几何平均，Prel为候选总似然归一后的百分比。同时输出MAE、RMSE、回归斜率截距。高Prel与小绝对误差应一起解释。

STS原文的几何和自由能权重来自B3LYP-D3(BJ)/TZVP(IEFPCM)，屏蔽为wB97X-D/6-31G(d)(IEFPCM)。Experimental Methods明确Gaussian09 D01或Gaussian16加g09defaults。当前共同气相小基组几何和ORCA/Shermo权重不同，保留自定义标记。[STS原文](https://doi.org/10.1021/acs.joc.0c01451)

## 6. 为什么区分PCM、SMD与Gaussian默认值

DP4+溶液参数与本版STS CDCl3表对应PCM/IEFPCM条件；SMD不是同一模型。Gaussian16改变了包括积分/响应网格、IEFPCM实现默认项和物理常数在内的默认值。G09Defaults用于恢复相关旧默认值，不保证所有版本计算逐位相同。[Gaussian官方说明](https://gaussian.com/relnotes/)

DP4+论文没有直接规定G16关键词，因为其使用Gaussian09；本项目选择G09Defaults保持旧方案兼容。STS明确提出G16加该关键词。三个项目选项是可追溯对照，不能混用同一批已计算数据后重新贴标签。

## 7. 如何解释候选评估

概率依赖候选集合、实验归属、训练参数和计算近似；不等于异构体含量，不是独立于模型的绝对可信度。普通非手性环境NMR通常不能区分严格对映体。存在错误键连、构型翻转、缺失低能构象或归属不确定时，高概率也可能误导。最终应结合单晶、二维NMR、耦合、NOE及其他独立证据。
