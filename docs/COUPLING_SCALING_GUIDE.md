# 耦合常数后处理指南（V2.2）

本页说明 `coupling` 模块导出什么、同碳 CH₂ 的环境如何判定、以及后处理时应记录什么。
流程本身**不再自动决定** CH₂ 的标度类别——它只给出原始 FC 与候选公式，判定留给后处理。

## 1. 两条分支各算什么

| 分支 | 作业 | 基组 | 输出 | 标度 |
|---|---|---|---|---|
| 常规 H–H | `coupling`（`08_coupling_fc.gjf`） | 非氢元素 6-31G(d,p)；H 为 8s1p（4 个解收缩 6-31G + 4 个紧 s + 1 个 p=1.1） | 完整 FC 矩阵（Hz，带符号） | `J = 0.9155 × FC` |
| 同碳 CH₂ | `coupling_du4`（`08_coupling_du4.gjf`） | C 4-31G；N/O/F/Li/Be/B 6-31G；Si/S/P/Cl/As/Se/Br/I 3-21G\*；H 为 DU4 的 8s1p | 完整 FC 矩阵 | 按环境取斜率/截距，**不自动应用** |

路由：`#p B3LYP/Gen 6D NMR=(FCOnly) Int=NoBasisTransform G09Defaults`

- 依据：Bally–Rablen 2011（[10.1021/jo200513q](https://doi.org/10.1021/jo200513q)）SI Table S4a 方法代码 `H F H D N`；
  Kutateladze–Mukhina 2014（[10.1021/jo501781b](https://doi.org/10.1021/jo501781b)）SI S2–S4。
- 两者都是 **FC-only 经验估计**，不是 Ramsey 四项全算。

## 2. 导出的文件

| 文件 | 内容 |
|---|---|
| `<候选>_J_result.csv` | 常规分支的加权 FC 矩阵（H×H，`X` 为对角） |
| `per_conformer_J.csv` | 逐构象常规 FC（未标度） |
| `CH2_geminal_FC_raw.csv` | 逐构象的同碳 CH₂ FC（DU4 基组，未标度） |
| `CH2_scaling_worksheet.csv` | **人工/AI 决策表**：每个 CH₂ 的加权 FC、环大小、π 邻接、四类公式各自算出的 J，以及三个空列待填 |
| `coupling_protocol.json` | 方法、基组、四类公式、CH₂ 清单与建议 |
| `dJ_DP4_tables.xlsx` | 原始矩阵、dJ 经验标度页、Bally 参考页 |

`CH2_scaling_worksheet.csv` 的待填列：`chosen_class`、`chosen_J_hz`、`decision_note`。

## 3. DU4 四类公式（SI S4 偕位面板）

| 类别 | 判据 | 公式 | 原文样本 |
|---|---|---|---|
| `sp3_general` | 普通 sp3 CH₂（含 5 元及以上环） | J = 1.240·FC + 0.57 | 23 点，0.13 Hz |
| `sp3_adjacent_sp2` | CH₂ 直接邻接 sp2 碳或氧（烯丙位、α-羰基等） | J = 1.072·FC − 1.66 | 10 点，0.09 Hz |
| `small_ring` | 3 或 4 元环内的 CH₂ | J = 1.249·FC + 1.05 | 25 点，0.14 Hz |
| `sp2` | sp2 碳上的 CH₂（=CH₂） | J = 1.337·FC + 0.52 | 20 点，0.24 Hz |

**判定要点（容易判错的地方）**

1. **CH₂–O 不等于 CH₂–Cπ**。醚/醇氧与羰基氧不同：`CH₂–O` 属 `sp3_general`，
   只有真正邻接 C=C / C=O / 芳环才用 `sp3_adjacent_sp2`。
2. **小环只指 3/4 元环**。5 元环按 `sp3_general`。若 3/4 元环，DU4 原文没有该类拟合，
   只能用 `sp3_general`，**必须在 `decision_note` 标注 extrapolative**。
3. **一个 CH₂ 邻接两个 π 中心**（例如同时接 C=C 和 C=O）属于外推，原文拟合中无此样本。
4. **缩醛/缩酮型 O–CH₂–O** 没有专门类别，按 `sp3_general`，可在 note 里标 "acetal-type"。

## 4. 建议的后处理步骤

1. 打开 `CH2_scaling_worksheet.csv`，逐个 CH₂ 核对 `ring_sizes`、`adjacent_sp2_atoms`。
2. 按第 3 节判据选定 `chosen_class`；`chosen_J_hz` 取对应 `J_<class>_hz` 列的值
   （加权 FC 已按构象自由能权重算好）。
3. 在 `decision_note` 写一句理由，例如
   `"C15: CH2-O, not pi-adjacent -> sp3_general"`。
4. 常规 H–H 对直接用 `0.9155 × FC`；**两套结果分别列、分别出表，不合并、不平均**。

## 5. 必须随结果一起报告的限制

- FC-only 经验拟合，不是从头算 J 的全部四项。
- 几何与权重是气相自定义条件（B3LYP/6-31G(d) 几何 + ORCA/Shermo 自由基布居），
  与原文 MMFF/电子能权重不同，属 exploratory。
- 气相拟合未包含 CDCl₃ 溶剂对 J 的影响。
- 基组覆盖：DU4 的 Si 用 3-21G\*（SI 明确）；Bally 原文未定义非 C/H/O 元素，
  按作者脚本一律 6-31G(d,p)，导出时会记录这一点。
- 若分子含其他元素（P、S、卤素），DU4 有定义（3-21G\*），Bally 仍按 6-31G(d,p)。

## 6. 使用标度法的需求（先看这一节）

标度法只是把**对应基组算出的 FC** 换算到经验 J 尺度，它不能补救错的输入。下面几条是硬需求：

| 需求 | 说明 | 不满足会怎样 |
|---|---|---|
| 基组与系数成对 | `0.9155` 只配 Bally 8s1p 的 FC；DU4 四类式只配 DU4 的 FC | 拿普通 6-31G(d,p) 的 FC 代进去，数值系统性偏差，且偏差随分子变 |
| 保留符号 | 一律对带符号的 FC 做乘法/仿射变换，不要先取绝对值 | W 型 ⁴J、端基 ²J 这类负值或近零值会被符号算错 |
| 同碳 CH₂ 逐个选类 | 按第 3 节判据人工选 `chosen_class` 并写 `decision_note` | 选错类，²JHH 可差 2–3 Hz |
| 两轨分别出表 | 常规 H–H 与同碳 CH₂ 各自成表，不合并、不取平均 | 破坏可追溯性，也无法判断哪一轨偏离实验 |
| 需要实验 J 才能评估准确度 | 没有实验值时只能报告预测值，不能声称"验证" | 容易被当成已经对过实验 |
| 做 dJ-DP4 需四件套 | `experiment.csv`、`assignments.csv`、`experiment_j.csv`、`assignments_j.csv`（见 README 第 5 节） | 流程不会给出相对概率 |
| 溶剂/温度变化需重新标度 | 两套系数都来自气相/特定条件，未含 CDCl₃ 溶剂与振动平均修正 | 跨条件直接套用会把系统误差算进结论 |

简言之：**标度法解决"同一基组下 FC→J 的系统偏差"，不解决几何、溶剂、振动平均和作业本身的误差。**
