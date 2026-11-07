# Auto ConfoNMR V2.2 工具说明

## 1. 定位

Auto ConfoNMR 是一个面向有机小分子的构象加权计算核磁流程。输入为候选三维结构（可多个）与该结构对应的分子式、电荷与多重度；输出为各构象的坐标与自由能权重、三套量子化学方法下的 ¹H/¹³C 屏蔽常数及其加权值、DP4+/STS 比较结果、以及耦合常数矩阵。

流程本身不实现量子化学方法，也不重新拟合任何经验标度；它负责的是调用顺序、构象去重与加权、输入输出约定、以及结果的可复现性。

## 2. 所依据的方法

| 环节 | 采用的方法 | 来源 |
|---|---|---|
| 构象搜索 | GFN1-xTB / GFN2-xTB + isostat 能量窗去重 | Grimme 等，xTB / CREST |
| 几何优化与频率 | B3LYP-D3(BJ)/6-31G(d)（气相） | — |
| 构象权重 | ORCA 单点 + Shermo 准 RRHO 自由能；或 GFN2 电子能 | Grimme / Li 等 |
| NMR 屏蔽 | 01 B3LYP/6-31G(d)；02 mPW1PW91/6-311+G(d,p)；03 wB97X-D/6-31G(d)，均 IEFPCM + G09Defaults | — |
| 概率比较 | DP4+（经典/scaled 分支）；STS 碳类型分段标度 | Grimblat–Zanardi–Sarotti；Li–Liu–Wang, *J. Org. Chem.* 2020, 85, 11350（[10.1021/acs.joc.0c01451](https://doi.org/10.1021/acs.joc.0c01451)） |
| 一般 H–H 耦合 | Bally–Rablen 协议：B3LYP/Gen，H 为 u+1s（8s1p），非氢元素 6-31G(d,p)，FC-only，`Int=NoBasisTransform`，整体乘有符号 0.9155 | Bally–Rablen, *J. Org. Chem.* 2011（[10.1021/jo200513q](https://doi.org/10.1021/jo200513q)）SI Table S4a |
| 同碳 CH₂ 耦合 | DU4 基组（C 4-31G；O 6-31G；H 8s1p）+ SI S4 四类环境式，人工选类并记录 | Kutateladze–Mukhina, *J. Org. Chem.* 2014（[10.1021/jo501781b](https://doi.org/10.1021/jo501781b)）SI S2–S4 |
| ECD | 三套 TDDFT 跃迁，玻尔兹曼加权 | — |

两条耦合分支的标度系数与基组严格对应，不交叉使用，结果分别出表、不取平均。

## 3. 计算流程

建项时写入溶剂、温度、压力、核数与资源上限，并冻结一份流程快照；此后流程升级不影响已建项目。执行顺序为：结构准备与去重 → 构象搜索 → 几何优化与频率 → 构象聚类 → 自由能加权 → 单点能 → 屏蔽与耦合计算 → 后处理导出。每一步产出结构化文件与收据（输入哈希、版本、参数），失败可定位到具体阶段。

资源调度当前提供 PBS 单节点提交；提交前执行 `--check-only` 预检，校验基组覆盖、泛函与溶剂组合、收据完整性与队列状态。

## 4. 数据与可复现性

- 每个项目独立目录，过程文件、导出文件与日志分层存放。
- Gaussian 输入统一以 UTF-8/LF 写入，收据哈希在 Windows 与 Linux 上一致。
- 导出的工作簿沿用上游 DP4+/STS 原版模板的公式与 TMS 常数，加权结果填入指定单元格，不重算公式。
- 后处理记录每个构象的权重、每个 CH₂ 的选类理由、每个耦合分支所用的标度式。

## 5. 验证

- 回归测试：Linux + 项目 venv 下 105 项通过、2 项跳过（含碘兼容、VMD 轨迹元素归一化、流程日志等新增用例）。
- 实际案例：C₁₁H₁₀O₃（四个手性中心，环戊烷/环丙烷/γ-内酯/环戊酮并环），GFN2 聚类得到单一构象（350 K 布居 100%，覆盖率 100%），三套屏蔽与两轨耦合均正常导出；16 核单节点全流程约 20 分钟（不含排队时间）。

## 6. 适用条件与限制

1. 面向 Linux 集群，当前提供 PBS 单节点提交方式；Gaussian、ORCA、GROMACS、VMD、xtb、crest、Shermo 需自行配置。
2. 位移与耦合均基于 GIAO 屏蔽与 FC-only 经验标度，不是 Ramsey 四项全算；标度系数来自特定基组与训练集，不可跨基组代用。
3. 几何为气相优化，溶剂效应以 IEFPCM/G09Defaults 单点近似，未含显式溶剂与振动平均修正。
4. DP4+/STS 给出的概率是在给定候选集合内的相对比较，不是绝对正确率；缺少实验归属时只能报告预测值。
5. 含碘体系仅做基组兼容（def2-SVP/ECP 参与几何与频率，全电子 DGDZVP 参与旧 NMR 分支），不含自旋–轨道效应。
6. 标度法的具体使用要求见 `docs/COUPLING_SCALING_GUIDE.md` 第 6 节。

## 7. 使用与许可

```bash
git clone https://github.com/Zoanthamine/Auto-ConfoNMR.git Auto_ConfoNMR
cd Auto_ConfoNMR && cp config/cluster.example.json config/cluster.json
bash 00_new_project.sh    # 建项与结构准备
bash 01_submit.sh         # 提交队列
bash 03_compare.sh        # 实验表就位后生成 DP4+/STS 比较
```

代码以 MIT 许可发布；随包的第三方二进制与工作簿保留各自许可，见 `LICENSE_NOTICE.md`。发行包见 [releases](https://github.com/Zoanthamine/Auto-ConfoNMR/releases/tag/V2.2)。