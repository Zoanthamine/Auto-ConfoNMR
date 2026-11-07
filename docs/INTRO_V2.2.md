# Auto ConfoNMR V2.2 发布介绍

## 一句话

给一个三维候选结构，自动跑完构象搜索 → DFT 精修 → 自由能加权 → 多套 NMR 屏蔽 → 与实验比对，直接得到 DP4+/STS 概率与可交付的耦合常数表。

## 它解决什么

算一个分子的计算核磁，过去是一串手工活：手挑构象、逐个写 Gaussian 输入、拼权重、换算 TMS、再往 DP4+/STS 的 Excel 里一格格粘。任何一步手抖，结论就得重来，而且很难复现。Auto ConfoNMR 把这条链路固化成一个项目目录和四个日常脚本，每一步都有收据、快照和日志。

## 现在能做什么

| 能力 | 具体内容 |
|---|---|
| 构象与布居 | GFN1/GFN2 搜索 + isostat 去重，可切 `fine`/`gfn_coarse` 等聚类档；快速档用电子能、严谨档用 ORCA+Shermo 准 RRHO 自由能做权重，给出覆盖率 |
| 三套 NMR 屏蔽 | 01 B3LYP/6-31G(d)、02 mPW1PW91/6-311+G(d,p)、03 wB97X-D/6-31G(d)，IEFPCM + G09Defaults，逐构象计算后按自由能加权 |
| DP4+ / STS | 自动把加权结果填进上游原版工作簿（公式与 TMS 常数保留），比较多个候选给出概率；STS 走碳类型分段标度，不需要实验值也能出预测 |
| 耦合常数（V2.2 双轨） | 常规 H–H 用真实 Bally–Rablen 协议（8s1p，FC-only，×0.9155）；CH₂ 同碳 ²J 用 DU4 基组 + 四类环境式，两轨分别出表、不平均 |
| 多候选与多样品 | 一次建项放多个候选结构；实验表按 `experiment_<样品名>.csv` 命名即可同项目对比多组数据 |
| ECD | 三套 TDDFT 跃迁、玻尔兹曼加权曲线、PNG/SVG，可叠加实验曲线 |
| 工程化 | 项目内冻结 workflow 快照（老项目不被新版本改写）、提交与队列双重防重、总流程日志、`--check-only` 预检、跨平台一致的收据哈希 |

## V2.2 新增

- **耦合双轨**：一般项换真实 Bally 基组（脱离 6-31G(d,p) 近似）；新增 DU4 同碳分支与四类环境式，选类理由留痕。
- **含碘兼容**：I 用 def2-SVP/ECP 参与几何与频率、全电子 DGDZVP 参与旧 NMR 分支；只做基组兼容，不含自旋–轨道效应，自动统计暂停。
- **标度法需求成文**：`docs/COUPLING_SCALING_GUIDE.md` 第 6 节把"基组与系数成对、保留符号、逐个选类、两轨不平均、评估准确度需实验 J、跨溶剂/温度重标"写成硬需求。
- 轨迹按 GRO 原子名归一化 VMD XYZ 元素列；Gaussian 输入统一 UTF-8/LF；新增总流程日志。

## 边界（写在前面）

- 面向 Linux 集群、当前提供 PBS 单节点提交；Gaussian、ORCA、GROMACS、VMD 等需自备。
- 耦合与位移的标度都是 FC/GIAO 经验拟合，不是四项全算；气相几何加溶剂模型不等于溶液中的真实振动平均。
- DP4+/STS 的"概率"只在给定候选集合内比较，不是绝对正确率。
- 已在多种子、多候选、多样品与真实有机分子上跑通回归（105 项测试通过）；作者自己的案例里，16 核、单构象的小分子全流程约 20 分钟（不含排队），复杂分子按构象数与基组放大。

## 三步上手

```bash
git clone https://github.com/Zoanthamine/Auto-ConfoNMR.git Auto_ConfoNMR
cd Auto_ConfoNMR && cp config/cluster.example.json config/cluster.json   # 填 CHANGE_ME
bash 00_new_project.sh      # 建项：选基础模块、填溶剂/温度/核数、放候选结构
bash 01_submit.sh           # 提交队列
bash 03_compare.sh          # 填好实验表后出 DP4+/STS
```

## 链接

- 仓库：https://github.com/Zoanthamine/Auto-ConfoNMR
- 下载与校验：https://github.com/Zoanthamine/Auto-ConfoNMR/releases/tag/V2.2
- 许可：MIT（第三方二进制与工作簿保留各自许可，见 `LICENSE_NOTICE.md`）