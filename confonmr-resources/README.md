# Linux 辅助资源包

本目录直接包含在 Git 仓库和 Auto ConfoNMR V1.2 运行包中。文件是实际程序和数据，不是软链接。适用于 **Linux x86-64**；不适用于 Windows、macOS 或 ARM。二进制来自已验证环境，未升级为官网当前最新版。

## 包含什么

|目录或文件|版本与用途|作者来源|
|---|---|---|
|`bin/crest`|CREST 2.12；GFN 构象精修|[官方 V2.12](https://github.com/crest-lab/crest/releases/tag/v2.12)|
|`bin/xtb`、`share/xtb/`|xTB 6.6.1；可执行文件及配套参数|[官方 V6.6.1](https://github.com/grimme-lab/xtb/releases/tag/v6.6.1)|
|`bin/Shermo`|Shermo 2.4.1；热化学与构象分布|[Tian Lu，Shermo](http://sobereva.com/soft/shermo/)|
|`multiwfn/`|Multiwfn 3.8(dev)，2023-12-08；无图形程序、设置及许可，用于电荷拟合|[Tian Lu，Multiwfn](http://sobereva.com/multiwfn/)|
|`isostat`、`mol2_xyz`、`xyz2QC`|Molclus 工具；聚类和格式转换|[Tian Lu，Molclus](http://www.keinsci.com/research/molclus.html)|
|`optional/molclus/molclus`|Molclus 主程序；供有需要者另行使用，本流程不调用它|[Molclus 使用说明](http://www.keinsci.com/research/molclus.html)|
|`sobtop`、`sobtop.ini`、`atomtype`、7 份 `.dat/.DEF/.DAT` 文件|Sobtop 1.0(dev3.1) 与现有配套原子类型、力场参数；拓扑生成|[Tian Lu，Sobtop](http://sobereva.com/soft/Sobtop/)|
|`upstream/DP4plus.xlsx`|DP4+ 工作簿接口；由上游 custom-dp4plus.xlsx 通过提取工具处理得到|[固定来源文件](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/blob/095dbc16db04b40371edfaa61a03aa581e5c7d93/DP4plus/tool/custom-dp4plus.xlsx)|
|`upstream/STS.xlsx`|STS 原始接口工作簿|[固定来源文件](https://github.com/chemicalchange/EasyCalcNMR_for_OrgSynLab/blob/095dbc16db04b40371edfaa61a03aa581e5c7d93/STS/tool/STS.xlsx)|
|`licenses/`、`notices/`、`sources/`|原始许可、引用要求、CREST/xTB 对应版本源码和依赖参考源码|各包原作者；归属保持不变|

Sobtop 二进制实际标识为 dev3.1，不能根据原安装目录名把它写成 dev3.2。参数文件保留现有流程版本；`sobtop.ini` 仅去除了站点路径，建项目时由配置重新填写。Multiwfn 的 `lastfile` 历史记录已清空，其他计算设置保留。程序版本、文件哈希见 `SOURCES.json` 和 `SHA256SUMS.txt`。

## 怎么接入

克隆仓库或解压运行包后，本目录已经位于程序根目录，无需另外下载或解压资源。可先浏览 `bin/`、`multiwfn/` 和 `upstream/`，确认程序和工作簿。目录结构为：

```text
Auto_ConfoNMR/
├── README.md
├── cluster_config.py
├── config/cluster.json
└── confonmr-resources/
    ├── bin/crest、xtb、Shermo
    ├── multiwfn/Multiwfn_noGUI、settings.ini、LICENSE.txt
    ├── isostat、sobtop、sobtop.ini、atomtype、mol2_xyz、xyz2QC
    ├── ATOMTYPE_AMBER.DEF、ATOMTYPE_GFF.DEF、CORR_NAME_TYPE.DAT
    ├── LJ_param.dat、assign_AT.dat、bondcrit.dat、bonded_param.dat
    ├── upstream/DP4plus.xlsx、STS.xlsx
    ├── share/xtb/
    └── licenses/、notices/、sources/
```

保留目录层级。在 `config/cluster.json` 中填入本机实际绝对路径：

```json
"resource_bundle": "/path/to/Auto_ConfoNMR/confonmr-resources"
```

`software.crest` 填 `<资源目录>/bin/crest`，`software.shermo` 填 `<资源目录>/bin/Shermo`，`software.multiwfn` 填 `<资源目录>/multiwfn/Multiwfn_noGUI`；将 `<资源目录>/bin` 加入 `path_prepend`。项目环境会自动设置本目录的 xTB 参数搜索路径；Sobtop 的两个外部程序路径会从 `software.multiwfn` 和 `software.obabel` 写入项目快照，不必逐项手改。

```bash
cd /path/to/Auto_ConfoNMR/confonmr-resources
sha256sum -c SHA256SUMS.txt
chmod +x bin/crest bin/xtb bin/Shermo multiwfn/Multiwfn_noGUI \
  isostat sobtop atomtype mol2_xyz xyz2QC optional/molclus/molclus
bin/crest --version
bin/xtb --version
```

CREST/xTB 为官方静态编译版，已与官方发行二进制核对 SHA256。其他程序依赖 Linux 系统运行库，请在实际计算节点用 `ldd` 检查；出现 `not found` 时先补充运行库或改用目标系统重新安装的兼容版本。不要从别处复制整套 ORCA 的库来凑依赖。

## 仍需自行准备

本包不是完整计算环境。**Gaussian/formchk、ORCA 及匹配的 OpenMPI、GROMACS、VMD、Open Babel、Python 与 requirements.txt 中的依赖、PBS** 仍由使用者或集群管理员配置。GROMACS、Open Babel 虽可免费使用，也需要其库和数据目录；VMD 需要完整安装，不能只拿一个执行文件。它们的安装位置与 module 名称写在 `config/cluster.json`。

Multiwfn 这里只提供 V1.2 电荷拟合所用的无图形程序和设置，并非带所有示例、手册和附加原子密度库的完整发行包。Molclus 主程序为附送工具，其独立使用所需的 settings.ini、计算模板和说明请按作者教程准备，不影响本流程调用 isostat。

## 许可与引用

资源保持原作者权利，不随 Auto ConfoNMR 重新授权。CREST/xTB 保留上游 GPL/LGPL 许可文本和对应版本源码；源码中的构建说明决定重编译所需的编译器、BLAS/LAPACK 及依赖。`sources/` 保留原始源码压缩档，供许可核查和重新编译使用，日常运行无需解压；实际辅助程序与参数均已展开存放。`sources/xtb-dependencies/` 保存上游 wrap 中固定版本的参考源包，不能据此宣称复现了官方二进制的全部构建环境。

Shermo 和 Multiwfn 的原始许可允许分发，并要求引用原论文。Molclus 附带引用要求；Sobtop 官网说明其学术和商业使用免费，使用时应正确引用。后两者与工作簿没有在这里被改为开源软件，也没有额外授予再许可。原作者更新版本或许可时，以对应版本的原始说明为准。

- CREST：Pracht, P.; Bohle, F.; Grimme, S. *Phys. Chem. Chem. Phys.* **2020**, 22, 7169–7192. [DOI](https://doi.org/10.1039/C9CP06869D)。
- GFN1-xTB：Grimme, S.; Bannwarth, C.; Shushkov, P. *J. Chem. Theory Comput.* **2017**, 13, 1989–2009. [DOI](https://doi.org/10.1021/acs.jctc.7b00118)。
- GFN2-xTB：Bannwarth, C.; Ehlert, S.; Grimme, S. *J. Chem. Theory Comput.* **2019**, 15, 1652–1671. [DOI](https://doi.org/10.1021/acs.jctc.8b01176)。
- Shermo：Lu, T.; Chen, Q. *Comput. Theor. Chem.* **2021**, 1200, 113249. [DOI](https://doi.org/10.1016/j.comptc.2021.113249)。
- Multiwfn：Lu, T.; Chen, F. *J. Comput. Chem.* **2012**, 33, 580–592. [DOI](https://doi.org/10.1002/jcc.22885)。
- Molclus、Sobtop：引用 Tian Lu 及上述官方主页，并记录实际版本和访问日期。
- DP4+/STS：见主 [README](../README.md#参考与许可) 的原始论文及 EasyCalcNMR 来源记录。
