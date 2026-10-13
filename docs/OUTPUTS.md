> 默认 NMR 基础参考；V2.1 全模块操作以 [README](../README.md) 和 [简明 SOP](QUICKSTART_V2.md) 为准。

# 输出路径、文件命名与日志

适用于V1.1起由00创建的新项目。已运行项目不能通过修改output_layout强行转换。

```text
项目/
  过程文件/
    001RR/
      sampling/resp_gas/Gau_RESP_gas_001RR.out
      sampling/resp_solvent/Gau_RESP_solvent_001RR.fchk
      conformers/c0001/
        opt1/Gau_opt_001RR_C0001_A01.gjf
        opt1/Gau_opt_001RR_C0001_A01.log
        opt1/Gau_opt_001RR_C0001_A01.chk
        freq/Gau_freq_001RR_C0001.log
        sp/Orca_sp_001RR_C0001.inp
        sp/Orca_sp_001RR_C0001.log
        sp/Orca_sp_001RR_C0001.gbw
        nmr01/Gau_nmr01_001RR_C0001.log
      shermo/Shermo_weight_001RR_C0001.log
      conformers.json
      state.json
    exports/
  运行日志/run_时间_作业号.log
  文件索引.csv
  构象布居与混合比例.csv
```

Gau表示Gaussian，Orca表示ORCA，opt/freq/sp分别为优化/频率/单点。C0001是候选内部构象编号，A01/A02/A03是优化尝试编号。nmr01为B3LYP→DP4+，nmr02为mPW1PW91→DP4+，nmr03为wB97XD→STS。输入、输出和检查点采用相同主文件名，ORCA生成的gbw等附属文件随输入名生成。RESP在构象搜索前进行，只有候选编号，没有构象编号。

软件固定生成的内部文件名（如CREST的crest_ensemble.xyz）、机器状态文件state.json/receipt.json以及MD原生文件名保留，避免破坏软件调用和断点恢复。它们保存在各候选独立目录，通过完整路径和文件索引辨认。优化失败尝试的目录和输出保留，不覆盖。

总日志保留MD搜索、Sobtop、CREST、聚类及其他工具的现有输出，也保留Gaussian/ORCA的起止、路径和返回值；其冗长计算详情只写入独立文件。软件失败时先读总日志定位路径，再打开对应独立日志。独立日志在运行中会增长，04不替换它们。

新项目无需额外scratch原始目录或查看软链接；大检查点和轨迹也是真实文件。00创建时固定output=过程文件，之后不要手工修改。项目根目录必须是计算节点可读写的共享存储，迁移时备份整个项目。04更新构型/权重/索引，并复制最终结果汇总；不删除原始数据。

历史正在运行项目继续使用原固定代码。只有作业完全停止、验收和比较结束后，才能另外生成真实文件整理副本和精简日志；原始数据与原总日志必须保留以便追溯。
