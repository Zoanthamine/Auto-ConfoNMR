# 集群配置文件填写说明

首次配置请读 [README 的迁移部分](../README.md#迁移配置)，按目标集群填写外部软件路径、module和队列。本页用于快速查字段；完整module示例和独立安装方式见迁移说明。

复制 cluster.example.json 为 cluster.json 后编辑真实值。cluster.json 不进入公共发布包。JSON 必须使用双引号，不加中文注释，不留尾逗号。完整安装步骤见 [迁移配置](../docs/MIGRATION.md)。

路径设置、软件设置和队列设置属于安装环境；泛函、溶剂、聚类和布居属于项目科学配置，分别保存，不要混用。

|配置|填写示例或解释|
|---|---|
|python|`/home/user/envs/confonmr/bin/python`|
|projects_root|`/home/user/projects`，不要放在 ConfoNMR 源码内部|
|scratch_root|兼容历史配置字段；新项目输出在项目/过程文件；请为projects_root选择足够容量的共享存储|
|resource_bundle|`/home/user/confonmr-resources`，含13个辅助资源及upstream的两个表|
|software.g16/formchk/orca/crest/shermo/multiwfn/obabel/vmd/gmx|实际可执行文件；建议绝对路径|
|software.openmpi|MPI安装根目录，含bin和lib，必须匹配ORCA构建|
|modules|例如 `["gcc/your_version","gromacs/your_version"]`；示例名称不能直接照抄|
|gaussian_profile|真实g16.profile；若由module完整初始化，可填空字符串并实测|
|path_prepend|必要bin目录数组；没有填 `[]`|
|ld_library_path_prepend|必要动态库目录数组；没有填 `[]`|
|scheduler.type|`pbs`|
|scheduler.queue|自己有权限的队列|
|scheduler.walltime|队列允许的最长时间|
|scheduler.node_resource|Torque：`nodes=1:ppn={cores}`；PBSPro：`select=1:ncpus={cores}:mem={memory_gb}gb`|
|scheduler.extra_qsub_args|如 `["-A","account"]`，没有填 `[]`|

配置完成先运行 `python cluster_config.py check`。该检查只验证格式和资源文件完整性；软件版本/许可证/批节点环境需要进一步验证。
