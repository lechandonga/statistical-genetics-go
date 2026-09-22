# statistical-genetics-go（sggo）

一套**纯本地文件**的遗传统计分析流程：用 TSV 文件模拟基因型矩阵、
连锁不平衡（LD）分数与多队列全基因组汇总统计（GWAS summary stats），
在**等位基因方向不一致、数据缺失、样本重叠**等非理想条件下，给出
可解释、可复现的结果。

纯 Go 标准库实现（无第三方数值依赖），数值层通过 `numwrap` 对后端
版本与错误类型做归一化封装。

## 能力一览

- **统一变体模型**：汇总统计 / LD 分数 / 基因型矩阵解析后归一到
  `variant.Variant`；字段缺失、非法数值、重复变体、无法识别的等位
  基因都带**可区分原因代码**拒绝，绝不静默跳过（严格模式默认）。
- **等位基因方向统一**：自动识别直接匹配、等位基因互换（效应取反）、
  链翻转（符号保留）及其组合；A/T、C/G 等链方向模糊位点默认拒绝，
  也可显式标记保留。
- **遗传力估计**：LD 分数回归（LDSC 风格），输出点估计、SE、95% CI；
  对负值/越界/低功效区间有 report / clamp / reject 三种明确策略。
- **多队列元分析**：固定效应 IVW、DerSimonian–Laird 随机效应、
  Cochran Q / I² 异质性；基于样本 ID 硬证据检测重叠与完全重复队列，
  支持 flag / drop / reject，不把重叠当作独立证据叠加。
- **个体风险评分（PRS）**：正确处理缺失基因型（期望剂量插补/排除）
  与效应位点缺失（拒绝/交集口径显式化），固定位点集合、符号约定与
  评分尺度的审计记录。
- **并发可复现**：确定性分片 + 固定顺序 Kahan 归约，任意并行度结果
  与串行**逐位一致**；`-race` 安全；失败不返回半更新结果。
- **可配置上限**：变体数、样本数、行/文件字节、墙钟耗时；超限可拒绝
  或截断，不残留部分统计状态。

## 包结构

```
variant/           统一变体模型、等位基因对齐、拒绝原因代码
ingest/            TSV 解析（汇总统计/LD 分数/基因型矩阵）、上限与报告
numwrap/           数值后端归一化（stdlib；gonum 预留为显式不支持）
stats/             h2.go 遗传力 / meta.go 元分析 / prs.go 风险评分
                   parallel.go 确定性分片并行
internal/simulate/ 固定种子的确定性测试数据生成器
pipeline/          端到端编排
cmd/sggo/          命令行入口
docs/              model.md（模型/对齐规则）、compatibility.md、verification.md
```

## 快速开始

```sh
go build -o sggo ./cmd/sggo
./sggo simulate -out testdata -variants 500 -samples 300 -overlap 80
./sggo h2   -sum testdata/cohort1_sumstats.tsv -ld testdata/ldscores.tsv -n 280
./sggo meta -in c1=testdata/cohort1_sumstats.tsv -in c2=testdata/cohort2_sumstats.tsv
./sggo prs  -sum testdata/cohort1_sumstats.tsv -gt testdata/genotypes.tsv \
            -n 280 -missing-sites ignore -scale std -out testdata/prs.tsv
```

## 测试

```sh
go vet ./...
go test ./... -race -count=1
```

测试以边界与异常数据为主（链翻转、模糊位点、负遗传力、队列重叠、
缺失字段、极端数值），日志打印输入摘要与判定依据。详见
[docs/verification.md](docs/verification.md)。

## 文档

- [docs/model.md](docs/model.md)：统计模型假设、变量定义、等位基因
  对齐规则、负值/重叠/缺失的处理语义；
- [docs/compatibility.md](docs/compatibility.md)：输入格式、列名别名、
  数值后端版本与错误归一化、Go 版本兼容范围；
- [docs/verification.md](docs/verification.md)：验证方法与确定性证据。

## 范围与限制（明确不做）

- 不联网、不读取压缩文件（gzip）/Parquet/VCF；仅本地未压缩 TSV。
- 不实现 LD 分数本身的估计（r² 矩阵计算），LD 分数为输入。
- 随机效应仅 DL；REML、协变量调整的 PRS 评估、多表型留待后续，
  当前请求这些选项会被显式拒绝。
- 异质性 P 值用 Wilson–Hilferty 近似（文档中已标注，非精确卡方）。
