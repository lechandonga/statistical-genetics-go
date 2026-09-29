# statistical-genetics-go（sggo）

一套**纯本地文件**的遗传统计分析流程：用 TSV 文件模拟基因型矩阵、
连锁不平衡（LD）分数与多队列全基因组汇总统计（GWAS summary stats），
在**等位基因方向不一致、数据缺失、样本重叠**等非理想条件下，给出
可解释、可复现的结果。

数值层 `numwrap` 提供两个后端：自包含的 Go 标准库实现（默认）与
版本锁定并 vendor 的成熟线代库 `gonum.org/v1/gonum v0.17.0`；
可按每次调用选择后端，或用 `crosscheck` 在一次调用内跑两边交叉验证
（数值容差内一致、判定逐字段一致，对不上明确报错）。断网可复现构建。

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
- **双后端交叉验证**：`stdlib` 与 vendored gonum；同后端逐位一致，
  两后端数值在 `rel=1e-9/abs=1e-12` 内一致，拒绝/标记/纳入等判定必须
  完全一致；不一致即 `backend_mismatch` 报错，不挑结果。
- **可配置输入预算**：汇总统计/LD 分数/基因型三类输入分别配置变体数、
  样本数、行/文件字节、墙钟耗时；超限拒绝或在完整行边界截断，
  不残留半截结果；每次运行输出字节/行/样本/耗时用量小结，并发互不串。
- **混合来源元分析**：只有 Z+N（+EAF）的队列按统一口径
  `SE=1/sqrt(2Np(1-p))`、`beta=Z*SE` 真正参与合并；缺 N、无 EAF
  （或 EAF=0/1）分别以可区分原因排除；逐队列权重、方向依据与排除原因
  全部进入审计输出。

## 包结构

```
variant/           统一变体模型、等位基因对齐、拒绝原因代码
ingest/            TSV 解析（汇总统计/LD 分数/基因型矩阵）、上限与报告
numwrap/           数值后端归一化（stdlib + vendored gonum + crosscheck）
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
# 双后端交叉验证 / Z-only 队列：
# ./sggo h2 ... -backend crosscheck
# ./sggo meta -in c1=beta.tsv -in c2=z.tsv -z c2 -audit-out audit.tsv
./sggo prs  -sum testdata/cohort1_sumstats.tsv -gt testdata/genotypes.tsv \
            -n 280 -missing-sites ignore -scale std -out testdata/prs.tsv
```

## 测试

```sh
GOFLAGS=-mod=vendor GOPROXY=off go vet ./...
GOFLAGS=-mod=vendor GOPROXY=off go test ./... -race -count=1
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

- 运行期不联网、不读取压缩文件（gzip）/Parquet/VCF；仅本地未压缩 TSV。
  gonum 源码随仓库 vendor，断网可构建/测试。
- 不实现 LD 分数本身的估计（r² 矩阵计算），LD 分数为输入。
- 随机效应仅 DL；REML、协变量调整的 PRS 评估、多表型留待后续，
  当前请求这些选项会被显式拒绝。
- 异质性 P 值用 Wilson–Hilferty 近似（文档中已标注，非精确卡方）。
