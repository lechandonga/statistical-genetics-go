# 验证方法

## 一键运行

```sh
go vet ./...
go test ./...            # 全部单测 + 端到端测试
go test ./... -race      # 竞态检测（并行分片路径）
go test -count=2 ./...   # 重复运行，确认无不确定行为
# 断网（依赖已 vendor、版本已锁定）
GOPROXY=off GOFLAGS=-mod=vendor go test ./... -race -count=1
```

所有测试都遵循“边界与异常数据优先”，并在日志中打印输入摘要与判定
依据（`[input] ... [judgement] ...`）。

## 覆盖的关键场景

| 能力 | 测试 | 关注点 |
| --- | --- | --- |
| 归一化/拒绝原因 | `variant/model_test.go` | 空等位基因、I/D 码、同等位基因、EAF 越界、反链、chrM |
| 对齐 | `variant/align_test.go` | 直接/互换/翻转/翻转+互换、A/T 与 C/G 模糊、禁止翻转、不匹配 |
| 汇总统计解析 | `ingest/summary_test.go` | 列别名、OR→logOR、Z-only、缺 N、重复变体、严格/宽松、截断 |
| LD/基因型解析 | `ingest/genotype_test.go` | 负 L2、剂量越界、ACGT 解码、缺失标记、样本上限、重复行 |
| 数值封装 | `numwrap/regression_test.go` | WLS 精确解、固定截距、NaN/奇异/负权重、不支持后端、分位数 |
| 双后端 | `numwrap/backend_test.go` | WLS/IVW/DL/矩/分位数两后端一致、同后端逐位、并发切换竞态、未知后端拒绝、gonum 错误翻译、版本三处锁定 |
| 遗传力 | `stats/h2_test.go` | 真值恢复、负 h2 三策略、越界、模糊剔除、workers=1..100 逐位一致 |
| 元分析 | `stats/meta_test.go`、`stats/backends_test.go` | IVW、Q/I²、DL、对齐、模糊、重叠、REML 拒绝；双后端一致、混合 Z-only 纳入/排除/审计、并发切换竞态 |
| PRS | `stats/prs_test.go`、`stats/backends_test.go` | 方向、缺失两策略、三尺度、阈值、评估回归、双后端一致 |
| 输入预算 | `ingest/limits_test.go` | 用量小结、极小 MaxFileBytes/MaxLineBytes/MaxElapsed 生效、拒绝无部分结果、变体/样本截断、并发用量隔离 |
| 端到端 | `pipeline/pipeline_test.go`、`pipeline/backends_test.go` | 注入行拒绝、重叠标记/拒绝、串行/4 worker 一致、CrossCheck 双后端通过、未知后端拒绝 |

## 确定性（可复现性）验证要点

1. `TestEstimateH2ParallelEqualsSerial`：同一输入在 workers ∈
   {1,2,3,4,8,16,100} 下 h2/SE/截距**逐位相等**。
2. `TestMetaRandomDLAndParallel`、`TestPRSScalesAndParallel`：
   元分析与标准化评分在多并行度下逐位相等。
3. `internal/simulate` 使用固定种子 xorshift64*（非全局 `rand`），
   两次 `sggo simulate -seed 42` 产出字节一致的数据文件。
4. 所有浮点归约采用固定顺序 Kahan 求和，累加顺序不随调度变化。

## CLI 手工验证

```sh
go build -o /tmp/sggo ./cmd/sggo
/tmp/sggo simulate -out /tmp/demo -variants 500 -samples 300 -overlap 80
/tmp/sggo h2   -sum /tmp/demo/cohort1_sumstats.tsv -ld /tmp/demo/ldscores.tsv -n 280
/tmp/sggo meta -in c1=/tmp/demo/cohort1_sumstats.tsv -in c2=/tmp/demo/cohort2_sumstats.tsv -model random
/tmp/sggo prs  -sum /tmp/demo/cohort1_sumstats.tsv -gt /tmp/demo/genotypes.tsv \
               -n 280 -missing-sites ignore -scale std -out /tmp/demo/prs.tsv
```

带重叠硬证据的元分析（为队列提供每行一个样本 ID 的文件）：

```sh
/tmp/sggo meta -in c1=sum1.tsv:ids1.txt -in c2=sum2.tsv:ids2.txt -overlap reject
```

`sggo version` 打印数值后端及其版本/接口描述。

## 手工核对两个后端一致

```sh
go build -o /tmp/sggo ./cmd/sggo
/tmp/sggo version            # 列出实际编译后端与锁定版本、容差
# 任一分析阶段加 -crosscheck：两后端各算一遍，不一致即非零退出
/tmp/sggo h2   -sum /tmp/demo/cohort1_sumstats.tsv \
              -ld  /tmp/demo/ldscores.tsv -n 280 -crosscheck
/tmp/sggo meta -in c1=/tmp/demo/cohort1_sumstats.tsv \
              -in c2=/tmp/demo/cohort2_sumstats.tsv -crosscheck -out /tmp/m.tsv
/tmp/sggo prs  -sum /tmp/demo/cohort1_sumstats.tsv \
              -gt  /tmp/demo/genotypes.tsv -n 280 -scale std \
              -missing-sites ignore -out /tmp/p.tsv -crosscheck
# 单独指定后端（未编译/拼错的名字会直接报错，不会回退）
/tmp/sggo meta -in c1=a.tsv -in c2=b.tsv -backend gonum
/tmp/sggo meta -in c1=a.tsv -in c2=b.tsv -backend openblas   # 显式失败
```

## 手工验证各类上限真的生效

```sh
# 文件字节：设得远小于文件 => 整体拒绝，无输出文件/部分结果
/tmp/sggo h2 -sum a.tsv -ld b.tsv -n 280 -max-file-bytes 200
# 行数：截断策略保留前 3 个变体，报告 truncated 与命中的上限
/tmp/sggo meta -in c1=a.tsv -in c2=b.tsv \
              -max-variants 3 -limit-policy truncate -out /tmp/m3.tsv
# 样本列：PRS 基因型只保留前 5 个样本
/tmp/sggo prs -sum a.tsv -gt g.tsv -n 280 -out /tmp/p5.tsv \
              -max-samples 5 -limit-policy truncate
# 墙钟预算：1ns 立即拒绝（无部分结果）
/tmp/sggo h2 -sum a.tsv -ld b.tsv -n 280 -max-elapsed 1ns
```

每次运行 stderr 打印 `usage summary: bytes=.. lines=.. variants=..
samples=.. elapsed=.. [limit_hit=..]`，据此回调下次上限。

## 手工验证混合来源队列

```sh
# c2 只有 Z 与 N：文件列为 SNP A1 A2 Z N；CLI 用第 3 段声明 z
# （name=path[:样本ID文件[:beta|or|z]]；留空则按表头自动识别）
/tmp/sggo meta -in c1=beta_se.tsv -in c2=z_only.tsv::z \
              -crosscheck -out /tmp/mixed.tsv
# TSV 末三列：included_cohorts / per_cohort_detail / excluded_detail
# 能看到 z_only 的权重数值、换算口径（p=2*Phi(-|Z|), N=..）、
# 对齐方向依据，以及 missing_sample_size / z_undirected 等排除原因。
```

## 预期行为速查

- 链模糊位点（A/T、C/G）：默认拒绝；不会“按符号相加”改变方向。
- 负遗传力：默认保留并打标志；`-boundary clamp|reject` 可截断或拒绝；
  截断不改变置信区间口径。
- 队列重叠：默认仅标记（含重复队列标志）；可配置剔除或整体拒绝。
- 缺失效应位点：默认拒绝整个 PRS 任务；`-missing-sites ignore` 时
  在交集上评分并显式列出口径差异。
- 超限：默认拒绝且无部分结果；解析层可配置为截断并记录。
