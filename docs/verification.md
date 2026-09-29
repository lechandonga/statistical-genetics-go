# 验证方法

## 一键运行

```sh
go vet ./...
go test ./...            # 全部单测 + 端到端测试
go test ./... -race      # 竞态检测（并行分片路径）
go test -count=2 ./...   # 重复运行，确认无不确定行为
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
| 双后端 | `numwrap/backend_crosscheck_test.go` | stdlib/gonum WLS/矩/分位数一致、同后端逐位一致、crosscheck 判定一致、未知后端拒绝、并发切换竞态安全 |
| 遗传力 | `stats/h2_test.go` | 真值恢复、负 h2 三策略、越界、模糊剔除、workers=1..100 逐位一致 |
| 元分析 | `stats/meta_test.go` | IVW 精确值、Q/I²、DL tau²、链翻转/互换、模糊、重叠三策略、REML 拒绝 |
| 混合来源元分析 | `stats/meta_mixed_test.go` | Z+N+EAF 按统一口径参与、缺 N/缺 EAF/EAF=0|1 可区分排除、互换取反、双后端判定一致、并发切换竞态安全 |
| 统计层双后端 | `stats/backend_crosscheck_test.go`、`stats/prs_backend_test.go` | 带噪 h2 的双后端数值+Flags+拒绝明细一致、PRS std 尺度与评估回归一致 |
| PRS | `stats/prs_test.go` | 方向、缺失 GT 两策略、缺失位点两策略、三尺度、阈值、评估回归 |
| 端到端 | `pipeline/pipeline_test.go` | 注入非法/重复行拒绝、共享样本重叠标记/拒绝、串行/4 worker 一致 |
| 上限与用量 | `ingest/usage_limits_test.go` | 字节/行/样本/耗时上限真实生效（含极小值）、reject 无半结果、truncate 完整行边界、三类输入用量、并发不串扰 |
| 新特性端到端 | `pipeline/backend_mixed_e2e_test.go` | crosscheck 跑通全流程、按类型上限拒绝/截断、混合 Z 队列纳入与审计原因 |

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
GOFLAGS=-mod=vendor GOPROXY=off go build -o /tmp/sggo ./cmd/sggo
/tmp/sggo simulate -out /tmp/demo -variants 500 -samples 300 -overlap 80
/tmp/sggo h2   -sum /tmp/demo/cohort1_sumstats.tsv -ld /tmp/demo/ldscores.tsv -n 280 -backend crosscheck
/tmp/sggo meta -in c1=/tmp/demo/cohort1_sumstats.tsv -in c2=/tmp/demo/cohort2_sumstats.tsv \
               -model random -backend gonum \
               -out /tmp/demo/meta.tsv -audit-out /tmp/demo/audit.tsv
/tmp/sggo prs  -sum /tmp/demo/cohort1_sumstats.tsv -gt /tmp/demo/genotypes.tsv \
               -n 280 -missing-sites ignore -scale std -backend crosscheck \
               -out /tmp/demo/prs.tsv
```

带重叠硬证据的元分析（为队列提供每行一个样本 ID 的文件）：

```sh
/tmp/sggo meta -in c1=sum1.tsv:ids1.txt -in c2=sum2.tsv:ids2.txt -overlap reject
```

`sggo version` 打印全部已编译后端及其锁定版本与 crosscheck 容差。

### 手工核对两个后端一致

1. 同一份输入分别跑 `-backend stdlib` 与 `-backend gonum`，比较
   `h2/se/CI/flags` 或 meta TSV 的数值列；
2. 直接跑 `-backend crosscheck`：它在一次调用内用两个后端分别求解，
   数值超过 `rel=1e-9/abs=1e-12` 或任何判定字段（拒绝、原因、
   flags、纳入/排除集合）不一致都会以 `backend_mismatch` 失败，
   成功时输出与 `-backend stdlib` 相同的基准值；
3. 重复跑同一后端两次（或改 `-workers`），报告值必须逐位相同。

### 手工验证各类上限真的生效

```sh
# 行数：截断到 50 行，stderr 用量小结显示 truncated=true
/tmp/sggo h2 -sum /tmp/demo/cohort1_sumstats.tsv -ld /tmp/demo/ldscores.tsv \
             -n 280 -sum-max-variants 50 -sum-limit-policy truncate
# 字节：上限设极小必须整体拒绝，且不产生结果文件
/tmp/sggo h2 -sum /tmp/demo/cohort1_sumstats.tsv -ld /tmp/demo/ldscores.tsv \
             -n 280 -sum-max-file-bytes 100
# 样本：基因型只保留前 10 列
/tmp/sggo prs -sum /tmp/demo/cohort1_sumstats.tsv -gt /tmp/demo/genotypes.tsv \
              -n 280 -missing-sites ignore -out /tmp/p.tsv \
              -gt-max-samples 10 -gt-limit-policy truncate
# 耗时：1ns 预算 + truncate 在完整行边界结束（结果仍为完整记录）
```

每次运行 stderr 的 `usage: bytes=... lines=... dataRows=...
samples=采用/提供 elapsed=... truncated=...` 即实际用量小结。

### 手工验证混合来源队列

```sh
# c1 为 beta+SE，c2 只有 Z/N/EAF：
/tmp/sggo meta -in c1=/path/beta.tsv -in c2=/path/z.tsv -z c2 \
               -backend crosscheck \
               -out /tmp/mix.tsv -audit-out /tmp/mix_audit.tsv
```

审计文件逐行给出每个队列是否 included、`basis`（`beta_se` 或
`z_n_eaf`）、权重、对齐 sign/relation、N/Z/EAF 与排除原因；
缺 N 行为 `z_missing_n`，缺 EAF（或 EAF=0/1）为
`z_unorientable_maf`，两种原因可区分。

## 预期行为速查

- 链模糊位点（A/T、C/G）：默认拒绝；不会“按符号相加”改变方向。
- 负遗传力：默认保留并打标志；`-boundary clamp|reject` 可截断或拒绝；
  截断不改变置信区间口径。
- 队列重叠：默认仅标记（含重复队列标志）；可配置剔除或整体拒绝。
- 缺失效应位点：默认拒绝整个 PRS 任务；`-missing-sites ignore` 时
  在交集上评分并显式列出口径差异。
- 超限：默认拒绝且无部分结果；解析层可配置为截断并记录。
