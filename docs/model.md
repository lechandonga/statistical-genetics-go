# 统计模型与统一变体模型

本文件说明 `sggo`（statistical-genetics-go）所实现流程的模型假设、
变量定义与可解释性约定。实现中的所有“拒绝/标记”都能追溯到这里的规则。

## 1. 统一变体模型（`variant.Variant`）

所有输入（汇总统计、LD 分数、基因型矩阵）都先归一化到同一结构：

| 字段 | 含义 | 归一化规则 |
| --- | --- | --- |
| `ID` | 变体标识 | 去空白、大写；作为首选匹配键 |
| `Chrom` | 染色体 | `chr1/1`→`1`；`23/24/25`→`X/Y/MT`；`chrM`→`MT` |
| `Pos` | 物理位置 | 1-based 正整数；与染色体必须同时出现 |
| `EA` | 效应等位基因 | 去空白、大写，仅允许 `A/C/G/T` |
| `NEA` | 非效应等位基因 | 同上；EA 与 NEA 不得相同 |
| `EAF` | 效应等位基因频率 | 缺省 `NaN`；若存在必须在 `[0,1]` |
| `Strand` | 链方向 | 反链声明在**构造时**立即互补到参考链 |

匹配键（`MatchKey`）：优先 `ID`，无 ID 时用 `CHROM:POS`。
匹配只在同键变体间进行；跨键从不猜测对应关系。

### 1.1 输入拒绝原因代码

拒绝不静默发生，每行错误都带 `Reason` 代码（见 `variant/errors.go`）：

`missing_id`、`missing_alleles`、`invalid_allele`、`invalid_chromosome`、
`invalid_position`、`invalid_frequency`、`invalid_numeric`、
`duplicate_variant`、`unresolvable_strand`、`ambiguous_strand`、
`flip_forbidden`、`allele_mismatch`、`missing_field`。

解析器默认**严格模式**：出现任何被拒绝行即整体失败，不返回部分结果；
设置 `Lenient=true` 时返回接受记录与完整 `ParseReport`（拒绝明细可审计）。

## 2. 等位基因对齐规则（`variant.Align`）

方向基准始终是**参考侧**：H2 用 LD 参考，元分析用每个变体第一个有效
队列，PRS 用个体基因型矩阵。设目标为 `(EA_t, NEA_t)`、参考为
`(EA_r, NEA_r)`，按以下顺序判定：

1. **直接匹配**：`EA_t=EA_r, NEA_t=NEA_r` → 效应符号 `+1`
2. **等位基因互换（swap）**：`EA_t=NEA_r, NEA_t=EA_r` → 效应符号 **`-1`**
3. **链翻转（flip）**：`comp(EA_t)=EA_r, comp(NEA_t)=NEA_r` → 符号 `+1`
4. **翻转+互换**：`comp(EA_t)=NEA_r, comp(NEA_t)=EA_r` → 符号 **`-1`****
5. 其余一律 `allele_mismatch` 拒绝。

链翻转为什么**不**改变 beta 符号：beta 定义为“每增加一个 EA 拷贝的
表型变化”。链翻转只重命名碱基（A↔T、C↔G），个体的 EA 拷贝数不变；
而互换改变了哪个等位基因被视为“效应”等位基因，拷贝数从 d 变为 2−d，
所以效应必须取反。**绝不允许在方向不明时按符号直接相加。**

### 2.1 链模糊位点（A/T、C/G）

对 SNV 而言，A/T 互补为 T/A，与“互换”在观测上不可区分；C/G 同理。
此时无法判断目标究竟做了链翻转还是等位基因互换，效应符号不可确定。

- `AmbiguousReject`（默认）：直接拒绝，原因 `ambiguous_strand`；
- `AmbiguousKeep`：仅接受无需判链的直接匹配/互换（并置 `Ambiguous=true`），
  任何需要借助链翻转才能成立的对齐仍拒绝；
- `FlipForbidden`：情况 3/4 一律以 `flip_forbidden` 拒绝。

## 3. 遗传力估计（`stats.EstimateH2`）

LDSC 风格的 LD 分数回归：

```
E[chi2_j] = intercept + slope * l2_j,   chi2_j = (beta_j / se_j)^2
slope     = Nbar * h2 / M
h2        = slope * M / Nbar
```

- `M` 为通过内连接与对齐的 SNP 数；N 优先取逐 SNP N（加权平均），
  否则使用 `H2Options.NBar`，两者皆无则拒绝（`missing_sample_size`）。
- 截距：`FreeIntercept=true` 自由估计（校正膨胀）；否则固定为 1，
  结果置 `intercept_fixed` 标志。
- 求解统一走 `numwrap.WLS`（加权最小二乘），权重可选 LDSC 简化形式
  `1/max(2*l2, eps)`；不暴露文档未支持的迭代选项。

### 3.1 负值、越界与不可解释区间

点估计不做任何静默修正。`H2BoundaryPolicy`：

- `H2ReportRaw`（默认）：保留原始值，用位标志 `negative_raw` /
  `above_max_raw` 表达；研究场景据此可审计。
- `H2Clamp`：报告值截断到 `[0, h2Max]`（默认 h2Max=1）并置 `clamped`，
  但**置信区间仍基于未截断分布**给出，不制造虚假的精确性。
- `H2Reject`：负值或 > h2Max 直接返回可区分错误。
- 当 95% CI 半宽超过 h2Max，置 `low_power`（区间宽到不可解释）。
- SNP 数低于 `MinVariants`（默认 20）返回 `insufficient_variants`。

不确定度：SE 来自 WLS 斜率方差乘以 `M/Nbar`；CI 用
`beta ± z(0.975)*SE`（`z` 由标准库 erfcinv 计算，跨平台一致）。

## 4. 多队列元分析（`stats.RunMeta`）

### 4.1 固定效应（IVW）

```
w_i = 1/se_i^2
beta_F = sum(w_i*beta_i) / sum(w_i)
se_F   = 1/sqrt(sum(w_i))
```

beta 在合并前先对齐到该变体第一个有效队列的方向；对齐失败的队列以
`PerCohort[i].Excluded=true` + 原因保留，**不静默丢弃**。

### 4.5 混合来源队列：Z-only 队列的统一口径

队列可以只提供 Z 分数与样本量 N（CLI `-z <cohort>`，解析效应模式
`ingest.EffectZ`），没有 beta 与 SE。这类队列**默认真正参与合并**，
按下列统一口径把 Z 换算到与 beta/SE 队列相同的逆方差尺度：

```
Var(X) = 2*N*p*(1-p)          # 加性剂量 0/1/2，p=效应等位基因频率
SE     = 1 / sqrt(2*N*p*(1-p))
beta   = Z * SE               # 方向直接取 Z 的符号
```

- 对齐到参考方向时 `Z_aligned = Z * AlignResult.Sign`，链互换
  （swap）取负、链翻转不变号，规则与 beta 队列完全相同；
- 审计字段 `MetaCohortLine.Basis` 标明口径：`beta_se`（直接提供）或
  `z_n_eaf`（按上式换算），`BasisDetail` 给出 N/p/SE 的完整代入文本；
- 只有以下情况才排除，并给出可区分原因：

| 原因 | 触发条件 |
| --- | --- |
| `z_missing_n` | Z-only 队列连换算所需的样本量 N 也缺失 |
| `z_unorientable_maf` | EAF 缺失，或 EAF 为 0/1（剂量方差为 0，beta 尺度无法定向） |
| `z_invalid` | Z 为 NaN/Inf |
| `non_finite_contribution` | 换算后 beta/SE/权重仍非有限（防御性兜底） |

变体层面参与队列数低于 `MinCohorts` 时仍标记为 `insufficient_cohorts`，
但每个队列各自的排除原因都在 `PerCohort` 与 `-audit-out` 文件中可见。
任何“参与了合并”的条目，其 beta/SE/权重都保证有限且带口径标记——
不允许出现参与了合并、结果不可解释、却没有任何标记的条目。

### 4.2 异质性

- Cochran Q = Σ w_i (beta_i − beta_F)²，自由度 k−1；
- I² = max(0, (Q − df)/Q)；
- P_Q 用 Wilson–Hilferty 卡方正态近似（明确为近似值，用于展示）。

### 4.3 随机效应（DerSimonian–Laird）

```
C     = Σw − Σw²/Σw
tau²  = max(0, (Q − df)/C)
w_i*  = 1/(se_i² + tau²)   再做 IVW
```

REML 后端未实现，请求 `MetaRandomREML` 会得到显式
`unsupported_model` 拒绝，而非退化成别的估计。

### 4.4 样本重叠与重复队列

重叠证据**只**接受硬证据：显式提供的样本 ID（或其稳定哈希）集合的
交集；没有提供 ID 的队列**绝不**基于效应相关性猜测重叠。

- 输出 `OverlapEvidence{共享数, Jaccard, Duplicate}`；
- `OverlapFlag`（默认）：置 `sample_overlap_detected` /
  `duplicate_cohort_detected` 标志，权重明细仍给出，调用方不得把
  被标记对当作相互独立证据；
- `OverlapDrop`：每组重叠保留先出现的队列，其余整体排除
  （`ExcludedCohorts` 稳定排序输出）；
- `OverlapReject`：检测到任意重叠即整体失败。

有效队列数低于 `MinCohorts`（默认 2）的变体被排除并说明原因。

## 5. 个体风险评分（`stats.BuildPRS`）

方向基准是个体基因型矩阵。每个模型位点对齐后：

```
raw score_j = Σ_i beta_i(aligned) * dosage_ij
```

- 缺失基因型：`GTMeanImpute` 用 `2*EAF`（EAF 缺省退化为 1，审计可见）；
  `GTExclude` 该位点不计入该样本，`NActive` 不增加。
- 效应位点缺失：`ModelStrict`（默认）整体拒绝（口径必须显式确认）；
  `ModelIgnore` 只在交集上评分，并在 `MissingModelSites` 中列出
  模型缺失与“矩阵有、模型无”的双向差异。
- P 值阈值过滤：无 P 值可用的位点在打开阈值时以 `pvalue_unavailable`
  记录，而非悄悄纳入。
- 尺度：`raw` / `perallele`（÷有效位点数）/ `std`（按本批样本标准化
  为均值 0、方差 1；零方差显式拒绝）。

`PRSResult` 固定记录位点集合、每个位点的符号/关系/缺失数与尺度，
位点集合变化必须显式改变口径（名称/审计差异可见），不会静默改变
符号约定或评分口径。

评估（`EvaluatePRS`）为简单线性回归 y~score（R²=corr²）；协变量调整
当前后端不支持，请求即 `unsupported_option` 拒绝。

## 6. 并发与可复现性

- 唯一并行原语是 `stats.mapShards`：分片边界只依赖 (n, workers)，
  worker 写入固定槽位，结果严格按分片顺序归并；fn 必须为纯函数。
- 所有跨分片合并量（遗传力充分统计量、IVW 权重/加权和、评分累加、
  标准化矩）都用固定顺序 Kahan 补偿求和，因此 workers=1 与任意
  workers 结果**逐位相同**；单测直接比对浮点位值。
- 任何分片错误/ctx 取消即整体失败，不返回部分结果。

## 7. 规模与耗时上限

`ingest.Limits`：`MaxVariants`、`MaxSamples`、`MaxLineBytes`、
`MaxFileBytes`、`MaxElapsed`，默认值见 `DefaultLimits()`。

- `LimitReject`（默认）：超限整体拒绝，不返回部分接受结果；
- `LimitTruncate`：截断并在报告中记录 `Truncated/TruncateAt`。

各项的真实语义（实现保证，CLI 按 `sum`/`ld`/`gt` 三类输入分别开放）：

- `MaxFileBytes`：流式读取，**读到超过上限的字节即触发**（不会先把
  整个文件读进内存再判断）；reject 立即失败且无部分结果，truncate
  在已读完整行边界结束并返回完整记录；
- `MaxVariants`：数据行计数（不含表头/空行，含被拒绝行），超限
  reject/truncate；
- `MaxSamples`：基因型表头样本列计数；reject 整体失败，truncate 按
  表头列顺序保留前 N 列，`Usage.Samples/SamplesOffered` 同时给出
  采用数与实际提供数；
- `MaxLineBytes`：单行畸形巨长，**无论 truncate 与否都拒绝**
  （截断一行中间没有安全语义）；
- `MaxElapsed`：每行读取前检查墙钟；truncate 时同样在完整行边界结束。

每次解析都在 `ParseReport.Usage` 中独立记录：`BytesRead`（实际从
输入读取的字节数，含表头与被拒绝行）、`LinesRead`（物理行）、
`DataRows`、`Samples`/`SamplesOffered`、`Elapsed`、`Truncated`；
计量不使用任何进程级共享状态，并发解析互不串扰。CLI 在每次运行后
把该小结打到 stderr。

流程级另有 `Config.Deadline` 墙钟上限（context 超时）。

## 8. 数值后端与交叉验证

见 [compatibility.md](compatibility.md) 的“数值/统计后端”一节：
`stdlib`（默认）与 vendored `gonum v0.17.0` 双后端，`crosscheck`
元后端按调用运行两边并比较——数值容差 `rel=1e-9, abs=1e-12`，
判定（拒绝/标记/纳入/排除/原因/符号）必须逐字段一致，对不上即
`backend_mismatch` 显式报错。后端选择是每次调用独立的，切换默认
后端不影响正在运行的任务，也不会产生算了一半的状态
（crosscheck 两边都完整跑完后才比较并返回）。
