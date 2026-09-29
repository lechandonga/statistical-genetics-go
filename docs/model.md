# 统计模型与统一变体模型

本文件说明 `sggo`（statistical-genetics-go）所实现流程的模型假设、
变量定义与可解释性约定。实现中的所有“拒绝/标记”都能追溯到这里的规则。

## 1. 统一变体模型（`variant.Variant`）

本版本新增三类能力，规则一并固化在本文件中：双数值后端与交叉验证
（见 §6）、可配置输入预算与用量小结（见 §7）、混合来源队列的
Z-only 换算口径（见 §4.1）。

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

### 4.5 混合来源队列：Z-only 的统一换算口径

完整队列（有 beta 与 SE）直接按 §4.1 参与；**仅有 Z 值与样本量 N 的
队列默认也真正参与**，按统一口径换算成等效 (beta, se)，不一律排除：

```
p_i      = 2 * Phi(-|Z_i|)                    # 双侧 P 值
I_i      = 2 * p_i * (1 - p_i) * N_i          # 等效 Fisher 信息量
se_i     = 1 / sqrt(I_i)
beta_i   = Z_i * se_i                         # Z 自带符号
w_i      = I_i = 1/se_i^2                     # 与完整队列同一逆方差权重
```

对齐方向与完整队列完全一致：来源 Z 的符号再乘 `AlignResult.Sign`
（互换取反、链翻转保留）。每个变体的审计行 `MetaCohortLine` 记录：
`SourceKind`（`beta_se`/`z_only`）、`Weight`、`WeightBasis`（含换算
公式与实际使用的 N）、`Sign`、对齐关系与 `DirectionNote`；CLI 元分析
TSV 的 `included_cohorts`/`per_cohort_detail`/`excluded_detail` 三列
直接展示这些内容。

只有以下情况才排除该队列在该变体的贡献，并给出**可区分**原因：

| 情况 | Reason 代码 |
| --- | --- |
| 缺换算所需样本量（N 列与 FallbackN 都没有，N<=0） | `missing_sample_size` |
| Z 为 NaN/Inf | `z_not_finite` |
| Z 恰好为 0：效应方向不存在 | `z_undirected` |
| 换算信息量非有限（极端 P/N） | `z_conversion_non_finite` |
| 对齐失败/链模糊/缺该变体等 | 沿用 §2 的既有原因 |

硬保证：凡出现在合并结果里的数值（beta/se/P/CI/Q/I²/权重）必须有限；
若后端数值异常则该变体以 `numerical_backend_error` 或
`non_finite_combined_result` 整体排除并标记，不存在“参与了合并、
结果却无法解释、又没有任何标记”的条目。参与队列数不足 `MinCohorts`
时，变体以 `insufficient_cohorts`（完全没有可对齐队列时为
`no_alignable_cohort`）排除。

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

三类输入（汇总统计、LD 分数、基因型矩阵）的上限均可在 CLI 按次设置：
`-max-variants`、`-max-samples`（仅基因型样本列）、`-max-line-bytes`、
`-max-file-bytes`、`-max-elapsed`，配合 `-limit-policy reject|truncate`。

- **真实生效**：字节上限在读入路径上累计并在到达上限时立即中止；
  行上限按物理数据行计数；样本列上限作用于表头识别出的样本列；
  墙钟上限在每次取行前检查。
- **reject（默认）**：整体失败，返回 nil 结果，不落任何半截输出文件
  （CLI 在解析成功后才创建分析输出）。
- **truncate**：确定性保留文件/表头顺序的前 N 行或前 N 个样本列，
  报告 `Truncated=true`、`TruncateAt=N` 与命中的上限名。
- **用量小结**：每次解析把实际 `BytesRead/Lines/Variants/Samples/
  Elapsed` 写进对应 `ParseReport.Usage`，CLI 在 stderr 打印；
  Usage 是每次解析的独立实例，并发任务之间不共享。

## 8. 双数值后端与交叉验证语义

- 后端：`stdlib`（默认，Go 标准库 math + Kahan）与 `gonum`
  （gonum/gonum 线代库，版本锁定见 `docs/compatibility.md`）。
- 选择粒度：进程默认（`numwrap.SetDefaultBackend`，加锁）与
  **按次显式选择**（各 Options.Backend / CLI `-backend`）。
  按次选择在入口解析一次后贯穿整个任务；切换默认后端不会打断、
  改变任何在跑的任务，也不存在半切换状态。
- 复现性：同一后端对同一输入（含相同并行度以外的任意并行度）
  **逐位一致**；两个后端之间只要求在声明容差内一致
  （rel `1e-9`、abs `1e-12`，见 `numwrap` 常量）。
- 判定一致性：交叉验证不仅比数值，还比所有离散结论——h2 flags、
  变体/队列纳入与排除原因、来源类型、PRS active 计数等必须完全相同。
- 失败语义：容差超限或判定不一致 => `BackendMismatchError` 硬失败；
  未知后端 => `KindUnsupportedOption`，绝不静默回退默认实现。
- 库错误隔离：gonum 的错误/panic 在 numwrap 内翻译为本包错误类型，
  不透出第三方错误文本（映射表见 `docs/compatibility.md`）。

- `LimitReject`（默认）：超限整体拒绝，不返回部分接受结果；
- `LimitTruncate`：截断并在报告中记录 `Truncated/TruncateAt`。
流程级另有 `Config.Deadline` 墙钟上限（context 超时）。
