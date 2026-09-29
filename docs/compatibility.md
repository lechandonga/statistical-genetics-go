# 兼容范围与第三方接口归一化

## 数值/统计后端（`numwrap`）

当前构建编译进两个相互独立的后端：

| 后端 | 版本/来源 | 状态 |
| --- | --- | --- |
| `stdlib` | Go 标准库 `math`（go1.25+ 工具链） | 默认，自包含，无第三方依赖 |
| `gonum` | `gonum.org/v1/gonum v0.15.1` | 已编译，经 `vendor/` 锁定离线可用 |

gonum 锁定版本同时出现在三处并由测试强制一致：

- `go.mod` / `go.sum`（校验和）
- `vendor/modules.txt` 与 `vendor/` 源码树
- `numwrap.GonumVersion` 常量（`numwrap.Version("gonum")` 输出它）

断网复现方式：

```sh
GOPROXY=off GOFLAGS=-mod=vendor go build ./...
GOPROXY=off GOFLAGS=-mod=vendor go test ./... -race
```

`numwrap.Version(backend)` 返回版本/接口描述；
`numwrap.AvailableBackends()` 列出实际编译的后端；
`SetDefaultBackend` 是进程级、加锁的切换，只影响之后新建的求解器。
**所有统计入口同时支持按次显式指定后端**（`H2Options.Backend`、
`MetaOptions.Backend`、`PRSOptions.Backend`/`EvalOptions.Backend`、
CLI `-backend`），入口处解析一次即贯穿整个任务，运行期间切换进程默认
后端不会影响在跑的任务，也不存在“算了一半换后端”的状态。
请求未编译/拼错的后端返回 `KindUnsupportedOption`，绝不静默回退默认。

对外 API 不返回任何第三方库类型；gonum 的错误或 panic 都被翻译成
下列本包错误类型（第三方错误文本不向外透出）。

### 双后端交叉验证

`CrossCheck=true`（CLI `-crosscheck`）时，同一输入分别用 `stdlib` 与
`gonum` 完整计算两遍：

- 同一后端重复运行：**逐位一致**（固定顺序求和/求解）。
- 两个后端之间：数值差异必须满足
  `|a-b| <= 1e-12 + 1e-9*max(|a|,|b|,1)`
  （`numwrap.CrossCheckRelTol` / `CrossCheckAbsTol`）。
- 离散判定（是否拒绝、是否标记、变体是否纳入、每个队列的纳入/排除
  原因与来源类型、h2 flags、PRS active 计数）必须**完全一致**。
- 任何不一致都返回 `*numwrap.BackendMismatchError`，任务整体失败；
  不存在“挑一侧结果返回”的路径。

### 错误分类（`numwrap.Error.Kind`）

| Kind | 触发情形 |
| --- | --- |
| `unsupported_option` | 请求未编译后端、REML、协变量调整等不支持选项 |
| `invalid_argument` | 长度/维度不一致、观测数不足、负权重、p∉(0,1) 等 |
| `singular` | 设计矩阵退化（x 零方差、权重和为 0） |
| `nan` / `inf` | 输入或中间结果出现 NaN / Inf |
| `no_convergence` | 预留给迭代算法（当前后端无迭代求解器） |

另有独立类型 `*numwrap.BackendMismatchError`：交叉验证发现两后端结论
不一致时作为硬错误返回（任务整体失败，不挑一侧结果）。

封装的数值原语：

- `WLS(x, y, opts)`：单自变量加权最小二乘，支持自由截距或固定截距；
- `KahanSum` / `MeanStd`：固定顺序补偿求和，确定性矩；
- `NormalCDFInv`：`-sqrt2 * erfcinv(2p)`，标准库实现、跨平台一致。
- `MetaCombine(betas, ses, opts)`：固定效应 IVW + 可选 DL 随机效应；
  stdlib 走 Kahan 求和，gonum 后端是独立的朴素顺序累加实现。

gonum 后端具体映射：

| 本包原语 | gonum 实现 |
| --- | --- |
| `WLS`（gonum） | `gonum/mat` QR 最小二乘 + `(XwᵀXw)⁻¹` 协方差 |
| `MetaCombine`（gonum） | 独立顺序累加的 IVW/DL（不共享 std 代码） |
| `MeanStd`（gonum） | `gonum/stat.Mean` / `stat.Variance` |
| `NormalCDFInv`（gonum） | `gonum/stat/distuv.UnitNormal.Quantile` |

### 错误翻译规则（gonum）

- 维度/长度非法、权重为负、NaN/Inf 输入：**在调用 gonum 之前**由本包
  校验拦截，直接给 `invalid_argument`/`nan`/`inf`；
- 设计矩阵秩亏/奇异（QR 求解或信息矩阵求逆失败）：`singular`；
- gonum 内部 panic：统一 recover 为 `singular`（线代）或
  `invalid_argument`（分位数），错误文本只保留失败类别，不带第三方堆栈。

第三方库的语法错误（例如选项拼写导致的不支持选项）与内部运行时异常
（奇异、NaN）因此被区分为不同错误，而不是以同一字符串透出。

## 输入文件格式

全部为**本地制表符分隔文本（TSV）**，UTF-8，首行表头。无网络/数据库
依赖。列名大小写不敏感，并忽略 `_`、`.`、空格、`-` 的差异
（`Effect_Allele` 与 `effectallele` 等价）。

### 汇总统计

必填：`SNP`（或 `CHROM+BP`）、`A1/EA/ALT`、`A2/NEA/REF`、效应列、
（beta/or 模式下）`SE`；N 列或调用方 `FallbackN` 至少其一。

| 语义 | 接受的列名（示例） |
| --- | --- |
| ID | SNP, SNPID, ID, VariantID, rsID |
| 染色体/位置 | CHR/CHROM/CHROMOSOME；BP/POS/POSITION |
| 等位基因 | A1/EA/EFFECT_ALLELE/ALT；A2/NEA/REF |
| 效应 | BETA/EFFECT/B；OR；Z（三种模式显式选择） |
| 其它 | SE, P/PVAL, N/NEFF, EAF/FRQ/MAF, STRAND |

- OR 模式内部取自然对数转为 logOR；OR<=0 拒绝。
- Z 模式为退化模式：直接采用 Z，不反演 SE（SE 记为 NaN）。
- STRAND 列接受 `+`/`-`；`-` 在构造变体时即互补到参考链。

### LD 分数

必填：`SNP`（或 CHR+BP）、`A1`、`A2`、`L2`（>=0）；
可选 `MAF∈[0,1]`、`CHISQ>=0`、`N>0`。

### 基因型矩阵

变体为行、样本为列。注释列默认识别
`SNP/CHR/CHROM/BP/POS/A1/A2/EA/NEA/STRAND`，其余列视为样本。

- `012`（默认）：EA 拷贝数剂量，取值仅允许 0/1/2；
- `acgt`：两碱基二联体（如 `AA`/`AG`），按该行 EA/NEA 数 EA 拷贝数；
  出现两等位基因之外的碱基拒绝；
- `space`：`A A` 形式（`GTSpaceAlleles`）。
- 缺失标记默认 `NA/NaN/./?/-/空`；样本数受 `MaxSamples` 限制。

## Go 版本与模块

- `go.mod` 声明 `go 1.25.0`；验证环境 `go1.26 linux/amd64`。
- 唯一直接第三方依赖为 `gonum.org/v1/gonum v0.15.1`，已 vendor；
  `go vet ./...` 与 `go test ./...` 在 `GOPROXY=off` 下离线可运行。
