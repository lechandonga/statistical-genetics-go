# 兼容范围与第三方接口归一化

## 数值/统计后端（`numwrap`）

当前构建只包含一个后端：

| 后端 | 版本/来源 | 状态 |
| --- | --- | --- |
| `stdlib` | Go 标准库 `math`（go1.25+ 工具链） | 默认，自包含，无第三方依赖 |
| `gonum` | gonum/gonum | **未编译**：请求即 `KindUnsupportedOption` |

`numwrap.Version(backend)` 返回后端版本/接口描述；
`SetDefaultBackend` 是进程级、加锁的切换，只影响之后新建的求解器。
对外 API 不返回任何第三方库类型；即使将来接入 gonum，其错误也会被
翻译成下列本包错误类型。

### 错误分类（`numwrap.Error.Kind`）

| Kind | 触发情形 |
| --- | --- |
| `unsupported_option` | 请求未编译后端、REML、协变量调整等不支持选项 |
| `invalid_argument` | 长度/维度不一致、观测数不足、负权重、p∉(0,1) 等 |
| `singular` | 设计矩阵退化（x 零方差、权重和为 0） |
| `nan` / `inf` | 输入或中间结果出现 NaN / Inf |
| `no_convergence` | 预留给迭代算法（当前后端无迭代求解器） |

封装的数值原语：

- `WLS(x, y, opts)`：单自变量加权最小二乘，支持自由截距或固定截距；
- `KahanSum` / `MeanStd`：固定顺序补偿求和，确定性矩；
- `NormalCDFInv`：`-sqrt2 * erfcinv(2p)`，标准库实现、跨平台一致。

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
- 纯标准库依赖，`go vet ./...` 与 `go test ./...` 离线可运行。
