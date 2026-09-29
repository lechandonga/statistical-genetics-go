# 兼容范围与第三方接口归一化

## 数值/统计后端（`numwrap`）

当前构建包含两个真实后端与一个交叉验证“元后端”：

| 后端 | 版本/来源 | 状态 |
| --- | --- | --- |
| `stdlib` | Go 标准库 `math`（go1.25+ 工具链） | 默认，自包含，无第三方依赖 |
| `gonum` | `gonum.org/v1/gonum v0.17.0`（只用 `gonum/mat`、`gonum/stat`、`gonum/stat/distuv`） | 已编译，版本锁定并随仓库 vendor |
| `crosscheck` | stdlib + gonum，每次调用分别求解后比较 | 已编译；不是新求解器 |

`numwrap.Version(backend)` 返回后端版本/接口描述；
`SetDefaultBackend` 是进程级、加锁的切换，只影响之后新建的求解器。
所有统计入口（`stats.EstimateH2`/`RunMeta`/`BuildPRS`/`EvaluatePRS`）
都接受按调用传入的 `Backend`，**不读取**进程级默认值：切换默认后端
不会影响任何正在执行或已显式指定后端的任务，并发任务的后端选择
互不串扰。

### 版本锁定与断网构建

- `go.mod` 要求 `gonum.org/v1/gonum v0.17.0`，`go.sum` 锁定校验和，
  完整源码随仓库放在 `vendor/`（见 `vendor/modules.txt`）。
- 断网环境使用 vendor 模式构建与测试，不访问任何模块代理：

```sh
GOFLAGS=-mod=vendor GOPROXY=off go build ./...
GOFLAGS=-mod=vendor GOPROXY=off go test ./... -race -count=1
GOFLAGS=-mod=vendor GOPROXY=off go vet ./...
```

升级 gonum 必须同时更新 `go.mod`、`go.sum`、`vendor/`、
`numwrap.backendVersions` 与本文档的版本表。

### 双后端一致性约定

- 同一后端、同一输入重复运行：报告值（系数/SE/协方差/评分）**逐位一致**；
- `stdlib` 与 `gonum` 之间：数值允许在声明容差内不同，
  `|a-b| <= max(1e-9*max(|a|,|b|), 1e-12)`
  （常量 `numwrap.CrosscheckRelTol/AbsTol`，`sggo version` 会打印）；
- **判定必须完全一致**：是否拒绝、拒绝原因代码、Flags、纳入/排除的
  队列与变体集合、权重口径、符号/对齐关系等逐字段比较；
- `crosscheck` 两边都成功且比较通过才返回 `stdlib` 基准值；任一边
  出错时要求两边错误分类一致，否则直接报 `backend_mismatch`
  （numwrap 层 `KindBackendMismatch` / stats 层
  `*stats.BackendMismatchError`），**绝不挑一边的结果返回**；
- 未知/未编译后端一律 `unsupported_option` 拒绝，不退回默认后端。

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
| `backend_mismatch` | 双后端数值超容差、判定字段不一致、或一边成功一边失败 |

封装的数值原语：

- `WLS(x, y, opts)`：单自变量加权最小二乘，支持自由截距或固定截距；
- `KahanSum` / `MeanStd`：固定顺序补偿求和，确定性矩；
- `NormalCDFInv`：stdlib 用 `-sqrt2 * erfcinv(2p)`，gonum 用
  `distuv.UnitNormal.Quantile`，两者在声明容差内一致。

gonum 内部的 panic（维度越界、矩阵不可逆、数值异常等）在 numwrap
边界统一 recover 并翻译：奇异→`singular`，非有限数→`nan`/`inf`，
其余→`unknown`；第三方错误类型与内部文本不会透出接口。输入校验
（长度、权重、p 范围、有限性）由两个后端共享同一函数完成，因此
同一份非法输入在两个后端得到完全相同的错误分类。

## 输入文件格式

全部为**本地制表符分隔文本（TSV）**，UTF-8，首行表头。无网络/数据库
依赖。列名大小写不敏感，并忽略 `_`、`.`、空格、`-` 的差异
（`Effect_Allele` 与 `effectallele` 等价）。

### 汇总统计

必填：`SNP`（或 `CHROM+BP`）、`A1/EA/ALT`、`A2/NEA/REF`、效应列、
（beta/or 模式下）`SE`；N 列或调用方 `FallbackN` 至少其一。

**例外（混合来源队列）**：`z` 效应模式允许缺 N——记录照常解析，
元分析在该变体上会以可区分原因 `z_missing_n` 排除该队列；beta/or
模式下缺 N 仍然直接拒绝（遗传力回归等下游直接依赖 N）。

| 语义 | 接受的列名（示例） |
| --- | --- |
| ID | SNP, SNPID, ID, VariantID, rsID |
| 染色体/位置 | CHR/CHROM/CHROMOSOME；BP/POS/POSITION |
| 等位基因 | A1/EA/EFFECT_ALLELE/ALT；A2/NEA/REF |
| 效应 | BETA/EFFECT/B；OR；Z（三种模式显式选择） |
| 其它 | SE, P/PVAL, N/NEFF, EAF/FRQ/MAF, STRAND |

- OR 模式内部取自然对数转为 logOR；OR<=0 拒绝。
- Z 模式直接采用 Z 作为方向，SE 记为 NaN；效应/SE 的统一换算发生在
  元分析层（见 [model.md](model.md) §4.5），解析层不猜测效应尺度。
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
- 唯一第三方依赖是 `gonum.org/v1/gonum v0.17.0`（已 vendor）；
  `go vet ./...` 与 `go test ./...` 在
  `GOFLAGS=-mod=vendor GOPROXY=off` 下离线可运行。

## 命令行新增参数

| 参数 | 适用子命令 | 含义 |
| --- | --- | --- |
| `-backend stdlib\|gonum\|crosscheck` | h2/meta/prs | 按本次调用选择数值后端；未知值显式报错 |
| `-z <cohort>`（可重复） | meta | 声明该队列为 Z-only（只有 Z/N，无 beta/SE） |
| `-audit-out <path>` | meta | 逐变体×逐队列审计 TSV（参与/排除、权重、口径、方向依据） |
| `-sum-*` / `-ld-*` / `-gt-*` | h2/meta/prs | 按输入类型独立设置上限（见下） |

按输入类型的上限参数（每类一套，互不覆盖）：

```
-<type>-max-variants N    数据行（变体）上限
-<type>-max-samples N     样本列上限（基因型）
-<type>-max-line-bytes N  单行字节上限（超限的畸形行始终拒绝）
-<type>-max-file-bytes N  文件总字节上限
-<type>-max-elapsed DUR   解析墙钟预算（如 500ms）
-<type>-limit-policy reject|truncate
```

`<type>` 为 `sum`（汇总统计）、`ld`（LD 分数）、`gt`（基因型矩阵）。
每次解析结束后在 stderr 打印一行用量小结：
实际读取字节数、物理行数、数据行数、样本数（采用/表头提供）、
墙钟耗时与是否截断。用量按每次解析独立计量（无进程级共享状态），
并发任务不会互相串。
