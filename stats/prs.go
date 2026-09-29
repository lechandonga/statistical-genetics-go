package stats

import (
	"context"
	"fmt"
	"math"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// MissingGTPolicy 控制个体缺失基因型的处理。
type MissingGTPolicy uint8

const (
	// GTMeanImpute 用 2*EAF 期望剂量插补（默认）；EAF 未知时用 1，
	// 并在审计中记录该位点没有真实观测。
	GTMeanImpute MissingGTPolicy = iota
	// GTExclude 该样本在该位点不贡献，NActive 不增加。
	GTExclude
)

// MissingModelPolicy 控制效应模型缺失位点的处理。
type MissingModelPolicy uint8

const (
	// ModelStrict 效应位点缺失即拒绝整个评分任务（默认：口径必须显式确认）。
	ModelStrict MissingModelPolicy = iota
	// ModelIgnore 只用交集位点评分，并在结果中列出缺失位点与口径变化。
	ModelIgnore
)

// PRSOptions 控制评分构建。
type PRSOptions struct {
	Name       string
	Align      variant.AlignConfig
	PThreshold float64 // <=0 表示不做 P 值过滤
	GT         MissingGTPolicy
	Model      MissingModelPolicy
	// ScoreScale:
	//
	//	"raw"       = sum(beta*dosage)
	//	"std"       = 按本批样本标准化（均值0、方差1，确定性单遍算法）
	//	"perallele" = 除以每样本实际计入的位点数
	ScoreScale string
	Workers    int
	// Backend 指定数值后端（空=默认；crosscheck=双后端验证）。
	Backend numwrap.Backend
}

// PRSVariantAudit 记录每个模型位点的使用情况。
type PRSVariantAudit struct {
	ID        string
	Used      bool
	Sign      float64
	Relation  variant.FlipRelation
	MissingGT int // 该位点缺失基因型的样本数
	Reason    variant.Reason
}

// PRSResult 是个体级评分结果。
type PRSResult struct {
	Name              string
	SampleIDs         []string
	Scores            []float64
	NActive           []int
	Audit             []PRSVariantAudit
	MissingModelSites []string
	Scale             string
	Diag              string
}

// prsError 是评分层的可区分错误。
type prsSite struct {
	gtI  int
	beta float64 // 已对齐到矩阵方向
	sign float64
	rel  variant.FlipRelation
	id   string
}

type prsError struct{ what, detail string }

func (e *prsError) Error() string { return "prs: " + e.what + ": " + e.detail }

// BuildPRS 按汇总统计效应与本地基因型矩阵构建个体风险评分。
//
// 等位基因方向以基因型矩阵（个体数据）为参考：模型 beta 对齐到矩阵的
// EA/NEA。链模糊位点在默认策略下被拒绝并记录，绝不通过“按符号相加”
// 猜测效应方向。位点集合（含/缺哪些位点）、尺度与缺失策略共同决定
// 评分口径，全部在 Audit/Diag 中可追溯。
func BuildPRS(ctx context.Context, model *ingest.SummaryTable, gt *ingest.GenotypeMatrix, opts PRSOptions) (*PRSResult, error) {
	if model == nil || gt == nil {
		return nil, &prsError{"nil_input", "model and genotype matrix required"}
	}
	if opts.Backend == numwrap.BackendCrosscheck {
		return crosscheckRun("prs",
			func(o PRSOptions, b numwrap.Backend) PRSOptions { o.Backend = b; return o },
			comparePRSResults,
			func(o PRSOptions) (*PRSResult, error) { return BuildPRS(ctx, model, gt, o) },
			opts)
	}
	if opts.Align == (variant.AlignConfig{}) {
		opts.Align = variant.DefaultAlignConfig()
	}
	scale := opts.ScoreScale
	if scale == "" {
		scale = "raw"
	}
	if scale != "raw" && scale != "std" && scale != "perallele" {
		return nil, &prsError{"unsupported_scale", "scale must be raw|std|perallele, got: " + scale}
	}
	if opts.Name == "" {
		opts.Name = "prs"
	}

	// 模型侧索引（重复键已在解析层拒绝）。
	modelByKey := make(map[variant.MatchKey]*ingest.SummaryRecord, len(model.Records))
	for i := range model.Records {
		modelByKey[model.Records[i].V.Key()] = &model.Records[i]
	}
	gtByKey := make(map[variant.MatchKey]int, len(gt.Variants))
	for i, v := range gt.Variants {
		gtByKey[v.Key()] = i
	}

	nS := gt.NSamples()
	scores := make([]float64, nS)
	active := make([]int, nS)
	audit := make([]PRSVariantAudit, 0, len(model.Records))
	var missingModel []string

	var sites []prsSite

	// 口径以“效应模型请求的位点集合”为基准遍历（模型顺序，确定性）。
	for mi := range model.Records {
		rec := &model.Records[mi]
		a := PRSVariantAudit{ID: rec.V.ID, Sign: 1}
		if opts.PThreshold > 0 {
			if math.IsNaN(rec.P) {
				a.Reason = variant.Reason("pvalue_unavailable")
				audit = append(audit, a)
				continue
			}
			if rec.P > opts.PThreshold {
				a.Reason = variant.Reason("pvalue_above_threshold")
				audit = append(audit, a)
				continue
			}
		}
		gi, ok := gtByKey[rec.V.Key()]
		if !ok {
			missingModel = append(missingModel, rec.V.ID)
			a.Reason = variant.Reason("effect_site_missing")
			audit = append(audit, a)
			continue
		}
		gv := gt.Variants[gi]
		al, err := variant.Align(gv, rec.V, opts.Align)
		if err != nil {
			r, _ := variant.AsReject(err)
			if r != nil {
				a.Reason = r.Reason
			} else {
				a.Reason = variant.ReasonAlleleMismatch
			}
			audit = append(audit, a)
			continue
		}
		a.Used = true
		a.Sign = al.Sign
		a.Relation = al.Relation
		sites = append(sites, prsSite{gtI: gi, beta: rec.Beta * al.Sign, sign: al.Sign, rel: al.Relation, id: rec.V.ID})
		audit = append(audit, a)
	}

	if len(missingModel) > 0 && opts.Model == ModelStrict {
		return nil, &prsError{"effect_sites_missing",
			fmt.Sprintf("%d requested effect site(s) absent in genotype matrix; first=%s; set ModelIgnore to score on intersection explicitly",
				len(missingModel), missingModel[0])}
	}
	if len(sites) == 0 {
		return nil, &prsError{"no_usable_sites", "no model sites aligned to genotype matrix"}
	}

	// 矩阵中有、模型中没有的位点也显式列出（集合差异必须对称可见）。
	for _, gv := range gt.Variants {
		if _, ok := modelByKey[gv.Key()]; !ok {
			missingModel = append(missingModel, gv.ID+"(gt-only)")
		}
	}

	// 逐样本累加。样本间相互独立可并行，但为保证“并行=串行”，
	// 统一在固定的 (site, sample) 字典序下做补偿求和；
	// Workers 仅用于把样本切分到 mapShards，归并顺序固定。
	compute := func(lo, hi int) ([]float64, []int, error) {
		sc := make([]float64, hi-lo)
		ac := make([]int, hi-lo)
		for _, st := range sites {
			for jLocal := lo; jLocal < hi; jLocal++ {
				dose, miss := gt.At(st.gtI, jLocal)
				k := jLocal - lo
				if miss {
					switch opts.GT {
					case GTMeanImpute:
						eaf := gt.Variants[st.gtI].EAF
						imp := 1.0
						if !math.IsNaN(eaf) {
							imp = 2 * eaf
						}
						sc[k] += st.beta * imp
						ac[k]++
					case GTExclude:
						// 不贡献
					}
					continue
				}
				sc[k] += st.beta * dose
				ac[k]++
			}
		}
		return sc, ac, nil
	}

	// 串行/并行走同一实现，仅外层分片；分片间不共享可写状态。
	workers := normalizeWorkers(opts.Workers, nS)
	if workers == 1 {
		sc, ac, err := compute(0, nS)
		if err != nil {
			return nil, err
		}
		scores, active = sc, ac
	} else {
		blocks, err := mapShards(ctx, nS, workers, func(start, end int) (scoreBlock, error) {
			sc, ac, err := compute(start, end)
			if err != nil {
				return scoreBlock{}, err
			}
			return scoreBlock{lo: start, hi: end, scores: sc, active: ac}, nil
		})
		if err != nil {
			return nil, err
		}
		for _, b := range blocks {
			copy(scores[b.lo:b.hi], b.scores)
			copy(active[b.lo:b.hi], b.active)
		}
	}
	// 位点缺失计数不参与评分数值；在选中位点上直接统计一次。
	for _, st := range sites {
		cnt := 0
		for j := 0; j < nS; j++ {
			if _, miss := gt.At(st.gtI, j); miss {
				cnt++
			}
		}
		audit = markAuditMissing(audit, st.id, cnt)
	}

	// 尺度口径（顺序确定）。
	switch scale {
	case "perallele":
		for j := 0; j < nS; j++ {
			if active[j] > 0 {
				scores[j] /= float64(active[j])
			}
		}
	case "std":
		mean, sd, err := numwrap.MeanStdBackend(scores, opts.Backend)
		if err != nil {
			return nil, &prsError{"standardize_failed", err.Error()}
		}
		if sd == 0 {
			return nil, &prsError{"zero_variance_scores", "cannot standardize zero-variance scores"}
		}
		for j := range scores {
			scores[j] = (scores[j] - mean) / sd
		}
	}

	nSignFlip := 0
	for _, a := range audit {
		if a.Used && a.Sign < 0 {
			nSignFlip++
		}
	}
	res := &PRSResult{
		Name: opts.Name, SampleIDs: append([]string(nil), gt.SampleIDs...),
		Scores: scores, NActive: active, Audit: audit,
		MissingModelSites: missingModel, Scale: scale,
	}
	res.Diag = fmt.Sprintf("sites_requested=%d sites_used=%d missing_sites=%d sign_flips=%d scale=%s gt_policy=%d samples=%d",
		len(model.Records), len(sites), len(missingModel), nSignFlip, scale, opts.GT, nS)
	return res, nil
}

type scoreBlock struct {
	lo, hi int
	scores []float64
	active []int
}

func markAuditMissing(a []PRSVariantAudit, id string, n int) []PRSVariantAudit {
	for i := range a {
		if a[i].ID == id {
			a[i].MissingGT = n
		}
	}
	return a
}

// EvalOptions 控制评分评估。
type EvalOptions struct {
	Phenotype  []float64   // 长度必须等于样本数，不允许 NaN/Inf
	Covariates [][]float64 // 预留：非空即显式拒绝
	Backend    numwrap.Backend
}

// EvalResult 是评分-表型关联评估（简单线性回归 y ~ score）。
type EvalResult struct {
	Beta float64
	SE   float64
	P    float64
	R2   float64
	N    int
}

// EvaluatePRS 用简单线性回归评估评分与表型的关联；
// 确定性：仅使用 numwrap 中的顺序累加与 WLS。
func EvaluatePRS(ctx context.Context, prs *PRSResult, opts EvalOptions) (*EvalResult, error) {
	if prs == nil {
		return nil, &prsError{"nil_prs", "nil PRSResult"}
	}
	n := len(prs.Scores)
	if len(opts.Covariates) > 0 {
		return nil, &prsError{"unsupported_option", "covariate adjustment not supported by current backend"}
	}
	if len(opts.Phenotype) != n {
		return nil, &prsError{"dimension_mismatch",
			fmt.Sprintf("phenotype length %d != samples %d", len(opts.Phenotype), n)}
	}
	for i := 0; i < n; i++ {
		if nonFinite(prs.Scores[i]) || nonFinite(opts.Phenotype[i]) {
			return nil, &prsError{"non_finite_value", fmt.Sprintf("at index %d", i)}
		}
	}
	fit, err := numwrap.WLS(prs.Scores, opts.Phenotype, numwrap.RegressionOptions{
		InterceptFree: true, Backend: opts.Backend})
	if err != nil {
		return nil, &prsError{"regression_failed", err.Error()}
	}
	beta := fit.Coefficients[1]
	se := fit.StdErr[1]
	z := beta / se
	p := 2 * normalSurvival(abs(z))

	// R2 = corr^2，用与 WLS 一致的确定性均值/方差。
	xm, xs, err := numwrap.MeanStdBackend(prs.Scores, opts.Backend)
	if err != nil {
		return nil, &prsError{"moments_failed", err.Error()}
	}
	ym, ys, err := numwrap.MeanStdBackend(opts.Phenotype, opts.Backend)
	if err != nil {
		return nil, &prsError{"moments_failed", err.Error()}
	}
	var sxy, comp float64
	for i := 0; i < n; i++ {
		v := (prs.Scores[i] - xm) * (opts.Phenotype[i] - ym)
		yv := v - comp
		t := sxy + yv
		comp = (t - sxy) - yv
		sxy = t
	}
	r2 := 0.0
	if xs > 0 && ys > 0 {
		corr := (sxy / float64(n-1)) / (xs * ys)
		r2 = corr * corr
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	return &EvalResult{Beta: beta, SE: se, P: p, R2: r2, N: n}, nil
}

func nonFinite(v float64) bool { return math.IsNaN(v) || math.IsInf(v, 0) }
func abs(v float64) float64 {
	if v < 0 {
		return -v
	}
	return v
}
