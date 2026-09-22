// Package simulate 生成确定性的本地测试数据：基因型矩阵、LD 分数与
// 多队列汇总统计。所有随机数来自固定种子的可复现流，便于端到端验证
// 链翻转、模糊位点、缺失字段与队列重叠的处理。
package simulate

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// Spec 描述一次合成数据集。
type Spec struct {
	Seed      int64
	NVariants int
	NSamples  int
	// FractionAmbiguous：A/T、C/G 位点比例（固定从第 0 个起按步长放置，
	// 不使用随机抽样，保证可精确断言）。
	FractionAmbiguous float64
	// FractionFlip / FractionSwap：汇总统计相对基因型参考的
	// 链翻转 / 等位基因互换比例（确定性交错）。
	FractionFlip float64
	FractionSwap float64
	// FractionMissingGT：每个样本每位点的缺失概率（按确定性 RNG）。
	FractionMissingGT float64
	Cohorts           []CohortSpec
	// DuplicateVariants：在每个汇总文件末尾追加的重复行数。
	DuplicateVariants int
	// InvalidRows：追加的非法行数（非法等位基因/数值）。
	InvalidRows int
}

// CohortSpec 描述合成队列。
type CohortSpec struct {
	Name       string
	N          int
	SharedWith string
	SharedN    int
	Duplicate  bool // 与 SharedWith 队列样本完全一致
}

// Dataset 是生成结果在内存中的表示。
type Dataset struct {
	SummaryFiles map[string][]byte
	LDScoreFile  []byte
	GenotypeFile []byte
	SampleIDs    []string
	Truth        Truth
}

// Truth 记录可断言的生成真值。
type Truth struct {
	H2           float64
	AmbiguousIDs []string
	FlippedIDs   []string
	SwappedIDs   []string
	OverlapPairs [][2]string
}

// rng 是确定性 xorshift64* PRNG（不依赖全局 rand，结果跨机器一致）。
type rng struct{ s uint64 }

func newRNG(seed int64) *rng {
	s := uint64(seed)
	if s == 0 {
		s = 0x9E3779B97F4A7C15
	}
	return &rng{s: s}
}

func (r *rng) next() uint64 {
	r.s ^= r.s >> 12
	r.s ^= r.s << 25
	r.s ^= r.s >> 27
	return r.s * 0x2545F4914F6CDD1D
}

func (r *rng) float64() float64 { return float64(r.next()>>11) / float64(uint64(1)<<53) }

func (r *rng) intn(n int) int {
	if n <= 0 {
		return 0
	}
	return int(r.next() % uint64(n))
}

// 非模糊位点等位基因组合（两碱基不同且非 A/T、C/G）。
var nonAmbiguous = [][2]string{
	{"A", "G"}, {"A", "C"}, {"C", "T"}, {"G", "T"}, {"C", "G"}, {"A", "T"} /* 末两个用于模糊 */}
var ambiguousPairs = [][2]string{{"A", "T"}, {"C", "G"}}

// Generate 按 Spec 确定性地生成完整数据集。
func Generate(spec Spec) (*Dataset, error) {
	if spec.NVariants <= 0 || spec.NSamples <= 0 {
		return nil, fmt.Errorf("simulate: NVariants/NSamples must be positive")
	}
	if len(spec.Cohorts) == 0 {
		return nil, fmt.Errorf("simulate: at least one cohort required")
	}
	r := newRNG(spec.Seed)

	type site struct {
		id        string
		ea, nea   string
		ambiguous bool
		l2        float64
	}
	sites := make([]site, spec.NVariants)
	nAmb := int(math.Round(float64(spec.NVariants) * spec.FractionAmbiguous))
	ambStep := 0
	if nAmb > 0 {
		ambStep = spec.NVariants / nAmb
		if ambStep < 1 {
			ambStep = 1
		}
	}
	truth := Truth{H2: 0.3}
	for i := 0; i < spec.NVariants; i++ {
		id := fmt.Sprintf("rs%06d", i+1)
		amb := ambStep > 0 && i%ambStep == 0 && len(truth.AmbiguousIDs) < nAmb
		var ea, nea string
		if amb {
			p := ambiguousPairs[i%2]
			ea, nea = p[0], p[1]
			truth.AmbiguousIDs = append(truth.AmbiguousIDs, id)
		} else {
			p := nonAmbiguous[r.intn(4)] // 只用前 4 个确定非模糊组合
			ea, nea = p[0], p[1]
		}
		sites[i] = site{id: id, ea: ea, nea: nea, ambiguous: amb,
			l2: 1 + float64(i%7) + r.float64()}
	}

	// 样本 ID
	sampleIDs := make([]string, spec.NSamples)
	for i := range sampleIDs {
		sampleIDs[i] = fmt.Sprintf("ind%06d", i+1)
	}

	// 基因型矩阵：剂量 0/1/2 + 确定性缺失
	var gt strings.Builder
	gt.WriteString("SNP\tCHR\tBP\tA1\tA2\t" + strings.Join(sampleIDs, "\t") + "\n")
	for i, st := range sites {
		row := []string{st.id, "1", fmt.Sprintf("%d", 1000+i), st.ea, st.nea}
		for j := 0; j < spec.NSamples; j++ {
			if r.float64() < spec.FractionMissingGT {
				row = append(row, "NA")
			} else {
				row = append(row, fmt.Sprintf("%d", r.intn(3)))
			}
		}
		gt.WriteString(strings.Join(row, "\t") + "\n")
	}

	// LD 分数
	var ld strings.Builder
	ld.WriteString("SNP\tCHR\tBP\tA1\tA2\tL2\tN\n")
	for i, st := range sites {
		fmt.Fprintf(&ld, "%s\t1\t%d\t%s\t%s\t%.6f\t%d\n",
			st.id, 1000+i, st.ea, st.nea, st.l2, spec.NSamples)
	}

	// 队列样本集合（确定性重叠）
	cohortSamples := map[string][]string{}
	for _, c := range spec.Cohorts {
		if c.Duplicate {
			base, ok := cohortSamples[c.SharedWith]
			if !ok {
				return nil, fmt.Errorf("simulate: duplicate cohort references unknown cohort %q", c.SharedWith)
			}
			cohortSamples[c.Name] = append([]string(nil), base...)
			continue
		}
		ids := make([]string, 0, c.N)
		if c.SharedN > 0 && c.SharedWith != "" {
			base, ok := cohortSamples[c.SharedWith]
			if !ok {
				return nil, fmt.Errorf("simulate: shared cohort %q unknown", c.SharedWith)
			}
			if c.SharedN > len(base) {
				return nil, fmt.Errorf("simulate: SharedN %d > cohort size %d", c.SharedN, len(base))
			}
			ids = append(ids, base[:c.SharedN]...)
		}
		// 其余样本从全局样本池确定性选取未用 ID
		used := map[string]bool{}
		for _, id := range ids {
			used[id] = true
		}
		next := 0
		for len(ids) < c.N && next < len(sampleIDs) {
			if !used[sampleIDs[next]] {
				ids = append(ids, sampleIDs[next])
				used[sampleIDs[next]] = true
			}
			next++
		}
		if len(ids) < c.N {
			return nil, fmt.Errorf("simulate: not enough unique samples for cohort %s (want %d, pool %d)",
				c.Name, c.N, len(sampleIDs))
		}
		cohortSamples[c.Name] = ids
	}
	// 重叠真值（两两交集）
	names := make([]string, 0, len(spec.Cohorts))
	for _, c := range spec.Cohorts {
		names = append(names, c.Name)
	}
	sort.Strings(names)
	for i := 0; i < len(names); i++ {
		for j := i + 1; j < len(names); j++ {
			a, b := cohortSamples[names[i]], cohortSamples[names[j]]
			set := map[string]bool{}
			for _, x := range a {
				set[x] = true
			}
			shared := 0
			for _, x := range b {
				if set[x] {
					shared++
				}
			}
			if shared > 0 {
				truth.OverlapPairs = append(truth.OverlapPairs, [2]string{names[i], names[j]})
			}
		}
	}

	// 汇总统计：效应与方向扰动（确定性交错）。
	//
	// 效应量按 LDSC 总体关系生成，使数值物理自洽：
	//   chi2 = 1 + (N*h2/M_eff)*l2 + 确定性扰动；beta = z*se
	// 模糊位点另行做互换扰动，仅用于对齐层拒绝测试（不纳入 h2 真值）。
	h2 := truth.H2
	mEff := float64(spec.NVariants - len(truth.AmbiguousIDs))
	flipEvery := 0
	if spec.FractionFlip > 0 {
		flipEvery = int(math.Round(1 / spec.FractionFlip))
	}
	swapEvery := 0
	if spec.FractionSwap > 0 {
		swapEvery = int(math.Round(1 / spec.FractionSwap))
	}
	files := map[string][]byte{}
	for ci, cspec := range spec.Cohorts {
		var b strings.Builder
		b.WriteString("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n")
		for i, st := range sites {
			if st.ambiguous {
				// 模糊位点在汇总侧做互换扰动，供对齐层拒绝测试
				ea, nea := st.nea, st.ea
				beta := 0.05 + 0.001*float64(i%5)
				fmt.Fprintf(&b, "%s\t1\t%d\t%s\t%s\t%.5f\t0.05\t1e-3\t%d\n",
					st.id, 1000+i, ea, nea, beta, cspec.N)
				continue
			}
			slope := float64(cspec.N) * h2 / mEff
			// 确定性“噪声”：由 i、ci 构成的有界交错扰动（无随机，可复现）。
			noise := 0.10 * math.Sin(float64(i*7+ci*13))
			chi2 := 1 + slope*st.l2 + noise
			if chi2 < 0.05 {
				chi2 = 0.05
			}
			z := math.Sqrt(chi2)
			se := 0.05
			beta := z * se
			ea, nea := st.ea, st.nea
			flipped, swapped := false, false
			if flipEvery > 0 && i%flipEvery == flipEvery/2 {
				ea = comp(ea)
				nea = comp(nea)
				flipped = true
			}
			if swapEvery > 0 && i%swapEvery == 0 {
				// 互换：仅改等位基因命名与 beta 符号（chi2 不变）
				ea, nea = nea, ea
				beta = -beta
				swapped = true
			}
			if flipped {
				truth.FlippedIDs = appendOnce(truth.FlippedIDs, st.id)
			}
			if swapped {
				truth.SwappedIDs = appendOnce(truth.SwappedIDs, st.id)
			}
			pv := 0.5 * math.Erfc(z/math.Sqrt2)
			fmt.Fprintf(&b, "%s\t1\t%d\t%s\t%s\t%.6f\t%.3f\t%.3g\t%d\n",
				st.id, 1000+i, ea, nea, beta, se, pv, cspec.N)
		}
		// 注入重复行（重复第一行）
		if spec.DuplicateVariants > 0 && len(sites) > 0 {
			st := sites[0]
			for k := 0; k < spec.DuplicateVariants; k++ {
				fmt.Fprintf(&b, "%s\t1\t%d\t%s\t%s\t0.05\t0.02\t1e-4\t%d\n",
					st.id, 1000, st.ea, st.nea, cspec.N)
			}
		}
		// 注入非法行
		for k := 0; k < spec.InvalidRows; k++ {
			switch k % 3 {
			case 0:
				fmt.Fprintf(&b, "rsBAD%d\t1\t%d\tI\tG\t0.05\t0.02\t1e-4\t%d\n", k, 9000+k, cspec.N)
			case 1:
				fmt.Fprintf(&b, "rsBAD%d\t1\t%d\tA\tG\tNaN\t0.02\t1e-4\t%d\n", k, 9000+k, cspec.N)
			case 2:
				fmt.Fprintf(&b, "rsBAD%d\t1\t%d\tA\tA\t0.05\t0.02\t1e-4\t%d\n", k, 9000+k, cspec.N)
			}
		}
		files[cspec.Name+"_sumstats.tsv"] = []byte(b.String())
	}

	return &Dataset{
		SummaryFiles: files,
		LDScoreFile:  []byte(ld.String()),
		GenotypeFile: []byte(gt.String()),
		SampleIDs:    sampleIDs,
		Truth:        truth,
	}, nil
}

func comp(s string) string {
	r := []byte(s)
	m := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	for i := range r {
		r[i] = m[r[i]]
	}
	return string(r)
}

func appendOnce(xs []string, x string) []string {
	for _, v := range xs {
		if v == x {
			return xs
		}
	}
	return append(xs, x)
}

// WriteToDir 将 Dataset 文件写入 dir（目录不存在则创建）。
func (d *Dataset) WriteToDir(dir string) (map[string]string, error) {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}
	paths := map[string]string{}
	write := func(name string, data []byte) error {
		p := filepath.Join(dir, name)
		paths[name] = p
		return os.WriteFile(p, data, 0o644)
	}
	if err := write("ldscores.tsv", d.LDScoreFile); err != nil {
		return nil, err
	}
	if err := write("genotypes.tsv", d.GenotypeFile); err != nil {
		return nil, err
	}
	names := make([]string, 0, len(d.SummaryFiles))
	for n := range d.SummaryFiles {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, n := range names {
		if err := write(n, d.SummaryFiles[n]); err != nil {
			return nil, err
		}
	}
	return paths, nil
}
