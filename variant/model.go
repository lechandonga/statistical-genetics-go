// Package variant 定义跨输入源（汇总统计、LD 分数、基因型矩阵）的统一变体模型。
//
// 设计约定见 docs/model.md。
package variant

import (
	"math"
	"strconv"
	"strings"
)

// Allele 表示经过归一化（大写）的单个等位基因字符串，例如 "A"、"ACGT"（indel）。
// 仅允许 A/C/G/T 字符；其它写法（"I"、"D"、缺失码）在解析层即被拒绝。
type Allele string

// Strand 表示 DNA 链方向。
type Strand uint8

const (
	StrandUnknown Strand = iota // 未给出链信息
	StrandForward               // 正链（参考链）
	StrandReverse               // 反链（互补链）
)

// Variant 是所有输入源对齐到的统一变体模型。
//
// 仅当 CHROM/POS 同时给出且合法时 Chrom/Pos 才有效；
// 仅有 SNPID + 等位基因时按 ID 与等位基因匹配。
type Variant struct {
	ID     string   // 变体标识（rsID 或来源主键），归一化后大写
	Chrom  string   // 染色体标签，归一化后如 "1".."22","X","Y","MT"
	Pos    int64    // 1-based 物理位置；未知时为 0
	EA     Allele   // effect allele（效应等位基因）
	NEA    Allele   // non-effect allele（非效应等位基因）
	Strand Strand   // 来源声明的链；未知为 StrandUnknown
	EAF    float64  // 效应等位基因频率；未知为 NaN
	Runes  []string // 规范化来源字段（审计用）：id, chrom, pos, ea, nea
}

// MatchKey 是去重/匹配使用的稳定键：优先 ID，否则 CHROM:POS。
type MatchKey string

var validBases = [256]bool{'A': true, 'C': true, 'G': true, 'T': true}

// NormalizeAllele 将原始等位基因字符串归一化：去除首尾空白、统一大写。
// 空串、含非 A/C/G/T 字符（含内部空白、I/D 缺口码）返回 ErrInvalidAllele。
func NormalizeAllele(raw string) (Allele, error) {
	s := strings.ToUpper(strings.TrimSpace(raw))
	if s == "" {
		return "", &RejectError{Reason: ReasonMissingAlleles, Detail: "empty allele"}
	}
	for i := 0; i < len(s); i++ {
		b := s[i]
		if !validBases[b] {
			return "", &RejectError{Reason: ReasonInvalidAllele,
				Detail: "allele contains non-ACGT character: " + raw}
		}
	}
	return Allele(s), nil
}

// NormalizeChrom 归一化染色体标签：
//
//	"chr1"/"1" -> "1"；"chrX" -> "X"；"chrM" -> "MT"；
//	"23"/"24"/"25" 分别映射为 X/Y/MT（部分旧式汇总统计约定）。
//
// 无法识别返回空串与 ReasonInvalidChromosome。
func NormalizeChrom(raw string) (string, error) {
	s := strings.ToUpper(strings.TrimSpace(raw))
	s = strings.TrimPrefix(s, "CHR")
	switch s {
	case "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
		"13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
		"X", "Y":
		return s, nil
	case "M", "MT":
		return "MT", nil
	case "23":
		return "X", nil
	case "24":
		return "Y", nil
	case "25":
		return "MT", nil
	default:
		return "", &RejectError{Reason: ReasonInvalidChrom,
			Detail: "unrecognized chromosome: " + raw}
	}
}

// ValidChromosomes 是本流程接受的染色体标签集合（归一化后）。
var ValidChromosomes = func() map[string]struct{} {
	m := map[string]struct{}{}
	for c := 1; c <= 22; c++ {
		m[strconv.Itoa(c)] = struct{}{}
	}
	for _, c := range []string{"X", "Y", "MT"} {
		m[c] = struct{}{}
	}
	return m
}()

// New 从原始字段构造归一化 Variant。
//
// 规则：
//   - ea/nea 必填且须为合法 ACGT 串，且二者不可相同；
//   - id 与 (chrom,pos) 至少提供其一；
//   - chrom 与 pos 必须同时出现，chrom 须可识别，pos >= 1；
//   - eaf 若非 NaN 必须落在 [0,1]；
//   - 若来源声明为反链，则在构造时立即翻转到参考链（见 docs/model.md），
//     调用方拿到的 Variant 永远是参考链方向。
func New(id, chrom string, pos int64, eaRaw, neaRaw string, eaf float64, strand Strand) (*Variant, error) {
	idN := strings.ToUpper(strings.TrimSpace(id))

	ea, err := NormalizeAllele(eaRaw)
	if err != nil {
		return nil, withID(idN, err)
	}
	nea, err := NormalizeAllele(neaRaw)
	if err != nil {
		return nil, withID(idN, err)
	}
	if ea == nea {
		return nil, &RejectError{ID: idN, Reason: ReasonInvalidAllele,
			Detail: "effect and non-effect alleles are identical: " + string(ea)}
	}

	var chromN string
	var posN int64
	chromRaw := strings.TrimSpace(chrom)
	switch {
	case chromRaw == "" && pos == 0:
		// 两者皆缺：合法，但此时必须有 id
	case chromRaw == "" && pos != 0:
		return nil, &RejectError{ID: idN, Reason: ReasonInvalidChrom,
			Detail: "position given without chromosome"}
	case chromRaw != "" && pos <= 0:
		return nil, &RejectError{ID: idN, Reason: ReasonInvalidPosition,
			Detail: "chromosome given without valid 1-based position"}
	default:
		c, err := NormalizeChrom(chromRaw)
		if err != nil {
			return nil, withID(idN, err)
		}
		chromN, posN = c, pos
	}

	if idN == "" && chromN == "" {
		return nil, &RejectError{Reason: ReasonMissingID,
			Detail: "variant requires ID or chromosome+position"}
	}

	if !math.IsNaN(eaf) {
		if math.IsInf(eaf, 0) || eaf < 0 || eaf > 1 {
			return nil, &RejectError{ID: idN, Reason: ReasonInvalidFrequency,
				Detail: "EAF out of [0,1] or non-finite"}
		}
	}

	// 来源声明反链：翻转到参考链。EA/NEA 同时互补，方向约定保持
	// （效应仍以参考链上的 EA 定义）。
	if strand == StrandReverse {
		ea = flipAllele(ea)
		nea = flipAllele(nea)
		strand = StrandForward
	}

	return &Variant{
		ID:     idN,
		Chrom:  chromN,
		Pos:    posN,
		EA:     ea,
		NEA:    nea,
		Strand: strand,
		EAF:    eaf,
		Runes:  []string{id, chrom, chromRawPos(pos), eaRaw, neaRaw},
	}, nil
}

func chromRawPos(pos int64) string {
	if pos == 0 {
		return ""
	}
	// 审计字段保留原始数字文本
	const digits = "0123456789"
	if pos < 0 {
		return ""
	}
	var b [20]byte
	i := len(b)
	for pos > 0 {
		i--
		b[i] = digits[pos%10]
		pos /= 10
	}
	return string(b[i:])
}

func withID(id string, err error) error {
	if r, ok := AsReject(err); ok && r.ID == "" {
		cp := *r
		cp.ID = id
		return &cp
	}
	return err
}

// flipAllele 返回等位基因的链互补；调用前须已确认全部字符为 ACGT。
func flipAllele(a Allele) Allele {
	b := []byte(a)
	for i := range b {
		if c, ok := ComplementBase(b[i]); ok {
			b[i] = c
		}
	}
	return Allele(b)
}

// Key 返回跨来源匹配的稳定键。
func (v *Variant) Key() MatchKey {
	if v.ID != "" {
		return MatchKey(v.ID)
	}
	return MatchKey(v.Chrom + ":" + chromRawPos(v.Pos))
}

// IsSNV 报告该变体是否为单碱基替换（两等位基因均为单个 A/C/G/T）。
func (v *Variant) IsSNV() bool { return len(v.EA) == 1 && len(v.NEA) == 1 }

// IsAmbiguous 报告该位点是否为 A/T 或 C/G 型链方向模糊位点。
//
// 这类 SNV 互补后与互换不可区分：A/T 互补为 T/A，
// 无法判断来源是“链翻转”还是“等位基因互换”，因此效应符号不可确定。
func (v *Variant) IsAmbiguous() bool {
	if !v.IsSNV() {
		return false
	}
	ea, nea := v.EA[0], v.NEA[0]
	return (ea == 'A' && nea == 'T') || (ea == 'T' && nea == 'A') ||
		(ea == 'C' && nea == 'G') || (ea == 'G' && nea == 'C')
}
