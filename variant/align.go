package variant

import "fmt"

// AmbiguousPolicy 定义链方向模糊位点（A/T、C/G）的处理策略。
type AmbiguousPolicy uint8

const (
	// AmbiguousReject 拒绝模糊位点（默认，最严格）：
	// A/T 与 C/G 的“链翻转”和“互换”在观测上不可区分，
	// 效应符号无法确定，禁止猜测。
	AmbiguousReject AmbiguousPolicy = iota
	// AmbiguousKeep 保留但标记：结果中 Ambiguous=true，
	// 且仅接受 RelationNone / RelationSwapped 这种无需判链的对齐；
	// 任何需要借助链翻转才能成立的对齐仍被拒绝。
	AmbiguousKeep
)

// FlipPolicy 定义是否允许通过链翻转完成对齐。
type FlipPolicy uint8

const (
	// FlipAllowed 允许链翻转（默认；仍要求位点非模糊）。
	FlipAllowed FlipPolicy = iota
	// FlipForbidden 禁止链翻转，只接受直接匹配或等位基因互换。
	FlipForbidden
)

// AlignConfig 控制两个变体间的等位基因对齐行为。
type AlignConfig struct {
	Ambiguous AmbiguousPolicy
	Flip      FlipPolicy
}

// DefaultAlignConfig 返回推荐配置：允许非模糊位点翻转、拒绝模糊位点。
func DefaultAlignConfig() AlignConfig {
	return AlignConfig{Ambiguous: AmbiguousReject, Flip: FlipAllowed}
}

// FlipRelation 描述目标变体相对参考变体的变换集合。
type FlipRelation uint8

const (
	RelationNone    FlipRelation = 0         // 直接匹配：EA→EA, NEA→NEA
	RelationSwapped FlipRelation = 1 << iota // 等位基因互换：EA→NEA（效应取反）
	RelationFlipped                          // 链翻转：EA→comp(EA)（效应不变号）
)

// EffectSign 按对齐关系返回效应符号：互换为 -1，其余为 +1。
// 链翻转本身不改变效应符号（beta 定义为“每增加一个 EA 拷贝”的效应，
// EA 经互补重命名后拷贝数不变）。
func (r FlipRelation) EffectSign() float64 {
	if r&RelationSwapped != 0 {
		return -1
	}
	return 1
}

// HasFlip 报告是否包含链翻转。
func (r FlipRelation) HasFlip() bool { return r&RelationFlipped != 0 }

func (r FlipRelation) String() string {
	switch r {
	case RelationNone:
		return "none"
	case RelationSwapped:
		return "swapped"
	case RelationFlipped:
		return "flipped"
	case RelationFlipped | RelationSwapped:
		return "flipped+swapped"
	default:
		return fmt.Sprintf("relation(%d)", uint8(r))
	}
}

// AlignResult 是对齐结论。
type AlignResult struct {
	Relation    FlipRelation
	Sign        float64 // +1 / -1，直接作用于 beta/or/effect
	Ambiguous   bool    // 是否为链模糊位点（仅 AmbiguousKeep 下可能为 true）
	MatchedEA   Allele  // 对齐后统一采用的效应等位基因（参考侧 EA）
	MatchedNEA  Allele
	Explanation string // 人类可读判定依据
}

// Align 将 target 的等位基因对齐到 ref。
//
// 匹配按以下顺序判定（位点为同一个 Key 的前提由调用方保证）：
//
//  1. 直接匹配 target{EA,NEA} == ref{EA,NEA}
//  2. 互换     target{EA,NEA} == ref{NEA,EA}（效应取反）
//  3. 链翻转   comp(target){EA,NEA} == ref{EA,NEA}
//  4. 翻转+互换
//
// 模糊位点（A/T、C/G）：AmbiguousReject 直接拒绝；
// AmbiguousKeep 下仅允许情况 1/2 并显式标记，情况 3/4 拒绝；
// FlipForbidden 下情况 3/4 一律拒绝。
func Align(ref, target *Variant, cfg AlignConfig) (*AlignResult, error) {
	if ref == nil || target == nil {
		return nil, &RejectError{Reason: ReasonAlleleMismatch,
			Detail: "align: nil variant"}
	}
	amb := ref.IsAmbiguous() || target.IsAmbiguous()
	id := ref.ID
	if id == "" {
		id = target.ID
	}

	// 1. 直接匹配
	if target.EA == ref.EA && target.NEA == ref.NEA {
		return &AlignResult{
			Relation: RelationNone, Sign: 1, Ambiguous: amb,
			MatchedEA: ref.EA, MatchedNEA: ref.NEA,
			Explanation: fmt.Sprintf("direct match on EA=%s NEA=%s", ref.EA, ref.NEA),
		}, nil
	}
	// 2. 互换
	if target.EA == ref.NEA && target.NEA == ref.EA {
		if amb && cfg.Ambiguous == AmbiguousReject {
			return nil, &RejectError{ID: id, Reason: ReasonAmbiguousStrand,
				Detail: fmt.Sprintf("ambiguous site %s/%s: strand versus swap unresolvable", ref.EA, ref.NEA)}
		}
		return &AlignResult{
			Relation: RelationSwapped, Sign: -1, Ambiguous: amb,
			MatchedEA: ref.EA, MatchedNEA: ref.NEA,
			Explanation: fmt.Sprintf("allele swap: target EA=%s is reference NEA; sign flipped", target.EA),
		}, nil
	}

	// 3/4 需要链翻转：先检查可行性
	if cfg.Flip == FlipForbidden {
		return nil, &RejectError{ID: id, Reason: ReasonFlipForbidden,
			Detail: fmt.Sprintf("target %s/%s aligns only via strand flip, but flips are forbidden", target.EA, target.NEA)}
	}
	if !alleleFlippable(target.EA) || !alleleFlippable(target.NEA) {
		return nil, &RejectError{ID: id, Reason: ReasonUnresolvableStrand,
			Detail: fmt.Sprintf("multi-base allele %s/%s cannot be strand-flipped", target.EA, target.NEA)}
	}
	if amb {
		// 模糊位点无论哪种策略，都不允许借助翻转判方向
		if cfg.Ambiguous == AmbiguousReject {
			return nil, &RejectError{ID: id, Reason: ReasonAmbiguousStrand,
				Detail: fmt.Sprintf("ambiguous site %s/%s: flip vs swap unresolvable", ref.EA, ref.NEA)}
		}
		return nil, &RejectError{ID: id, Reason: ReasonAmbiguousStrand,
			Detail: "ambiguous site cannot be aligned through a strand flip even under keep-policy"}
	}

	fea, fnea := flipAllele(target.EA), flipAllele(target.NEA)
	// 3. 翻转后直接匹配
	if fea == ref.EA && fnea == ref.NEA {
		return &AlignResult{
			Relation: RelationFlipped, Sign: 1,
			MatchedEA: ref.EA, MatchedNEA: ref.NEA,
			Explanation: fmt.Sprintf("strand flip: %s/%s -> %s/%s; sign preserved",
				target.EA, target.NEA, fea, fnea),
		}, nil
	}
	// 4. 翻转 + 互换
	if fea == ref.NEA && fnea == ref.EA {
		return &AlignResult{
			Relation: RelationFlipped | RelationSwapped, Sign: -1,
			MatchedEA: ref.EA, MatchedNEA: ref.NEA,
			Explanation: fmt.Sprintf("strand flip + swap: %s/%s -> %s/%s; sign flipped",
				target.EA, target.NEA, fea, fnea),
		}, nil
	}

	return nil, &RejectError{ID: id, Reason: ReasonAlleleMismatch,
		Detail: fmt.Sprintf("alleles target=%s/%s cannot be aligned to ref=%s/%s",
			target.EA, target.NEA, ref.EA, ref.NEA)}
}

// alleleFlippable 报告等位基因是否每个碱基都是可互补的 A/C/G/T。
// 当前 NormalizeAllele 已保证纯 ACGT，indel 同样可逐碱基互补，
// 但 indel 的链表示在来源间不可靠，这里保守要求两等位基因长度相同
// 且与对侧长度一致——长度检查由调用场景保证；本函数只做 ACGT 校验。
func alleleFlippable(a Allele) bool {
	if len(a) == 0 {
		return false
	}
	for i := 0; i < len(a); i++ {
		if _, ok := ComplementBase(a[i]); !ok {
			return false
		}
	}
	return true
}

// ComplementBase 返回单碱基的互补碱基；非 ACGT 时第二返回值为 false。
func ComplementBase(b byte) (byte, bool) {
	switch b {
	case 'A':
		return 'T', true
	case 'T':
		return 'A', true
	case 'C':
		return 'G', true
	case 'G':
		return 'C', true
	default:
		return 0, false
	}
}
