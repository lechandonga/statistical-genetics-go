package variant

import "fmt"

// Reason 是可区分的输入拒绝原因代码。
type Reason string

const (
	ReasonMissingID          Reason = "missing_id"
	ReasonMissingAlleles     Reason = "missing_alleles"
	ReasonInvalidAllele      Reason = "invalid_allele"
	ReasonInvalidChrom       Reason = "invalid_chromosome"
	ReasonInvalidPosition    Reason = "invalid_position"
	ReasonInvalidFrequency   Reason = "invalid_frequency"
	ReasonInvalidNumeric     Reason = "invalid_numeric"
	ReasonDuplicateVariant   Reason = "duplicate_variant"
	ReasonUnresolvableStrand Reason = "unresolvable_strand"
	ReasonAmbiguousStrand    Reason = "ambiguous_strand"
	ReasonFlipForbidden      Reason = "flip_forbidden"
	ReasonAlleleMismatch     Reason = "allele_mismatch"
	ReasonMissingField       Reason = "missing_field"
)

// RejectError 携带行号、变体标识、原因代码与判定依据。
type RejectError struct {
	Source string // 输入源名称（文件名/逻辑名）
	Line   int    // 1-based 行号（含表头时为物理行号）
	ID     string // 尽力解析出的变体标识
	Reason Reason
	Detail string
}

func (e *RejectError) Error() string {
	return fmt.Sprintf("reject: source=%s line=%d id=%q reason=%s detail=%s",
		e.Source, e.Line, e.ID, e.Reason, e.Detail)
}

// AsReject 从任意错误中提取 *RejectError。
func AsReject(err error) (*RejectError, bool) {
	for err != nil {
		if r, ok := err.(*RejectError); ok {
			return r, true
		}
		type unwrapper interface{ Unwrap() error }
		u, ok := err.(unwrapper)
		if !ok {
			return nil, false
		}
		err = u.Unwrap()
	}
	return nil, false
}
