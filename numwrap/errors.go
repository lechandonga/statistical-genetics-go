// Package numwrap 对第三方数值/统计库（以及标准库数学函数）做归一化封装。
//
// 目标：
//   - 固定后端版本语义，对外只暴露本包类型；
//   - 库的语法错误（如不支持选项）与内部运行时异常（NaN/Inf、矩阵奇异、
//     维度不匹配）映射为可区分的错误类型，不向接口外透出第三方错误；
//   - 不支持的选项通过 ErrUnsupportedOption 显式拒绝。
package numwrap

import "fmt"

// ErrorKind 区分数值层错误类别。
type ErrorKind uint8

const (
	KindUnknown           ErrorKind = iota
	KindUnsupportedOption           // 请求了后端不支持的选项
	KindInvalidArgument             // 维度/长度/参数非法（调用方编程错误）
	KindSingular                    // 矩阵奇异，无法求解
	KindNaN                         // 输入或结果出现 NaN
	KindInf                         // 输入或结果出现 Inf
	KindNoConvergence               // 迭代未收敛
	KindBackendMismatch             // 双后端交叉验证结果不一致（含数值与判定字段）
)

// Error 是数值层唯一对外错误类型。
type Error struct {
	Kind   ErrorKind
	Op     string // 发生错误的操作名
	Detail string
}

func (e *Error) Error() string {
	return fmt.Sprintf("numwrap: op=%s kind=%s: %s", e.Op, kindName(e.Kind), e.Detail)
}

func kindName(k ErrorKind) string {
	switch k {
	case KindUnsupportedOption:
		return "unsupported_option"
	case KindInvalidArgument:
		return "invalid_argument"
	case KindSingular:
		return "singular"
	case KindNaN:
		return "nan"
	case KindInf:
		return "inf"
	case KindNoConvergence:
		return "no_convergence"
	case KindBackendMismatch:
		return "backend_mismatch"
	default:
		return "unknown"
	}
}

// AsError 从任意错误提取 *Error。
func AsError(err error) (*Error, bool) {
	if e, ok := err.(*Error); ok {
		return e, true
	}
	return nil, false
}

// Sentinel errors for direct errors.Is comparison.
type typedSentinel struct {
	kind ErrorKind
	msg  string
}
