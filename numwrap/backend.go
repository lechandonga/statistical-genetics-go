package numwrap

import (
	"fmt"
	"sync"
)

// Backend 标识数值后端实现。
type Backend string

const (
	// BackendStd 使用 Go 标准库 math 实现的自包含后端（默认、无第三方依赖）。
	BackendStd Backend = "stdlib"
	// BackendGonum 基于 gonum/gonum 线代库的第二实现，用于双后端交叉验证。
	BackendGonum Backend = "gonum"
)

// GonumVersion 是本构建锁定并实际编译的 gonum 版本；
// 必须与 go.mod / vendor/modules.txt 中的 gonum.org/v1/gonum 版本一致，
// TestGonumVersionPinned 会校验三者，避免文档/代码/依赖漂移。
const GonumVersion = "v0.15.1"

// backendVersions 记录各后端的版本/来源（见 docs/compatibility.md）。
var backendVersions = map[Backend]string{
	BackendStd:   "go1.25+ (stdlib math, build-tag fixed)",
	BackendGonum: "gonum.org/v1/gonum " + GonumVersion + " (mat/gonum/stat/distuv)",
}

// AvailableBackends 返回当前构建实际编译可用的全部后端，顺序固定。
func AvailableBackends() []Backend {
	return []Backend{BackendStd, BackendGonum}
}

// Available 报告后端是否真实编译进本构建；不存在的后端返回 false，
// 调用方据此显式拒绝，绝不静默回退默认实现。
func Available(b Backend) bool {
	_, ok := backendVersions[b]
	return ok
}

// resolve 把空后端解析为当前默认，并显式拒绝未编译的后端。
func resolve(b Backend, op string) (Backend, error) {
	if b == "" {
		b = DefaultBackend()
	}
	if !Available(b) {
		return "", &Error{Kind: KindUnsupportedOption, Op: op,
			Detail: "backend not compiled into this build: " + string(b)}
	}
	return b, nil
}

var (
	mu      sync.RWMutex
	current = BackendStd
)

// Version 返回后端的版本/接口描述字符串。
func Version(b Backend) (string, error) {
	if v, ok := backendVersions[b]; ok {
		return v, nil
	}
	return "", &Error{Kind: KindUnsupportedOption, Op: "Version",
		Detail: "backend not compiled into this build: " + string(b)}
}

// SetDefaultBackend 设置进程级默认后端（并发安全）；
// 只影响之后新创建的求解器，正在执行的计算不被打断。
func SetDefaultBackend(b Backend) error {
	if !Available(b) {
		return &Error{Kind: KindUnsupportedOption, Op: "SetDefaultBackend",
			Detail: "backend unavailable: " + string(b)}
	}
	mu.Lock()
	current = b
	mu.Unlock()
	return nil
}

// DefaultBackend 返回当前默认后端。
func DefaultBackend() Backend {
	mu.RLock()
	b := current
	mu.RUnlock()
	return b
}

// CrossCheckTol 是两个后端之间数值一致性的声明容差（相对 + 绝对）。
// 同一后端重复运行必须逐位一致；跨后端只要求落在该容差内，
// 但所有离散判定（是否拒绝/标记/纳入）必须完全一致。
//
// 选 1e-9：两侧均为 float64 直接公式求解（无迭代），典型差异仅在
// 最后 1-2 个 ulp 量级；1e-9 对本流程全部下游结论（含 P 值判定）足够宽。
const (
	CrossCheckRelTol = 1e-9
	CrossCheckAbsTol = 1e-12
)

// BackendMismatchError 表示两个后端对同一份输入给出了超出声明容差的结果。
// 它是硬错误：交叉验证模式下宁可整体失败，也绝不挑选其中一侧的结果。
type BackendMismatchError struct {
	Op     string
	Detail string
}

func (e *BackendMismatchError) Error() string {
	return fmt.Sprintf("numwrap: backend cross-check mismatch: op=%s: %s", e.Op, e.Detail)
}

// equalWithinTol 按声明的相对/绝对容差比较两个标量。
func equalWithinTol(a, b float64) bool {
	d := a - b
	if d < 0 {
		d = -d
	}
	scale := a
	if b < 0 {
		scale = -b
	}
	if scale < 1 {
		scale = 1
	}
	return d <= CrossCheckAbsTol+CrossCheckRelTol*scale
}

// crossCheckSlices 比较两侧浮点切片并在不一致时返回 BackendMismatchError。
func crossCheckSlices(op, field string, a, b []float64) error {
	if len(a) != len(b) {
		return &BackendMismatchError{Op: op,
			Detail: fmt.Sprintf("%s length %d != %d", field, len(a), len(b))}
	}
	for i := range a {
		if !equalWithinTol(a[i], b[i]) {
			return &BackendMismatchError{Op: op,
				Detail: fmt.Sprintf("%s[%d] stdlib=%.17g gonum=%.17g diff=%.3g",
					field, i, a[i], b[i], absd(a[i], b[i]))}
		}
	}
	return nil
}

// CrossCheckScalar 比较两侧标量；导出供 stats 层交叉验证复用同一容差。
func CrossCheckScalar(op, field string, a, b float64) error {
	if !equalWithinTol(a, b) {
		return &BackendMismatchError{Op: op,
			Detail: fmt.Sprintf("%s stdlib=%.17g gonum=%.17g diff=%.3g",
				field, a, b, absd(a, b))}
	}
	return nil
}

// CrossCheckInt 比较两侧整数（离散量必须完全相等）。
func CrossCheckInt(op, field string, a, b int) error {
	if a != b {
		return &BackendMismatchError{Op: op,
			Detail: fmt.Sprintf("%s %d != %d", field, a, b)}
	}
	return nil
}

// CrossCheckBool 比较两侧布尔判定（是否拒绝/标记必须完全一致）。
func CrossCheckBool(op, field string, a, b bool) error {
	if a != b {
		return &BackendMismatchError{Op: op,
			Detail: fmt.Sprintf("%s %v != %v", field, a, b)}
	}
	return nil
}

func absd(a, b float64) float64 {
	d := a - b
	if d < 0 {
		return -d
	}
	return d
}
