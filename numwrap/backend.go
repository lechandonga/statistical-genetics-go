package numwrap

import "sync"

// Backend 标识数值后端实现。
type Backend string

const (
	// BackendStd 使用 Go 标准库 math 实现的自包含后端（默认、无第三方依赖）。
	BackendStd Backend = "stdlib"
	// BackendGonum 使用 gonum/gonum（mat/stat）实现的第三方交叉验证后端。
	// 版本在 go.mod / vendor/modules.txt 中锁死（见 docs/compatibility.md）。
	BackendGonum Backend = "gonum"
	// BackendCrosscheck 不是一个求解器：它在一次调用内分别用 stdlib 与
	// gonum 求解，并按 CrosscheckTol 比较；数值不一致或判定相关字段
	// 不一致时返回 KindBackendMismatch 错误，绝不挑一边的结果返回。
	BackendCrosscheck Backend = "crosscheck"
)

// backendVersions 记录各后端的兼容版本范围（见 docs/compatibility.md）。
var backendVersions = map[Backend]string{
	BackendStd:        "go1.25+ (stdlib math, build-tag fixed)",
	BackendGonum:      "gonum.org/v1/gonum v0.17.0 (vendored)",
	BackendCrosscheck: "stdlib v0 + gonum.org/v1/gonum v0.17.0 (vendored), compared per call",
}

// AvailableBackends 返回当前构建实际可用的后端（顺序固定）。
// crosscheck 始终可用（它是两种真实后端的组合）。
func AvailableBackends() []Backend {
	return []Backend{BackendStd, BackendGonum, BackendCrosscheck}
}

// ResolveBackend 把字符串解析为已编译后端；未知后端返回
// KindUnsupportedOption，绝不静默回退默认后端。
func ResolveBackend(s string) (Backend, error) {
	b := Backend(s)
	if _, ok := backendVersions[b]; ok {
		return b, nil
	}
	return "", &Error{Kind: KindUnsupportedOption, Op: "ResolveBackend",
		Detail: "unknown or uncompiled backend: " + s}
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
	if _, ok := backendVersions[b]; !ok {
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
