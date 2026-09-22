package numwrap

import "sync"

// Backend 标识数值后端实现。
type Backend string

const (
	// BackendStd 使用 Go 标准库 math 实现的自包含后端（默认、无第三方依赖）。
	BackendStd Backend = "stdlib"
	// BackendGonum 为预留后端；未编译进当前构建时请求会得到
	// KindUnsupportedOption 错误，不透出任何第三方包错误。
	BackendGonum Backend = "gonum"
)

// backendVersions 记录各后端的兼容版本范围（见 docs/compatibility.md）。
var backendVersions = map[Backend]string{
	BackendStd: "go1.25+ (stdlib math, build-tag fixed)",
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
