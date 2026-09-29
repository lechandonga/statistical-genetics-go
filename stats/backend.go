package stats

import (
	"fmt"

	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

// 本文件实现统计层的双后端交叉验证。
//
// 约定（见 docs/compatibility.md / docs/verification.md）：
//   - 数值字段（beta/se/CI/h2/tau2/Q/...）在 numwrap.CrosscheckRelTol/
//     AbsTol 内一致即视为通过；
//   - 判定字段（Flags、Excluded、ExcludeReason、NCohorts、逐队列
//     Excluded/Reason/Sign/Relation、样本数、位点集合）必须逐字段一致；
//   - 任一边出错：错误分类必须一致，否则报 backend_mismatch；
//   - 对不上一律返回 backend_mismatch 错误，绝不挑选任一边的结果。

// BackendMismatchError 表示双后端在同一输入上结论不一致。
// 它不包含任何第三方库内部文本。
type BackendMismatchError struct {
	Stage  string
	Detail string
}

func (e *BackendMismatchError) Error() string {
	return fmt.Sprintf("backend crosscheck failed at %s: %s", e.Stage, e.Detail)
}

func floatsAgree(a, b float64) bool { return numwrap.WithinTol(a, b) }

func floatSlicesAgree(a, b []float64) error {
	if len(a) != len(b) {
		return fmt.Errorf("slice length %d != %d", len(a), len(b))
	}
	for i := range a {
		if !floatsAgree(a[i], b[i]) {
			return fmt.Errorf("element %d: %.17g vs %.17g", i, a[i], b[i])
		}
	}
	return nil
}

// crosscheckRun 用两个后端各跑一遍同一分析并比较。
// run 必须是纯函数：只依赖其入参（后端在 b 中指定），不得读取
// numwrap 的进程级默认后端，以保证并发任务互不串扰。
func crosscheckRun[O, R any](stage string, setBackend func(O, numwrap.Backend) O,
	compare func(a, b R) error, run func(O) (R, error), opts O) (R, error) {
	var zero R
	o1 := setBackend(opts, numwrap.BackendStd)
	r1, err1 := run(o1)
	o2 := setBackend(opts, numwrap.BackendGonum)
	r2, err2 := run(o2)
	if err1 != nil || err2 != nil {
		if errorsSameKind(err1, err2) {
			return zero, err1
		}
		return zero, &BackendMismatchError{Stage: stage,
			Detail: fmt.Sprintf("one backend errored and the other did not or kinds differ: stdlib=%v gonum=%v", err1, err2)}
	}
	if err := compare(r1, r2); err != nil {
		return zero, &BackendMismatchError{Stage: stage, Detail: err.Error()}
	}
	return r1, nil
}

// errorsSameKind 比较两个错误是否属于同一“判定类别”：
// numwrap 错误按 Kind；stats 层错误按类型/原因字符串；nil==nil。
func errorsSameKind(a, b error) bool {
	if a == nil || b == nil {
		return a == nil && b == nil
	}
	if ea, ok := numwrap.AsError(a); ok {
		if eb, ok2 := numwrap.AsError(b); ok2 {
			return ea.Kind == eb.Kind
		}
		return false
	}
	return typeAndReason(a) == typeAndReason(b)
}

func typeAndReason(e error) string {
	return fmt.Sprintf("%T", e)
}
