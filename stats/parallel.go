package stats

import (
	"context"
	"runtime"
	"sync"
)

// mapShards 把长度为 n 的输入按确定性边界切成至多 workers 个分片，
// 对每个分片调用 fn(start,end)，结果按分片顺序写回返回切片。
//
// 不变量（满足“并行结果与串行一致”的要求）：
//   - 分片边界只依赖 n 与 workers，与调度无关；
//   - fn 必须是纯函数：不得读写共享可变状态，不得依赖 goroutine ID；
//   - 返回结果严格按分片 0..k-1 的顺序排列；
//   - 任意一个分片返回错误（或 ctx 取消）即整体失败，
//     且不会返回任何部分结果（无半更新输出）。
//
// workers<=1 时直接串行执行（同样的分片边界与顺序）。
func mapShards[T any](ctx context.Context, n, workers int, fn func(start, end int) (T, error)) ([]T, error) {
	if workers < 1 {
		workers = 1
	}
	if workers > n {
		workers = n
	}
	bounds := shardBounds(n, workers)
	out := make([]T, len(bounds))
	if len(bounds) == 1 || workers == 1 {
		for i, b := range bounds {
			if err := ctx.Err(); err != nil {
				return nil, err
			}
			v, err := fn(b[0], b[1])
			if err != nil {
				return nil, err
			}
			out[i] = v
		}
		return out, nil
	}

	errCh := make(chan error, len(bounds))
	var wg sync.WaitGroup
	// 用固定槽位写入，避免结果顺序受调度影响。
	for i, b := range bounds {
		wg.Add(1)
		go func(i, start, end int) {
			defer wg.Done()
			if err := ctx.Err(); err != nil {
				errCh <- err
				return
			}
			v, err := fn(start, end)
			if err != nil {
				errCh <- err
				return
			}
			out[i] = v
		}(i, b[0], b[1])
	}
	wg.Wait()
	close(errCh)
	for err := range errCh {
		if err != nil {
			return nil, err
		}
	}
	return out, nil
}

// shardBounds 返回 [start,end) 边界；除最后一片外每片大小为 ceil(n/k)
// 或按余数均匀分配——这里采用“除最后一片外等长”的稳定划分。
func shardBounds(n, k int) [][2]int {
	if n == 0 {
		return nil
	}
	if k > n {
		k = n
	}
	base := n / k
	rem := n % k
	bounds := make([][2]int, 0, k)
	start := 0
	for i := 0; i < k; i++ {
		size := base
		if i < rem {
			size++
		}
		if size == 0 {
			continue
		}
		bounds = append(bounds, [2]int{start, start + size})
		start += size
	}
	return bounds
}

// normalizeWorkers 收敛并行度；0 => GOMAXPROCS（受上限 n 约束）。
func normalizeWorkers(w, n int) int {
	if w < 0 {
		w = 0
	}
	if w == 0 {
		w = runtime.GOMAXPROCS(0)
	}
	if w > n && n > 0 {
		w = n
	}
	if w < 1 {
		w = 1
	}
	return w
}
