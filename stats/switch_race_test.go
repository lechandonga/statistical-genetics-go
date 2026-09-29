package stats

import (
	"context"
	"sync"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

// TestMetaBackendSwitchMidFlightRace：一批元分析在显式后端下运行的同时，
// 另一个 goroutine 反复切换进程默认后端。已在跑的任务必须按各自启动时
// 选定的后端完成，结论不随全局切换改变（-race 下检验竞态安全）。
func TestMetaBackendSwitchMidFlightRace(t *testing.T) {
	var recs1, recs2 []ingest.SummaryRecord
	for i := 0; i < 60; i++ {
		recs1 = append(recs1, metaRec(idFor(i), "A", "G", 0.1, 0.1, 100))
		recs2 = append(recs2, metaRec(idFor(i), "A", "G", 0.2, 0.1, 100))
	}
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: recs1}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: recs2}}

	var wg sync.WaitGroup
	stop := make(chan struct{})
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			select {
			case <-stop:
				return
			default:
				_ = numwrap.SetDefaultBackend(numwrap.BackendGonum)
				_ = numwrap.SetDefaultBackend(numwrap.BackendStd)
			}
		}
	}()

	for g := 0; g < 8; g++ {
		wg.Add(1)
		go func(g int) {
			defer wg.Done()
			b := numwrap.BackendStd
			if g%2 == 0 {
				b = numwrap.BackendGonum
			}
			for k := 0; k < 20; k++ {
				res, err := RunMeta(context.Background(), []Cohort{c1, c2},
					MetaOptions{MinCohorts: 2, Backend: b})
				if err != nil {
					t.Error(err)
					return
				}
				if res.Backend != b {
					t.Errorf("task backend drifted: want %s got %s", b, res.Backend)
					return
				}
				if err := numwrap.CrossCheckScalar("test", "beta",
					res.Variants[0].Beta, 0.15); err != nil {
					t.Errorf("beta drifted under switch: %v", res.Variants[0].Beta)
					return
				}
			}
		}(g)
	}
	close(stop)
	wg.Wait()
	t.Log("[judgement] in-flight meta tasks unaffected by concurrent default-backend switching")
}
