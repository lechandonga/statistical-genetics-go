// Package parseutil contains the small deterministic helpers shared by the
// delimited-file parsers (sumstats, ldscore, genotype): splitting, numeric
// validation with NaN/Inf rejection, and column lookup.
package parseutil

import (
	"math"
	"strconv"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Split splits one line by delimiter. Delim "whitespace" (or "") splits on
// runs of spaces/tabs. Single-byte delimiters are matched exactly.
func Split(line, delim string) []string {
	if delim == "" || delim == "whitespace" {
		return strings.Fields(line)
	}
	if len(delim) == 1 {
		return strings.Split(line, delim)
	}
	return strings.Split(line, delim)
}

// IsBlank reports an empty/whitespace-only line.
func IsBlank(line string) bool { return strings.TrimSpace(line) == "" }

// Float parses a float and refuses NaN and Inf explicitly.
func Float(field, col string) (float64, error) {
	v, err := strconv.ParseFloat(strings.TrimSpace(field), 64)
	if err != nil {
		return 0, stgerr.New(stgerr.KindInvalidValue, "parseutil",
			"column "+col+": unparsable number "+strconv.Quote(field))
	}
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, stgerr.New(stgerr.KindInvalidValue, "parseutil",
			"column "+col+": NaN/Inf not permitted")
	}
	return v, nil
}

// OptionalFloat parses an empty field as absent (ok=false, no error).
func OptionalFloat(field, col string) (float64, bool, error) {
	if strings.TrimSpace(field) == "" || strings.EqualFold(strings.TrimSpace(field), "NA") {
		return 0, false, nil
	}
	v, err := Float(field, col)
	if err != nil {
		return 0, false, err
	}
	return v, true, nil
}

// Int parses a positive integer position.
func Int(field, col string) (int, error) {
	v, err := Float(field, col)
	if err != nil {
		return 0, err
	}
	if v != math.Trunc(v) {
		return 0, stgerr.New(stgerr.KindInvalidValue, "parseutil",
			"column "+col+": expected integer, got "+field)
	}
	return int(v), nil
}

// IndexOf finds a column name in a header (case-insensitive).
func IndexOf(header []string, names ...string) int {
	want := make(map[string]bool, len(names))
	for _, n := range names {
		want[strings.ToUpper(n)] = true
	}
	for i, h := range header {
		if want[strings.ToUpper(strings.TrimSpace(h))] {
			return i
		}
	}
	return -1
}

// Cell returns row[idx] with bounds checking.
func Cell(row []string, idx int, col string) (string, error) {
	if idx < 0 || idx >= len(row) {
		return "", stgerr.New(stgerr.KindMissingField, "parseutil",
			"column "+col+" absent from row")
	}
	return row[idx], nil
}
