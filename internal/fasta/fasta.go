// Package fasta parses FASTA files and computes pairwise distance matrices.
package fasta

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"seq-align/internal/align"
	"seq-align/internal/scoring"
)

// Mode selects the alignment algorithm used for distance computation.
type Mode int

const (
	ModeGlobal Mode = iota
	ModeLocal
)

// ParseMode converts a string ("global"/"local") to a Mode value.
func ParseMode(s string) (Mode, error) {
	switch strings.ToLower(s) {
	case "global", "":
		return ModeGlobal, nil
	case "local":
		return ModeLocal, nil
	default:
		return 0, fmt.Errorf("fasta: unknown mode %q (want global|local)", s)
	}
}

// Record is a single FASTA record: ID is the first whitespace-delimited token after '>'.
type Record struct {
	ID  string
	Seq string
}

// Parse reads a FASTA file: '>id description' header + multi-line sequence concatenation.
// Tolerates UTF-8 BOM on the first line.
func Parse(path string) ([]Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("fasta: open %s: %w", path, err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	var recs []Record
	var cur *Record
	first := true
	for sc.Scan() {
		line := sc.Text()
		if first {
			line = strings.TrimPrefix(line, "\ufeff")
			first = false
		}
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			head := strings.TrimSpace(strings.TrimPrefix(line, ">"))
			if head == "" {
				return nil, fmt.Errorf("fasta: empty header in %s", path)
			}
			id := strings.Fields(head)[0]
			recs = append(recs, Record{ID: id})
			cur = &recs[len(recs)-1]
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("fasta: sequence data before any header in %s", path)
		}
		cur.Seq += line
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("fasta: read %s: %w", path, err)
	}
	if len(recs) == 0 {
		return nil, fmt.Errorf("fasta: no records in %s", path)
	}
	for _, r := range recs {
		if r.Seq == "" {
			return nil, fmt.Errorf("fasta: record %s has empty sequence", r.ID)
		}
	}
	return recs, nil
}

// DistanceMatrix computes all-pairs alignment, returning a full symmetric matrix
// and ID list. Distance = 100 - identity. Requires at least 2 records.
// The mode parameter selects global or local alignment.
func DistanceMatrix(rs []Record, s scoring.Scheme, mode Mode) ([][]float64, []string, error) {
	if len(rs) < 2 {
		return nil, nil, fmt.Errorf("fasta: need at least 2 records for distance matrix")
	}
	if err := scoring.Validate(s); err != nil {
		return nil, nil, err
	}
	ids := make([]string, len(rs))
	m := make([][]float64, len(rs))
	for i := range rs {
		ids[i] = rs[i].ID
		m[i] = make([]float64, len(rs))
	}
	for i := 0; i < len(rs); i++ {
		for j := i + 1; j < len(rs); j++ {
			var res align.Result
			var err error
			switch mode {
			case ModeLocal:
				res, err = align.Local(rs[i].Seq, rs[j].Seq, s)
			default:
				res, err = align.Global(rs[i].Seq, rs[j].Seq, s)
			}
			if err != nil {
				return nil, nil, fmt.Errorf("fasta: align %s vs %s: %w", rs[i].ID, rs[j].ID, err)
			}
			d := 100 - res.Identity
			if d < 0 {
				d = 0
			}
			// Symmetry guarantee: d(i,j) == d(j,i)
			m[i][j] = d
			m[j][i] = d
		}
	}
	return m, ids, nil
}
