// Package align provides Needleman-Wunsch global and Smith-Waterman local
// sequence alignment with support for both linear and affine gap penalties.
package align

import (
	"fmt"

	"seq-align/internal/scoring"
)

const negInf = -(1 << 30)

// Result holds the outcome of a pairwise alignment.
type Result struct {
	Score      int
	A, B       string  // aligned sequences (with gaps '-')
	Identity   float64 // percent identity over aligned length
	AlignedLen int
	StartA     int // 0-based start position in original seq A (local only)
	StartB     int // 0-based start position in original seq B (local only)
}

// Global performs Needleman-Wunsch global alignment.
// Supports both linear gap (single Gap field) and affine gap (GapOpen/GapExtend).
func Global(a, b string, s scoring.Scheme) (Result, error) {
	if err := scoring.Validate(s); err != nil {
		return Result{}, err
	}
	if len(a) == 0 || len(b) == 0 {
		return Result{}, fmt.Errorf("align: global alignment requires non-empty sequences")
	}
	if s.IsAffine() {
		return globalAffine(a, b, s)
	}
	return globalLinear(a, b, s)
}

// Local performs Smith-Waterman local alignment.
// Supports both linear gap and affine gap penalties.
func Local(a, b string, s scoring.Scheme) (Result, error) {
	if err := scoring.Validate(s); err != nil {
		return Result{}, err
	}
	if len(a) == 0 || len(b) == 0 {
		return Result{}, fmt.Errorf("align: local alignment requires non-empty sequences")
	}
	if s.IsAffine() {
		return localAffine(a, b, s)
	}
	return localLinear(a, b, s)
}

// --- Linear gap implementation (original Needleman-Wunsch / Smith-Waterman) ---

type cell struct {
	score int
	prev  byte // 'd' diagonal, 'u' up, 'l' left, 0 stop (local)
}

func globalLinear(a, b string, s scoring.Scheme) (Result, error) {
	rows, cols := len(a)+1, len(b)+1
	dp := make([][]cell, rows)
	for i := range dp {
		dp[i] = make([]cell, cols)
	}
	for j := 1; j < cols; j++ {
		dp[0][j] = cell{score: dp[0][j-1].score + s.Gap, prev: 'l'}
	}
	for i := 1; i < rows; i++ {
		dp[i][0] = cell{score: dp[i-1][0].score + s.Gap, prev: 'u'}
		for j := 1; j < cols; j++ {
			diag := dp[i-1][j-1].score
			if a[i-1] == b[j-1] {
				diag += s.Match
			} else {
				diag += s.Mismatch
			}
			up := dp[i-1][j].score + s.Gap
			left := dp[i][j-1].score + s.Gap
			best, prev := diag, byte('d')
			if up > best {
				best, prev = up, 'u'
			}
			if left > best {
				best, prev = left, 'l'
			}
			dp[i][j] = cell{score: best, prev: prev}
		}
	}
	return tracebackLinear(dp, rows-1, cols-1, a, b), nil
}

func localLinear(a, b string, s scoring.Scheme) (Result, error) {
	rows, cols := len(a)+1, len(b)+1
	dp := make([][]cell, rows)
	for i := range dp {
		dp[i] = make([]cell, cols)
	}
	maxI, maxJ, maxScore := 0, 0, 0
	for i := 1; i < rows; i++ {
		for j := 1; j < cols; j++ {
			diag := dp[i-1][j-1].score
			if a[i-1] == b[j-1] {
				diag += s.Match
			} else {
				diag += s.Mismatch
			}
			up := dp[i-1][j].score + s.Gap
			left := dp[i][j-1].score + s.Gap
			best, prev := 0, byte(0)
			if diag > best {
				best, prev = diag, 'd'
			}
			if up > best {
				best, prev = up, 'u'
			}
			if left > best {
				best, prev = left, 'l'
			}
			dp[i][j] = cell{score: best, prev: prev}
			if best > maxScore {
				maxScore, maxI, maxJ = best, i, j
			}
		}
	}
	return tracebackLinear(dp, maxI, maxJ, a, b), nil
}

func tracebackLinear(dp [][]cell, i, j int, a, b string) Result {
	startI, startJ := i, j
	var ra, rb []byte
	for i > 0 || j > 0 {
		c := dp[i][j]
		if c.prev == 0 {
			break
		}
		switch c.prev {
		case 'd':
			ra = append(ra, a[i-1])
			rb = append(rb, b[j-1])
			i--
			j--
		case 'u':
			ra = append(ra, a[i-1])
			rb = append(rb, '-')
			i--
		case 'l':
			ra = append(ra, '-')
			rb = append(rb, b[j-1])
			j--
		}
	}
	reverse(ra)
	reverse(rb)
	same := 0
	for k := range ra {
		if ra[k] == rb[k] {
			same++
		}
	}
	n := len(ra)
	var identity float64
	if n > 0 {
		identity = float64(same) / float64(n) * 100
	}
	return Result{
		Score:      dp[startI][startJ].score,
		A:          string(ra),
		B:          string(rb),
		Identity:   identity,
		AlignedLen: n,
		StartA:     i, // position where traceback stopped
		StartB:     j,
	}
}

// --- Affine gap implementation (Gotoh three-matrix) ---
// M[i][j]: best score ending with a match/mismatch at (i,j)
// Ix[i][j]: best score ending with a gap in sequence B (deletion from A)
// Iy[i][j]: best score ending with a gap in sequence A (insertion into A)

type affineCell struct {
	m, ix, iy int
}

// traceDir encodes which matrix we came from at each cell
type traceEntry struct {
	mFrom  byte // 'd'=M, 'x'=Ix, 'y'=Iy
	ixFrom byte // 'o'=open from M, 'e'=extend from Ix
	iyFrom byte // 'o'=open from M, 'e'=extend from Iy
}

func globalAffine(a, b string, s scoring.Scheme) (Result, error) {
	rows, cols := len(a)+1, len(b)+1
	dp := make([][]affineCell, rows)
	tr := make([][]traceEntry, rows)
	for i := range dp {
		dp[i] = make([]affineCell, cols)
		tr[i] = make([]traceEntry, cols)
	}

	// Initialize boundaries
	dp[0][0] = affineCell{m: 0, ix: negInf, iy: negInf}
	for j := 1; j < cols; j++ {
		dp[0][j] = affineCell{
			m:  negInf,
			ix: negInf,
			iy: s.GapOpen + (j-1)*s.GapExtend,
		}
		tr[0][j] = traceEntry{iyFrom: 'e'}
		if j == 1 {
			tr[0][j] = traceEntry{iyFrom: 'o'}
		}
	}
	for i := 1; i < rows; i++ {
		dp[i][0] = affineCell{
			m:  negInf,
			ix: s.GapOpen + (i-1)*s.GapExtend,
			iy: negInf,
		}
		tr[i][0] = traceEntry{ixFrom: 'e'}
		if i == 1 {
			tr[i][0] = traceEntry{ixFrom: 'o'}
		}
	}

	for i := 1; i < rows; i++ {
		for j := 1; j < cols; j++ {
			// Ix: gap in B (consuming A[i])
			ixFromM := dp[i-1][j].m + s.GapOpen
			ixFromIx := dp[i-1][j].ix + s.GapExtend
			if ixFromM >= ixFromIx {
				dp[i][j].ix = ixFromM
				tr[i][j].ixFrom = 'o'
			} else {
				dp[i][j].ix = ixFromIx
				tr[i][j].ixFrom = 'e'
			}

			// Iy: gap in A (consuming B[j])
			iyFromM := dp[i][j-1].m + s.GapOpen
			iyFromIy := dp[i][j-1].iy + s.GapExtend
			if iyFromM >= iyFromIy {
				dp[i][j].iy = iyFromM
				tr[i][j].iyFrom = 'o'
			} else {
				dp[i][j].iy = iyFromIy
				tr[i][j].iyFrom = 'e'
			}

			// M: match/mismatch at (i,j)
			sub := s.Mismatch
			if a[i-1] == b[j-1] {
				sub = s.Match
			}
			mFromM := dp[i-1][j-1].m + sub
			mFromIx := dp[i-1][j-1].ix + sub
			mFromIy := dp[i-1][j-1].iy + sub
			best := mFromM
			from := byte('d')
			if mFromIx > best {
				best = mFromIx
				from = 'x'
			}
			if mFromIy > best {
				best = mFromIy
				from = 'y'
			}
			dp[i][j].m = best
			tr[i][j].mFrom = from
		}
	}

	// Find terminal best among M, Ix, Iy at (rows-1, cols-1)
	last := dp[rows-1][cols-1]
	best := last.m
	startMatrix := byte('m')
	if last.ix > best {
		best = last.ix
		startMatrix = 'x'
	}
	if last.iy > best {
		best = last.iy
		startMatrix = 'y'
	}

	ra, rb := tracebackAffine(tr, dp, rows-1, cols-1, startMatrix, a, b, false)
	return buildResult(ra, rb, best, 0, 0), nil
}

func localAffine(a, b string, s scoring.Scheme) (Result, error) {
	rows, cols := len(a)+1, len(b)+1
	dp := make([][]affineCell, rows)
	tr := make([][]traceEntry, rows)
	for i := range dp {
		dp[i] = make([]affineCell, cols)
		tr[i] = make([]traceEntry, cols)
		for j := range dp[i] {
			dp[i][j] = affineCell{m: 0, ix: 0, iy: 0}
		}
	}

	maxScore := 0
	maxI, maxJ := 0, 0
	maxMatrix := byte('m')

	for i := 1; i < rows; i++ {
		for j := 1; j < cols; j++ {
			// Ix
			ixFromM := dp[i-1][j].m + s.GapOpen
			ixFromIx := dp[i-1][j].ix + s.GapExtend
			if ixFromM >= ixFromIx && ixFromM > 0 {
				dp[i][j].ix = ixFromM
				tr[i][j].ixFrom = 'o'
			} else if ixFromIx > 0 {
				dp[i][j].ix = ixFromIx
				tr[i][j].ixFrom = 'e'
			}
			// else stays 0

			// Iy
			iyFromM := dp[i][j-1].m + s.GapOpen
			iyFromIy := dp[i][j-1].iy + s.GapExtend
			if iyFromM >= iyFromIy && iyFromM > 0 {
				dp[i][j].iy = iyFromM
				tr[i][j].iyFrom = 'o'
			} else if iyFromIy > 0 {
				dp[i][j].iy = iyFromIy
				tr[i][j].iyFrom = 'e'
			}

			// M
			sub := s.Mismatch
			if a[i-1] == b[j-1] {
				sub = s.Match
			}
			mFromM := dp[i-1][j-1].m + sub
			mFromIx := dp[i-1][j-1].ix + sub
			mFromIy := dp[i-1][j-1].iy + sub
			best := 0
			from := byte(0)
			if mFromM > best {
				best = mFromM
				from = 'd'
			}
			if mFromIx > best {
				best = mFromIx
				from = 'x'
			}
			if mFromIy > best {
				best = mFromIy
				from = 'y'
			}
			dp[i][j].m = best
			tr[i][j].mFrom = from

			// Track maximum
			if dp[i][j].m > maxScore {
				maxScore = dp[i][j].m
				maxI, maxJ = i, j
				maxMatrix = 'm'
			}
			if dp[i][j].ix > maxScore {
				maxScore = dp[i][j].ix
				maxI, maxJ = i, j
				maxMatrix = 'x'
			}
			if dp[i][j].iy > maxScore {
				maxScore = dp[i][j].iy
				maxI, maxJ = i, j
				maxMatrix = 'y'
			}
		}
	}

	if maxScore == 0 {
		return Result{}, nil
	}

	ra, rb := tracebackAffine(tr, dp, maxI, maxJ, maxMatrix, a, b, true)
	// startA/startB: where the local alignment begins in original sequences
	startA := maxI - countNonGap(ra)
	startB := maxJ - countNonGap(rb)
	return buildResult(ra, rb, maxScore, startA, startB), nil
}

func tracebackAffine(tr [][]traceEntry, dp [][]affineCell, i, j int, matrix byte, a, b string, local bool) ([]byte, []byte) {
	var ra, rb []byte
	for i > 0 || j > 0 {
		if local {
			switch matrix {
			case 'm':
				if dp[i][j].m <= 0 {
					return ra, rb
				}
			case 'x':
				if dp[i][j].ix <= 0 {
					return ra, rb
				}
			case 'y':
				if dp[i][j].iy <= 0 {
					return ra, rb
				}
			}
		}

		switch matrix {
		case 'm':
			ra = append(ra, a[i-1])
			rb = append(rb, b[j-1])
			from := tr[i][j].mFrom
			i--
			j--
			switch from {
			case 'd':
				matrix = 'm'
			case 'x':
				matrix = 'x'
			case 'y':
				matrix = 'y'
			default:
				// reached origin
				return ra, rb
			}
		case 'x':
			ra = append(ra, a[i-1])
			rb = append(rb, '-')
			from := tr[i][j].ixFrom
			i--
			switch from {
			case 'o':
				matrix = 'm'
			case 'e':
				matrix = 'x'
			default:
				return ra, rb
			}
		case 'y':
			ra = append(ra, '-')
			rb = append(rb, b[j-1])
			from := tr[i][j].iyFrom
			j--
			switch from {
			case 'o':
				matrix = 'm'
			case 'e':
				matrix = 'y'
			default:
				return ra, rb
			}
		default:
			return ra, rb
		}
	}
	return ra, rb
}

func buildResult(ra, rb []byte, score, startA, startB int) Result {
	reverse(ra)
	reverse(rb)
	same := 0
	for k := range ra {
		if ra[k] == rb[k] {
			same++
		}
	}
	n := len(ra)
	var identity float64
	if n > 0 {
		identity = float64(same) / float64(n) * 100
	}
	return Result{
		Score:      score,
		A:          string(ra),
		B:          string(rb),
		Identity:   identity,
		AlignedLen: n,
		StartA:     startA,
		StartB:     startB,
	}
}

func countNonGap(s []byte) int {
	n := 0
	for _, c := range s {
		if c != '-' {
			n++
		}
	}
	return n
}

func reverse(b []byte) {
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
}
