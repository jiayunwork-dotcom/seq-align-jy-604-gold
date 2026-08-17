// Package scoring defines alignment scoring schemes and provides per-position scoring.
package scoring

import "fmt"

// Scheme is a scoring scheme for sequence alignment.
// For affine gap scoring, set GapOpen and GapExtend (both negative).
// When GapOpen == 0 and GapExtend == 0, the linear Gap field is used instead.
type Scheme struct {
	Match     int
	Mismatch  int
	Gap       int // linear gap penalty (used when GapOpen/GapExtend are zero)
	GapOpen   int // affine: cost to open a new gap (negative)
	GapExtend int // affine: cost to extend an existing gap (negative)
}

// IsAffine reports whether the scheme uses affine gap penalties.
func (s Scheme) IsAffine() bool {
	return s.GapOpen != 0 || s.GapExtend != 0
}

// Default returns the default linear scoring scheme {2, -1, -2}.
func Default() Scheme {
	return Scheme{Match: 2, Mismatch: -1, Gap: -2}
}

// DefaultAffine returns a default affine gap scoring scheme:
// match=2, mismatch=-1, gapOpen=-5, gapExtend=-1.
func DefaultAffine() Scheme {
	return Scheme{Match: 2, Mismatch: -1, GapOpen: -5, GapExtend: -1}
}

// Validate checks that a scheme has valid parameters.
// For linear mode: match>0, mismatch<0, gap<0.
// For affine mode: match>0, mismatch<0, gapOpen<0, gapExtend<0.
func Validate(s Scheme) error {
	if s.Match <= 0 {
		return fmt.Errorf("scoring: match must be > 0, got %d", s.Match)
	}
	if s.Mismatch >= 0 {
		return fmt.Errorf("scoring: mismatch must be < 0, got %d", s.Mismatch)
	}
	if s.IsAffine() {
		if s.GapOpen >= 0 {
			return fmt.Errorf("scoring: gap_open must be < 0, got %d", s.GapOpen)
		}
		if s.GapExtend >= 0 {
			return fmt.Errorf("scoring: gap_extend must be < 0, got %d", s.GapExtend)
		}
	} else {
		if s.Gap >= 0 {
			return fmt.Errorf("scoring: gap must be < 0, got %d", s.Gap)
		}
	}
	return nil
}

// GapCost returns the cost of a gap of length n in this scheme.
// For linear: n * Gap. For affine: GapOpen + (n-1) * GapExtend.
// Returns 0 if n <= 0.
func (s Scheme) GapCost(n int) int {
	if n <= 0 {
		return 0
	}
	if s.IsAffine() {
		return s.GapOpen + (n-1)*s.GapExtend
	}
	return n * s.Gap
}

// ScoreAlignment scores two aligned sequences position-by-position.
// Requires equal length; only A-Za-z and '-' are valid characters.
// Neither sequence may be entirely gaps. Supports both linear and affine gap scoring.
func ScoreAlignment(a, b string, s Scheme) (int, error) {
	if err := Validate(s); err != nil {
		return 0, err
	}
	if len(a) != len(b) {
		return 0, fmt.Errorf("scoring: aligned strings differ in length (%d vs %d)", len(a), len(b))
	}
	ok := func(c byte) bool { return c == '-' || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') }
	total := 0
	aAllGap, bAllGap := true, true
	aInGap, bInGap := false, false
	for i := 0; i < len(a); i++ {
		ca, cb := a[i], b[i]
		if !ok(ca) || !ok(cb) {
			return 0, fmt.Errorf("scoring: invalid character %q/%q at position %d", ca, cb, i)
		}
		if ca != '-' {
			aAllGap = false
		}
		if cb != '-' {
			bAllGap = false
		}
		switch {
		case ca == '-':
			if s.IsAffine() {
				if !aInGap {
					total += s.GapOpen
					aInGap = true
				} else {
					total += s.GapExtend
				}
			} else {
				total += s.Gap
			}
			bInGap = false
		case cb == '-':
			if s.IsAffine() {
				if !bInGap {
					total += s.GapOpen
					bInGap = true
				} else {
					total += s.GapExtend
				}
			} else {
				total += s.Gap
			}
			aInGap = false
		case ca == cb:
			total += s.Match
			aInGap = false
			bInGap = false
		default:
			total += s.Mismatch
			aInGap = false
			bInGap = false
		}
	}
	if aAllGap || bAllGap {
		return 0, fmt.Errorf("scoring: one sequence is entirely gaps")
	}
	return total, nil
}
