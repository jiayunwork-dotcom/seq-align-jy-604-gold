package align

import (
	"testing"

	"seq-align/internal/scoring"
)

// TestGlobalScoreAtLeastLocalScore verifies the invariant that for the same
// input sequences and scoring scheme, the global alignment score is always
// <= the local alignment score (local can cherry-pick the best region).
// This is a cross-mode contract: modifying one algorithm must not violate this.
func TestGlobalScoreAtLeastLocalScore(t *testing.T) {
	cases := []struct {
		name string
		a, b string
		s    scoring.Scheme
	}{
		{
			name: "short_mismatch_flanks",
			a:    "TTACGTACGTTT",
			b:    "AAACGTACGTAA",
			s:    scoring.Default(),
		},
		{
			name: "affine_long_gap",
			a:    "ACGTACGTACGT",
			b:    "ACGTACGT",
			s:    scoring.DefaultAffine(),
		},
		{
			name: "divergent_flanks_affine",
			a:    "GGGGACGTACGTCCCC",
			b:    "TTTTACGTACGTAAAA",
			s:    scoring.Scheme{Match: 3, Mismatch: -2, GapOpen: -6, GapExtend: -1},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			gRes, err := Global(tc.a, tc.b, tc.s)
			if err != nil {
				t.Fatalf("Global: %v", err)
			}
			lRes, err := Local(tc.a, tc.b, tc.s)
			if err != nil {
				t.Fatalf("Local: %v", err)
			}
			// Local best score >= Global score (local can skip bad flanks)
			if lRes.Score < gRes.Score {
				t.Errorf("contract violation: local score %d < global score %d", lRes.Score, gRes.Score)
			}
		})
	}
}

// TestLocalAlignmentBoundedBySequenceLength ensures local alignment coordinates
// are valid: StartA + non-gap chars in A <= len(original A), same for B.
func TestLocalAlignmentBoundedBySequenceLength(t *testing.T) {
	cases := []struct {
		a, b string
		s    scoring.Scheme
	}{
		{"ACGTACGTACGT", "ACGT", scoring.Default()},
		{"XXXACGTYYY", "ZZZACGTZZZ", scoring.DefaultAffine()},
		{"AAAAACGTACGTAAAA", "CGTACGT", scoring.Scheme{Match: 2, Mismatch: -3, GapOpen: -4, GapExtend: -1}},
	}
	for _, tc := range cases {
		res, err := Local(tc.a, tc.b, tc.s)
		if err != nil {
			t.Fatalf("Local(%q,%q): %v", tc.a, tc.b, err)
		}
		if res.AlignedLen == 0 {
			continue
		}
		nonGapA := countNonGap([]byte(res.A))
		nonGapB := countNonGap([]byte(res.B))
		if res.StartA+nonGapA > len(tc.a) {
			t.Errorf("StartA(%d)+consumed(%d) > len(a)(%d)", res.StartA, nonGapA, len(tc.a))
		}
		if res.StartB+nonGapB > len(tc.b) {
			t.Errorf("StartB(%d)+consumed(%d) > len(b)(%d)", res.StartB, nonGapB, len(tc.b))
		}
	}
}

// TestAffineGapOpenExtendInvariant verifies that for affine scoring,
// extending a gap is always cheaper than opening a new one. This contract
// means re-scoring an alignment with a single long gap must cost less than
// the same total gap length split into multiple blocks.
func TestAffineGapOpenExtendInvariant(t *testing.T) {
	s := scoring.Scheme{Match: 2, Mismatch: -1, GapOpen: -5, GapExtend: -1}
	// Single gap of length 4: open + 3*extend = -5 + 3*(-1) = -8
	// Four gaps of length 1: 4 * open = 4*(-5) = -20
	singleGapCost := s.GapCost(4)
	multiGapCost := 4 * s.GapCost(1)
	// Costs are negative; "cheaper" means closer to zero (higher value).
	if singleGapCost <= multiGapCost {
		t.Errorf("single gap cost (%d) should be cheaper (higher) than multi gap cost (%d)",
			singleGapCost, multiGapCost)
	}

	// Verify alignment behavior: align sequences where optimal has one long gap
	a := "ACGTXXXXACGT"
	b := "ACGTACGT"
	res, err := Global(a, b, s)
	if err != nil {
		t.Fatal(err)
	}
	// Should have exactly 1 gap block in B (the XXXX deletion)
	blocks := 0
	inGap := false
	for _, c := range res.B {
		if c == '-' {
			if !inGap {
				blocks++
				inGap = true
			}
		} else {
			inGap = false
		}
	}
	if blocks != 1 {
		t.Errorf("expected 1 gap block, got %d (A=%s B=%s)", blocks, res.A, res.B)
	}
}

// TestGlobalAffineLocalModeBoundary is the HIGH-DIFFICULTY test surface.
// It verifies that the boundary condition where local and global alignment
// produce different traceback paths on the SAME scoring scheme is handled
// correctly. Specifically: a sequence pair where global must use mismatching
// flanks (penalized) but local clips them, yielding different alignment strings
// but both internally consistent (rescore == reported score).
func TestGlobalAffineLocalModeBoundary(t *testing.T) {
	// Sequences with matching core and mismatching flanks
	a := "CCCCACGTACGTGGGG"
	b := "TTTTACGTACGTAAAA"
	s := scoring.Scheme{Match: 3, Mismatch: -2, GapOpen: -6, GapExtend: -1}

	gRes, err := Global(a, b, s)
	if err != nil {
		t.Fatalf("Global: %v", err)
	}
	lRes, err := Local(a, b, s)
	if err != nil {
		t.Fatalf("Local: %v", err)
	}

	// Global must align the full sequences (aligned length covers all positions)
	gNonGapA := countNonGap([]byte(gRes.A))
	gNonGapB := countNonGap([]byte(gRes.B))
	if gNonGapA != len(a) {
		t.Errorf("global should consume all of A: got %d, want %d", gNonGapA, len(a))
	}
	if gNonGapB != len(b) {
		t.Errorf("global should consume all of B: got %d, want %d", gNonGapB, len(b))
	}

	// Local should find only the matching core (8 chars ACGTACGT)
	if lRes.AlignedLen < 8 {
		t.Errorf("local should find at least 8-char match, got %d", lRes.AlignedLen)
	}
	// Local identity should be higher than global (it skips mismatches)
	if lRes.Identity <= gRes.Identity {
		t.Errorf("local identity (%.1f) should exceed global identity (%.1f)",
			lRes.Identity, gRes.Identity)
	}

	// Both must rescore correctly
	gSc, err := scoring.ScoreAlignment(gRes.A, gRes.B, s)
	if err != nil {
		t.Fatalf("rescore global: %v", err)
	}
	if gSc != gRes.Score {
		t.Errorf("global rescore mismatch: %d vs %d", gSc, gRes.Score)
	}
	lSc, err := scoring.ScoreAlignment(lRes.A, lRes.B, s)
	if err != nil {
		t.Fatalf("rescore local: %v", err)
	}
	if lSc != lRes.Score {
		t.Errorf("local rescore mismatch: %d vs %d", lSc, lRes.Score)
	}

	// Contract: local score >= global score
	if lRes.Score < gRes.Score {
		t.Errorf("contract: local score %d < global score %d", lRes.Score, gRes.Score)
	}
}
