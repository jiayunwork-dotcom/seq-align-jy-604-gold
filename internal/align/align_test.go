package align

import (
	"strings"
	"testing"

	"seq-align/internal/scoring"
)

func TestGlobalShort(t *testing.T) {
	res, err := Global("ACGT", "AGT", scoring.Default())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.Score != 4 {
		t.Fatalf("want score 4, got %d", res.Score)
	}
	if strings.Count(res.A, "-")+strings.Count(res.B, "-") != 1 {
		t.Fatalf("want exactly one gap, got %s/%s", res.A, res.B)
	}
	if res.AlignedLen != 4 {
		t.Fatalf("want aligned len 4, got %d", res.AlignedLen)
	}
	if sc, err := scoring.ScoreAlignment(res.A, res.B, scoring.Default()); err != nil || sc != res.Score {
		t.Fatalf("rescore mismatch: %d vs %d (err=%v)", sc, res.Score, err)
	}
}

func TestGlobalIdentical(t *testing.T) {
	res, err := Global("ACGTACGT", "ACGTACGT", scoring.Default())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.Identity != 100 {
		t.Fatalf("want identity 100, got %v", res.Identity)
	}
	if res.Score != 16 || res.A != "ACGTACGT" || res.B != "ACGTACGT" {
		t.Fatalf("unexpected result: %+v", res)
	}
}

func TestLocalFindsIsland(t *testing.T) {
	res, err := Local("AAAAAACGGCAAAAAA", "TTTTTCGGCTTTT", scoring.Default())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.A != "CGGC" || res.B != "CGGC" {
		t.Fatalf("want CGGC/CGGC, got %s/%s", res.A, res.B)
	}
	if res.Identity != 100 || res.Score != 8 || res.AlignedLen != 4 {
		t.Fatalf("unexpected result: %+v", res)
	}
}

func TestGlobalErrors(t *testing.T) {
	if _, err := Global("", "ACG", scoring.Default()); err == nil {
		t.Error("empty a: want error")
	}
	if _, err := Global("ACG", "", scoring.Default()); err == nil {
		t.Error("empty b: want error")
	}
	if _, err := Local("", "", scoring.Default()); err == nil {
		t.Error("both empty: want error")
	}
	if _, err := Global("ACG", "ACG", scoring.Scheme{Match: 0, Mismatch: -1, Gap: -1}); err == nil {
		t.Error("invalid scheme: want error")
	}
}

// --- Affine gap tests ---

func TestGlobalAffinePrefersSingleLongGap(t *testing.T) {
	// With affine gap, a single long gap is cheaper than multiple short gaps.
	// seq A: ACGTTTACGT (10)
	// seq B: ACGACGT    (7) -> needs 3 gaps somewhere
	s := scoring.Scheme{Match: 2, Mismatch: -3, GapOpen: -5, GapExtend: -1}
	res, err := Global("ACGTTTACGT", "ACGACGT", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Affine gap should produce one contiguous gap block of length 3
	// cost = open(-5) + 2*extend(-1) = -7 total gap penalty
	gapsA := countGapBlocks(res.A)
	gapsB := countGapBlocks(res.B)
	totalBlocks := gapsA + gapsB
	if totalBlocks != 1 {
		t.Errorf("affine should produce 1 gap block, got %d (A=%s, B=%s)", totalBlocks, res.A, res.B)
	}
	// Verify rescore consistency
	sc, err := scoring.ScoreAlignment(res.A, res.B, s)
	if err != nil {
		t.Fatalf("rescore error: %v", err)
	}
	if sc != res.Score {
		t.Errorf("rescore mismatch: alignment score %d vs rescored %d", res.Score, sc)
	}
}

func TestGlobalAffineTracebackConsistency(t *testing.T) {
	// Traceback-produced alignment must rescore to the same value as Result.Score.
	s := scoring.DefaultAffine()
	pairs := [][2]string{
		{"AACCGGTT", "ACGT"},
		{"GATTACA", "GCATGCU"},
		{"HEAGAWGHEE", "PAWHEAE"},
	}
	for _, p := range pairs {
		res, err := Global(p[0], p[1], s)
		if err != nil {
			t.Fatalf("Global(%q,%q): %v", p[0], p[1], err)
		}
		sc, err := scoring.ScoreAlignment(res.A, res.B, s)
		if err != nil {
			t.Fatalf("rescore error: %v", err)
		}
		if sc != res.Score {
			t.Errorf("Global(%q,%q): score=%d but rescore=%d (A=%s B=%s)",
				p[0], p[1], res.Score, sc, res.A, res.B)
		}
	}
}

func TestLocalAffineFindsIsland(t *testing.T) {
	s := scoring.Scheme{Match: 3, Mismatch: -3, GapOpen: -5, GapExtend: -1}
	res, err := Local("XXXACGTACGTXXX", "YYYACGTACGTYYY", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// The local alignment should find the common ACGTACGT region
	if !strings.Contains(res.A, "ACGTACGT") {
		t.Errorf("expected ACGTACGT in alignment A, got %s", res.A)
	}
	if res.Score <= 0 {
		t.Errorf("expected positive score, got %d", res.Score)
	}
	// Verify coordinates
	if res.StartA < 0 || res.StartB < 0 {
		t.Errorf("start coordinates should be non-negative: StartA=%d StartB=%d", res.StartA, res.StartB)
	}
}

func TestGlobalAffineVsLinearDifferentResults(t *testing.T) {
	// For sequences where gap placement matters, affine and linear should differ.
	a := "ACGTACGTACGT"
	b := "ACGTACGT" // 4 chars shorter
	linear := scoring.Scheme{Match: 2, Mismatch: -1, Gap: -2}
	affine := scoring.Scheme{Match: 2, Mismatch: -1, GapOpen: -5, GapExtend: -1}

	resLin, err := Global(a, b, linear)
	if err != nil {
		t.Fatal(err)
	}
	resAff, err := Global(a, b, affine)
	if err != nil {
		t.Fatal(err)
	}
	// Both should produce valid alignments
	if resLin.AlignedLen == 0 || resAff.AlignedLen == 0 {
		t.Fatal("both alignments should produce non-empty results")
	}
	// Affine should prefer fewer gap blocks
	linBlocks := countGapBlocks(resLin.A) + countGapBlocks(resLin.B)
	affBlocks := countGapBlocks(resAff.A) + countGapBlocks(resAff.B)
	if affBlocks > linBlocks && linBlocks > 1 {
		t.Errorf("affine should not have more gap blocks than linear: aff=%d lin=%d", affBlocks, linBlocks)
	}
}

func TestLocalAffineEmptyResult(t *testing.T) {
	// When all matches score negative (impossible with positive match), result is empty.
	// But with very harsh penalties on all-mismatch sequences, local can return empty.
	s := scoring.Scheme{Match: 1, Mismatch: -5, GapOpen: -10, GapExtend: -5}
	res, err := Local("AAAA", "CCCC", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.Score != 0 && res.AlignedLen != 0 {
		t.Errorf("expected empty local alignment for all-mismatch, got score=%d len=%d", res.Score, res.AlignedLen)
	}
}

// --- helpers ---

func countGapBlocks(s string) int {
	blocks := 0
	inGap := false
	for _, c := range s {
		if c == '-' {
			if !inGap {
				blocks++
				inGap = true
			}
		} else {
			inGap = false
		}
	}
	return blocks
}
