package scoring

import "testing"

func TestScoreAlignment(t *testing.T) {
	s := Default()
	// ACGT vs A-GT: match(2) gap(-2) match(2) match(2) => 4
	got, err := ScoreAlignment("ACGT", "A-GT", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != 4 {
		t.Fatalf("want 4, got %d", got)
	}
	// identical sequences
	got, err = ScoreAlignment("ACGT", "ACGT", s)
	if err != nil || got != 8 {
		t.Fatalf("identical: want 8, err=%v", err)
	}
}

func TestScoreAlignmentAffine(t *testing.T) {
	s := Scheme{Match: 2, Mismatch: -1, GapOpen: -5, GapExtend: -1}
	// AC--GT vs ACTTGT: pos0 match(2), pos1 match(2), pos2 gap-open(-5), pos3 gap-ext(-1), pos4 match(2), pos5 match(2)
	// total = 2+2-5-1+2+2 = 2
	got, err := ScoreAlignment("AC--GT", "ACTTGT", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != 2 {
		t.Fatalf("want 2, got %d", got)
	}
	// single gap: open only
	// A-G vs ATG: match(2) gap-open(-5) match(2) = -1
	got, err = ScoreAlignment("A-G", "ATG", s)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != -1 {
		t.Fatalf("want -1, got %d", got)
	}
}

func TestScoreAlignmentErrors(t *testing.T) {
	s := Default()
	if _, err := ScoreAlignment("ACG", "ACGT", s); err == nil {
		t.Error("length mismatch: want error")
	}
	if _, err := ScoreAlignment("AC1", "AC2", s); err == nil {
		t.Error("invalid char: want error")
	}
	if _, err := ScoreAlignment("---", "ACG", s); err == nil {
		t.Error("all-gap: want error")
	}
	if _, err := ScoreAlignment("ACG", "ACG", Scheme{Match: 0, Mismatch: -1, Gap: -1}); err == nil {
		t.Error("invalid scheme: want error")
	}
}

func TestValidate(t *testing.T) {
	if err := Validate(Default()); err != nil {
		t.Fatalf("Default should be valid: %v", err)
	}
	if err := Validate(DefaultAffine()); err != nil {
		t.Fatalf("DefaultAffine should be valid: %v", err)
	}
	if err := Validate(Scheme{Match: -1, Mismatch: -1, Gap: -1}); err == nil {
		t.Error("match<=0: want error")
	}
	if err := Validate(Scheme{Match: 2, Mismatch: 1, Gap: -1}); err == nil {
		t.Error("mismatch>=0: want error")
	}
	if err := Validate(Scheme{Match: 2, Mismatch: -1, Gap: 0}); err == nil {
		t.Error("gap>=0: want error")
	}
	// affine validation
	if err := Validate(Scheme{Match: 2, Mismatch: -1, GapOpen: 0, GapExtend: -1}); err == nil {
		t.Error("affine gap_open>=0: want error")
	}
	if err := Validate(Scheme{Match: 2, Mismatch: -1, GapOpen: -5, GapExtend: 0}); err == nil {
		t.Error("affine gap_extend>=0: want error")
	}
}

func TestGapCost(t *testing.T) {
	linear := Default()
	if got := linear.GapCost(3); got != -6 {
		t.Errorf("linear gap cost(3): want -6, got %d", got)
	}
	if got := linear.GapCost(0); got != 0 {
		t.Errorf("linear gap cost(0): want 0, got %d", got)
	}
	affine := DefaultAffine()
	// open(-5) + 2*extend(-1) = -7
	if got := affine.GapCost(3); got != -7 {
		t.Errorf("affine gap cost(3): want -7, got %d", got)
	}
	if got := affine.GapCost(1); got != -5 {
		t.Errorf("affine gap cost(1): want -5, got %d", got)
	}
}

func TestIsAffine(t *testing.T) {
	if Default().IsAffine() {
		t.Error("Default should not be affine")
	}
	if !DefaultAffine().IsAffine() {
		t.Error("DefaultAffine should be affine")
	}
}
