package fasta

import (
	"math"
	"os"
	"path/filepath"
	"testing"

	"seq-align/internal/scoring"
)

func writeFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestFastaParse(t *testing.T) {
	p := writeFile(t, "in.fasta", ">seqA first record\nACGTAC\nGTAC\n>seqB\nACGTACGG\n")
	recs, err := Parse(p)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("want 2 records, got %d", len(recs))
	}
	if recs[0].ID != "seqA" || recs[0].Seq != "ACGTACGTAC" {
		t.Fatalf("record 0: %+v", recs[0])
	}
	if recs[1].ID != "seqB" || recs[1].Seq != "ACGTACGG" {
		t.Fatalf("record 1: %+v", recs[1])
	}
}

func TestFastaParseErrors(t *testing.T) {
	if _, err := Parse(writeFile(t, "nohead.fasta", "ACGTACGT\n")); err == nil {
		t.Error("no header: want error")
	}
	if _, err := Parse(writeFile(t, "empty.fasta", ">seqA\n\n")); err == nil {
		t.Error("empty seq: want error")
	}
	if _, err := Parse(writeFile(t, "bom.fasta", "\ufeff>seqA\nACGT\n")); err != nil {
		t.Errorf("BOM should be tolerated: %v", err)
	}
	if _, err := Parse("/nonexistent.fasta"); err == nil {
		t.Error("missing file: want error")
	}
}

func TestDistanceMatrix(t *testing.T) {
	p := writeFile(t, "in.fasta", ">a\nACGTACGT\n>b\nACGTACGT\n>c\nTTTTTTTT\n")
	recs, err := Parse(p)
	if err != nil {
		t.Fatal(err)
	}
	m, ids, err := DistanceMatrix(recs, scoring.Default(), ModeGlobal)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(ids) != 3 || len(m) != 3 || len(m[0]) != 3 {
		t.Fatalf("want 3x3 matrix, got %v", m)
	}
	for i := 0; i < 3; i++ {
		if m[i][i] != 0 {
			t.Errorf("diagonal [%d][%d] = %v, want 0", i, i, m[i][i])
		}
	}
	if m[0][1] != 0 || m[1][0] != 0 {
		t.Errorf("a vs b identical, want distance 0, got %v/%v", m[0][1], m[1][0])
	}
	if m[0][2] <= 0 || m[0][2] != m[2][0] {
		t.Errorf("a vs c want positive symmetric distance, got %v/%v", m[0][2], m[2][0])
	}
	if _, _, err := DistanceMatrix(recs[:1], scoring.Default(), ModeGlobal); err == nil {
		t.Error("<2 records: want error")
	}
}

func TestDistanceMatrixSymmetry(t *testing.T) {
	p := writeFile(t, "sym.fasta", ">x\nACGTACGTAA\n>y\nACGTTTTTAA\n>z\nGGGGACGTAA\n")
	recs, err := Parse(p)
	if err != nil {
		t.Fatal(err)
	}
	m, _, err := DistanceMatrix(recs, scoring.Default(), ModeGlobal)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < len(m); i++ {
		for j := 0; j < len(m); j++ {
			if math.Abs(m[i][j]-m[j][i]) > 1e-9 {
				t.Errorf("not symmetric: m[%d][%d]=%.6f != m[%d][%d]=%.6f", i, j, m[i][j], j, i, m[j][i])
			}
		}
	}
}

func TestDistanceMatrixLocalMode(t *testing.T) {
	// With local mode, identical subsequences should yield distance 0 even if flanking differs.
	p := writeFile(t, "local.fasta", ">a\nXXXACGTACGTXXX\n>b\nYYYACGTACGTYYY\n")
	recs, err := Parse(p)
	if err != nil {
		t.Fatal(err)
	}
	mLocal, _, err := DistanceMatrix(recs, scoring.Default(), ModeLocal)
	if err != nil {
		t.Fatal(err)
	}
	mGlobal, _, err := DistanceMatrix(recs, scoring.Default(), ModeGlobal)
	if err != nil {
		t.Fatal(err)
	}
	// Local should find perfect match on ACGTACGT region -> distance close to 0
	// Global will have mismatches on flanking X/Y regions -> larger distance
	if mLocal[0][1] >= mGlobal[0][1] {
		t.Errorf("local distance (%v) should be less than global distance (%v) for flanked sequences",
			mLocal[0][1], mGlobal[0][1])
	}
}

func TestParseMode(t *testing.T) {
	if m, err := ParseMode("global"); err != nil || m != ModeGlobal {
		t.Errorf("ParseMode(global): %v, %v", m, err)
	}
	if m, err := ParseMode("local"); err != nil || m != ModeLocal {
		t.Errorf("ParseMode(local): %v, %v", m, err)
	}
	if m, err := ParseMode(""); err != nil || m != ModeGlobal {
		t.Errorf("ParseMode(empty): %v, %v", m, err)
	}
	if _, err := ParseMode("invalid"); err == nil {
		t.Error("ParseMode(invalid): want error")
	}
}
