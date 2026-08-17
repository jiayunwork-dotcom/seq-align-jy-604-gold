package matrix

import (
	"math"
	"strings"
	"testing"
)

func sampleMatrix() *Matrix {
	return &Matrix{
		IDs: []string{"a", "b", "c"},
		Data: [][]float64{
			{0, 10.5, 25.0},
			{10.5, 0, 18.3},
			{25.0, 18.3, 0},
		},
	}
}

func TestNewValid(t *testing.T) {
	ids := []string{"x", "y"}
	data := [][]float64{{0, 5}, {5, 0}}
	m, err := New(ids, data)
	if err != nil {
		t.Fatal(err)
	}
	if m.Dim() != 2 {
		t.Errorf("Dim = %d, want 2", m.Dim())
	}
}

func TestNewErrors(t *testing.T) {
	if _, err := New([]string{"a"}, [][]float64{{0, 1}, {1, 0}}); err == nil {
		t.Error("ID count mismatch: want error")
	}
	if _, err := New([]string{"a", "b"}, [][]float64{{0, 1}, {1}}); err == nil {
		t.Error("non-square row: want error")
	}
}

func TestValidateSymmetric(t *testing.T) {
	m := sampleMatrix()
	if err := m.Validate(1e-9); err != nil {
		t.Errorf("valid matrix failed Validate: %v", err)
	}
}

func TestValidateNotSymmetric(t *testing.T) {
	m := sampleMatrix()
	m.Data[0][1] = 99.0 // break symmetry
	err := m.Validate(1e-9)
	if err == nil {
		t.Error("expected ErrNotSymmetric")
	}
	if !strings.Contains(err.Error(), "not symmetric") {
		t.Errorf("unexpected error: %v", err)
	}
}

func TestValidateDiagonalNonZero(t *testing.T) {
	m := sampleMatrix()
	m.Data[1][1] = 0.5
	err := m.Validate(1e-9)
	if err == nil {
		t.Error("expected ErrDiagonalNonZero")
	}
}

func TestClosest(t *testing.T) {
	m := sampleMatrix()
	i, j, d := m.Closest()
	if i != 0 || j != 1 || d != 10.5 {
		t.Errorf("Closest = (%d,%d,%.1f), want (0,1,10.5)", i, j, d)
	}
}

func TestFarthest(t *testing.T) {
	m := sampleMatrix()
	i, j, d := m.Farthest()
	if i != 0 || j != 2 || d != 25.0 {
		t.Errorf("Farthest = (%d,%d,%.1f), want (0,2,25.0)", i, j, d)
	}
}

func TestMeanDistance(t *testing.T) {
	m := sampleMatrix()
	mean := m.MeanDistance()
	expected := (10.5 + 25.0 + 18.3) / 3.0
	if math.Abs(mean-expected) > 1e-9 {
		t.Errorf("MeanDistance = %v, want %v", mean, expected)
	}
}

func TestRenderTSV(t *testing.T) {
	m := sampleMatrix()
	out, err := m.Render(FormatTSV, 2)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(out), "\n")
	if len(lines) != 4 { // header + 3 rows
		t.Fatalf("TSV should have 4 lines, got %d", len(lines))
	}
	if !strings.Contains(lines[0], "a\tb\tc") {
		t.Errorf("header missing IDs: %s", lines[0])
	}
}

func TestRenderCSV(t *testing.T) {
	m := sampleMatrix()
	out, err := m.Render(FormatCSV, 1)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out, ",") {
		t.Error("CSV output should contain commas")
	}
}

func TestRenderJSON(t *testing.T) {
	m := sampleMatrix()
	out, err := m.Render(FormatJSON, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out, `"ids"`) || !strings.Contains(out, `"data"`) {
		t.Errorf("JSON missing expected fields: %s", out[:100])
	}
}

func TestRenderPHYLIP(t *testing.T) {
	m := sampleMatrix()
	out, err := m.Render(FormatPHYLIP, 2)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(out), "\n")
	if lines[0] != "3" {
		t.Errorf("PHYLIP first line should be dimension: got %q", lines[0])
	}
	// First data line (row 0) has no distances (lower triangle)
	if !strings.HasPrefix(lines[1], "a") {
		t.Errorf("PHYLIP row 0: %q", lines[1])
	}
}

func TestParseJSON(t *testing.T) {
	m := sampleMatrix()
	out, _ := m.Render(FormatJSON, 2)
	parsed, err := ParseJSON([]byte(out))
	if err != nil {
		t.Fatal(err)
	}
	if parsed.Dim() != 3 {
		t.Errorf("Dim = %d, want 3", parsed.Dim())
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(parsed.Data[i][j]-m.Data[i][j]) > 0.01 {
				t.Errorf("[%d][%d] = %v, want %v", i, j, parsed.Data[i][j], m.Data[i][j])
			}
		}
	}
}

func TestParseTSVRoundTrip(t *testing.T) {
	m := sampleMatrix()
	out, _ := m.Render(FormatTSV, 4)
	parsed, err := ParseTSV(out)
	if err != nil {
		t.Fatal(err)
	}
	if parsed.Dim() != m.Dim() {
		t.Fatalf("Dim mismatch: %d vs %d", parsed.Dim(), m.Dim())
	}
	for i := 0; i < m.Dim(); i++ {
		if parsed.IDs[i] != m.IDs[i] {
			t.Errorf("ID[%d] = %q, want %q", i, parsed.IDs[i], m.IDs[i])
		}
		for j := 0; j < m.Dim(); j++ {
			if math.Abs(parsed.Data[i][j]-m.Data[i][j]) > 0.001 {
				t.Errorf("[%d][%d] = %v, want %v", i, j, parsed.Data[i][j], m.Data[i][j])
			}
		}
	}
}

func TestParseCSVRoundTrip(t *testing.T) {
	m := sampleMatrix()
	out, _ := m.Render(FormatCSV, 3)
	parsed, err := ParseCSV(out)
	if err != nil {
		t.Fatal(err)
	}
	if err := parsed.Validate(0.01); err != nil {
		t.Errorf("parsed CSV not valid: %v", err)
	}
}

func TestParseErrors(t *testing.T) {
	if _, err := ParseTSV(""); err == nil {
		t.Error("empty: want error")
	}
	if _, err := ParseTSV("\ta\tb\na\t0\t1\n"); err == nil {
		t.Error("missing row: want error")
	}
	if _, err := ParseJSON([]byte(`{"ids":[],"data":[]}`)); err == nil {
		t.Error("empty ids: want error")
	}
	if _, err := ParseJSON([]byte(`not json`)); err == nil {
		t.Error("invalid json: want error")
	}
}

func TestSmallMatrix(t *testing.T) {
	m := &Matrix{IDs: []string{"x"}, Data: [][]float64{{0}}}
	i, j, _ := m.Closest()
	if i != -1 || j != -1 {
		t.Error("Closest on 1x1 should return -1,-1")
	}
	if m.MeanDistance() != 0 {
		t.Error("MeanDistance on 1x1 should be 0")
	}
}
