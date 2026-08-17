// Package matrix provides distance matrix serialization, deserialization,
// validation, and multi-format output. It operates on square distance matrices
// produced by pairwise sequence alignment.
package matrix

import (
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"strings"
)

// ErrNotSquare is returned when the matrix dimensions are inconsistent.
var ErrNotSquare = errors.New("matrix: not a square matrix")

// ErrIDCountMismatch is returned when len(IDs) != matrix dimension.
var ErrIDCountMismatch = errors.New("matrix: ID count does not match dimension")

// ErrNotSymmetric is returned when the matrix fails the symmetry check.
var ErrNotSymmetric = errors.New("matrix: not symmetric")

// ErrDiagonalNonZero is returned when a diagonal entry is not zero.
var ErrDiagonalNonZero = errors.New("matrix: diagonal entry is not zero")

// Matrix holds a square distance matrix with associated sequence IDs.
type Matrix struct {
	IDs  []string    `json:"ids"`
	Data [][]float64 `json:"data"`
}

// New creates a Matrix from IDs and a 2D distance slice.
// It validates that Data is square and len(IDs) == dimension.
func New(ids []string, data [][]float64) (*Matrix, error) {
	n := len(data)
	if len(ids) != n {
		return nil, ErrIDCountMismatch
	}
	for i, row := range data {
		if len(row) != n {
			return nil, fmt.Errorf("%w: row %d has length %d, want %d", ErrNotSquare, i, len(row), n)
		}
	}
	return &Matrix{IDs: ids, Data: data}, nil
}

// Dim returns the dimension (number of sequences) of the matrix.
func (m *Matrix) Dim() int {
	return len(m.IDs)
}

// Get returns the distance between sequence i and j.
func (m *Matrix) Get(i, j int) float64 {
	return m.Data[i][j]
}

// Validate checks structural invariants:
//   - Square shape
//   - ID count matches dimension
//   - Diagonal entries are zero
//   - Symmetric within tolerance
func (m *Matrix) Validate(tolerance float64) error {
	n := m.Dim()
	if len(m.Data) != n {
		return ErrNotSquare
	}
	for i, row := range m.Data {
		if len(row) != n {
			return fmt.Errorf("%w: row %d has length %d, want %d", ErrNotSquare, i, len(row), n)
		}
	}
	for i := 0; i < n; i++ {
		if math.Abs(m.Data[i][i]) > tolerance {
			return fmt.Errorf("%w: m[%d][%d] = %v", ErrDiagonalNonZero, i, i, m.Data[i][i])
		}
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if math.Abs(m.Data[i][j]-m.Data[j][i]) > tolerance {
				return fmt.Errorf("%w: m[%d][%d]=%.6f != m[%d][%d]=%.6f",
					ErrNotSymmetric, i, j, m.Data[i][j], j, i, m.Data[j][i])
			}
		}
	}
	return nil
}

// Closest returns the pair (i, j) with the smallest non-zero distance.
// If the matrix has fewer than 2 entries, returns (-1,-1, 0).
func (m *Matrix) Closest() (int, int, float64) {
	n := m.Dim()
	if n < 2 {
		return -1, -1, 0
	}
	bestI, bestJ := 0, 1
	bestD := m.Data[0][1]
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if m.Data[i][j] < bestD {
				bestD = m.Data[i][j]
				bestI, bestJ = i, j
			}
		}
	}
	return bestI, bestJ, bestD
}

// Farthest returns the pair (i, j) with the largest distance.
func (m *Matrix) Farthest() (int, int, float64) {
	n := m.Dim()
	if n < 2 {
		return -1, -1, 0
	}
	bestI, bestJ := 0, 1
	bestD := m.Data[0][1]
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if m.Data[i][j] > bestD {
				bestD = m.Data[i][j]
				bestI, bestJ = i, j
			}
		}
	}
	return bestI, bestJ, bestD
}

// MeanDistance returns the average of all upper-triangle distances.
// Returns 0 if n < 2.
func (m *Matrix) MeanDistance() float64 {
	n := m.Dim()
	if n < 2 {
		return 0
	}
	sum := 0.0
	count := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			sum += m.Data[i][j]
			count++
		}
	}
	return sum / float64(count)
}

// --- Serialization ---

// Format selects the output format for matrix rendering.
type Format int

const (
	FormatTSV  Format = iota // tab-separated values
	FormatCSV                // comma-separated values
	FormatJSON               // JSON object
	FormatPHYLIP             // PHYLIP lower-triangle format
)

// Render produces a string representation of the matrix in the given format.
func (m *Matrix) Render(f Format, precision int) (string, error) {
	switch f {
	case FormatTSV:
		return m.renderDelimited("\t", precision), nil
	case FormatCSV:
		return m.renderDelimited(",", precision), nil
	case FormatJSON:
		return m.renderJSON(precision)
	case FormatPHYLIP:
		return m.renderPHYLIP(precision), nil
	default:
		return "", fmt.Errorf("matrix: unknown format %d", f)
	}
}

func (m *Matrix) renderDelimited(sep string, prec int) string {
	var sb strings.Builder
	n := m.Dim()
	// Header row
	sb.WriteString(sep)
	sb.WriteString(strings.Join(m.IDs, sep))
	sb.WriteByte('\n')
	// Data rows
	for i := 0; i < n; i++ {
		sb.WriteString(m.IDs[i])
		for j := 0; j < n; j++ {
			sb.WriteString(sep)
			sb.WriteString(fmt.Sprintf("%.*f", prec, m.Data[i][j]))
		}
		sb.WriteByte('\n')
	}
	return sb.String()
}

func (m *Matrix) renderJSON(prec int) (string, error) {
	type jsonMatrix struct {
		IDs  []string    `json:"ids"`
		Data [][]float64 `json:"data"`
		Dim  int         `json:"dim"`
	}
	// Round data to specified precision for clean output
	rounded := make([][]float64, m.Dim())
	factor := math.Pow(10, float64(prec))
	for i, row := range m.Data {
		rounded[i] = make([]float64, len(row))
		for j, v := range row {
			rounded[i][j] = math.Round(v*factor) / factor
		}
	}
	out := jsonMatrix{IDs: m.IDs, Data: rounded, Dim: m.Dim()}
	b, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		return "", err
	}
	return string(b), nil
}

func (m *Matrix) renderPHYLIP(prec int) string {
	var sb strings.Builder
	n := m.Dim()
	sb.WriteString(fmt.Sprintf("%d\n", n))
	// PHYLIP lower-triangle format
	for i := 0; i < n; i++ {
		// Pad ID to 10 chars (PHYLIP convention)
		id := m.IDs[i]
		if len(id) > 10 {
			id = id[:10]
		}
		sb.WriteString(fmt.Sprintf("%-10s", id))
		for j := 0; j < i; j++ {
			sb.WriteString(fmt.Sprintf("  %.*f", prec, m.Data[i][j]))
		}
		sb.WriteByte('\n')
	}
	return sb.String()
}

// --- Deserialization ---

// ParseJSON reads a Matrix from a JSON byte slice.
func ParseJSON(data []byte) (*Matrix, error) {
	var m Matrix
	if err := json.Unmarshal(data, &m); err != nil {
		return nil, fmt.Errorf("matrix: parse JSON: %w", err)
	}
	if len(m.IDs) == 0 {
		return nil, fmt.Errorf("matrix: empty IDs in JSON")
	}
	if len(m.Data) != len(m.IDs) {
		return nil, ErrIDCountMismatch
	}
	for i, row := range m.Data {
		if len(row) != len(m.IDs) {
			return nil, fmt.Errorf("%w: row %d", ErrNotSquare, i)
		}
	}
	return &m, nil
}

// ParseTSV reads a Matrix from TSV format (header row + data rows).
func ParseTSV(text string) (*Matrix, error) {
	return parseDelimited(text, "\t")
}

// ParseCSV reads a Matrix from CSV format (header row + data rows).
func ParseCSV(text string) (*Matrix, error) {
	return parseDelimited(text, ",")
}

func parseDelimited(text, sep string) (*Matrix, error) {
	lines := strings.Split(strings.TrimRight(text, "\n\r"), "\n")
	if len(lines) < 2 {
		return nil, fmt.Errorf("matrix: need at least header + 1 data row")
	}
	// Header: first cell empty (or label), rest are IDs
	headerParts := strings.Split(lines[0], sep)
	if len(headerParts) < 2 {
		return nil, fmt.Errorf("matrix: header has too few columns")
	}
	ids := headerParts[1:] // skip first cell (empty corner)
	// Trim whitespace from IDs
	for i := range ids {
		ids[i] = strings.TrimSpace(ids[i])
	}
	n := len(ids)
	if len(lines)-1 != n {
		return nil, fmt.Errorf("matrix: expected %d data rows, got %d", n, len(lines)-1)
	}

	data := make([][]float64, n)
	for i := 0; i < n; i++ {
		parts := strings.Split(lines[i+1], sep)
		if len(parts) != n+1 {
			return nil, fmt.Errorf("matrix: row %d has %d columns, want %d", i, len(parts), n+1)
		}
		data[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			var v float64
			cell := strings.TrimSpace(parts[j+1])
			if _, err := fmt.Sscanf(cell, "%f", &v); err != nil {
				return nil, fmt.Errorf("matrix: parse cell [%d][%d] %q: %w", i, j, cell, err)
			}
			data[i][j] = v
		}
	}
	return &Matrix{IDs: ids, Data: data}, nil
}
