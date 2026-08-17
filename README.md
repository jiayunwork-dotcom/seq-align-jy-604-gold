# seq-align

Sequence alignment and scoring CLI supporting Needleman-Wunsch global alignment
and Smith-Waterman local alignment with both linear and affine gap penalties
(Gotoh three-matrix algorithm). Outputs aligned sequences (with `-` gap markers),
raw score, percent identity, and aligned length.

For multi-sequence FASTA files, computes all-pairs distance matrices
(distance = 100 - identity%) in either global or local mode.

## Build

```sh
go build -o seq-align .
```

## Usage

```sh
# Pairwise alignment (default: global, linear gap)
seq-align pair ACGT AGT
seq-align pair ACGT AGT -mode local

# Custom linear scoring
seq-align pair ACGT AGT -m 3 -mm -2 -g -1

# Affine gap scoring (gap open + gap extend)
seq-align pair AACCGGTT ACGT -go -5 -ge -1

# FASTA all-pairs distance matrix (global mode)
seq-align matrix example/seqs.fasta

# FASTA distance matrix in local mode
seq-align matrix example/seqs.fasta -mode local

# Affine gap with distance matrix
seq-align matrix example/seqs.fasta -go -5 -ge -1 -mode global
```

## Flags

Flags apply to both `pair` and `matrix` commands:

| flag | meaning | default |
|------|---------|---------|
| `-mode` | `global` or `local` | `global` |
| `-m` | match score (positive) | `2` |
| `-mm` | mismatch penalty (negative) | `-1` |
| `-g` | linear gap penalty (negative) | `-2` |
| `-go` | affine gap open penalty (negative, 0 = disabled) | `0` |
| `-ge` | affine gap extend penalty (negative, 0 = disabled) | `0` |

When `-go` and `-ge` are both zero the linear `-g` penalty is used.
When either is non-zero, affine gap mode is activated and `-g` is ignored.

Exit codes: missing/unknown subcommand → 2; parameter or file error → 1.

## Package Structure

- `internal/scoring` — Scoring scheme (`Scheme`) with linear and affine gap
  support, validation, per-position rescoring, and gap cost calculation.
- `internal/align` — `Global` (Needleman-Wunsch) and `Local` (Smith-Waterman)
  with automatic selection between linear DP and Gotoh three-matrix DP based
  on the scoring scheme. Returns aligned strings, score, identity, aligned
  length, and start coordinates.
- `internal/fasta` — FASTA parsing (BOM-tolerant, multi-line sequences) and
  all-pairs distance matrix computation with configurable alignment mode.

## Testing

```sh
go test ./...
```

## License

MIT
