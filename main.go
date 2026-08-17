// seq-align: sequence alignment and scoring CLI.
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"seq-align/internal/align"
	"seq-align/internal/fasta"
	"seq-align/internal/scoring"
)

const usageText = `usage: seq-align <command> [flags] [args]

commands:
  pair <seqA> <seqB>   align two sequences
  matrix <fasta>       all-pairs alignment distance matrix

flags:
  -mode global|local   alignment mode (default global)
  -m int               match score (default 2)
  -mm int              mismatch score (default -1)
  -g int               gap penalty, negative (default -2, linear mode)
  -go int              gap open penalty (default 0, enables affine when set)
  -ge int              gap extend penalty (default 0, enables affine when set)
`

func reorder(args []string) []string {
	var flags, rest []string
	i := 0
	for i < len(args) {
		if strings.HasPrefix(args[i], "-") && args[i] != "-" && args[i] != "--" {
			flags = append(flags, args[i])
			if i+1 < len(args) && !strings.HasPrefix(args[i+1], "-") &&
				!strings.Contains(args[i], "=") &&
				needsValue(args[i]) {
				flags = append(flags, args[i+1])
				i += 2
				continue
			}
			i++
			continue
		}
		rest = append(rest, args[i])
		i++
	}
	return append(flags, rest...)
}

func needsValue(name string) bool {
	n := strings.TrimLeft(name, "-")
	return n == "m" || n == "mm" || n == "g" || n == "go" || n == "ge" || n == "mode"
}

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, usageText)
		os.Exit(2)
	}
	cmd := os.Args[1]
	var err error
	switch cmd {
	case "pair":
		err = runPair(reorder(os.Args[2:]))
	case "matrix":
		err = runMatrix(reorder(os.Args[2:]))
	case "-h", "--help", "help":
		fmt.Print(usageText)
	default:
		fmt.Fprintf(os.Stderr, "seq-align: unknown command %q\n\n%s", cmd, usageText)
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "seq-align: %v\n", err)
		os.Exit(1)
	}
}

func parseScheme(fs *flag.FlagSet) scoring.Scheme {
	s := scoring.Scheme{Match: 2, Mismatch: -1, Gap: -2}
	fs.IntVar(&s.Match, "m", s.Match, "match score")
	fs.IntVar(&s.Mismatch, "mm", s.Mismatch, "mismatch score")
	fs.IntVar(&s.Gap, "g", s.Gap, "gap penalty (negative, linear mode)")
	fs.IntVar(&s.GapOpen, "go", 0, "gap open penalty (negative, enables affine)")
	fs.IntVar(&s.GapExtend, "ge", 0, "gap extend penalty (negative, enables affine)")
	return s
}

func runPair(args []string) error {
	fs := flag.NewFlagSet("pair", flag.ContinueOnError)
	s := parseScheme(fs)
	mode := fs.String("mode", "global", "alignment mode: global|local")
	if err := fs.Parse(args); err != nil {
		return err
	}
	// Re-read values after Parse (flag pointers point to s fields)
	s.GapOpen = lookupInt(fs, "go")
	s.GapExtend = lookupInt(fs, "ge")
	if err := scoring.Validate(s); err != nil {
		return err
	}
	m, err := fasta.ParseMode(*mode)
	if err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 2 {
		return fmt.Errorf("pair requires exactly 2 sequences, got %d", len(rest))
	}
	var res align.Result
	switch m {
	case fasta.ModeLocal:
		res, err = align.Local(rest[0], rest[1], s)
	default:
		res, err = align.Global(rest[0], rest[1], s)
	}
	if err != nil {
		return err
	}
	fmt.Println(res.A)
	fmt.Println(res.B)
	fmt.Printf("score=%d identity=%.1f%% aligned_len=%d mode=%s\n", res.Score, res.Identity, res.AlignedLen, *mode)
	return nil
}

func runMatrix(args []string) error {
	fs := flag.NewFlagSet("matrix", flag.ContinueOnError)
	s := parseScheme(fs)
	mode := fs.String("mode", "global", "alignment mode: global|local")
	if err := fs.Parse(args); err != nil {
		return err
	}
	s.GapOpen = lookupInt(fs, "go")
	s.GapExtend = lookupInt(fs, "ge")
	if err := scoring.Validate(s); err != nil {
		return err
	}
	m, err := fasta.ParseMode(*mode)
	if err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 1 {
		return fmt.Errorf("matrix requires exactly 1 FASTA file, got %d", len(rest))
	}
	recs, err := fasta.Parse(rest[0])
	if err != nil {
		return err
	}
	mat, ids, err := fasta.DistanceMatrix(recs, s, m)
	if err != nil {
		return err
	}
	fmt.Print("\t")
	fmt.Println(strings.Join(ids, "\t"))
	for i, row := range mat {
		cells := make([]string, len(row))
		for j, v := range row {
			cells[j] = fmt.Sprintf("%.2f", v)
		}
		fmt.Printf("%s\t%s\n", ids[i], strings.Join(cells, "\t"))
	}
	return nil
}

func lookupInt(fs *flag.FlagSet, name string) int {
	f := fs.Lookup(name)
	if f == nil {
		return 0
	}
	if g, ok := f.Value.(flag.Getter); ok {
		if v, ok := g.Get().(int); ok {
			return v
		}
	}
	return 0
}
