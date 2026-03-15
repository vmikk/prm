package prm

import (
	"fmt"
	"io"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/shenwei356/xopen"
)

func WriteHumanSummary(w io.Writer, summary Summary) error {
	if _, err := fmt.Fprintf(w,
		"prm summary\ninput1:\t%s\ninput2:\t%s\nforward:\t%s\nforward_rc:\t%s\nreverse:\t%s\nreverse_rc:\t%s\nmismatches:\t%d\nhead:\t%s\nthreads:\t%d\nconcrete_variants:\t%d\nFWD_variants:\t%d\nFWD_RC_variants:\t%d\nREV_variants:\t%d\nREV_RC_variants:\t%d\n\n",
		summary.Config.Input1,
		emptyToDash(summary.Config.Input2),
		summary.Forward,
		summary.ForwardRC,
		summary.Reverse,
		summary.ReverseRC,
		summary.Config.Mismatches,
		formatHeadValue(summary.Config.Head),
		summary.Config.Threads,
		summary.ConcreteVariants,
		summary.VariantCounts["FWD"],
		summary.VariantCounts["FWD_RC"],
		summary.VariantCounts["REV"],
		summary.VariantCounts["REV_RC"],
	); err != nil {
		return err
	}

	return writeMatchTable(w, summary)
}

func writeMatchTable(w io.Writer, summary Summary) error {
	if _, err := fmt.Fprintln(w, "primer matches"); err != nil {
		return err
	}
	tabWriter := tabwriter.NewWriter(w, 0, 0, 2, ' ', 0)
	if _, err := fmt.Fprintln(tabWriter, "Read\tFwd\tFwdRC\tRev\tRevRC"); err != nil {
		return err
	}

	for _, row := range matchRows(summary) {
		if _, err := fmt.Fprintf(
			tabWriter,
			"%s\t%d\t%d\t%d\t%d\n",
			row.Label,
			row.Values[0],
			row.Values[1],
			row.Values[2],
			row.Values[3],
		); err != nil {
			return err
		}
	}
	return tabWriter.Flush()
}

type orientationRow struct {
	Label  string
	Values [4]uint64
}

func matchRows(summary Summary) []orientationRow {
	if !summary.Paired {
		return []orientationRow{
			{Label: "SE", Values: readMatchCounts(summary.Read1Counts)},
		}
	}

	return []orientationRow{
		{Label: "R1", Values: readMatchCounts(summary.Read1Counts)},
		{Label: "R2", Values: readMatchCounts(summary.Read2Counts)},
	}
}

func readMatchCounts(states [16]uint64) [4]uint64 {
	return [4]uint64{
		sumStatesWithBit(states, bitFwd),
		sumStatesWithBit(states, bitFwdRC),
		sumStatesWithBit(states, bitRev),
		sumStatesWithBit(states, bitRevRC),
	}
}

func sumStatesWithBit(states [16]uint64, bit uint8) uint64 {
	var total uint64
	for state, count := range states {
		if uint8(state)&bit != 0 {
			total += count
		}
	}
	return total
}

func WriteTSV(path string, summary Summary) error {
	outFile, err := xopen.Wopen(path)
	if err != nil {
		return err
	}
	defer outFile.Close()

	headers, rows := tsvRows(summary)
	if _, err := outFile.WriteString(strings.Join(headers, "\t") + "\n"); err != nil {
		return err
	}
	for _, row := range rows {
		if _, err := outFile.WriteString(strings.Join(row, "\t") + "\n"); err != nil {
			return err
		}
	}
	return nil
}

func tsvRows(summary Summary) ([]string, [][]string) {
	headers := []string{
		"read",
		"file",
		"reads",
		"fwd",
		"fwd_rc",
		"rev",
		"rev_rc",
		"forward",
		"forward_rc",
		"reverse",
		"reverse_rc",
		"mismatches",
	}

	rows := [][]string{
		tsvCountRow(
			"SE",
			summary.Config.Input1,
			summary.Read1Total,
			readMatchCounts(summary.Read1Counts),
			summary,
		),
	}
	if summary.Paired {
		rows = [][]string{
			tsvCountRow(
				"R1",
				summary.Config.Input1,
				summary.Read1Total,
				readMatchCounts(summary.Read1Counts),
				summary,
			),
			tsvCountRow(
				"R2",
				summary.Config.Input2,
				summary.Read2Total,
				readMatchCounts(summary.Read2Counts),
				summary,
			),
		}
	}
	return headers, rows
}

func tsvCountRow(label string, file string, total uint64, counts [4]uint64, summary Summary) []string {
	return []string{
		label,
		file,
		strconv.FormatUint(total, 10),
		strconv.FormatUint(counts[0], 10),
		strconv.FormatUint(counts[1], 10),
		strconv.FormatUint(counts[2], 10),
		strconv.FormatUint(counts[3], 10),
		summary.Forward,
		summary.ForwardRC,
		summary.Reverse,
		summary.ReverseRC,
		strconv.Itoa(summary.Config.Mismatches),
	}
}

func formatHeadValue(value int) string {
	if value == 0 {
		return "all"
	}
	return strconv.Itoa(value)
}

func emptyToDash(value string) string {
	if value == "" {
		return "-"
	}
	return value
}
