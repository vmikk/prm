package prm

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/shenwei356/xopen"
)

func WriteHumanSummary(w io.Writer, summary Summary) error {
	palette := newReportPalette(isTerminalWriter(w))
	for _, line := range humanSummaryLines(summary, palette) {
		if _, err := fmt.Fprintln(w, line); err != nil {
			return err
		}
	}
	return nil
}

type summaryLine struct {
	Label string
	Value string
}

type reportPalette struct {
	enabled bool
}

type rgb struct {
	r int
	g int
	b int
}

var (
	labelColor   = rgb{135, 148, 166}
	valueColor   = rgb{120, 176, 220}
	headerColor  = rgb{148, 163, 184}
	zeroColor    = rgb{184, 194, 204}
	nonZeroColor = rgb{236, 164, 122}
	sectionColor = rgb{125, 211, 192}
	iupacColors  = map[byte]rgb{
		'R': {244, 162, 97},
		'Y': {255, 209, 102},
		'M': {130, 205, 255},
		'K': {255, 173, 173},
		'S': {162, 210, 255},
		'W': {202, 240, 248},
		'H': {255, 200, 221},
		'B': {186, 230, 126},
		'V': {189, 178, 255},
		'D': {255, 214, 165},
		'N': {201, 214, 255},
		'I': {255, 199, 206},
	}
)

func newReportPalette(enabled bool) reportPalette {
	return reportPalette{enabled: enabled}
}

func (p reportPalette) color(text string, shade rgb) string {
	if !p.enabled || text == "" {
		return text
	}
	return fmt.Sprintf("\x1b[38;2;%d;%d;%dm%s\x1b[0m", shade.r, shade.g, shade.b, text)
}

func (p reportPalette) label(text string) string {
	return p.color(text, labelColor)
}

func (p reportPalette) value(text string) string {
	return p.color(text, valueColor)
}

func (p reportPalette) header(text string) string {
	return p.color(text, headerColor)
}

func (p reportPalette) section(text string) string {
	return p.color(text, sectionColor)
}

func (p reportPalette) count(text string, value uint64) string {
	if value == 0 {
		return p.color(text, zeroColor)
	}
	return p.color(text, nonZeroColor)
}

func (p reportPalette) primer(seq string) string {
	if !p.enabled {
		return seq
	}
	var out strings.Builder
	out.Grow(len(seq) * 12)
	for i := 0; i < len(seq); i++ {
		base := seq[i]
		shade, ok := iupacColors[base]
		if !ok {
			out.WriteByte(base)
			continue
		}
		out.WriteString(p.color(string(base), shade))
	}
	return out.String()
}

func (p reportPalette) primerPair(seq string, rc string) string {
	return p.primer(seq) + p.label(" / ") + p.primer(rc)
}

func isTerminalWriter(w io.Writer) bool {
	file, ok := w.(*os.File)
	if !ok {
		return false
	}
	info, err := file.Stat()
	if err != nil {
		return false
	}
	term := os.Getenv("TERM")
	return info.Mode()&os.ModeCharDevice != 0 && term != "" && term != "dumb"
}

func humanSummaryLines(summary Summary, palette reportPalette) []string {
	lines := renderSummaryBlock(summaryLines(summary, palette), palette)
	lines = append(lines, "")
	lines = append(lines, renderMatchTable(matchRows(summary), palette)...)
	return lines
}

func summaryLines(summary Summary, palette reportPalette) []summaryLine {
	lines := []summaryLine{
		{Label: "input1", Value: palette.value(summary.Config.Input1)},
	}
	if summary.Paired {
		lines = append(lines, summaryLine{Label: "input2", Value: palette.value(summary.Config.Input2)})
	}
	lines = append(lines,
		summaryLine{Label: "forward / rc ", Value: palette.primerPair(summary.Forward, summary.ForwardRC)},
		summaryLine{Label: "reverse / rc", Value: palette.primerPair(summary.Reverse, summary.ReverseRC)},
		summaryLine{Label: "mismatches", Value: palette.value(strconv.Itoa(summary.Config.Mismatches))},
	)
	if shouldShowVariants(summary) {
		lines = append(lines,
			summaryLine{Label: "primer_variants", Value: palette.value(strconv.Itoa(summary.ConcreteVariants))},
			summaryLine{Label: "fwd_variants", Value: palette.value(strconv.Itoa(summary.VariantCounts["FWD"]))},
			summaryLine{Label: "fwd_rc_variants", Value: palette.value(strconv.Itoa(summary.VariantCounts["FWD_RC"]))},
			summaryLine{Label: "rev_variants", Value: palette.value(strconv.Itoa(summary.VariantCounts["REV"]))},
			summaryLine{Label: "rev_rc_variants", Value: palette.value(strconv.Itoa(summary.VariantCounts["REV_RC"]))},
		)
	}
	return lines
}

func renderSummaryBlock(lines []summaryLine, palette reportPalette) []string {
	width := 0
	for _, line := range lines {
		if len(line.Label) > width {
			width = len(line.Label)
		}
	}

	rendered := make([]string, 0, len(lines))
	for _, line := range lines {
		label := fmt.Sprintf("%-*s", width+1, line.Label+":")
		rendered = append(rendered, fmt.Sprintf("%s %s", palette.label(label), line.Value))
	}
	return rendered
}

func renderMatchTable(rows []orientationRow, palette reportPalette) []string {
	headers := []string{"Read", "Fwd", "FwdRC", "Rev", "RevRC"}
	widths := []int{
		len(headers[0]),
		len(headers[1]),
		len(headers[2]),
		len(headers[3]),
		len(headers[4]),
	}
	for _, row := range rows {
		if len(row.Label) > widths[0] {
			widths[0] = len(row.Label)
		}
		for idx, value := range row.Values {
			if n := len(strconv.FormatUint(value, 10)); n > widths[idx+1] {
				widths[idx+1] = n
			}
		}
	}

	rendered := make([]string, 0, len(rows)+2)
	rendered = append(rendered, palette.section("primer matches"))
	rendered = append(rendered,
		fmt.Sprintf(
			"%-*s  %*s  %*s  %*s  %*s",
			widths[0], palette.header(headers[0]),
			widths[1], palette.header(headers[1]),
			widths[2], palette.header(headers[2]),
			widths[3], palette.header(headers[3]),
			widths[4], palette.header(headers[4]),
		),
	)
	for _, row := range rows {
		rendered = append(rendered,
			fmt.Sprintf(
				"%-*s  %*s  %*s  %*s  %*s",
				widths[0], palette.label(row.Label),
				widths[1], palette.count(strconv.FormatUint(row.Values[0], 10), row.Values[0]),
				widths[2], palette.count(strconv.FormatUint(row.Values[1], 10), row.Values[1]),
				widths[3], palette.count(strconv.FormatUint(row.Values[2], 10), row.Values[2]),
				widths[4], palette.count(strconv.FormatUint(row.Values[3], 10), row.Values[3]),
			),
		)
	}
	return rendered
}

func shouldShowVariants(summary Summary) bool {
	return summary.Config.Mismatches > 0 || summary.HasIUPAC
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
