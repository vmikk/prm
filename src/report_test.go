package prm

import (
	"bytes"
	"regexp"
	"strings"
	"testing"
)

func TestOrientationRowsSingle(t *testing.T) {
	summary := Summary{}
	summary.Read1Counts[int(bitFwd)] = 2
	summary.Read1Counts[int(bitFwdRC)] = 3
	summary.Read1Counts[int(bitRev|bitRevRC)] = 5
	summary.Read1Counts[int(bitFwd|bitFwdRC|bitRev|bitRevRC)] = 7

	rows := matchRows(summary)
	if len(rows) != 1 {
		t.Fatalf("expected 1 row, got %d", len(rows))
	}

	assertOrientationRow(t, rows[0], "SE", [4]uint64{9, 10, 12, 12})
}

func TestOrientationRowsPaired(t *testing.T) {
	summary := Summary{Paired: true}
	summary.Read1Counts[int(bitFwd)] = 2
	summary.Read1Counts[int(bitFwdRC)] = 1
	summary.Read1Counts[int(bitFwd|bitFwdRC|bitRev|bitRevRC)] = 1
	summary.Read2Counts[int(bitRev)] = 3
	summary.Read2Counts[int(bitRevRC)] = 4
	summary.Read2Counts[int(bitFwd|bitFwdRC|bitRev|bitRevRC)] = 2

	rows := matchRows(summary)
	if len(rows) != 2 {
		t.Fatalf("expected 2 rows, got %d", len(rows))
	}

	assertOrientationRow(t, rows[0], "R1", [4]uint64{3, 2, 1, 1})
	assertOrientationRow(t, rows[1], "R2", [4]uint64{2, 2, 5, 6})
}

func TestWriteHumanSummaryPlainLayout(t *testing.T) {
	summary := Summary{
		Config: Config{
			Input1:     "reads_R1.fq.gz",
			Input2:     "reads_R2.fq.gz",
			Mismatches: 0,
			Head:       5,
			Threads:    8,
		},
		Paired:    true,
		Forward:   "ACGT",
		ForwardRC: "ACGT",
		Reverse:   "TGCA",
		ReverseRC: "TGCA",
	}
	summary.Read1Counts[int(bitFwd)] = 11
	summary.Read1Counts[int(bitFwdRC)] = 22
	summary.Read1Counts[int(bitRev)] = 33
	summary.Read1Counts[int(bitRevRC)] = 44
	summary.Read2Counts[int(bitFwd)] = 55
	summary.Read2Counts[int(bitFwdRC)] = 66
	summary.Read2Counts[int(bitRev)] = 77
	summary.Read2Counts[int(bitRevRC)] = 88

	var out bytes.Buffer
	if err := WriteHumanSummary(&out, summary); err != nil {
		t.Fatalf("write summary: %v", err)
	}

	text := out.String()
	if strings.Contains(text, "prm summary") {
		t.Fatalf("unexpected legacy header in summary:\n%s", text)
	}
	if strings.Contains(text, "head:") {
		t.Fatalf("unexpected head line in summary:\n%s", text)
	}
	if strings.Contains(text, "threads:") {
		t.Fatalf("unexpected threads line in summary:\n%s", text)
	}
	if strings.Contains(text, "primer_variants:") {
		t.Fatalf("unexpected variants in exact-match summary:\n%s", text)
	}
	if strings.Contains(text, "\x1b[") {
		t.Fatalf("unexpected ANSI escapes in plain writer output:\n%s", text)
	}

	lines := strings.Split(strings.TrimRight(text, "\n"), "\n")
	valueColumn := valueColumnIndex(t, lines, "input1:", "reads_R1.fq.gz")
	for _, item := range []struct {
		label string
		value string
	}{
		{label: "input2:", value: "reads_R2.fq.gz"},
		{label: "forward / rc:", value: "ACGT / ACGT"},
		{label: "reverse / rc:", value: "TGCA / TGCA"},
		{label: "mismatches:", value: "0"},
	} {
		if got := valueColumnIndex(t, lines, item.label, item.value); got != valueColumn {
			t.Fatalf("value column misaligned for %s: got %d want %d", item.label, got, valueColumn)
		}
	}
	if strings.Contains(text, "forward_rc:") || strings.Contains(text, "reverse_rc:") {
		t.Fatalf("unexpected separate rc lines in summary:\n%s", text)
	}
	assertTableTokenColumn(t, lines, "R1", "11", "R2", "55")
	assertTableTokenColumn(t, lines, "R1", "22", "R2", "66")
	assertTableTokenColumn(t, lines, "R1", "33", "R2", "77")
	assertTableTokenColumn(t, lines, "R1", "44", "R2", "88")
}

func TestWriteHumanSummaryShowsVariantsForMismatches(t *testing.T) {
	summary := Summary{
		Config: Config{
			Input1:     "reads.fq.gz",
			Mismatches: 1,
		},
		Forward:          "ACGT",
		ForwardRC:        "ACGT",
		Reverse:          "TGCA",
		ReverseRC:        "TGCA",
		ConcreteVariants: 12,
		VariantCounts: map[string]int{
			"FWD":    3,
			"FWD_RC": 3,
			"REV":    3,
			"REV_RC": 3,
		},
	}

	var out bytes.Buffer
	if err := WriteHumanSummary(&out, summary); err != nil {
		t.Fatalf("write summary: %v", err)
	}

	text := out.String()
	if !strings.Contains(text, "primer variants (fwd + rev): 12 = (3 + 3) * 2") {
		t.Fatalf("missing combined variant summary:\n%s", text)
	}
	for _, token := range []string{
		"primer_variants:",
		"fwd_variants:",
		"fwd_rc_variants:",
		"rev_variants:",
		"rev_rc_variants:",
	} {
		if strings.Contains(text, token) {
			t.Fatalf("unexpected legacy variant token %s in summary:\n%s", token, text)
		}
	}
}

func TestWriteHumanSummaryShowsVariantsForIUPAC(t *testing.T) {
	summary := Summary{
		Config: Config{
			Input1:     "reads.fq.gz",
			Mismatches: 0,
		},
		Forward:          "ATGR",
		ForwardRC:        "YCAT",
		Reverse:          "CCCY",
		ReverseRC:        "RGGG",
		HasIUPAC:         true,
		ConcreteVariants: 8,
		VariantCounts: map[string]int{
			"FWD":    2,
			"FWD_RC": 2,
			"REV":    2,
			"REV_RC": 2,
		},
	}

	var out bytes.Buffer
	if err := WriteHumanSummary(&out, summary); err != nil {
		t.Fatalf("write summary: %v", err)
	}

	text := out.String()
	if !strings.Contains(text, "primer variants (fwd + rev): 8 = (2 + 2) * 2") {
		t.Fatalf("expected combined variants for IUPAC summary:\n%s", text)
	}
}

func TestPrimerStylingColorsCanonicalAndAmbiguousBases(t *testing.T) {
	rendered := newReportPalette(true).primer("ATGRYN")
	if !strings.HasPrefix(rendered, "\x1b[") {
		t.Fatalf("expected canonical bases to be colored, got %q", rendered)
	}
	palette := newReportPalette(true)
	for _, base := range []string{"A", "T", "G"} {
		if !strings.Contains(rendered, palette.color(base, valueColor)) {
			t.Fatalf("expected canonical base %s to use value color, got %q", base, rendered)
		}
	}
	if got := strings.Count(rendered, "\x1b[0m"); got != 6 {
		t.Fatalf("expected 6 colored bases, got %d in %q", got, rendered)
	}
}

func TestPrimerPairStylingCombinesPrimerAndRC(t *testing.T) {
	rendered := newReportPalette(false).primerPair("ATGR", "YCAT")
	if rendered != "ATGR / YCAT" {
		t.Fatalf("unexpected primer pair rendering: %q", rendered)
	}
}

func TestRenderMatchTableAlignsColoredOutput(t *testing.T) {
	rows := []orientationRow{
		{Label: "R1", Values: [4]uint64{1234, 22, 33333, 4444}},
		{Label: "R2", Values: [4]uint64{5678, 90, 77777, 8888}},
	}

	lines := renderMatchTable(rows, newReportPalette(true))
	text := strings.Join(lines, "\n")
	if !strings.Contains(text, "\x1b[") {
		t.Fatalf("expected ANSI escapes in colored table output:\n%s", text)
	}

	plain := stripANSI(text)
	plainLines := strings.Split(plain, "\n")
	assertTableTokenColumn(t, plainLines, "R1", "1234", "R2", "5678")
	assertTableTokenColumn(t, plainLines, "R1", "22", "R2", "90")
	assertTableTokenColumn(t, plainLines, "R1", "33333", "R2", "77777")
	assertTableTokenColumn(t, plainLines, "R1", "4444", "R2", "8888")
}

func assertOrientationRow(t *testing.T, row orientationRow, wantLabel string, want [4]uint64) {
	t.Helper()
	if row.Label != wantLabel {
		t.Fatalf("unexpected row label: got %s want %s", row.Label, wantLabel)
	}
	if row.Values != want {
		t.Fatalf("unexpected row values for %s: got %v want %v", wantLabel, row.Values, want)
	}
}

func valueColumnIndex(t *testing.T, lines []string, label string, value string) int {
	t.Helper()
	for _, line := range lines {
		if !strings.Contains(line, label) {
			continue
		}
		idx := strings.Index(line, value)
		if idx < 0 {
			t.Fatalf("missing value %q in line %q", value, line)
		}
		return idx
	}
	t.Fatalf("missing label %q", label)
	return -1
}

func assertTableTokenColumn(t *testing.T, lines []string, row1 string, value1 string, row2 string, value2 string) {
	t.Helper()
	var first, second int = -1, -1
	for _, line := range lines {
		if strings.HasPrefix(line, row1) {
			first = strings.Index(line, value1)
		}
		if strings.HasPrefix(line, row2) {
			second = strings.Index(line, value2)
		}
	}
	if first < 0 || second < 0 {
		t.Fatalf("missing table tokens %q/%q", value1, value2)
	}
	if first != second {
		t.Fatalf("table column misaligned: %q at %d, %q at %d", value1, first, value2, second)
	}
}

var ansiPattern = regexp.MustCompile(`\x1b\[[0-9;]*m`)

func stripANSI(text string) string {
	return ansiPattern.ReplaceAllString(text, "")
}
