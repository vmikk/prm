package main

import (
	"bytes"
	"compress/gzip"
	"encoding/csv"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSingleEndTSVAndSummary(t *testing.T) {
	dir := t.TempDir()
	outTSV := filepath.Join(dir, "single.tsv")

	stdout, _ := executeCommand(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", filepath.Join("testdata", "single_states.fa"),
		"--out-tsv", outTSV,
	)

	assertSummaryLineTokens(t, stdout, []string{"Read", "Fwd", "FwdRC", "Rev", "RevRC"})
	assertSummaryLineTokens(t, stdout, []string{"SE", "8", "8", "8", "8"})
	rows := readTSVRows(t, outTSV)
	if len(rows) != 1 {
		t.Fatalf("expected 1 TSV row, got %d", len(rows))
	}
	row := rows[0]
	assertRowValue(t, row, "read", "SE")
	assertRowValue(t, row, "file", filepath.Join("testdata", "single_states.fa"))
	assertRowValue(t, row, "reads", "16")
	assertRowValue(t, row, "fwd", "8")
	assertRowValue(t, row, "fwd_rc", "8")
	assertRowValue(t, row, "rev", "8")
	assertRowValue(t, row, "rev_rc", "8")
}

func TestSingleEndCompressedFASTQ(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "single_states.fq.gz")
	gzipFixture(t, filepath.Join("testdata", "single_states.fq"), in)

	outTSV := filepath.Join(dir, "single.tsv.gz")
	_, _ = executeCommand(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", in,
		"--out-tsv", outTSV,
	)

	rows := readTSVRows(t, outTSV)
	if len(rows) != 1 {
		t.Fatalf("expected 1 TSV row, got %d", len(rows))
	}
	assertRowValue(t, rows[0], "fwd", "8")
	assertRowValue(t, rows[0], "rev_rc", "8")
}

func TestPairedEndTSV(t *testing.T) {
	dir := t.TempDir()
	outTSV := filepath.Join(dir, "paired.tsv")

	stdout, _ := executeCommand(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", filepath.Join("testdata", "paired_r1.fq"),
		"--input2", filepath.Join("testdata", "paired_r2.fq"),
		"--out-tsv", outTSV,
	)

	assertSummaryLineTokens(t, stdout, []string{"Read", "Fwd", "FwdRC", "Rev", "RevRC"})
	assertSummaryLineTokens(t, stdout, []string{"R1", "2", "2", "2", "2"})
	assertSummaryLineTokens(t, stdout, []string{"R2", "1", "1", "1", "1"})

	rows := readTSVRows(t, outTSV)
	if len(rows) != 2 {
		t.Fatalf("expected 2 TSV rows, got %d", len(rows))
	}
	assertRowValue(t, rows[0], "read", "R1")
	assertRowValue(t, rows[0], "file", filepath.Join("testdata", "paired_r1.fq"))
	assertRowValue(t, rows[0], "reads", "4")
	assertRowValue(t, rows[0], "fwd", "2")
	assertRowValue(t, rows[0], "fwd_rc", "2")
	assertRowValue(t, rows[0], "rev", "2")
	assertRowValue(t, rows[0], "rev_rc", "2")
	assertRowValue(t, rows[1], "read", "R2")
	assertRowValue(t, rows[1], "file", filepath.Join("testdata", "paired_r2.fq"))
	assertRowValue(t, rows[1], "reads", "4")
	assertRowValue(t, rows[1], "fwd", "1")
	assertRowValue(t, rows[1], "fwd_rc", "1")
	assertRowValue(t, rows[1], "rev", "1")
	assertRowValue(t, rows[1], "rev_rc", "1")
}

func TestHeadSingleAndPaired(t *testing.T) {
	dir := t.TempDir()
	singleTSV := filepath.Join(dir, "head_single.tsv")
	_, _ = executeCommand(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", filepath.Join("testdata", "single_states.fa"),
		"--head", "3",
		"--out-tsv", singleTSV,
	)
	rows := readTSVRows(t, singleTSV)
	if len(rows) != 1 {
		t.Fatalf("expected 1 TSV row, got %d", len(rows))
	}
	assertRowValue(t, rows[0], "reads", "3")

	pairedTSV := filepath.Join(dir, "head_paired.tsv")
	_, _ = executeCommand(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", filepath.Join("testdata", "paired_r1.fq"),
		"--input2", filepath.Join("testdata", "paired_r2.fq"),
		"--head", "2",
		"--out-tsv", pairedTSV,
	)
	rows = readTSVRows(t, pairedTSV)
	if len(rows) != 2 {
		t.Fatalf("expected 2 TSV rows, got %d", len(rows))
	}
	assertRowValue(t, rows[0], "reads", "2")
	assertRowValue(t, rows[0], "fwd", "1")
	assertRowValue(t, rows[0], "fwd_rc", "1")
	assertRowValue(t, rows[1], "reads", "2")
	assertRowValue(t, rows[1], "rev", "1")
	assertRowValue(t, rows[1], "rev_rc", "1")
}

func TestPairedUnevenRecordsFails(t *testing.T) {
	dir := t.TempDir()
	short := filepath.Join(dir, "short_r2.fq")
	data, err := os.ReadFile(filepath.Join("testdata", "paired_r2.fq"))
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	if len(lines) < 12 {
		t.Fatalf("fixture too short")
	}
	if err := os.WriteFile(short, []byte(strings.Join(lines[:12], "\n")+"\n"), 0o644); err != nil {
		t.Fatalf("write short fixture: %v", err)
	}

	_, stderr := executeCommandExpectError(t,
		"--forward", "AAAA",
		"--reverse", "CCCC",
		"--input", filepath.Join("testdata", "paired_r1.fq"),
		"--input2", short,
	)
	if !strings.Contains(stderr, "different record counts") {
		t.Fatalf("unexpected stderr: %s", stderr)
	}
}

func BenchmarkPipelineSingle(b *testing.B) {
	for i := 0; i < b.N; i++ {
		if _, _, err := executeCommandRaw(
			"--forward", "AAAA",
			"--reverse", "CCCC",
			"--input", filepath.Join("testdata", "single_states.fa"),
			"--head", "16",
		); err != nil {
			b.Fatalf("execute failed: %v", err)
		}
	}
}

func executeCommand(t *testing.T, args ...string) (string, string) {
	t.Helper()
	stdout, stderr, err := executeCommandRaw(args...)
	if err != nil {
		t.Fatalf("execute failed: %v, stderr: %s", err, stderr)
	}
	return stdout, stderr
}

func executeCommandRaw(args ...string) (string, string, error) {
	cmd := rootCommand()
	var stdout bytes.Buffer
	var stderr bytes.Buffer
	cmd.SetOut(&stdout)
	cmd.SetErr(&stderr)
	cmd.SetArgs(args)
	err := cmd.Execute()
	return stdout.String(), stderr.String(), err
}

func executeCommandExpectError(t *testing.T, args ...string) (string, string) {
	t.Helper()
	stdout, _, err := executeCommandRaw(args...)
	if err == nil {
		t.Fatalf("expected command failure")
	}
	return stdout, err.Error()
}

func readTSVRows(t *testing.T, path string) []map[string]string {
	t.Helper()
	var data []byte
	var err error
	if strings.HasSuffix(path, ".gz") {
		file, err := os.Open(path)
		if err != nil {
			t.Fatalf("open gz tsv: %v", err)
		}
		defer file.Close()
		gzipReader, err := gzip.NewReader(file)
		if err != nil {
			t.Fatalf("gzip reader: %v", err)
		}
		defer gzipReader.Close()
		data, err = io.ReadAll(gzipReader)
		if err != nil {
			t.Fatalf("read gz tsv: %v", err)
		}
	} else {
		data, err = os.ReadFile(path)
		if err != nil {
			t.Fatalf("read tsv: %v", err)
		}
	}

	reader := csv.NewReader(strings.NewReader(string(data)))
	reader.Comma = '\t'
	records, err := reader.ReadAll()
	if err != nil {
		t.Fatalf("parse tsv: %v", err)
	}
	if len(records) < 2 {
		t.Fatalf("expected at least 2 TSV rows, got %d", len(records))
	}
	rows := make([]map[string]string, 0, len(records)-1)
	for _, record := range records[1:] {
		row := make(map[string]string, len(records[0]))
		for idx, header := range records[0] {
			row[header] = record[idx]
		}
		rows = append(rows, row)
	}
	return rows
}

func gzipFixture(t *testing.T, src string, dst string) {
	t.Helper()
	data, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	file, err := os.Create(dst)
	if err != nil {
		t.Fatalf("create gz fixture: %v", err)
	}
	defer file.Close()
	gzipWriter := gzip.NewWriter(file)
	if _, err := gzipWriter.Write(data); err != nil {
		t.Fatalf("write gz fixture: %v", err)
	}
	if err := gzipWriter.Close(); err != nil {
		t.Fatalf("close gz fixture: %v", err)
	}
}

func assertRowValue(t *testing.T, row map[string]string, key string, want string) {
	t.Helper()
	if got := row[key]; got != want {
		t.Fatalf("unexpected %s: got %s want %s", key, got, want)
	}
}

func assertSummaryLineTokens(t *testing.T, stdout string, want []string) {
	t.Helper()
	for _, line := range strings.Split(stdout, "\n") {
		if fields := strings.Fields(line); len(fields) == len(want) {
			match := true
			for i := range want {
				if fields[i] != want[i] {
					match = false
					break
				}
			}
			if match {
				return
			}
		}
	}
	t.Fatalf("stdout missing line %v:\n%s", want, stdout)
}
