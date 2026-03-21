package prm

import "testing"

func TestNewMatcherNormalizesAndDeduplicates(t *testing.T) {
	m, err := NewMatcher("augr", "cccy", 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if m.Forward != "ATGR" {
		t.Fatalf("unexpected forward primer: %s", m.Forward)
	}
	if m.ForwardRC != "YCAT" {
		t.Fatalf("unexpected forward reverse complement: %s", m.ForwardRC)
	}
	if m.Reverse != "CCCY" {
		t.Fatalf("unexpected reverse primer: %s", m.Reverse)
	}
	if m.ReverseRC != "RGGG" {
		t.Fatalf("unexpected reverse reverse complement: %s", m.ReverseRC)
	}
	if m.VariantCounts["FWD"] != 2 {
		t.Fatalf("expected 2 FWD variants, got %d", m.VariantCounts["FWD"])
	}
	if !m.HasIUPAC {
		t.Fatal("expected matcher to record IUPAC usage")
	}
}

func TestNewMatcherSupportsInosine(t *testing.T) {
	m, err := NewMatcher("AIG", "CCI", 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if m.Forward != "AIG" {
		t.Fatalf("unexpected forward primer: %s", m.Forward)
	}
	if m.ForwardRC != "CIT" {
		t.Fatalf("unexpected forward reverse complement: %s", m.ForwardRC)
	}
	if m.Reverse != "CCI" {
		t.Fatalf("unexpected reverse primer: %s", m.Reverse)
	}
	if m.ReverseRC != "IGG" {
		t.Fatalf("unexpected reverse reverse complement: %s", m.ReverseRC)
	}
	if m.VariantCounts["FWD"] != 4 {
		t.Fatalf("expected 4 FWD variants, got %d", m.VariantCounts["FWD"])
	}
	if m.VariantCounts["REV"] != 4 {
		t.Fatalf("expected 4 REV variants, got %d", m.VariantCounts["REV"])
	}
}

func TestNormalizePrimerRejectsUnsupportedBase(t *testing.T) {
	_, err := normalizePrimer("ABZ")
	if err == nil {
		t.Fatal("expected unsupported base error")
	}
	if got := err.Error(); got != `unsupported primer base 'Z' at position 3` {
		t.Fatalf("unexpected error: %v", err)
	}
}

func TestHasAmbiguousBases(t *testing.T) {
	if hasAmbiguousBases([]byte("ACGT")) {
		t.Fatal("did not expect exact DNA to be ambiguous")
	}
	if !hasAmbiguousBases([]byte("ATGN")) {
		t.Fatal("expected N to count as ambiguous")
	}
	if !hasAmbiguousBases([]byte("ATGI")) {
		t.Fatal("expected I to count as ambiguous")
	}
	if hasAmbiguousBases([]byte("ATGT")) {
		t.Fatal("did not expect normalized U->T sequence to be ambiguous")
	}
}

func TestMatcherPreservesSemanticBitsForPalindromes(t *testing.T) {
	m, err := NewMatcher("ATAT", "CCCC", 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	state := m.Match([]byte("GGGATATCCC"))
	if state != (bitFwd | bitFwdRC) {
		t.Fatalf("unexpected state: %s", StateLabel(state))
	}
}

func TestMatcherExactAndMismatch(t *testing.T) {
	m, err := NewMatcher("AAAA", "CCCC", 1)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if state := m.Match([]byte("CAAATGCA")); state != bitFwd {
		t.Fatalf("expected FWD hit with one mismatch, got %s", StateLabel(state))
	}
	if state := m.Match([]byte("AANNGTCA")); state != 0 {
		t.Fatalf("expected no hit with Ns counted as mismatches, got %s", StateLabel(state))
	}
}

func TestMatcherNoHit(t *testing.T) {
	m, err := NewMatcher("AAAA", "CCCC", 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if state := m.Match([]byte("ACGTACGT")); state != 0 {
		t.Fatalf("expected no hit, got %s", StateLabel(state))
	}
}

func TestExpandVariantsLimit(t *testing.T) {
	if _, err := expandVariants([]byte("NNNNNNN"), 4096); err == nil {
		t.Fatal("expected expansion limit error")
	}
}

func BenchmarkMatcherExact(b *testing.B) {
	m, err := NewMatcher("AAAA", "CCCC", 0)
	if err != nil {
		b.Fatalf("unexpected error: %v", err)
	}
	seq := []byte("NNNNAAAANNNNGGGGNNNNTTTTNNNNCCCCNNNN")
	for i := 0; i < b.N; i++ {
		buf := append([]byte(nil), seq...)
		_ = m.Match(buf)
	}
}

func BenchmarkMatcherMismatch(b *testing.B) {
	m, err := NewMatcher("AAAA", "CCCC", 1)
	if err != nil {
		b.Fatalf("unexpected error: %v", err)
	}
	seq := []byte("NNNNAAATNNNNGGGGNNNNTTTTNNNNCCCTNNNN")
	for i := 0; i < b.N; i++ {
		buf := append([]byte(nil), seq...)
		_ = m.Match(buf)
	}
}
