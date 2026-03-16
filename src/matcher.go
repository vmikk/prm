package prm

import (
	"bytes"
	"fmt"
	"sort"
	"strings"
)

// Orientation bits represent which primer variant matched
// A match state is a uint8 combining these bits: FWD|FWD_RC|REV|REV_RC
// Example: state 0b1010 (s1010) means FWD and REV primers were found
const (
	bitFwd   uint8 = 1 << 3 // forward primer (original)
	bitFwdRC uint8 = 1 << 2 // forward primer reverse complement
	bitRev   uint8 = 1 << 1 // reverse primer (original)
	bitRevRC uint8 = 1 << 0 // reverse primer reverse complement

	maxExpandedVariants = 16384
)

var orientationBits = []struct {
	Name string
	Bit  uint8
}{
	{Name: "FWD", Bit: bitFwd},
	{Name: "FWD_RC", Bit: bitFwdRC},
	{Name: "REV", Bit: bitRev},
	{Name: "REV_RC", Bit: bitRevRC},
}

var dnaExpansion = map[byte]string{
	'A': "A",
	'C': "C",
	'G': "G",
	'T': "T",
	'U': "T",
	'R': "AG",
	'Y': "CT",
	'M': "AC",
	'K': "GT",
	'S': "CG",
	'W': "AT",
	'H': "ACT",
	'B': "CGT",
	'V': "ACG",
	'D': "AGT",
	'N': "ACGT",
	'I': "ACGT",
}

var dnaComplement = map[byte]byte{
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A',
	'R': 'Y',
	'Y': 'R',
	'M': 'K',
	'K': 'M',
	'S': 'S',
	'W': 'W',
	'H': 'D',
	'B': 'V',
	'V': 'B',
	'D': 'H',
	'N': 'N',
	'I': 'I',
}

type Matcher struct {
	Forward          string
	ForwardRC        string
	Reverse          string
	ReverseRC        string
	Mismatches       int
	HasIUPAC         bool
	VariantCounts    map[string]int
	ConcreteVariants int
	variants         []variant
}

type variant struct {
	pattern         []byte
	orientationBits uint8
}

// Convert DNA sequence to uppercase and change U to T
func normalizeDNA(seq []byte) {
	for i, b := range seq {
		switch {
		case b >= 'a' && b <= 'z':
			b -= 'a' - 'A'
		}
		if b == 'U' {
			b = 'T'
		}
		seq[i] = b
	}
}

func normalizePrimer(value string) ([]byte, error) {
	value = strings.TrimSpace(value)
	if value == "" {
		return nil, fmt.Errorf("blank primer is not allowed")
	}
	b := []byte(value)
	normalizeDNA(b)
	if err := validatePrimerBases(b); err != nil {
		return nil, err
	}
	return b, nil
}

func validatePrimerBases(seq []byte) error {
	for i, base := range seq {
		if _, ok := dnaExpansion[base]; ok {
			continue
		}
		return fmt.Errorf("unsupported primer base %q at position %d", base, i+1)
	}
	return nil
}

func reverseComplement(src []byte) ([]byte, error) {
	out := make([]byte, len(src))
	for i := range src {
		base := src[len(src)-1-i]
		comp, ok := dnaComplement[base]
		if !ok {
			return nil, fmt.Errorf("unsupported primer base %q", base)
		}
		out[i] = comp
	}
	return out, nil
}

func hasAmbiguousBases(seq []byte) bool {
	for _, base := range seq {
		switch base {
		case 'A', 'C', 'G', 'T':
			continue
		default:
			return true
		}
	}
	return false
}

func expandVariants(primer []byte, limit int) ([][]byte, error) {
	variants := [][]byte{{}}
	for _, base := range primer {
		letters, ok := dnaExpansion[base]
		if !ok {
			return nil, fmt.Errorf("unsupported base %q", base)
		}

		next := make([][]byte, 0, len(variants)*len(letters))
		seen := make(map[string]struct{}, len(variants)*len(letters))
		for _, prefix := range variants {
			for i := 0; i < len(letters); i++ {
				item := append(append([]byte(nil), prefix...), letters[i])
				key := string(item)
				if _, exists := seen[key]; exists {
					continue
				}
				seen[key] = struct{}{}
				next = append(next, item)
				if len(next) > limit {
					return nil, fmt.Errorf("expanded to more than %d concrete variants", limit)
				}
			}
		}
		variants = next
	}
	sort.Slice(variants, func(i, j int) bool {
		return bytes.Compare(variants[i], variants[j]) < 0
	})
	return variants, nil
}

func (m *Matcher) Match(seq []byte) uint8 {
	normalizeDNA(seq)
	var state uint8
	for _, pattern := range m.variants {
		// Skip variants for orientations already found
		if state&pattern.orientationBits == pattern.orientationBits {
			continue
		}
		if m.Mismatches == 0 {
			if bytes.Index(seq, pattern.pattern) >= 0 {
				state |= pattern.orientationBits
			}
		} else if containsWithMismatches(seq, pattern.pattern, m.Mismatches) {
			state |= pattern.orientationBits
		}
		// 0x0F = all 4 orientation bits set (FWD|FWD_RC|REV|REV_RC)
		if state == 0x0F {
			break
		}
	}
	return state
}

func containsWithMismatches(seq []byte, pattern []byte, mismatches int) bool {
	if len(pattern) == 0 {
		return true
	}
	if len(seq) < len(pattern) {
		return false
	}
	last := len(seq) - len(pattern)
	for start := 0; start <= last; start++ {
		mismatchCount := 0
		for i := range pattern {
			if seq[start+i] == pattern[i] {
				continue
			}
			mismatchCount++
			if mismatchCount > mismatches {
				break
			}
		}
		if mismatchCount <= mismatches {
			return true
		}
	}
	return false
}

// StateLabel formats a match state as sXXXX (e.g., s1000, s0011).
// Each bit indicates if the corresponding orientation was found in the sequence.
func StateLabel(state uint8) string {
	return fmt.Sprintf("s%04b", state)
}

func AllStateLabels() []string {
	labels := make([]string, 16)
	for i := range 16 {
		labels[i] = StateLabel(uint8(i))
	}
	return labels
}

// NewMatcher creates a Matcher with all primer variants expanded for matching.
// It normalizes primers, computes reverse complements, and expands IUPAC codes.
func NewMatcher(forward, reverse string, mismatches int) (*Matcher, error) {
	fwd, err := normalizePrimer(forward)
	if err != nil {
		return nil, fmt.Errorf("forward primer: %w", err)
	}
	rev, err := normalizePrimer(reverse)
	if err != nil {
		return nil, fmt.Errorf("reverse primer: %w", err)
	}
	if mismatches > len(fwd) {
		return nil, fmt.Errorf("forward primer length %d is shorter than mismatches %d", len(fwd), mismatches)
	}
	if mismatches > len(rev) {
		return nil, fmt.Errorf("reverse primer length %d is shorter than mismatches %d", len(rev), mismatches)
	}

	fwdRC, err := reverseComplement(fwd)
	if err != nil {
		return nil, fmt.Errorf("forward primer reverse complement: %w", err)
	}
	revRC, err := reverseComplement(rev)
	if err != nil {
		return nil, fmt.Errorf("reverse primer reverse complement: %w", err)
	}

	orientationSeqs := []struct {
		name string
		bit  uint8
		seq  []byte
	}{
		{name: "FWD", bit: bitFwd, seq: fwd},
		{name: "FWD_RC", bit: bitFwdRC, seq: fwdRC},
		{name: "REV", bit: bitRev, seq: rev},
		{name: "REV_RC", bit: bitRevRC, seq: revRC},
	}

	merged := make(map[string]uint8, 16)
	variantCounts := make(map[string]int, len(orientationSeqs))
	for _, item := range orientationSeqs {
		variants, err := expandVariants(item.seq, maxExpandedVariants)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", item.name, err)
		}
		variantCounts[item.name] = len(variants)
		for _, v := range variants {
			key := string(v)
			merged[key] |= item.bit
		}
	}

	patterns := make([]string, 0, len(merged))
	for pattern := range merged {
		patterns = append(patterns, pattern)
	}
	sort.Strings(patterns)

	variants := make([]variant, 0, len(patterns))
	for _, pattern := range patterns {
		variants = append(variants, variant{
			pattern:         []byte(pattern),
			orientationBits: merged[pattern],
		})
	}

	return &Matcher{
		Forward:          string(fwd),
		ForwardRC:        string(fwdRC),
		Reverse:          string(rev),
		ReverseRC:        string(revRC),
		Mismatches:       mismatches,
		HasIUPAC:         hasAmbiguousBases(fwd) || hasAmbiguousBases(rev),
		VariantCounts:    variantCounts,
		ConcreteVariants: len(variants),
		variants:         variants,
	}, nil
}
