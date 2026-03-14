package prm

import (
	"bytes"
	"fmt"
	"sort"
	"strings"

	"github.com/shenwei356/bio/seq"
)

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
	if err := seq.DNAredundant.IsValid(b); err != nil {
		return nil, err
	}
	return b, nil
}

func reverseComplement(src []byte) ([]byte, error) {
	cp := append([]byte(nil), src...)
	s, err := seq.NewSeq(seq.DNAredundant, cp)
	if err != nil {
		return nil, err
	}
	out := append([]byte(nil), s.RevCom().Seq...)
	normalizeDNA(out)
	return out, nil
}

