package prm

import (
	"context"
	"fmt"
	"io"
	"runtime"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

const batchSize = 512

type partialCounts struct {
	read1 [16]uint64
	read2 [16]uint64
	pairs [256]uint64
}

type recordBatch struct {
	read1 [][]byte
	read2 [][]byte
}
