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

type Summary struct {
	Config           Config
	Paired           bool
	Forward          string
	ForwardRC        string
	Reverse          string
	ReverseRC        string
	HasIUPAC         bool
	VariantCounts    map[string]int
	ConcreteVariants int
	Read1Counts      [16]uint64
	Read2Counts      [16]uint64
	PairCounts       [256]uint64
	Read1Total       uint64
	Read2Total       uint64
	PairTotal        uint64
}

type ProgressUpdate struct {
	Processed uint64
}

type RunOptions struct {
	Progress func(ProgressUpdate)
}

type partialCounts struct {
	read1 [16]uint64
	read2 [16]uint64
	pairs [256]uint64
}

type recordBatch struct {
	read1 [][]byte
	read2 [][]byte
}

func Run(ctx context.Context, cfg Config) (Summary, error) {
	return RunWithOptions(ctx, cfg, RunOptions{})
}

func RunWithOptions(ctx context.Context, cfg Config, opts RunOptions) (Summary, error) {
	runtime.GOMAXPROCS(cfg.Threads)

	matcher, err := NewMatcher(cfg.Forward, cfg.Reverse, cfg.Mismatches)
	if err != nil {
		return Summary{}, err
	}

	summary := Summary{
		Config:           cfg,
		Paired:           cfg.Input2 != "",
		Forward:          matcher.Forward,
		ForwardRC:        matcher.ForwardRC,
		Reverse:          matcher.Reverse,
		ReverseRC:        matcher.ReverseRC,
		HasIUPAC:         matcher.HasIUPAC,
		VariantCounts:    matcher.VariantCounts,
		ConcreteVariants: matcher.ConcreteVariants,
	}

	parts, err := runPipeline(ctx, cfg, matcher, opts)
	if err != nil {
		return Summary{}, err
	}
	for _, part := range parts {
		mergeReadCounts(&summary.Read1Counts, &part.read1)
		mergeReadCounts(&summary.Read2Counts, &part.read2)
		mergePairCounts(&summary.PairCounts, &part.pairs)
	}

	summary.Read1Total = sum16(summary.Read1Counts)
	summary.Read2Total = sum16(summary.Read2Counts)
	summary.PairTotal = sum256(summary.PairCounts)
	return summary, nil
}

func runPipeline(ctx context.Context, cfg Config, matcher *Matcher, opts RunOptions) ([]partialCounts, error) {
	ctx, cancel := context.WithCancel(ctx)
	defer cancel()

	jobs := make(chan recordBatch, cfg.Threads*2)
	results := make(chan partialCounts, cfg.Threads*2)
	producerErr := make(chan error, 1)

	go func() {
		var err error
		if cfg.Input2 == "" {
			err = produceSingle(ctx, cfg, jobs, opts)
		} else {
			err = producePaired(ctx, cfg, jobs, opts)
		}
		producerErr <- err
		close(jobs)
	}()

	var workerWG sync.WaitGroup
	for i := 0; i < cfg.Threads; i++ {
		workerWG.Add(1)
		go func() {
			defer workerWG.Done()
			for batch := range jobs {
				select {
				case <-ctx.Done():
					return
				default:
				}
				part := partialCounts{}
				for idx, seq1 := range batch.read1 {
					state1 := matcher.Match(seq1)
					part.read1[state1]++
					if len(batch.read2) == 0 {
						continue
					}
					state2 := matcher.Match(batch.read2[idx])
					part.read2[state2]++
					part.pairs[(int(state1)<<4)|int(state2)]++
				}
				select {
				case results <- part:
				case <-ctx.Done():
					return
				}
			}
		}()
	}

	go func() {
		workerWG.Wait()
		close(results)
	}()

	parts := make([]partialCounts, 0, 16)
	for part := range results {
		parts = append(parts, part)
	}

	if err := <-producerErr; err != nil {
		cancel()
		return nil, err
	}
	return parts, nil
}

func produceSingle(ctx context.Context, cfg Config, jobs chan<- recordBatch, opts RunOptions) error {
	reader, err := fastx.NewReader(nil, cfg.Input1, "")
	if err != nil {
		return err
	}
	defer reader.Close()

	var processed uint64
	batch := recordBatch{read1: make([][]byte, 0, batchSize)}
	checkedAlphabet := false

	for {
		if cfg.Head > 0 && processed >= uint64(cfg.Head) {
			break
		}
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return err
		}
		if !checkedAlphabet {
			if err := rejectProtein(reader.Alphabet(), cfg.Input1); err != nil {
				return err
			}
			checkedAlphabet = true
		}
		batch.read1 = append(batch.read1, append([]byte(nil), record.Seq.Seq...))
		processed++
		if len(batch.read1) == batchSize {
			if err := sendBatch(ctx, jobs, batch); err != nil {
				return err
			}
			emitProgress(opts, processed)
			batch = recordBatch{read1: make([][]byte, 0, batchSize)}
		}
	}
	if len(batch.read1) > 0 {
		if err := sendBatch(ctx, jobs, batch); err != nil {
			return err
		}
		emitProgress(opts, processed)
	}
	return nil
}

func producePaired(ctx context.Context, cfg Config, jobs chan<- recordBatch, opts RunOptions) error {
	reader1, err := fastx.NewReader(nil, cfg.Input1, "")
	if err != nil {
		return err
	}
	defer reader1.Close()

	reader2, err := fastx.NewReader(nil, cfg.Input2, "")
	if err != nil {
		return err
	}
	defer reader2.Close()

	var processed uint64
	batch := recordBatch{
		read1: make([][]byte, 0, batchSize),
		read2: make([][]byte, 0, batchSize),
	}
	checkedAlpha1 := false
	checkedAlpha2 := false

	for {
		if cfg.Head > 0 && processed >= uint64(cfg.Head) {
			break
		}

		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()

		if err1 == io.EOF || err2 == io.EOF {
			if err1 == io.EOF && err2 == io.EOF {
				break
			}
			return fmt.Errorf("paired inputs have different record counts")
		}
		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}

		if !checkedAlpha1 {
			if err := rejectProtein(reader1.Alphabet(), cfg.Input1); err != nil {
				return err
			}
			checkedAlpha1 = true
		}
		if !checkedAlpha2 {
			if err := rejectProtein(reader2.Alphabet(), cfg.Input2); err != nil {
				return err
			}
			checkedAlpha2 = true
		}

		batch.read1 = append(batch.read1, append([]byte(nil), record1.Seq.Seq...))
		batch.read2 = append(batch.read2, append([]byte(nil), record2.Seq.Seq...))
		processed++
		if len(batch.read1) == batchSize {
			if err := sendBatch(ctx, jobs, batch); err != nil {
				return err
			}
			emitProgress(opts, processed)
			batch = recordBatch{
				read1: make([][]byte, 0, batchSize),
				read2: make([][]byte, 0, batchSize),
			}
		}
	}

	if len(batch.read1) > 0 {
		if err := sendBatch(ctx, jobs, batch); err != nil {
			return err
		}
		emitProgress(opts, processed)
	}
	return nil
}

func emitProgress(opts RunOptions, processed uint64) {
	if opts.Progress == nil {
		return
	}
	opts.Progress(ProgressUpdate{Processed: processed})
}

func sendBatch(ctx context.Context, jobs chan<- recordBatch, batch recordBatch) error {
	select {
	case jobs <- batch:
		return nil
	case <-ctx.Done():
		return ctx.Err()
	}
}

func rejectProtein(alphabet *seq.Alphabet, path string) error {
	if alphabet == seq.Protein {
		return fmt.Errorf("protein input is not supported: %s", path)
	}
	return nil
}

func mergeReadCounts(dst *[16]uint64, src *[16]uint64) {
	for i := range dst {
		dst[i] += src[i]
	}
}

func mergePairCounts(dst *[256]uint64, src *[256]uint64) {
	for i := range dst {
		dst[i] += src[i]
	}
}

func sum16(values [16]uint64) uint64 {
	var total uint64
	for _, value := range values {
		total += value
	}
	return total
}

func sum256(values [256]uint64) uint64 {
	var total uint64
	for _, value := range values {
		total += value
	}
	return total
}
