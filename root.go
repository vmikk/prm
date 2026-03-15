package main

import (
	"errors"
	"fmt"
	"runtime"

	"prm/internal/prm"

	"github.com/spf13/cobra"
)

func newRootCommand() *cobra.Command {
	cfg := prm.Config{}

	cmd := &cobra.Command{
		Use:           "prm",
		Short:         "Count primer orientation combinations in FASTA/FASTQ reads",
		Args:          cobra.NoArgs,
		SilenceUsage:  true,
		SilenceErrors: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			if err := validateConfig(cfg); err != nil {
				return err
			}
			summary, err := prm.Run(cmd.Context(), cfg)
			if err != nil {
				return err
			}
			if err := prm.WriteHumanSummary(cmd.OutOrStdout(), summary); err != nil {
				return err
			}
			if cfg.OutTSV == "" {
				return nil
			}
			if err := prm.WriteTSV(cfg.OutTSV, summary); err != nil {
				return err
			}
			return nil
		},
	}

	cmd.Flags().StringVarP(&cfg.Forward, "forward", "f", "", "forward primer sequence")
	cmd.Flags().StringVarP(&cfg.Reverse, "reverse", "r", "", "reverse primer sequence")
	cmd.Flags().StringVarP(&cfg.Input1, "input", "i", "", `input FASTA/FASTQ file ("-" allowed only in single-end mode)`)
	cmd.Flags().StringVarP(&cfg.Input2, "input2", "I", "", "optional second FASTA/FASTQ file for paired-end reads")
	cmd.Flags().IntVarP(&cfg.Mismatches, "mismatches", "m", 0, "allowed substitution mismatches")
	cmd.Flags().IntVarP(&cfg.Head, "head", "n", 0, "number of records to process (or read pairs in paired mode, 0 for all)")
	cmd.Flags().IntVarP(&cfg.Threads, "threads", "j", runtime.NumCPU(), "number of worker threads")
	cmd.Flags().StringVarP(&cfg.OutTSV, "out-tsv", "o", "", `optional TSV output path (suffix ".gz" writes gzipped output)`)

	_ = cmd.MarkFlagRequired("forward")
	_ = cmd.MarkFlagRequired("reverse")
	_ = cmd.MarkFlagRequired("input")

	return cmd
}

func validateConfig(cfg prm.Config) error {
	if cfg.Threads < 1 {
		return errors.New("--threads must be greater than 0")
	}
	if cfg.Mismatches < 0 {
		return errors.New("--mismatches must be greater than or equal to 0")
	}
	if cfg.Head < 0 {
		return errors.New("--head must be greater than or equal to 0")
	}
	if cfg.Input1 == "" {
		return errors.New("--input is required")
	}
	if cfg.Input2 != "" && cfg.Input1 == "-" {
		return errors.New(`--input "-" is not supported with --input2`)
	}
	if cfg.Input2 == "-" {
		return errors.New(`--input2 "-" is not supported`)
	}
	if cfg.Forward == "" || cfg.Reverse == "" {
		return errors.New("--forward and --reverse are required")
	}
	if cfg.Mismatches > 0 {
		if len(cfg.Forward) < cfg.Mismatches || len(cfg.Reverse) < cfg.Mismatches {
			return fmt.Errorf("--mismatches (%d) must not exceed primer lengths", cfg.Mismatches)
		}
	}
	return nil
}
