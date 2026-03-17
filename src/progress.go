package prm

import (
	"fmt"
	"io"
	"strconv"
	"sync"
	"time"
)

var spinnerFrames = []string{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"}

type ProgressSpinner struct {
	writer   io.Writer
	unit     string
	delay    time.Duration
	interval time.Duration

	mu        sync.Mutex
	started   bool
	shown     bool
	success   bool
	frame     int
	processed uint64

	done      chan struct{}
	finished  chan struct{}
	startOnce sync.Once
	stopOnce  sync.Once
}

func NewProgressSpinner(w io.Writer, unit string) *ProgressSpinner {
	return newProgressSpinner(w, unit, 200*time.Millisecond, 100*time.Millisecond)
}

func newProgressSpinner(w io.Writer, unit string, delay time.Duration, interval time.Duration) *ProgressSpinner {
	return &ProgressSpinner{
		writer:   w,
		unit:     unit,
		delay:    delay,
		interval: interval,
		done:     make(chan struct{}),
		finished: make(chan struct{}),
	}
}

func (s *ProgressSpinner) Start() {
	if s == nil {
		return
	}
	s.startOnce.Do(func() {
		s.mu.Lock()
		s.started = true
		s.mu.Unlock()
		go s.run()
	})
}

func (s *ProgressSpinner) Update(update ProgressUpdate) {
	if s == nil {
		return
	}
	s.mu.Lock()
	s.processed = update.Processed
	s.mu.Unlock()
}

func (s *ProgressSpinner) Stop(success bool) {
	if s == nil {
		return
	}
	s.mu.Lock()
	started := s.started
	s.success = success
	s.mu.Unlock()
	if !started {
		return
	}
	s.stopOnce.Do(func() {
		close(s.done)
		<-s.finished
	})
}

func (s *ProgressSpinner) run() {
	defer close(s.finished)

	delayTimer := time.NewTimer(s.delay)
	defer delayTimer.Stop()

	var ticker *time.Ticker
	for {
		if ticker == nil {
			select {
			case <-s.done:
				return
			case <-delayTimer.C:
				s.mu.Lock()
				s.shown = true
				line := s.processingLineLocked()
				s.mu.Unlock()
				s.write("\r\x1b[2K" + line)
				ticker = time.NewTicker(s.interval)
			}
			continue
		}

		select {
		case <-s.done:
			ticker.Stop()
			s.finish()
			return
		case <-ticker.C:
			s.mu.Lock()
			line := s.processingLineLocked()
			s.mu.Unlock()
			s.write("\r\x1b[2K" + line)
		}
	}
}

func (s *ProgressSpinner) finish() {
	s.mu.Lock()
	shown := s.shown
	success := s.success
	processed := s.processed
	unit := s.unit
	s.mu.Unlock()

	if !shown {
		return
	}
	if success {
		s.write("\r\x1b[2K" + fmt.Sprintf("Processed %s %s\n", formatGroupedUint(processed), unit))
		return
	}
	s.write("\r\x1b[2K")
}

func (s *ProgressSpinner) processingLineLocked() string {
	frame := spinnerFrames[s.frame%len(spinnerFrames)]
	s.frame++
	return fmt.Sprintf("%s Processing %s %s", frame, formatGroupedUint(s.processed), s.unit)
}

func (s *ProgressSpinner) write(text string) {
	_, _ = io.WriteString(s.writer, text)
}

func formatGroupedUint(value uint64) string {
	text := strconv.FormatUint(value, 10)
	if len(text) <= 3 {
		return text
	}
	out := make([]byte, 0, len(text)+(len(text)-1)/3)
	prefix := len(text) % 3
	if prefix == 0 {
		prefix = 3
	}
	out = append(out, text[:prefix]...)
	for i := prefix; i < len(text); i += 3 {
		out = append(out, ',')
		out = append(out, text[i:i+3]...)
	}
	return string(out)
}
