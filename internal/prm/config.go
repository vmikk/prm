package prm

type Config struct {
	Forward    string
	Reverse    string
	Input1     string
	Input2     string
	Mismatches int
	Head       int
	Threads    int
	OutTSV     string
}
