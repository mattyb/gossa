package gossa

import (
	. "math"
	"testing"
)

var populationTests = []struct {
	initPop     []int
	rxnVectors  [][]int
	rxns        []int
	expectedPop []int
}{
	{
		[]int{5, 0},
		[][]int{{-1, 1}},
		[]int{5},
		[]int{0, 5},
	},
	{
		[]int{3, 5, 0},
		[][]int{{-1, 1, 0},
			{0, -1, 1}},
		[]int{1, 3},
		[]int{2, 3, 3},
	},
}

func TestGetPopulation(t *testing.T) {
	for i, test := range populationTests {
		pop := getPopulation(test.initPop, test.rxnVectors, test.rxns)
		for j, p := range pop {
			if p != test.expectedPop[j] {
				t.Errorf("Wrong output for getPopulation test %d", i)
			}
		}
	}
}

var leapTests = []struct {
	name        string
	pop         []int
	rxnVectors  [][]int
	rxnsK       []float64
	r1          float64
	r2          float64
	expectedTau float64
	expectedRxn int
}{
	{
		"Basic Test",
		[]int{5, 0},
		[][]int{{-1, 1}},
		[]float64{.2},
		1.0 / E,
		0.49,
		1,
		0,
	},
	{
		"Pick first if r2 low",
		[]int{5, 5},
		[][]int{{-1, 1},
			{1, -1}},
		[]float64{.2, .2},
		1.0 / E,
		0.25,
		0.5,
		0,
	},
	{
		"Pick second if r2 high",
		[]int{5, 5},
		[][]int{{-1, 1},
			{1, -1}},
		[]float64{.2, .2},
		1.0 / E,
		0.75,
		0.5,
		1,
	},
	{
		"Pick second if first impossible",
		[]int{0, 5},
		[][]int{{-1, 1},
			{1, -1}},
		[]float64{.2, .2},
		1.0 / E,
		0.25,
		1,
		1,
	},
}

func TestClassicLeap(t *testing.T) {
	for i, test := range leapTests {
		seq := &DeterministicSequence{}
		seq.Sequence = []float64{test.r1, test.r2}
		seq.Reset()

		tau, rxnsDelta := classicLeap(test.rxnVectors, test.rxnsK, test.pop, seq)

		//Should have the right time leap
		if tau != test.expectedTau {
			t.Errorf("Failed leap test %d. tau was %f, expected %f", i, tau, test.expectedTau)
		}

		//Should have the right reactions
		for j, v := range rxnsDelta {
			if (j == test.expectedRxn && v == 0) || (j != test.expectedRxn && v == 1) {
				t.Errorf("Failed leap test (%s). Rxn Delta %v does not equal that of expected reaction %d (%v).", test.name, rxnsDelta, test.expectedRxn, test.rxnVectors[test.expectedRxn])
			}
		}
	}
}

type DeterministicSequence struct {
	Sequence []float64
	index    int
}

func (s *DeterministicSequence) Float64() (val float64) {
	val = s.Sequence[s.index]
	s.index += 1
	if s.index >= len(s.Sequence) {
		s.index = 0
	}
	return
}

/* not sure how to make this function work
func (s DeterministicSequence) SetSequence(seq []float64) {
	//s.sequence = make([]float64, len(seq))
	//copy(s.sequence, seq)
	for _, v := range seq {
		s.sequence = append(s.sequence, v)
	}
}
*/

func (s *DeterministicSequence) Reset() {
	s.index = 0
}
