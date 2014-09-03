package gossa

import (
	"encoding/csv"
	"log"
	"math"
	"math/rand"
	"strconv"
	"time"
)

func Ssa(leapMethod string, leapOption float64, rxnVectors *[][]int, rxnsK *[]float64, initPop *[]int, dur float64, outputDeltaT float64, outputFile *csv.Writer) (population []int, rxns []int) {
	//TODO: handle io.Writer as well as csv.Writer
	rng := rand.New(rand.NewSource(time.Now().UTC().UnixNano()))

	if leapMethod != "classic" {
		log.Fatalf("Leap method	%s not supported", leapMethod)
	}
	t := float64(0)
	tau := float64(0.0)
	lastOutputTime := float64(0)
	currentPop := make([]int, len(*initPop))
	rxns = make([]int, len(*rxnVectors))
	rxnsDelta := make([]int, len(*rxnVectors))
	copy(currentPop, *initPop)

	if outputFile != nil {
		writeOutput(outputFile, t, &currentPop)
	}

	//TODO: pseudorandom seed
	rand.Seed(5)
	for t <= dur {
		switch leapMethod {
		case "classic":
			tau, rxnsDelta = classicLeap(*rxnVectors, *rxnsK, currentPop, rng)
		}
		if tau+t >= lastOutputTime+outputDeltaT {
			lastOutputTime += outputDeltaT
			if lastOutputTime > dur {
				lastOutputTime = dur
			}
			writeOutput(outputFile, lastOutputTime, &currentPop)
		}
		t += tau
		if t <= dur {
			for rxnIndex, delta := range rxnsDelta {
				rxns[rxnIndex] += delta
			}
			currentPop = getPopulation(*initPop, *rxnVectors, rxns)
		}
	}
	return
}

func writeOutput(f *csv.Writer, t float64, pop *[]int) {
	vals := make([]string, len(*pop)+1)
	//TODO: automatically set t format
	vals[0] = strconv.FormatFloat(t, 'f', 5, 64)
	for i, v := range *pop {
		vals[i+1] = strconv.Itoa(v)
	}
	f.Write(vals)
}

func getPopulation(initPop []int, rxnVectors [][]int, rxns []int) []int {
	pop := make([]int, len(initPop))
	copy(pop, initPop)

	for i := range rxnVectors {
		if rxns[i] > 0 {
			for j := range pop {
				pop[j] += rxns[i] * rxnVectors[i][j]
			}
		}
	}
	for i, v := range pop {
		if v < 0 {
			pop[i] = 0
		}
	}
	return pop
}

type randomSource interface {
	Float64() float64
}

func classicLeap(rxnVectors [][]int, rxnsK []float64, pop []int, rng randomSource) (tau float64, rxnsDelta []int) {
	rxnsDelta = make([]int, len(rxnsK))
	propensities := make([]float64, len(rxnsK))
	for i := range propensities {
		for j, v := range rxnVectors[i] {
			if v == -1 {
				propensities[i] += rxnsK[i] * float64(pop[j])
			} else if v == -2 && pop[j] > 1 {
				propensities[i] = rxnsK[i] * float64(pop[j]) * float64(pop[j]-1)
			}
		}
	}
	propSum := 0.0
	for _, v := range propensities {
		propSum += v
	}

	// Draw two random numbers
	r1 := rng.Float64()
	r2 := rng.Float64()

	// Determine tau
	tau = 1.0 / propSum * math.Log(1.0/r1)

	// Determine reaction
	rxn := len(propensities)
	var propCum float64
	for i, v := range propensities {
		propCum += v
		if r2*propSum <= propCum {
			rxn = i
			break
		}
	}
	rxnsDelta = rxnVectors[rxn]
	return
}
