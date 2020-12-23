// Copyright © 2020 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package cmd

import (
	"math"
)

// CalcSignatureSize is from https://github.com/bingmann/cobs/blob/master/cobs/util/calc_signature_size.cpp .
// but we roundup to 2^n.
func CalcSignatureSize(numElements uint64, numHashes int, falsePositiveRate float64) uint64 {
	ratio := float64(-numHashes) / (math.Log(1 - math.Pow(falsePositiveRate, 1/float64(numHashes))))
	return roundup64(uint64(math.Ceil(float64(numElements) * ratio)))
}

// get the two basic hash function values for data.
// Based on early version of https://github.com/willf/bloom/blob/master/bloom.go .
func baseHashes(hash uint64) (uint32, uint32) {
	return uint32(hash >> 32), uint32(hash)
}

// return locations in bitset for a hash
func hashLocations(hash uint64, numHashes int, numSigs uint64) []int {
	if numHashes < 1 {
		return nil
	}

	locs := make([]int, numHashes)
	if numHashes == 1 {
		locs[0] = int(hash % numSigs)
		return locs
	}

	a, b := baseHashes(hash)
	for i := uint32(0); i < uint32(numHashes); i++ {
		locs[i] = int(uint64(a+b*i) % numSigs)
	}
	return locs
}

// return locations in bitset for a hash, faster with AND operation
func hashLocationsFaster(hash uint64, numHashes int, numSigsM1 uint64) []int {
	if numHashes < 1 {
		return nil
	}

	locs := make([]int, numHashes)
	if numHashes == 1 {
		locs[0] = int(hash & numSigsM1)
		return locs
	}

	a, b := baseHashes(hash)
	for i := uint32(0); i < uint32(numHashes); i++ {
		locs[i] = int(uint64(a+b*i) & numSigsM1)
	}
	return locs
}

// return hashes for a hash
func hashValues(hash uint64, numHashes int) []uint64 {
	if numHashes < 1 {
		return nil
	}

	hashes := make([]uint64, numHashes)

	if numHashes == 1 {
		hashes[0] = hash
		return hashes
	}

	a, b := baseHashes(hash)
	for i := uint32(0); i < uint32(numHashes); i++ {
		hashes[i] = uint64(a + b*i)
	}
	return hashes
}

// https://gist.github.com/badboy/6267743 .
// version with mask: https://gist.github.com/lh3/974ced188be2f90422cc .
func hash64(key uint64) uint64 {
	key = (^key) + (key << 21) // key = (key << 21) - key - 1
	key = key ^ (key >> 24)
	key = (key + (key << 3)) + (key << 8) // key * 265
	key = key ^ (key >> 14)
	key = (key + (key << 2)) + (key << 4) // key * 21
	key = key ^ (key >> 28)
	key = key + (key << 31)
	return key
}
