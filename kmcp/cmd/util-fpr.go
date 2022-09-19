// Copyright © 2020-2022 Wei Shen <shenwei356@gmail.com>
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
	"math/big"
)

// Theorem 2 in doi:10.1038/nbt.3442
// n: the number of query k-mers.
// k: the number of matched k-mers.
// fpr: the false positive of a single k-mer.
func QueryFPR(n int, k int, fpr float64) float64 {
	var r float64 = 1
	var coeff float64

	for i := 0; i <= k; i++ {
		coeff = BinomialCoeff(n, i)
		if math.IsInf(coeff, 1) {
			return 0
		}

		r -= coeff * math.Pow(fpr, float64(i)) * math.Pow(1-fpr, float64(n-i))

		if r < 0 {
			return 0
		}
	}

	return r
}

// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
// It's slow. You can call BinomialCoeffWithCache to create a same function with cache.
func BinomialCoeff(n int, k int) float64 {
	res := big.NewFloat(1)

	// Since C(n, k) = C(n, n-k)
	if k > n-k {
		k = n - k
	}

	// Calculate value of
	// [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	for i := 0; i < k; i++ {
		res = res.Mul(res, big.NewFloat(float64(n-i)))
		res = res.Quo(res, big.NewFloat(float64(i+1)))
	}

	coeff, _ := res.Float64()
	return coeff
}

// BinomialCoeffWithCache returns a BinomialCoeff function with cache.
// When n < bufSize, it returns cached values.
func BinomialCoeffWithCache(bufSize int) func(n int, k int) float64 {
	if bufSize < 1 {
		panic("buffer size should be > 0")
	}

	w := bufSize + 1
	h := (bufSize+1)/2 + 1 // a half is enough, because C(n, k) = C(n, n-k)
	buf := make([]float64, w*h)

	return func(n int, k int) float64 {
		if n > bufSize {
			return BinomialCoeff(n, k)
		}

		var idx int
		if k > n-k {
			idx = n*h + n - k
		} else {
			idx = n*h + k
		}

		coef := buf[idx]
		if coef == 0 {
			coef = BinomialCoeff(n, k)
			buf[idx] = coef
		}

		return coef
	}
}

// QueryFPRWithCache returns a QueryFPR function with cache.
// Since fpr, being a float, is difficult to use a slice to cache the values.
// But it could be achieved by splitting the maxFPR into bins.
// When n < bufSize and fpr <= maxFPR, it returns cached values.
func QueryFPRWithCache(bufSize int, maxFPR float64, bins int) func(n int, k int, fpr float64) float64 {
	if bufSize < 1 {
		panic("buffer size should be > 0")
	}
	if maxFPR <= 0 || maxFPR >= 1 {
		panic("maxFPR should be in range of (0, 1)")
	}
	if bins < 10 {
		panic("bins should be >= 10")
	}

	w := bufSize + 1
	h := (bufSize+1)/2 + 1 // a half is enough, because C(n, k) = C(n, n-k)
	wh := w * h
	// bins * w * h, using a slice is faster than a 3-dimension slice.
	buf := make([]float64, (bins+1)*wh)
	for i := 0; i < len(buf); i++ {
		buf[i] = -1 // fpr may be 0, so we initialize it to -1
	}

	binsFloat := float64(bins)

	bc := BinomialCoeffWithCache(bufSize)

	return func(n int, k int, fpr float64) float64 {
		if n > bufSize || fpr > maxFPR {
			return QueryFPR(n, k, fpr)
		}

		bin := int(fpr / maxFPR * binsFloat)

		// r := buf[bin][n][k]
		var idx int
		if k > n-k {
			idx = bin*wh + n*h + n - k
		} else {
			idx = bin*wh + n*h + k
		}

		r := buf[idx]

		if r > -1 { // cached
			// r0 := QueryFPR(n, k, fpr)
			// if r0 != r {
			//     fmt.Println(r0, r)
			// }
			return r
		}

		r = 1
		var coeff float64
		for i := 0; i <= k; i++ {
			coeff = bc(n, i)
			if math.IsInf(coeff, 1) {
				r = 0
				break
			}

			r -= coeff * math.Pow(fpr, float64(i)) * math.Pow(1-fpr, float64(n-i))

			if r < 0 {
				r = 0
				break
			}
		}

		buf[idx] = r
		return r
	}
}

/*
Solomon and Kingsford also apply a Chernoff bound and show that the false positive probability
for a query to be detected in a document is ≤ exp(−l(K − p)2 /(2(1 − p))).

p, fpr of single bloom filter.
k, theshold of query coverage.
l, number of k-mers.

import math
fpr = lambda p,k,l: math.exp(-l * (k - p) * (k - p) / (2 * (1 - p)))

fpr(0.3, 0.5, 70)
*/
func maxFPR(p float64, k float64, l int) float64 {
	return math.Exp(-float64(l) * (k - p) * (k - p) / (2 * (1 - p)))
}

func maxFPRf(p float64, k float64, l float64) float64 {
	return math.Exp(-l * (k - p) * (k - p) / (2 * (1 - p)))
}
