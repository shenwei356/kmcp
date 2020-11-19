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
	"fmt"

	"github.com/shenwei356/unikmer"
)

const extDataFile = ".unik"

func checkCompatibility(reader0 *unikmer.Reader, reader *unikmer.Reader, file string) {
	if reader0.K != reader.K {
		checkError(fmt.Errorf(`k-mer length not consistent (%d != %d), please check with "unikmer stats": %s`, reader0.K, reader.K, file))
	}
	if reader0.IsCanonical() != reader.IsCanonical() {
		checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer stats": %s`, file))
	}
	if reader0.IsHashed() != reader.IsHashed() {
		checkError(fmt.Errorf(`'hashed' flags not consistent, please check with "unikmer stats": %s`, file))
	}
	if reader0.IsScaled() != reader.IsScaled() {
		checkError(fmt.Errorf(`'scaled' flags not consistent, please check with "unikmer stats": %s`, file))
	}
}
