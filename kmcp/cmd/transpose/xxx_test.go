package transpose

import (
	"bytes"
	"math/rand"
	"testing"
)

var data [64]*[]byte
var ncolumn = 8

func init() {
	for i := 0; i < 64; i++ {
		tmp := make([]byte, ncolumn)
		for j := 0; j < ncolumn; j++ {
			if rand.Float64() < 0.5 {
				tmp[j] = byte(rand.Uint32() & 0b11)
			}
		}
		data[i] = &tmp

		// fmt.Println(tmp)
	}

}

func TestTranspose(t *testing.T) {
	var buf [64]byte
	var buf0 [64]byte

	for j := 0; j < 1; j++ {
		for i := 0; i < 64; i++ {
			buf0[i] = (*data[i])[j]
		}

		// fmt.Println(buf0)
		TransposeOneColumn(&buf, &data, j)
		// fmt.Println(buf)

		if !bytes.Equal(buf[:], buf0[:]) {
			t.Error("error")
		}
	}
}
