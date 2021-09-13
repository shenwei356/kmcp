//go:build ignore
// +build ignore

package main

import (
	"fmt"

	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
)

func main() {
	TEXT("TransposeOneColumn", NOSPLIT|NOPTR, "func(buf *[64]byte, pointers *[64]*[]byte, pos int)")

	Comment("pointer of buf")
	buf := Mem{Base: Load(Param("buf"), GP64())}

	Comment("pointer of pointers")
	pointers := Mem{Base: Load(Param("pointers"), GP64())}
	// pointers := Load(Param("pointers"), GP64())

	Comment("position pos")
	pos := Load(Param("pos"), GP64())

	Comment("//////////////////////////////////////////////////////////////////////")

	tmp := GP64()
	b := GP8()

	for j := 0; j < 8; j++ {
		Comment(fmt.Sprintf("=======================[group %d/8]=======================", j+1))
		Comment("reset the m64")
		q := GP64()
		XORQ(q, q)

		for i := 0; i < 8; i++ {
			Comment(fmt.Sprintf("------------------[byte %d/8]------------------", i+1))
			Comment("address of []byte, not just *[]byte")
			x := pointers.Offset(j*64 + i*8)
			MOVQ(x, tmp)

			Comment("address of slice data")
			MOVQ(Mem{Base: tmp}, tmp)

			Comment("read byte")
			ADDQ(pos, tmp)
			MOVB(Mem{Base: tmp}, b)

			Comment("add to m64 and ror")
			MOVBQZX(b, tmp)

			ADDQ(tmp, q)
			RORQ(Imm(8), q)

			// MOVB(b, buf.Offset(j*8+i))
		}

		Comment("-----------------------------------------")
		Comment("write 8 bytes")
		MOVQ(q, buf.Offset(j*8))
	}

	Comment("//////////////////////////////////////////////////////////////////////")

	// Store(tmp, ReturnIndex(0))
	RET()

	Generate()
}
