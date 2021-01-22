package jubjub

import (
	"math/big"
	"math/bits"

	"github.com/pkg/errors"
)

var (
	ErrScalarOutOfRange = errors.New("scalar was not in the correct range")
)

// Scalar is an element of an internally-specified subgroup of specific curve.
type Scalar struct {
	n          *big.Int
	fieldOrder *big.Int
}

// newScalar returns n mod order. It additionally returns an out-of-range error
// if the value needed to be reduced.
func newScalar(n, order *big.Int) (*Scalar, error) {
	if n == nil {
		n = new(big.Int)
	}

	if n.Cmp(order) == 1 || n.Cmp(big.NewInt(0)) == -1 {
		n.Mod(n, order)
		return &Scalar{n, order}, ErrScalarOutOfRange
	}

	return &Scalar{n, order}, nil
}

// ScalarFromBytes reads a scalar value from little-endian bytes and returns
// it. If the value of the Int is outside the order of the subgroup, ScalarFromBytes
// reduces it.
func (curve *Jubjub) ScalarFromBytes(in []byte) (*Scalar, error) {
	sc, _ := newScalar(nil, curve.subgroupOrder)
	return sc.fromBytes(in)
}

// ScalarFromBig converts a big.Int into a Scalar value in the correct range.
// If the value of the Int is outside the order of the subgroup, ScalarFromBig
// additionally returns an error indicating this was the case.
func (curve *Jubjub) ScalarFromBig(n *big.Int) (*Scalar, error) {
	return newScalar(n, curve.subgroupOrder)
}

// FromBytes sets the scalar to the value of a little-endian bytestring and returns the scalar.
// If the value would be outside the order of the subgroup, FromBytes reduces it.
// FromBytes additionally returns an error if the value needed to be reduced.
// It is unexported because it requires knowledge of the represented field and should be used from a concrete parameter set.
func (sc *Scalar) fromBytes(in []byte) (*Scalar, error) {
	words := make([]big.Word, (sc.fieldOrder.BitLen()+7)/bits.UintSize)
	for n := range words {
		for i := 0; i < bits.UintSize; i += 8 {
			if len(in) == 0 {
				break
			}
			words[n] |= big.Word(in[0]) << big.Word(i)
			in = in[1:]
		}
	}

	sc.n.SetBits(words)

	if sc.n.Cmp(sc.fieldOrder) == 1 || sc.n.Cmp(big.NewInt(0)) == -1 {
		sc.n.Mod(sc.n, sc.fieldOrder)
		return sc, ErrScalarOutOfRange
	}

	return sc, nil
}

// ToBytes reduces then converts the scalar to a little-endian bytestring.
func (sc Scalar) ToBytes() []byte {
	sc.n.Mod(sc.n, sc.fieldOrder)

	buf := make([]byte, 0, (sc.fieldOrder.BitLen()+7)/8)
	for _, word := range sc.n.Bits() {
		for i := 0; i < bits.UintSize; i += 8 {
			if len(buf) >= cap(buf) {
				break
			}
			buf = append(buf, byte(word))
			word >>= 8
		}
	}

	return buf[:(sc.fieldOrder.BitLen()+7)/8]
}
