package jubjub

import (
	"math/big"
	"math/bits"
)

// FieldElement is an element of an arbitrary integer field.
type FieldElement struct {
	n          *big.Int
	fieldOrder *big.Int
}

// newFieldElement returns a new element of the curve's base field initialized to the value of `n`.
func (curve *Jubjub) newFieldElement(n *big.Int) *FieldElement {
	return newFieldElement(n, curve.fieldOrder)
}

// FeFromBytes reads a field element from little-endian bytes and returns it.
// If the value is larger than the size of the field, FeFromBytes will return a reduced value.
func (curve *Jubjub) FeFromBytes(in []byte) *FieldElement {
	fe := newFieldElement(nil, curve.fieldOrder)
	return fe.fromBytes(in)
}

// newFieldElement returns a new element of the field defined by `order` initialized to the value of `n`.
func newFieldElement(n, order *big.Int) *FieldElement {
	if n == nil {
		n = new(big.Int)
	}
	if n.Cmp(order) == 1 {
		// XXX timing
		n.Mod(n, order)
	}
	return &FieldElement{n, order}
}

// Set sets z to x and returns z.
func (z *FieldElement) Set(x *FieldElement) *FieldElement {
	z.n.Set(x.n)
	return z
}

// Equals compares two field elements and returns true if they are equal.
func (z *FieldElement) Equals(x *FieldElement) bool {
	return z.n.Cmp(x.n) == 0
}

// Exp sets z = x**y mod |m| (i.e. the sign of m is ignored), and returns z.
// If m == nil or m == 0, z = x**y unless y <= 0 then z = 1. If m > 0, y < 0,
// and x and n are not relatively prime, z is unchanged and nil is returned.
// m is always the order of the field.
func (z *FieldElement) Exp(x, y *FieldElement) *FieldElement {
	z.n.Exp(x.n, y.n, z.fieldOrder)
	return z
}

// Mul sets z to the product x*y, reducing the result by the field order, and returns z.
func (z *FieldElement) Mul(x, y *FieldElement) *FieldElement {
	z.n.Mul(x.n, y.n).Mod(z.n, z.fieldOrder)
	return z
}

// Sub sets z to the difference x-y, reducing the result by the field order, and returns z.
func (z *FieldElement) Sub(x, y *FieldElement) *FieldElement {
	z.n.Sub(x.n, y.n).Mod(z.n, z.fieldOrder)
	return z
}

// Add sets z to the sum x+y, reducing the result by the field order, and returns z.
func (z *FieldElement) Add(x, y *FieldElement) *FieldElement {
	z.n.Add(x.n, y.n).Mod(z.n, z.fieldOrder)
	return z
}

// Cmp compares x and y and returns
//     -1 if x < y
//      0 if x == y
//     +1 if x > y
func (z *FieldElement) Cmp(x *FieldElement) int {
	return z.n.Cmp(x.n)
}

// Neg sets z to -x and returns z.
func (z *FieldElement) Neg(x *FieldElement) *FieldElement {
	z.n.Neg(x.n).Mod(z.n, z.fieldOrder)
	return z
}

// ModInverse sets z to the multiplicative inverse of x in the field and returns z.
func (z *FieldElement) ModInverse(x *FieldElement) *FieldElement {
	z.n.ModInverse(x.n, z.fieldOrder)
	return z
}

// ModSqrt sets z to a square root of x in the field if such a square root exists,
// and returns z. If x is not a square in the field, ModSqrt leaves z unchanged
// and returns nil.
func (z *FieldElement) ModSqrt(x *FieldElement) *FieldElement {
	exists := z.n.ModSqrt(x.n, z.fieldOrder)
	if exists == nil {
		return nil
	}
	return z
}

// FromBytes decodes a little endian bytestring as a field element, sets z to that value, and returns z.
// It is unexported because it requires knowledge of the represented field and should be used from a concrete parameter set.
func (z *FieldElement) fromBytes(ser []byte) *FieldElement {
	in := ser[:]

	words := make([]big.Word, (z.fieldOrder.BitLen()+7)/bits.UintSize)
	for n := range words {
		for i := 0; i < bits.UintSize; i += 8 {
			if len(in) == 0 {
				break
			}
			words[n] |= big.Word(in[0]) << big.Word(i)
			in = in[1:]
		}
	}
	z.n.SetBits(words)
	z.n.Mod(z.n, z.fieldOrder)
	return z
}

// ToBytes converts z to a little-endian bytestring and returns the bytes.
func (z *FieldElement) ToBytes() []byte {
	z.n.Mod(z.n, z.fieldOrder)

	buf := make([]byte, 0, (z.fieldOrder.BitLen()+7)/8)
	for _, word := range z.n.Bits() {
		for i := 0; i < bits.UintSize; i += 8 {
			if len(buf) >= cap(buf) {
				break
			}
			buf = append(buf, byte(word))
			word >>= 8
		}
	}

	return buf[:(z.fieldOrder.BitLen()+7)/8]
}
