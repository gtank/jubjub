// Package jubjub provides an implementation of the Jubjub elliptic curve used in Zcash.
package jubjub

import (
	"math/big"

	"github.com/pkg/errors"
)

var (
	ErrInvalidPoint error = errors.New("not a valid jubjub point")
	ErrIdentity           = errors.New("point was in the h-torsion")
)

// Jubjub provides a context for working with the Jubjub elliptic curve.
type Jubjub struct {
	fieldOrder    *big.Int
	subgroupOrder *big.Int
	generatorY    *FieldElement
	d             *FieldElement
	cofactor      *Scalar
	fieldZero     *FieldElement
	fieldOne      *FieldElement
	fieldTwo      *FieldElement
}

// Curve initializes a bunch of values needed for working with the Jubjub curve and returns a handle to that context.
func Curve() *Jubjub {
	fieldOrder, _ := new(big.Int).SetString("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10)
	subgroupOrder, _ := new(big.Int).SetString("6554484396890773809930967563523245729705921265872317281365359162392183254199", 10)

	// d = -10240/10241
	d := big.NewInt(10241)
	d.ModInverse(d, fieldOrder)
	d.Mul(d, big.NewInt(-10240))

	// Used in formulas and comparisons, nice to have cached
	feZero := newFieldElement(big.NewInt(0), fieldOrder)
	feOne := newFieldElement(big.NewInt(1), fieldOrder)

	h, _ := newScalar(big.NewInt(8), subgroupOrder)

	jubjub := &Jubjub{
		fieldOrder:    fieldOrder,
		subgroupOrder: subgroupOrder,
		generatorY:    newFieldElement(big.NewInt(11), fieldOrder),
		d:             newFieldElement(d, fieldOrder),
		cofactor:      h,
		fieldZero:     feZero,
		fieldOne:      feOne,
	}
	return jubjub
}

// Identity returns the curve's identity point
func (curve *Jubjub) Identity() *Point {
	return &Point{
		curve,
		curve.newFieldElement(big.NewInt(0)),
		curve.newFieldElement(big.NewInt(1)),
	}
}

// Generator returns a generator for the full 8*q group on Jubjub, the positive point with y-value 11.
func (curve *Jubjub) Generator() *Point {
	g, _ := curve.Decompress(curve.generatorY.ToBytes())
	return g
}

// SubgroupGenerator returns a generator for the prime-order subgroup of Jubjub.
func (curve *Jubjub) SubgroupGenerator() *Point {
	return curve.Generator().MulByCofactor()
}

// ScalarMult multiplies the point by the scalar and returns a newly allocated result point.
// It returns an error if the point is not on the curve.
func (curve *Jubjub) ScalarMult(scalar *Scalar, point *Point) (*Point, error) {
	if !point.IsOnCurve() {
		// TODO: is it worth having this check here instead of at callsites?
		return nil, ErrInvalidPoint
	}

	r0, r1 := curve.Identity(), point.Clone()
	for i := scalar.n.BitLen() - 1; i >= 0; i-- {
		if scalar.n.Bit(i) == 0 {
			r1.Add(r0, r1)
			r0.Double(r0)
		} else {
			r0.Add(r0, r1)
			r1.Double(r1)
		}
	}

	return r0, nil
}

// Decompress reads a compressed Edwards point and returns that point or an error if it is invalid.
func (curve *Jubjub) Decompress(compressed []byte) (*Point, error) {
	p := &Point{curve, curve.newFieldElement(nil), curve.newFieldElement(nil)}
	err := p.UnmarshalBinary(compressed)

	if err != nil {
		return nil, err
	}
	return p, nil
}

// Point is a point on Jubjub.
type Point struct {
	curve *Jubjub
	x, y  *FieldElement
}

// newPoint returns a newly allocated point for the given curve, set to the identity.
func newPoint(curve *Jubjub) *Point {
	return curve.Identity()
}

// Clone returns a newly allocated copy of p.
func (p *Point) Clone() *Point {
	newPoint := &Point{p.curve, p.curve.newFieldElement(nil), p.curve.newFieldElement(nil)}
	newPoint.x.Set(p.x)
	newPoint.y.Set(p.y)
	return newPoint
}

// Equals returns true if p == q and false if they are not.
func (p *Point) Equals(q *Point) bool {
	return p.x.Cmp(q.x) == 0 && p.y.Cmp(q.y) == 0
}

// IsOnCurve returns true if the point is on the curve and false if not.
func (p *Point) IsOnCurve() bool {
	// a*x^2+y^2 = 1 + d*x^2*y^2, a = -1
	//   => -x^2 + y^2 - 1 - d*x^2*y^2 = 0

	xx := newFieldElement(nil, p.curve.fieldOrder).Mul(p.x, p.x)
	yy := newFieldElement(nil, p.curve.fieldOrder).Mul(p.y, p.y)

	// 1 + d*x^2*y^2
	dxxyy := newFieldElement(nil, p.curve.fieldOrder).Mul(xx, yy)
	dxxyy.Mul(dxxyy, p.curve.d)
	dxxyy.Add(dxxyy, p.curve.fieldOne)

	// -1 - d*x^2*y^2
	dxxyy.Neg(dxxyy)

	// -x^2 + y^2 + (-1 - d*x^2*y^2) == 0
	result := newFieldElement(nil, p.curve.fieldOrder).Neg(xx)
	result.Add(result, yy)
	result.Add(result, dxxyy)

	return result.Cmp(p.curve.fieldZero) == 0
}

// IsIdentity returns true if the point is the identity point, and false if not.
func (p *Point) IsIdentity() bool {
	return p.x.Equals(p.curve.fieldZero) && p.y.Equals(p.curve.fieldOne)
}

// Compress returns a representation of the point in compressed Edwards y format,
// ignoring whether or not the point is valid. If you are not confident in the
// provenance of your point, use MarshalBinary directly to receive the error from the check.
func (p *Point) Compress() []byte {
	repr, _ := p.MarshalBinary()
	return repr
}

// MarshalBinary returns the point in "compressed Edwards y" format.
func (p *Point) MarshalBinary() ([]byte, error) {
	feX := p.x.ToBytes()
	feY := p.y.ToBytes()

	// TODO fixed length
	feY[31] |= (feX[0] & 1) << 7

	return feY, nil
}

// UnmarshalBinary reads a Jubjub point in compressed Edwards y format and attempts to decompress it.
func (p *Point) UnmarshalBinary(compressed []byte) error {
	// TODO fixed length
	// Recall that Jubjub is a slightly smaller curve, fits in 32
	if len(compressed) != 32 {
		return ErrInvalidPoint
	}

	fieldOne := p.curve.fieldOne
	fieldOrder := p.curve.fieldOrder

	// Extract & clear sign bit
	// TODO fixed length
	in := make([]byte, 32)
	copy(in, compressed)

	sign := in[31] >> 7
	in[31] &= 0x7F

	// We want to know sqrt((y^2 - 1) / (dy^2 + 1))

	y := newFieldElement(nil, fieldOrder).fromBytes(in)
	yy := newFieldElement(nil, fieldOrder).Mul(y, y)
	v := newFieldElement(nil, fieldOrder)
	u := newFieldElement(nil, fieldOrder).Sub(yy, fieldOne) // u = y^2 - 1
	v.Mul(yy, p.curve.d).Add(v, fieldOne)                   // v = d*y^2 + 1
	v.ModInverse(v)                                         // 1 / d*y^2 + 1
	u.Mul(u, v)                                             // y^2 - 1 / d*y^2 + 1

	// 5.4.8.3 Jubjub
	// When computing square roots in Fq in order to decompress a point encoding,
	// the implementation MUST NOT assume that the square root exists, or that
	// the encoding represents a point on the curve
	if u.ModSqrt(u) == nil {
		return ErrInvalidPoint
	}

	decompressed := u.ToBytes()[0] & 1
	if sign != decompressed {
		u.Neg(u)
	}

	p.x.Set(u)
	p.y.Set(y)

	if !p.IsOnCurve() {
		return ErrInvalidPoint
	}

	return nil
}

// Neg sets p to the negated form of q and returns p.
func (p *Point) Neg(q *Point) *Point {
	p.x.Neg(q.x)
	p.y.Set(q.y)
	return p
}

// MulByCofactor sets p to the value of h*p and returns p.
func (p *Point) MulByCofactor() *Point {
	res, _ := p.curve.ScalarMult(p.curve.cofactor, p)

	p.x.Set(res.x)
	p.y.Set(res.y)

	return p
}

// Add adds p1+p2 and returns a newly allocated result point.
func (curve *Jubjub) Add(p1 *Point, p2 *Point) *Point {
	// Affine addition formulas: (x1,y1) + (x2,y2) = (x3,y3) where
	//  x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
	//  y3 = (y1*y2 - a*x1*x2) / (1 - d*x1*x2*y1*y2)
	// Recall a = -1

	x1y2 := newFieldElement(nil, curve.fieldOrder).Mul(p1.x, p2.y)
	x2y1 := newFieldElement(nil, curve.fieldOrder).Mul(p2.x, p1.y)

	y1y2 := newFieldElement(nil, curve.fieldOrder).Mul(p1.y, p2.y)
	x1x2 := newFieldElement(nil, curve.fieldOrder).Mul(p1.x, p2.x)

	// d*x1*x2*y1*y2
	commonTerm := newFieldElement(nil, curve.fieldOrder).Mul(x1x2, y1y2)
	commonTerm.Mul(commonTerm, curve.d)

	// 1 / (1 + d*x1*x2*y1*y2)
	tmp1 := newFieldElement(nil, curve.fieldOrder).Add(curve.fieldOne, commonTerm)
	tmp1.ModInverse(tmp1)

	// 1 / (1 - d*x1*x2*y1*y1)
	tmp2 := newFieldElement(nil, curve.fieldOrder).Sub(curve.fieldOne, commonTerm)
	tmp2.ModInverse(tmp2)

	// x3 = (x1*y2 + x2*y1) / (1 + d*x1*x2*y1*y1)
	x3 := newFieldElement(nil, curve.fieldOrder)
	x3.Add(x1y2, x2y1).Mul(x3, tmp1)

	// y3 = (y1*y2 + x1*x2) / (1 - d*x1*x2*y1*y1)
	y3 := newFieldElement(nil, curve.fieldOrder)
	y3.Add(y1y2, x1x2).Mul(y3, tmp2)

	return &Point{curve, x3, y3}
}

// Double adds p1+p1 and returns a newly allocated result point.
func (curve *Jubjub) Double(p1 *Point) *Point {
	// Affine doubling formulas: 2(x1,y1) = (x3,y3) where
	// x3 = (x1*y1 + y1*x1) / (1 + d*x1*x1*y1*y1)
	// y3 = (y1*y1 - a*x1*x1) / (1 - d*x1*x1*y1*y1)
	// Recall a = -1

	x1x1 := newFieldElement(nil, curve.fieldOrder).Mul(p1.x, p1.x)
	y1y1 := newFieldElement(nil, curve.fieldOrder).Mul(p1.y, p1.y)
	x1y1 := newFieldElement(nil, curve.fieldOrder).Mul(p1.x, p1.y)

	// d*x1*x1*y1*y1
	commonTerm := newFieldElement(nil, curve.fieldOrder).Mul(x1x1, y1y1)
	commonTerm.Mul(commonTerm, curve.d)

	// 1 / (1 + d*x1*x1*y1*y1)
	tmp := newFieldElement(nil, curve.fieldOrder).Add(curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// x3 = (x1*y1 + y1*x1) / (1 + d*x1*x1*y1*y1)
	x3 := newFieldElement(nil, curve.fieldOrder)
	x3.Add(x1y1, x1y1).Mul(x3, tmp)

	// 1 / (1 - d*x1*x1*y1*y1)
	tmp.Sub(curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// y3 = (y1*y1 - a*x1*x1) / (1 - d*x1*x1*y1*y1)
	y3 := newFieldElement(nil, curve.fieldOrder)
	y3.Add(y1y1, x1x1).Mul(y3, tmp)

	return &Point{curve, x3, y3}
}

// Add sets p to the sum p1+p2 and returns p.
func (p *Point) Add(p1 *Point, p2 *Point) *Point {
	// Affine addition formulas: (x1,y1) + (x2,y2) = (x3,y3) where
	//  x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
	//  y3 = (y1*y2 - a*x1*x2) / (1 - d*x1*x2*y1*y2)
	// Recall a = -1

	x1y2 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.x, p2.y)
	x2y1 := newFieldElement(nil, p.curve.fieldOrder).Mul(p2.x, p1.y)

	y1y2 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.y, p2.y)
	x1x2 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.x, p2.x)

	// d*x1*x2*y1*y2
	commonTerm := newFieldElement(nil, p.curve.fieldOrder).Mul(x1x2, y1y2)
	commonTerm.Mul(commonTerm, p.curve.d)

	// 1 / (1 + d*x1*x2*y1*y2)
	tmp := newFieldElement(nil, p.curve.fieldOrder).Add(p.curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// x3 = (x1*y2 + x2*y1) / (1 + d*x1*x2*y1*y1)
	p.x.Add(x1y2, x2y1).Mul(p.x, tmp)

	// 1 / (1 - d*x1*x2*y1*y1)
	tmp.Sub(p.curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// y3 = (y1*y2 + x1*x2) / (1 - d*x1*x2*y1*y1)
	p.y.Add(y1y2, x1x2).Mul(p.y, tmp)

	return p
}

// Double sets p to the sum p1+p1 and returns p.
func (p *Point) Double(p1 *Point) *Point {
	// Affine doubling formulas: 2(x1,y1) = (x3,y3) where
	// x3 = (x1*y1 + y1*x1) / (1 + d*x1*x1*y1*y1)
	// y3 = (y1*y1 - a*x1*x1) / (1 - d*x1*x1*y1*y1)
	// Recall a = -1

	x1x1 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.x, p1.x)
	y1y1 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.y, p1.y)
	x1y1 := newFieldElement(nil, p.curve.fieldOrder).Mul(p1.x, p1.y)

	// d*x1*x1*y1*y1
	commonTerm := newFieldElement(nil, p.curve.fieldOrder).Mul(x1x1, y1y1)
	commonTerm.Mul(commonTerm, p.curve.d)

	// 1 / (1 + d*x1*x1*y1*y1)
	tmp := newFieldElement(nil, p.curve.fieldOrder).Add(p.curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// x3 = (x1*y1 + y1*x1) / (1 + d*x1*x1*y1*y1)
	p.x.Add(x1y1, x1y1).Mul(p.x, tmp)

	// 1 / (1 - d*x1*x1*y1*y1)
	tmp.Sub(p.curve.fieldOne, commonTerm)
	tmp.ModInverse(tmp)

	// y3 = (y1*y1 - a*x1*x1) / (1 - d*x1*x1*y1*y1)
	p.y.Add(y1y1, x1x1).Mul(p.y, tmp)

	return p
}
