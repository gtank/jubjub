package jubjub

import (
	"bytes"
	"encoding/hex"
	"math/big"
	"testing"
	"testing/quick"
)

var compressedPoints = []string{
	// From https://github.com/zcash-hackworks/zcash-test-vectors/blob/master/sapling_note_encryption.py
	"db4cd2b0aac4f7eb8ca131f16567c445a9555126d3c29f14e3d776e841ae7415",
	"a6b13ea336ddb7a67bb09a0e68e9d3cfb39210831ea3a296ba09a922060fd38b",
	"66141739514b28f05def8a18eeee5eed4d44c6225c3c65d88dd9907708012f5a",
	"25eb55fccf761fc64e85a588efe6ead7832fb1f0f7a83165895bdff942925f5c",
	"8b2a337f03622c24ff381d4c546f6977f90522e92fde44c9d1bb099714b9db2b",
	"6b27daccb5a8207f532d10ca238f9786648a11b5966e51a2f7d89e15d29b8fdf",
	"d11da01f0b43bdd5288d32385b8771d223493c69802544043f77cf1d71c1cb8c",
	"32cb2806b882f1368b0d4a898f72c4c8f728132cc12456946e7f4cb0fb058da9",
	"9e64174b4ab981405c323b5e12475945a46d4fedf8060828041cd20e62fd2cef",
	"b68e9ee0c0678d7b3036931c831a25255f7ee487385a30316e15f6482b874fda",
}

// quickCheckConfig will make each quickcheck test run (1024 * -quickchecks)
// times. The default value of -quickchecks is 100, indicated by 0.
var quickCheckConfig = &quick.Config{MaxCountScale: 16} // 1024 / 16 = 64 quickchecks

func TestCompressionRoundtrip(t *testing.T) {
	curve := Curve()
	for i, s := range compressedPoints {
		compressed, _ := hex.DecodeString(s)
		point, err := curve.Decompress(compressed)
		if err != nil {
			t.Error(err)
			continue
		}

		rt, _ := point.MarshalBinary()

		if !bytes.Equal(rt, compressed) {
			t.Errorf("Incorrect result for test %d:\nWant: %x\nHave: %x", i, compressed, rt)
		}
	}
}

func TestGeneratorBehavior(t *testing.T) {
	curve := Curve()

	// This is the affine generator from the jubjub Rust crate.
	gX, _ := hex.DecodeString("feada7f15dd3b3e4af81bf291b5df5ca87810ad6dd030f8bc88737bfb8cbed62")
	gY := newFieldElement(big.NewInt(11), curve.fieldOrder).ToBytes()

	// We should get the expected x-coordinate when we decompress the encoded y-coordinate.
	group, err := curve.Decompress(gY)
	if err != nil {
		t.Fatal("Couldn't decompress generator")
	}

	if !bytes.Equal(group.x.ToBytes(), gX) {
		t.Fatal("Decompressed to different generator than expected")
	}

	subgroup := group.Clone().MulByCofactor()
	if err != nil {
		t.Error(err)
	}

	if !subgroup.Equals(curve.SubgroupGenerator()) {
		t.Fatal("Construction resulted in a different subgroup generator")
	}

	// Multiplying the generator by the order of the group should yield the identity point.
	subgroupOrder, _ := new(big.Int).SetString("6554484396890773809930967563523245729705921265872317281365359162392183254199", 10)

	scalar, err := curve.ScalarFromBig(subgroupOrder)
	if err != nil {
		t.Error("didn't like subgroup order scalar")
	}

	identity, _ := curve.ScalarMult(scalar, subgroup)
	if !identity.Equals(curve.Identity()) {
		t.Fatal("q*8*G != (0, 1)")
	}
}

func BenchmarkScalarMult(b *testing.B) {
	var p *Point

	curve := Curve()
	g := curve.Generator()

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		scalar, err := newScalar(big.NewInt(int64(i)), curve.subgroupOrder)
		if err != nil {
			b.Fatal("Wow! The iteration counter exceeded scalar range!")
		}
		p, _ = curve.ScalarMult(scalar, g)
	}
	if !p.IsOnCurve() {
		b.Fail()
	}
}

func TestPointClone(t *testing.T) {
	curve := Curve()
	g := curve.Generator()
	cloned := g.Clone()

	if !g.Equals(cloned) {
		t.Fatal("Clone is broken")
	}

	g.x.Set(curve.newFieldElement(big.NewInt(15)))

	if g.Equals(cloned) {
		t.Fatal("Clone is broken")
	}
}

func TestBasicProperties(t *testing.T) {
	curve := Curve()
	G := curve.Generator()

	scZero, _ := curve.ScalarFromBig(big.NewInt(0))
	scOne, _ := curve.ScalarFromBig(big.NewInt(1))
	scTwo, _ := curve.ScalarFromBig(big.NewInt(2))
	scThree, _ := curve.ScalarFromBig(big.NewInt(3))

	// Does a point equal itself?
	{
		if !G.Equals(G) {
			t.Error("Point equality is broken")
		}
	}

	// Does scalar multiplication by 0 return the identity?
	{
		zeroG, _ := curve.ScalarMult(scZero, G)
		if !zeroG.Equals(curve.Identity()) {
			t.Error("0*G != (0, 1)")
		}
	}

	// Does negation work?
	{
		identity := newPoint(curve)
		identity = identity.Neg(G).Add(identity, G)

		if !identity.Equals(curve.Identity()) {
			t.Error("Negation isn't working")
		}
	}

	if !testing.Short() {
		negationWorks := func(x uint64) bool {
			scalar, _ := curve.ScalarFromBig(big.NewInt(int64(x)))
			point, _ := curve.ScalarMult(scalar, curve.Generator())

			negPoint := point.Clone().Neg(point)
			identity := curve.Add(point, negPoint)

			return identity.Equals(curve.Identity())
		}

		if err := quick.Check(negationWorks, quickCheckConfig); err != nil {
			t.Error("quickcheck: negation doesn't work")
		}
	}

	// Do add and multiply agree?
	{
		addOneG := G.Clone()
		addTwoG := newPoint(curve).Add(G, G)

		addThreeG := G.Clone()
		addThreeG.Add(G, G).Add(addThreeG, G)

		mulOneG, _ := curve.ScalarMult(scOne, G)
		mulTwoG, _ := curve.ScalarMult(scTwo, G)
		mulThreeG, _ := curve.ScalarMult(scThree, G)

		if !addOneG.Equals(mulOneG) {
			t.Error("G != 1*G")
		}

		if !addTwoG.Equals(mulTwoG) {
			t.Error("G+G != 2*G")
		}

		if !addThreeG.Equals(mulThreeG) {
			t.Error("G+G+G != 3*G")
		}
	}

	if !testing.Short() {
		addMulAgree := func(x uint8) bool {
			accumulator := curve.Identity()
			point := curve.SubgroupGenerator()

			for i := 0; i < int(x); i++ {
				accumulator.Add(accumulator, point)
			}

			scalar, _ := curve.ScalarFromBig(big.NewInt(int64(x)))
			mulPoint, _ := curve.ScalarMult(scalar, point)

			return accumulator.Equals(mulPoint)
		}

		if err := quick.Check(addMulAgree, quickCheckConfig); err != nil {
			t.Error("quickcheck: add and multiply don't agree")
		}
	}

	// Do add and double agree?
	{
		addTwoG := newPoint(curve).Add(G, G)
		mulTwoG, _ := curve.ScalarMult(scTwo, G)
		doubledG := newPoint(curve).Double(G)

		if !addTwoG.Equals(doubledG) {
			t.Error("Adding twice != doubled")
		}

		if !doubledG.Equals(mulTwoG) {
			t.Error("Multiplying by two != doubled")
		}
	}

	if !testing.Short() {
		addDoubleAgree := func(x uint8) bool {
			accumulator := curve.SubgroupGenerator()
			doubler := curve.SubgroupGenerator()

			for i := 0; i < int(x); i++ {
				accumulator.Add(accumulator, accumulator)
				doubler.Double(doubler)
			}

			return accumulator.Equals(doubler)
		}

		if err := quick.Check(addDoubleAgree, quickCheckConfig); err != nil {
			t.Error("quickcheck: add and double don't agree")
		}
	}

	// Is scalar multiplication associative?
	{
		mulOneG, _ := curve.ScalarMult(scOne, G)
		mulTwoG, _ := curve.ScalarMult(scTwo, G)
		mulThreeG, _ := curve.ScalarMult(scThree, G)

		// 1G + 2G == (1+2)G
		addThreeG := newPoint(curve).Add(mulOneG, mulTwoG)
		if !addThreeG.Equals(mulThreeG) {
			t.Error("Scalar multiplication is not associative")
		}

		// 2*3G == 3G + 3G
		mulSixG, _ := curve.ScalarMult(scTwo, addThreeG)
		addSixG := newPoint(curve).Add(mulThreeG, mulThreeG)

		if !mulSixG.Equals(addSixG) {
			t.Error("Scalar multiplication is not associative")
		}
	}
}
