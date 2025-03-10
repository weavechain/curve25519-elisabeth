/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import org.hamcrest.Matchers;
import org.junit.*;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.is;

public class ConstantsTest {
    @Test
    public void checkEdwardsD() {
        assertThat(Constants.EDWARDS_D, is(FieldElement
                .fromByteArray(Utils.hexToBytes("a3785913ca4deb75abd841414d0a700098e879777940c78c73fe6f2bee6c0352"))));
    }

    @Test
    public void checkEdwards2D() {
        FieldElement two = FieldElement.ONE.add(FieldElement.ONE);
        assertThat(Constants.EDWARDS_2D, is(Constants.EDWARDS_D.multiply(two)));
    }

    @Test
    public void checkSqrtADMinusOne() {
        System.out.println(
                FieldElement.fromByteArray(Constants.SQRT_AD_MINUS_ONE.toByteArray()).printInternalRepresentation());
        assertThat(Constants.SQRT_AD_MINUS_ONE.square().add(FieldElement.ONE).negate(), Matchers.is(Constants.EDWARDS_D));
    }

    @Test
    public void checkInvSqrtAMinusD() {
        assertThat(Constants.INVSQRT_A_MINUS_D.invert().square().add(FieldElement.ONE).negate(),
                Matchers.is(Constants.EDWARDS_D));
    }

    @Test
    public void checkSqrtM1() {
        assertThat(Constants.SQRT_M1, is(FieldElement
                .fromByteArray(Utils.hexToBytes("b0a00e4a271beec478e42fad0618432fa7d7fb3d99004d2b0bdfc14f8024832b"))));
    }

    @Test
    public void checkEd25519Basepoint() throws InvalidEncodingException {
        CompressedEdwardsY encoded = new CompressedEdwardsY(
                Utils.hexToBytes("5866666666666666666666666666666666666666666666666666666666666666"));
        EdwardsPoint B = encoded.decompress();
        assertThat(Constants.ED25519_BASEPOINT.X, Matchers.is(B.X));
        assertThat(Constants.ED25519_BASEPOINT.Y, Matchers.is(B.Y));
        assertThat(Constants.ED25519_BASEPOINT.Z, Matchers.is(B.Z));
        assertThat(Constants.ED25519_BASEPOINT.T, Matchers.is(B.T));
    }
}
