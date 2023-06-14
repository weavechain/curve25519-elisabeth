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

public class EdwardsBasepointTableTest {
    @Test
    public void scalarMulVsEd25519py() {
        EdwardsBasepointTable Bt = new EdwardsBasepointTable(Constants.ED25519_BASEPOINT);
        EdwardsPoint aB = Bt.multiply(EdwardsPointTest.A_SCALAR);
        assertThat(aB.compress(), Matchers.is(EdwardsPointTest.A_TIMES_BASEPOINT));
    }
}
