/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class MultiscalarMulTest {

    static final String[] FROM_UNIFORM_BYTES_INPUTS = new String[] {
            "5d1be09e3d0c82fc538112490e35701979d99e06ca3e2b5b54bffe8b4dc772c14d98b696a1bbfb5ca32c436cc61c16563790306c79eaca7705668b47dffe5bb6",
            "f116b34b8f17ceb56e8732a60d913dd10cce47a6d53bee9204be8b44f6678b270102a56902e2488c46120e9276cfe54638286b9e4b3cdb470b542d46c2068d38",
            "8422e1bbdaab52938b81fd602effb6f89110e1e57208ad12d9ad767e2e25510c27140775f9337088b982d83d7fcf0b2fa1edffe51952cbe7365e95c86eaf325c",
            "ac22415129b61427bf464e17baee8db65940c233b98afce8d17c57beeb7876c2150d15af1cb1fb824bbd14955f2b57d08d388aab431a391cfc33d5bafb5dbbaf",
            "165d697a1ef3d5cf3c38565beefcf88c0f282b8e7dbd28544c483432f1cec7675debea8ebb4e5fe7d6f6e5db15f15587ac4d4d4a1de7191e0c1ca6664abcc413",
            "a836e6c9a9ca9f1e8d486273ad56a78c70cf18f0ce10abb1c7172ddd605d7fd2979854f47ae1ccf204a33102095b4200e5befc0465accc263175485f0e17ea5c",
            "2cdc11eaeb95daf01189417cdddbf95952993aa9cb9c640eb5058d09702c74622c9965a697a3b345ec24ee56335b556e677b30e6f90ac77d781064f866a3c982" };


    @Test
    public void testMul() {
        RistrettoElement p1 = RistrettoElement.fromUniformBytes(Utils.hexToBytes(FROM_UNIFORM_BYTES_INPUTS[0]));
        Scalar s1 = Scalar.fromBits(Utils.hexToBytes(FROM_UNIFORM_BYTES_INPUTS[0].substring(0, 64)));

        List<RistrettoElement> p2 = new ArrayList<>();
        List<Scalar> s2 = new ArrayList<>();

        List<RistrettoElement> p3 = new ArrayList<>();
        List<Scalar> s3 = new ArrayList<>();

        for (int i = 0; i < FROM_UNIFORM_BYTES_INPUTS.length; i++) {
            RistrettoElement p = RistrettoElement.fromUniformBytes(Utils.hexToBytes(FROM_UNIFORM_BYTES_INPUTS[i]));
            Scalar s = Scalar.fromBits(Utils.hexToBytes(FROM_UNIFORM_BYTES_INPUTS[0].substring(0, 64)));
            p2.add(p);
            s2.add(s);

            for (int j = 0; j < 100; j++) {
                p3.add(p);
                s3.add(s);
            }
        }

        //warmup
        MulUtils.multiscalarMul(s1, s2, s3, p1, p2, p3);
        MulUtils.multiscalarMulStraus(s1, s2, s3, p1, p2, p3);
        MulUtils.multiscalarMulPippenger(s1, s2, s3, p1, p2, p3);

        RistrettoElement r4 = MulUtils.mulStraus(s1, s2.get(0), p1, p2.get(0));
        RistrettoElement r5 = p1.multiply(s1).add(p2.get(0).multiply(s2.get(0)));
        assertTrue(r4.equals(r5));

        long t1 = System.currentTimeMillis();
        RistrettoElement r1 = MulUtils.multiscalarMul(s1, s2, s3, p1, p2, p3);
        long t2 = System.currentTimeMillis();
        RistrettoElement r2 = MulUtils.multiscalarMulStraus(s1, s2, s3, p1, p2, p3);
        long t3 = System.currentTimeMillis();
        RistrettoElement r3 = MulUtils.multiscalarMulPippenger(s1, s2, s3, p1, p2, p3);
        long t4 = System.currentTimeMillis();

        System.out.println((t2 - t1) + " vs " + (t3 - t2) + " vs " + (t4 - t3));
        assertTrue(r1.equals(r2));
        assertTrue(r1.equals(r3));

        int iter = 1000;
        for (int k = 5; k < 100; k += 10) {
            List<Scalar> st = s3.subList(0, k);
            List<RistrettoElement> pt = p3.subList(0, k);

            long m1 = System.currentTimeMillis();
            for (int i = 0; i < iter; i++) {
                MulUtils.mulStraus(st, pt);
            }
            long m2 = System.currentTimeMillis();
            for (int i = 0; i < iter; i++) {
                MulUtils.mulPippenger(st, pt);
            }
            long m3 = System.currentTimeMillis();
            for (int i = 0; i < iter; i++) {
                MulUtils.multiscalarMulOpt(st, pt);
            }
            long m4 = System.currentTimeMillis();
            System.out.println(k + ": " + (m2 - m1) + " vs " + (m3 - m2) + " | " + (m4 - m3));
        }
    }
}
