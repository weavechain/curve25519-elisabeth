/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import java.util.ArrayList;
import java.util.List;

public class MulUtils {

    public static RistrettoElement multiscalarMul(Scalar s1, List<Scalar> s2, List<Scalar> s3, RistrettoElement p1, List<RistrettoElement> p2, List<RistrettoElement> p3) {
        RistrettoElement res = p1.multiply(s1);

        if (s2 != null && p2 != null && s2.size() == p2.size()) {
            for (int i = 0; i < s2.size(); i++) {
                res = res.add(p2.get(i).multiply(s2.get(i)));
            }
        }

        if (s3 != null && p3 != null && s3.size() == p3.size()) {
            for (int i = 0; i < s3.size(); i++) {
                res = res.add(p3.get(i).multiply(s3.get(i)));
            }
        }

        return res;
    }

    public static RistrettoElement mulStraus(List<Scalar> s, List<RistrettoElement> p) {
        RistrettoElement res = RistrettoElement.IDENTITY;

        List<ProjectiveNielsPoint.LookupTable> lookups = new ArrayList<>();
        for (RistrettoElement it : p) {
            lookups.add(it.lookupTable());
        }

        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2(4);

            for (int j = 0; j < p.size(); j++) {
                final byte[] e = s.get(j).toRadix16();

                Q = Q.add(lookups.get(j).select(e[i])).toExtended();
            }
        }
        res = res.add(new RistrettoElement(Q));

        return res;
    }

    public static RistrettoElement mulStraus(List<Scalar> s1, List<RistrettoElement> p1, List<Scalar> s2, List<RistrettoElement> p2) {
        RistrettoElement res = RistrettoElement.IDENTITY;

        List<ProjectiveNielsPoint.LookupTable> lookups = new ArrayList<>();
        if (p1 != null) {
            for (RistrettoElement it : p1) {
                lookups.add(it.lookupTable());
            }
        }
        int offset = lookups.size();
        if (p2 != null) {
            for (RistrettoElement it : p2) {
                lookups.add(it.lookupTable());
            }
        }

        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2(4);

            if (p1 != null) {
                for (int j = 0; j < p1.size(); j++) {
                    final byte[] e = s1.get(j).toRadix16();
                    Q = Q.add(lookups.get(j).select(e[i])).toExtended();
                }
            }
            if (p2 != null) {
                for (int j = 0; j < p2.size(); j++) {
                    final byte[] e = s2.get(j).toRadix16();
                    Q = Q.add(lookups.get(j + offset).select(e[i])).toExtended();
                }
            }
        }
        res = res.add(new RistrettoElement(Q));

        return res;
    }

    public static RistrettoElement multiscalarMulStraus(Scalar s1, List<Scalar> s2, List<Scalar> s3, RistrettoElement p1, List<RistrettoElement> p2, List<RistrettoElement> p3) {
        RistrettoElement res = p1.multiply(s1);

        if (s2 != null && p2 != null && s2.size() == p2.size() || s3 != null && p3 != null && s3.size() == p3.size()) {
            res = res.add(mulStraus(s2, p2, s3, p3));
        }

        return res;
    }
}
