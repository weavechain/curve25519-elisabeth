/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MulUtils {

    private static final double LOG2 = Math.log(2);

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

    public static RistrettoElement mulStraus(Scalar s1, Scalar s2, RistrettoElement p1, RistrettoElement p2) {
        ProjectiveNielsPoint.LookupTable n1 = p1.lookupTable();
        ProjectiveNielsPoint.LookupTable n2 = p2.lookupTable();

        byte[] b1 = s1.toRadix16();
        byte[] b2 = s2.toRadix16();

        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2(4);

            Q = Q.add(n1.select(b1[i])).toExtended();
            Q = Q.add(n2.select(b2[i])).toExtended();
        }
        return new RistrettoElement(Q);
    }

    public static RistrettoElement mulStraus(List<Scalar> s, List<RistrettoElement> p) {
        List<ProjectiveNielsPoint.LookupTable> lookups = new ArrayList<>();
        for (RistrettoElement it : p) {
            lookups.add(it.lookupTable());
        }

        List<byte[]> scalars = new ArrayList<>();
        for (Scalar it : s) {
            scalars.add(it.toRadix16());
        }

        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2(4);

            for (int j = 0; j < p.size(); j++) {
                Q = Q.add(lookups.get(j).select(scalars.get(j)[i])).toExtended();
            }
        }
        return new RistrettoElement(Q);
    }

    public static RistrettoElement mulPippenger(List<Scalar> s, List<RistrettoElement> p) {
        //int l = (int)Math.max(1, Math.log(p.size() / LOG2));
        //int c = Math.min(8, Math.max(6, l - l / 3));
        int c = p.size() < 500 ? 6 : p.size() < 800 ? 7 : 8;

        List<ProjectiveNielsPoint> pts = new ArrayList<>();
        for (RistrettoElement it : p) {
            pts.add(it.repr.toProjectiveNiels());
        }
        List<byte[]> scalars = new ArrayList<>();
        for (Scalar it : s) {
            scalars.add(it.toRadix2w(c));
        }

        int bucketsCount = 1 << (c - 1);
        int digits = c == 8 ? (256 + c - 1) / c + 1 : (256 + c - 1) / c;

        final ArrayList<EdwardsPoint> identity = new ArrayList<>(Collections.nCopies(bucketsCount, EdwardsPoint.IDENTITY));

        EdwardsPoint Q = null;
        for (int k = digits - 1; k >= 0; k--) {
            final ArrayList<EdwardsPoint> buckets = new ArrayList<>(identity);
            for (int i = 0; i < p.size(); i++) {
                int d = scalars.get(i)[k];
                if (d != 0) {
                    int idx = d > 0 ? d - 1 : -d - 1;
                    EdwardsPoint pt = (d > 0 ? buckets.get(idx).add(pts.get(i)) : buckets.get(idx).subtract(pts.get(i))).toExtended();
                    buckets.set(idx, pt);
                }
            }

            EdwardsPoint sum = buckets.get(bucketsCount - 1);
            EdwardsPoint bsum = buckets.get(bucketsCount - 1);
            for (int i = bucketsCount - 2; i >= 0; i--) {
                sum = sum.add(buckets.get(i));
                bsum = bsum.add(sum);
            }

            if (Q == null) {
                Q = bsum;
            } else {
                Q = Q.multiplyByPow2(c).add(bsum);
            }
        }

        return new RistrettoElement(Q);
    }

    public static RistrettoElement mulStraus(List<Scalar> s1, List<RistrettoElement> p1, List<Scalar> s2, List<RistrettoElement> p2) {
        List<ProjectiveNielsPoint.LookupTable> lookups = new ArrayList<>();
        List<byte[]> scalars = new ArrayList<>();
        if (p1 != null) {
            for (RistrettoElement it : p1) {
                lookups.add(it.lookupTable());
            }
            for (Scalar it : s1) {
                scalars.add(it.toRadix16());
            }
        }
        int offset = lookups.size();
        if (p2 != null) {
            for (RistrettoElement it : p2) {
                lookups.add(it.lookupTable());
            }
            for (Scalar it : s2) {
                scalars.add(it.toRadix16());
            }
        }

        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2(4);

            if (p1 != null) {
                for (int j = 0; j < p1.size(); j++) {
                    final byte[] e = scalars.get(j);
                    Q = Q.add(lookups.get(j).select(e[i])).toExtended();
                }
            }
            if (p2 != null) {
                for (int j = 0; j < p2.size(); j++) {
                    final byte[] e = scalars.get(j + offset);
                    Q = Q.add(lookups.get(j + offset).select(e[i])).toExtended();
                }
            }
        }

        return new RistrettoElement(Q);
    }

    public static RistrettoElement multiscalarMulStraus(Scalar s1, List<Scalar> s2, List<Scalar> s3, RistrettoElement p1, List<RistrettoElement> p2, List<RistrettoElement> p3) {
        RistrettoElement res = p1.multiply(s1);

        if (s2 != null && p2 != null && s2.size() == p2.size() || s3 != null && p3 != null && s3.size() == p3.size()) {
            res = res.add(mulStraus(s2, p2, s3, p3));
        }

        return res;
    }

    public static RistrettoElement multiscalarMulPippenger(Scalar s1, List<Scalar> s2, List<Scalar> s3, RistrettoElement p1, List<RistrettoElement> p2, List<RistrettoElement> p3) {
        int len2 = s2 != null && p2 != null && s2.size() == p2.size() ? p2.size() : 0;
        int len3 = s3 != null && p3 != null && s3.size() == p3.size() ? p3.size() : 0;
        if (len2 > 0 || len3 > 0) {
            List<Scalar> ss = new ArrayList<>();
            List<RistrettoElement> pp = new ArrayList<>();
            ss.add(s1);
            pp.add(p1);

            if (s2 != null && p2 != null && s2.size() == p2.size()) {
                for (int i = 0; i < s2.size(); i++) {
                    ss.add(s2.get(i));
                    pp.add(p2.get(i));
                }
            }

            if (s3 != null && p3 != null && s3.size() == p3.size()) {
                for (int i = 0; i < s3.size(); i++) {
                    ss.add(s3.get(i));
                    pp.add(p3.get(i));
                }
            }

            return MulUtils.mulPippenger(ss, pp);
        } else {
            return p1.multiply(s1);
        }
    }

    public static RistrettoElement multiscalarMulOpt(Scalar s1, List<Scalar> s2, List<Scalar> s3, RistrettoElement p1, List<RistrettoElement> p2, List<RistrettoElement> p3) {
        int len2 = s2 != null && p2 != null && s2.size() == p2.size() ? p2.size() : 0;
        int len3 = s3 != null && p3 != null && s3.size() == p3.size() ? p3.size() : 0;
        if (len2 > 0 || len3 > 0) {
            List<Scalar> ss = new ArrayList<>();
            List<RistrettoElement> pp = new ArrayList<>();
            ss.add(s1);
            pp.add(p1);

            if (s2 != null && p2 != null && s2.size() == p2.size()) {
                for (int i = 0; i < s2.size(); i++) {
                    ss.add(s2.get(i));
                    pp.add(p2.get(i));
                }
            }

            if (s3 != null && p3 != null && s3.size() == p3.size()) {
                for (int i = 0; i < s3.size(); i++) {
                    ss.add(s3.get(i));
                    pp.add(p3.get(i));
                }
            }

            return len2 + len3 >= 30 ? MulUtils.mulPippenger(ss, pp) : MulUtils.mulStraus(ss, pp);
        } else {
            return p1.multiply(s1);
        }
    }

    public static RistrettoElement multiscalarMulOpt(List<Scalar> ss, List<RistrettoElement> pp) {
        return pp.size() >= 30 ? MulUtils.mulPippenger(ss, pp) : MulUtils.mulStraus(ss, pp);
    }
}
