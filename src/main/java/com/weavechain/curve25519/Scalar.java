/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import java.io.IOException;
import java.io.InvalidObjectException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

import com.weavechain.subtle.ConstantTime;

/**
 * An integer $s \lt 2^{255}$ which represents an element of the field
 * $\mathbb{Z} / \ell$.
 */
public class Scalar implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final Scalar ZERO = new Scalar(new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
    public static final Scalar ONE = new Scalar(new byte[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

    /**
     * The 32-byte little-endian encoding of an integer representing a scalar modulo
     * the group order.
     * <p>
     * Invariant: the highest bit must be zero ($s[31] \le 127$).
     */
    private transient byte[] s;

    Scalar(byte[] s) {
        if (s.length != 32 || (((s[31] >> 7) & 0x01) != 0)) {
            throw new IllegalArgumentException("Invalid scalar representation");
        }
        // Store a copy to prevent interior mutability
        this.s = Arrays.copyOf(s, s.length);
    }

    public boolean test(int bit) {
        int idx = bit / 8;
        return (s[idx] & (bit % 8)) != 0;
    }

    public int get(int idx) {
        return s[idx];
    }

    /**
     * Overrides class serialization to use the canonical encoded format.
     */
    private void writeObject(ObjectOutputStream out) throws IOException {
        out.write(this.toByteArray());
    }

    /**
     * Overrides class serialization to use the canonical encoded format.
     */
    private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
        byte[] scalar = new byte[32];
        in.readFully(scalar);
        if (((scalar[31] >> 7) & 0x01) != 0) {
            throw new InvalidObjectException("Invalid scalar representation");
        }
        this.s = scalar;
    }

    @SuppressWarnings("unused")
    private void readObjectNoData() throws ObjectStreamException {
        throw new InvalidObjectException("Cannot deserialize Scalar from no data");
    }

    /**
     * Construct a Scalar by reducing a 256-bit little-endian integer modulo the
     * group order $\ell$.
     */
    public static Scalar fromBytesModOrder(final byte[] input) {
        if (input.length != 32) {
            throw new IllegalArgumentException("Input must by 32 bytes");
        }

        return Scalar.reduce(input);
    }

    /**
     * Construct a Scalar by reducing a 512-bit little-endian integer modulo the
     * group order $\ell$.
     */
    public static Scalar fromBytesModOrderWide(byte[] input) {
        if (input.length != 64) {
            throw new IllegalArgumentException("Input must by 64 bytes");
        }

        // s0,..., s22 have 21 bits, s23 has 29 bits
        long s0 = 0x1FFFFF & FieldElement.load_3(input, 0);
        long s1 = 0x1FFFFF & (FieldElement.load_4(input, 2) >> 5);
        long s2 = 0x1FFFFF & (FieldElement.load_3(input, 5) >> 2);
        long s3 = 0x1FFFFF & (FieldElement.load_4(input, 7) >> 7);
        long s4 = 0x1FFFFF & (FieldElement.load_4(input, 10) >> 4);
        long s5 = 0x1FFFFF & (FieldElement.load_3(input, 13) >> 1);
        long s6 = 0x1FFFFF & (FieldElement.load_4(input, 15) >> 6);
        long s7 = 0x1FFFFF & (FieldElement.load_3(input, 18) >> 3);
        long s8 = 0x1FFFFF & FieldElement.load_3(input, 21);
        long s9 = 0x1FFFFF & (FieldElement.load_4(input, 23) >> 5);
        long s10 = 0x1FFFFF & (FieldElement.load_3(input, 26) >> 2);
        long s11 = 0x1FFFFF & (FieldElement.load_4(input, 28) >> 7);
        long s12 = 0x1FFFFF & (FieldElement.load_4(input, 31) >> 4);
        long s13 = 0x1FFFFF & (FieldElement.load_3(input, 34) >> 1);
        long s14 = 0x1FFFFF & (FieldElement.load_4(input, 36) >> 6);
        long s15 = 0x1FFFFF & (FieldElement.load_3(input, 39) >> 3);
        long s16 = 0x1FFFFF & FieldElement.load_3(input, 42);
        long s17 = 0x1FFFFF & (FieldElement.load_4(input, 44) >> 5);
        long s18 = 0x1FFFFF & (FieldElement.load_3(input, 47) >> 2);
        long s19 = 0x1FFFFF & (FieldElement.load_4(input, 49) >> 7);
        long s20 = 0x1FFFFF & (FieldElement.load_4(input, 52) >> 4);
        long s21 = 0x1FFFFF & (FieldElement.load_3(input, 55) >> 1);
        long s22 = 0x1FFFFF & (FieldElement.load_4(input, 57) >> 6);
        long s23 = (FieldElement.load_4(input, 60) >> 3);
        long carry0;
        long carry1;
        long carry2;
        long carry3;
        long carry4;
        long carry5;
        long carry6;
        long carry7;
        long carry8;
        long carry9;
        long carry10;
        long carry11;
        long carry12;
        long carry13;
        long carry14;
        long carry15;
        long carry16;

        // @formatter:off
        /**
         * Lots of magic numbers :)
         * To understand what's going on below, note that
         *
         * (1) q = 2^252 + q0 where q0 = 27742317777372353535851937790883648493.
         * (2) s11 is the coefficient of 2^(11*21), s23 is the coefficient of 2^(23*21) and 2^252 = 2^((23-11) * 21)).
         * (3) 2^252 congruent -q0 modulo q.
         * (4) -q0 = 666643 * 2^0 + 470296 * 2^21 + 654183 * 2^(2*21) - 997805 * 2^(3*21) + 136657 * 2^(4*21) - 683901 * 2^(5*21)
         *
         * Thus
         * s23 * 2^(23*11) = s23 * 2^(12*21) * 2^(11*21) = s23 * 2^252 * 2^(11*21) congruent
         * s23 * (666643 * 2^0 + 470296 * 2^21 + 654183 * 2^(2*21) - 997805 * 2^(3*21) + 136657 * 2^(4*21) - 683901 * 2^(5*21)) * 2^(11*21) modulo q =
         * s23 * (666643 * 2^(11*21) + 470296 * 2^(12*21) + 654183 * 2^(13*21) - 997805 * 2^(14*21) + 136657 * 2^(15*21) - 683901 * 2^(16*21)).
         *
         * The same procedure is then applied for s22,...,s18.
         */
        // @formatter:on
        s11 += s23 * 666643;
        s12 += s23 * 470296;
        s13 += s23 * 654183;
        s14 -= s23 * 997805;
        s15 += s23 * 136657;
        s16 -= s23 * 683901;

        s10 += s22 * 666643;
        s11 += s22 * 470296;
        s12 += s22 * 654183;
        s13 -= s22 * 997805;
        s14 += s22 * 136657;
        s15 -= s22 * 683901;

        s9 += s21 * 666643;
        s10 += s21 * 470296;
        s11 += s21 * 654183;
        s12 -= s21 * 997805;
        s13 += s21 * 136657;
        s14 -= s21 * 683901;

        s8 += s20 * 666643;
        s9 += s20 * 470296;
        s10 += s20 * 654183;
        s11 -= s20 * 997805;
        s12 += s20 * 136657;
        s13 -= s20 * 683901;

        s7 += s19 * 666643;
        s8 += s19 * 470296;
        s9 += s19 * 654183;
        s10 -= s19 * 997805;
        s11 += s19 * 136657;
        s12 -= s19 * 683901;

        s6 += s18 * 666643;
        s7 += s18 * 470296;
        s8 += s18 * 654183;
        s9 -= s18 * 997805;
        s10 += s18 * 136657;
        s11 -= s18 * 683901;

        /**
         * Time to reduce the coefficient in order not to get an overflow.
         */
        // @formatter:off
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry12 = (s12 + (1<<20)) >> 21; s13 += carry12; s12 -= carry12 << 21;
        carry14 = (s14 + (1<<20)) >> 21; s15 += carry14; s14 -= carry14 << 21;
        carry16 = (s16 + (1<<20)) >> 21; s17 += carry16; s16 -= carry16 << 21;

        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12 += carry11; s11 -= carry11 << 21;
        carry13 = (s13 + (1<<20)) >> 21; s14 += carry13; s13 -= carry13 << 21;
        carry15 = (s15 + (1<<20)) >> 21; s16 += carry15; s15 -= carry15 << 21;
        // @formatter:on

        /**
         * Continue with above procedure.
         */
        s5 += s17 * 666643;
        s6 += s17 * 470296;
        s7 += s17 * 654183;
        s8 -= s17 * 997805;
        s9 += s17 * 136657;
        s10 -= s17 * 683901;

        s4 += s16 * 666643;
        s5 += s16 * 470296;
        s6 += s16 * 654183;
        s7 -= s16 * 997805;
        s8 += s16 * 136657;
        s9 -= s16 * 683901;

        s3 += s15 * 666643;
        s4 += s15 * 470296;
        s5 += s15 * 654183;
        s6 -= s15 * 997805;
        s7 += s15 * 136657;
        s8 -= s15 * 683901;

        s2 += s14 * 666643;
        s3 += s14 * 470296;
        s4 += s14 * 654183;
        s5 -= s14 * 997805;
        s6 += s14 * 136657;
        s7 -= s14 * 683901;

        s1 += s13 * 666643;
        s2 += s13 * 470296;
        s3 += s13 * 654183;
        s4 -= s13 * 997805;
        s5 += s13 * 136657;
        s6 -= s13 * 683901;

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        /**
         * Reduce coefficients again.
         */
        // @formatter:off
        carry0  = (s0  + (1<<20)) >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry2  = (s2  + (1<<20)) >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry4  = (s4  + (1<<20)) >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;

        carry1  = (s1  + (1<<20)) >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry3  = (s3  + (1<<20)) >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry5  = (s5  + (1<<20)) >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry11 = s11 >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        // @formatter:on

        // s0, ..., s11 got 21 bits each.
        byte[] result = new byte[32];
        result[0] = (byte) s0;
        result[1] = (byte) (s0 >> 8);
        result[2] = (byte) ((s0 >> 16) | (s1 << 5));
        result[3] = (byte) (s1 >> 3);
        result[4] = (byte) (s1 >> 11);
        result[5] = (byte) ((s1 >> 19) | (s2 << 2));
        result[6] = (byte) (s2 >> 6);
        result[7] = (byte) ((s2 >> 14) | (s3 << 7));
        result[8] = (byte) (s3 >> 1);
        result[9] = (byte) (s3 >> 9);
        result[10] = (byte) ((s3 >> 17) | (s4 << 4));
        result[11] = (byte) (s4 >> 4);
        result[12] = (byte) (s4 >> 12);
        result[13] = (byte) ((s4 >> 20) | (s5 << 1));
        result[14] = (byte) (s5 >> 7);
        result[15] = (byte) ((s5 >> 15) | (s6 << 6));
        result[16] = (byte) (s6 >> 2);
        result[17] = (byte) (s6 >> 10);
        result[18] = (byte) ((s6 >> 18) | (s7 << 3));
        result[19] = (byte) (s7 >> 5);
        result[20] = (byte) (s7 >> 13);
        result[21] = (byte) s8;
        result[22] = (byte) (s8 >> 8);
        result[23] = (byte) ((s8 >> 16) | (s9 << 5));
        result[24] = (byte) (s9 >> 3);
        result[25] = (byte) (s9 >> 11);
        result[26] = (byte) ((s9 >> 19) | (s10 << 2));
        result[27] = (byte) (s10 >> 6);
        result[28] = (byte) ((s10 >> 14) | (s11 << 7));
        result[29] = (byte) (s11 >> 1);
        result[30] = (byte) (s11 >> 9);
        result[31] = (byte) (s11 >> 17);
        return new Scalar(result);
    }

    /**
     * Attempt to construct a Scalar from a canonical byte representation.
     *
     * @return the Scalar if the input was its canonical representation.
     */
    public static Scalar fromCanonicalBytes(byte[] input) {
        Scalar s = new Scalar(input);
        if (!s.isCanonical()) {
            throw new IllegalArgumentException("Non-canonical scalar representation");
        }
        return s;
    }

    /**
     * Construct a Scalar from the low 255 bits of a 256-bit integer.
     * <p>
     * This function is intended for applications like X25519 which require specific
     * bit-patterns when performing scalar multiplication.
     */
    public static Scalar fromBits(byte[] input) {
        // Ensure that s < 2^255 by masking the high bit
        input[31] &= 0x7f;
        return new Scalar(input);
    }

    /**
     * Convert this Scalar to its underlying sequence of bytes.
     *
     * @return the 32-byte little-endian encoding of this Scalar.
     */
    public byte[] toByteArray() {
        // Return a copy to prevent interior mutability
        return Arrays.copyOf(this.s, this.s.length);
    }

    /**
     * Constant-time equality check.
     * <p>
     * Compares the encodings of the two Scalars.
     *
     * @return 1 if self and other are equal, 0 otherwise.
     */
    public int ctEquals(Scalar other) {
        return ConstantTime.equal(s, other.s);
    }

    /**
     * Equality check overridden to be constant-time.
     * <p>
     * Fails fast if the objects are of different types.
     *
     * @return true if self and other are equal, false otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Scalar)) {
            return false;
        }

        Scalar other = (Scalar) obj;
        return ctEquals(other) == 1;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(s);
    }

    @Override
    public String toString() {
        return "Scalar(" + StrUtils.bytesToHex(this.s) + ")";
    }

    /**
     * Check whether this Scalar is the canonical representative mod $\ell$.
     */
    boolean isCanonical() {
        return this.ctEquals(Scalar.reduce(this.s)) == 1;
    }

    /**
     * Reduce this Scalar modulo $\ell$.
     *
     * @return the reduced Scalar.
     */
    public Scalar reduce() {
        return Scalar.reduce(this.s);
    }

    /**
     * Reduce the given scalar modulo $\ell$.
     *
     * @return the reduced Scalar.
     */
    static Scalar reduce(byte[] x) {
        long[] xR = UnpackedScalar.fromByteArray(x).mulInternal(Constants.R);
        return new Scalar(UnpackedScalar.montgomeryReduce(xR).toByteArray());
    }

    /**
     * Compute $a + b \bmod \ell$.
     * <p>
     * If $a$ and $b$ are both canonical Scalars, the result is guaranteed to be a
     * canonical Scalar. In all other cases, the result may be outside the range
     * $[0, \ell)$.
     *
     * @param b the Scalar to add to this.
     * @return $a + b \bmod \ell$
     */
    public Scalar add(final Scalar b) {
        return new Scalar(UnpackedScalar.fromByteArray(this.s).add(UnpackedScalar.fromByteArray(b.s)).toByteArray());
    }

    /**
     * Compute $a - b \bmod \ell$.
     * <p>
     * If $a$ and $b$ are both canonical Scalars, the result is guaranteed to be a
     * canonical Scalar. In all other cases, the result may be outside the range
     * $[0, \ell)$.
     *
     * @param b the Scalar to subtract from this.
     * @return $a - b \bmod \ell$
     */
    public Scalar subtract(final Scalar b) {
        return new Scalar(
                UnpackedScalar.fromByteArray(this.s).subtract(UnpackedScalar.fromByteArray(b.s)).toByteArray());
    }

    /**
     * Compute $a * b \bmod \ell$.
     *
     * @param b the Scalar to multiply with this.
     * @return $a * b \bmod \ell$
     */
    public Scalar multiply(final Scalar b) {
        return this.multiplyAndAdd(b, Scalar.ZERO);
    }

    /**
     * Compute $a * b + c \bmod \ell$.
     *
     * @param b the Scalar to multiply this by.
     * @param c the Scalar to add to the product.
     * @return $a * b + c \bmod \ell$
     */
    public Scalar multiplyAndAdd(Scalar b, Scalar c) {
        long a0 = 0x1FFFFF & FieldElement.load_3(this.s, 0);
        long a1 = 0x1FFFFF & (FieldElement.load_4(this.s, 2) >> 5);
        long a2 = 0x1FFFFF & (FieldElement.load_3(this.s, 5) >> 2);
        long a3 = 0x1FFFFF & (FieldElement.load_4(this.s, 7) >> 7);
        long a4 = 0x1FFFFF & (FieldElement.load_4(this.s, 10) >> 4);
        long a5 = 0x1FFFFF & (FieldElement.load_3(this.s, 13) >> 1);
        long a6 = 0x1FFFFF & (FieldElement.load_4(this.s, 15) >> 6);
        long a7 = 0x1FFFFF & (FieldElement.load_3(this.s, 18) >> 3);
        long a8 = 0x1FFFFF & FieldElement.load_3(this.s, 21);
        long a9 = 0x1FFFFF & (FieldElement.load_4(this.s, 23) >> 5);
        long a10 = 0x1FFFFF & (FieldElement.load_3(this.s, 26) >> 2);
        long a11 = (FieldElement.load_4(this.s, 28) >> 7);
        long b0 = 0x1FFFFF & FieldElement.load_3(b.s, 0);
        long b1 = 0x1FFFFF & (FieldElement.load_4(b.s, 2) >> 5);
        long b2 = 0x1FFFFF & (FieldElement.load_3(b.s, 5) >> 2);
        long b3 = 0x1FFFFF & (FieldElement.load_4(b.s, 7) >> 7);
        long b4 = 0x1FFFFF & (FieldElement.load_4(b.s, 10) >> 4);
        long b5 = 0x1FFFFF & (FieldElement.load_3(b.s, 13) >> 1);
        long b6 = 0x1FFFFF & (FieldElement.load_4(b.s, 15) >> 6);
        long b7 = 0x1FFFFF & (FieldElement.load_3(b.s, 18) >> 3);
        long b8 = 0x1FFFFF & FieldElement.load_3(b.s, 21);
        long b9 = 0x1FFFFF & (FieldElement.load_4(b.s, 23) >> 5);
        long b10 = 0x1FFFFF & (FieldElement.load_3(b.s, 26) >> 2);
        long b11 = (FieldElement.load_4(b.s, 28) >> 7);
        long c0 = 0x1FFFFF & FieldElement.load_3(c.s, 0);
        long c1 = 0x1FFFFF & (FieldElement.load_4(c.s, 2) >> 5);
        long c2 = 0x1FFFFF & (FieldElement.load_3(c.s, 5) >> 2);
        long c3 = 0x1FFFFF & (FieldElement.load_4(c.s, 7) >> 7);
        long c4 = 0x1FFFFF & (FieldElement.load_4(c.s, 10) >> 4);
        long c5 = 0x1FFFFF & (FieldElement.load_3(c.s, 13) >> 1);
        long c6 = 0x1FFFFF & (FieldElement.load_4(c.s, 15) >> 6);
        long c7 = 0x1FFFFF & (FieldElement.load_3(c.s, 18) >> 3);
        long c8 = 0x1FFFFF & FieldElement.load_3(c.s, 21);
        long c9 = 0x1FFFFF & (FieldElement.load_4(c.s, 23) >> 5);
        long c10 = 0x1FFFFF & (FieldElement.load_3(c.s, 26) >> 2);
        long c11 = (FieldElement.load_4(c.s, 28) >> 7);
        long s0;
        long s1;
        long s2;
        long s3;
        long s4;
        long s5;
        long s6;
        long s7;
        long s8;
        long s9;
        long s10;
        long s11;
        long s12;
        long s13;
        long s14;
        long s15;
        long s16;
        long s17;
        long s18;
        long s19;
        long s20;
        long s21;
        long s22;
        long s23;
        long carry0;
        long carry1;
        long carry2;
        long carry3;
        long carry4;
        long carry5;
        long carry6;
        long carry7;
        long carry8;
        long carry9;
        long carry10;
        long carry11;
        long carry12;
        long carry13;
        long carry14;
        long carry15;
        long carry16;
        long carry17;
        long carry18;
        long carry19;
        long carry20;
        long carry21;
        long carry22;

        // @formatter:off
        s0  = c0  + a0*b0;
        s1  = c1  + a0*b1  + a1*b0;
        s2  = c2  + a0*b2  + a1*b1  + a2*b0;
        s3  = c3  + a0*b3  + a1*b2  + a2*b1  + a3*b0;
        s4  = c4  + a0*b4  + a1*b3  + a2*b2  + a3*b1  + a4*b0;
        s5  = c5  + a0*b5  + a1*b4  + a2*b3  + a3*b2  + a4*b1  + a5*b0;
        s6  = c6  + a0*b6  + a1*b5  + a2*b4  + a3*b3  + a4*b2  + a5*b1  + a6*b0;
        s7  = c7  + a0*b7  + a1*b6  + a2*b5  + a3*b4  + a4*b3  + a5*b2  + a6*b1  + a7*b0;
        s8  = c8  + a0*b8  + a1*b7  + a2*b6  + a3*b5  + a4*b4  + a5*b3  + a6*b2  + a7*b1  + a8*b0;
        s9  = c9  + a0*b9  + a1*b8  + a2*b7  + a3*b6  + a4*b5  + a5*b4  + a6*b3  + a7*b2  + a8*b1  + a9*b0;
        s10 = c10 + a0*b10 + a1*b9  + a2*b8  + a3*b7  + a4*b6  + a5*b5  + a6*b4  + a7*b3  + a8*b2  + a9*b1  + a10*b0;
        s11 = c11 + a0*b11 + a1*b10 + a2*b9  + a3*b8  + a4*b7  + a5*b6  + a6*b5  + a7*b4  + a8*b3  + a9*b2  + a10*b1  + a11*b0;
        s12 =                a1*b11 + a2*b10 + a3*b9  + a4*b8  + a5*b7  + a6*b6  + a7*b5  + a8*b4  + a9*b3  + a10*b2  + a11*b1;
        s13 =                         a2*b11 + a3*b10 + a4*b9  + a5*b8  + a6*b7  + a7*b6  + a8*b5  + a9*b4  + a10*b3  + a11*b2;
        s14 =                                  a3*b11 + a4*b10 + a5*b9  + a6*b8  + a7*b7  + a8*b6  + a9*b5  + a10*b4  + a11*b3;
        s15 =                                           a4*b11 + a5*b10 + a6*b9  + a7*b8  + a8*b7  + a9*b6  + a10*b5  + a11*b4;
        s16 =                                                    a5*b11 + a6*b10 + a7*b9  + a8*b8  + a9*b7  + a10*b6  + a11*b5;
        s17 =                                                             a6*b11 + a7*b10 + a8*b9  + a9*b8  + a10*b7  + a11*b6;
        s18 =                                                                      a7*b11 + a8*b10 + a9*b9  + a10*b8  + a11*b7;
        s19 =                                                                               a8*b11 + a9*b10 + a10*b9  + a11*b8;
        s20 =                                                                                        a9*b11 + a10*b10 + a11*b9;
        s21 =                                                                                                 a10*b11 + a11*b10;
        s22 =                                                                                                           a11*b11;

        carry0  = (s0  + (1<<20)) >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry2  = (s2  + (1<<20)) >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry4  = (s4  + (1<<20)) >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry12 = (s12 + (1<<20)) >> 21; s13 += carry12; s12 -= carry12 << 21;
        carry14 = (s14 + (1<<20)) >> 21; s15 += carry14; s14 -= carry14 << 21;
        carry16 = (s16 + (1<<20)) >> 21; s17 += carry16; s16 -= carry16 << 21;
        carry18 = (s18 + (1<<20)) >> 21; s19 += carry18; s18 -= carry18 << 21;
        carry20 = (s20 + (1<<20)) >> 21; s21 += carry20; s20 -= carry20 << 21;
        carry22 = (s22 + (1<<20)) >> 21; s23  = carry22; s22 -= carry22 << 21;

        carry1  = (s1  + (1<<20)) >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry3  = (s3  + (1<<20)) >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry5  = (s5  + (1<<20)) >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12 += carry11; s11 -= carry11 << 21;
        carry13 = (s13 + (1<<20)) >> 21; s14 += carry13; s13 -= carry13 << 21;
        carry15 = (s15 + (1<<20)) >> 21; s16 += carry15; s15 -= carry15 << 21;
        carry17 = (s17 + (1<<20)) >> 21; s18 += carry17; s17 -= carry17 << 21;
        carry19 = (s19 + (1<<20)) >> 21; s20 += carry19; s19 -= carry19 << 21;
        carry21 = (s21 + (1<<20)) >> 21; s22 += carry21; s21 -= carry21 << 21;
        // @formatter:on

        s11 += s23 * 666643;
        s12 += s23 * 470296;
        s13 += s23 * 654183;
        s14 -= s23 * 997805;
        s15 += s23 * 136657;
        s16 -= s23 * 683901;

        s10 += s22 * 666643;
        s11 += s22 * 470296;
        s12 += s22 * 654183;
        s13 -= s22 * 997805;
        s14 += s22 * 136657;
        s15 -= s22 * 683901;

        s9 += s21 * 666643;
        s10 += s21 * 470296;
        s11 += s21 * 654183;
        s12 -= s21 * 997805;
        s13 += s21 * 136657;
        s14 -= s21 * 683901;

        s8 += s20 * 666643;
        s9 += s20 * 470296;
        s10 += s20 * 654183;
        s11 -= s20 * 997805;
        s12 += s20 * 136657;
        s13 -= s20 * 683901;

        s7 += s19 * 666643;
        s8 += s19 * 470296;
        s9 += s19 * 654183;
        s10 -= s19 * 997805;
        s11 += s19 * 136657;
        s12 -= s19 * 683901;

        s6 += s18 * 666643;
        s7 += s18 * 470296;
        s8 += s18 * 654183;
        s9 -= s18 * 997805;
        s10 += s18 * 136657;
        s11 -= s18 * 683901;

        // @formatter:off
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry12 = (s12 + (1<<20)) >> 21; s13 += carry12; s12 -= carry12 << 21;
        carry14 = (s14 + (1<<20)) >> 21; s15 += carry14; s14 -= carry14 << 21;
        carry16 = (s16 + (1<<20)) >> 21; s17 += carry16; s16 -= carry16 << 21;

        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12 += carry11; s11 -= carry11 << 21;
        carry13 = (s13 + (1<<20)) >> 21; s14 += carry13; s13 -= carry13 << 21;
        carry15 = (s15 + (1<<20)) >> 21; s16 += carry15; s15 -= carry15 << 21;
        // @formatter:on

        s5 += s17 * 666643;
        s6 += s17 * 470296;
        s7 += s17 * 654183;
        s8 -= s17 * 997805;
        s9 += s17 * 136657;
        s10 -= s17 * 683901;

        s4 += s16 * 666643;
        s5 += s16 * 470296;
        s6 += s16 * 654183;
        s7 -= s16 * 997805;
        s8 += s16 * 136657;
        s9 -= s16 * 683901;

        s3 += s15 * 666643;
        s4 += s15 * 470296;
        s5 += s15 * 654183;
        s6 -= s15 * 997805;
        s7 += s15 * 136657;
        s8 -= s15 * 683901;

        s2 += s14 * 666643;
        s3 += s14 * 470296;
        s4 += s14 * 654183;
        s5 -= s14 * 997805;
        s6 += s14 * 136657;
        s7 -= s14 * 683901;

        s1 += s13 * 666643;
        s2 += s13 * 470296;
        s3 += s13 * 654183;
        s4 -= s13 * 997805;
        s5 += s13 * 136657;
        s6 -= s13 * 683901;

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = (s0  + (1<<20)) >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry2  = (s2  + (1<<20)) >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry4  = (s4  + (1<<20)) >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;

        carry1  = (s1  + (1<<20)) >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry3  = (s3  + (1<<20)) >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry5  = (s5  + (1<<20)) >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry11 = s11 >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        // @formatter:on

        byte[] result = new byte[32];
        result[0] = (byte) s0;
        result[1] = (byte) (s0 >> 8);
        result[2] = (byte) ((s0 >> 16) | (s1 << 5));
        result[3] = (byte) (s1 >> 3);
        result[4] = (byte) (s1 >> 11);
        result[5] = (byte) ((s1 >> 19) | (s2 << 2));
        result[6] = (byte) (s2 >> 6);
        result[7] = (byte) ((s2 >> 14) | (s3 << 7));
        result[8] = (byte) (s3 >> 1);
        result[9] = (byte) (s3 >> 9);
        result[10] = (byte) ((s3 >> 17) | (s4 << 4));
        result[11] = (byte) (s4 >> 4);
        result[12] = (byte) (s4 >> 12);
        result[13] = (byte) ((s4 >> 20) | (s5 << 1));
        result[14] = (byte) (s5 >> 7);
        result[15] = (byte) ((s5 >> 15) | (s6 << 6));
        result[16] = (byte) (s6 >> 2);
        result[17] = (byte) (s6 >> 10);
        result[18] = (byte) ((s6 >> 18) | (s7 << 3));
        result[19] = (byte) (s7 >> 5);
        result[20] = (byte) (s7 >> 13);
        result[21] = (byte) s8;
        result[22] = (byte) (s8 >> 8);
        result[23] = (byte) ((s8 >> 16) | (s9 << 5));
        result[24] = (byte) (s9 >> 3);
        result[25] = (byte) (s9 >> 11);
        result[26] = (byte) ((s9 >> 19) | (s10 << 2));
        result[27] = (byte) (s10 >> 6);
        result[28] = (byte) ((s10 >> 14) | (s11 << 7));
        result[29] = (byte) (s11 >> 1);
        result[30] = (byte) (s11 >> 9);
        result[31] = (byte) (s11 >> 17);
        return new Scalar(result);
    }

    public Scalar square() {
        long a0 = 0x1FFFFF & FieldElement.load_3(this.s, 0);
        long a1 = 0x1FFFFF & (FieldElement.load_4(this.s, 2) >> 5);
        long a2 = 0x1FFFFF & (FieldElement.load_3(this.s, 5) >> 2);
        long a3 = 0x1FFFFF & (FieldElement.load_4(this.s, 7) >> 7);
        long a4 = 0x1FFFFF & (FieldElement.load_4(this.s, 10) >> 4);
        long a5 = 0x1FFFFF & (FieldElement.load_3(this.s, 13) >> 1);
        long a6 = 0x1FFFFF & (FieldElement.load_4(this.s, 15) >> 6);
        long a7 = 0x1FFFFF & (FieldElement.load_3(this.s, 18) >> 3);
        long a8 = 0x1FFFFF & FieldElement.load_3(this.s, 21);
        long a9 = 0x1FFFFF & (FieldElement.load_4(this.s, 23) >> 5);
        long a10 = 0x1FFFFF & (FieldElement.load_3(this.s, 26) >> 2);
        long a11 = (FieldElement.load_4(this.s, 28) >> 7);
        long s0;
        long s1;
        long s2;
        long s3;
        long s4;
        long s5;
        long s6;
        long s7;
        long s8;
        long s9;
        long s10;
        long s11;
        long s12;
        long s13;
        long s14;
        long s15;
        long s16;
        long s17;
        long s18;
        long s19;
        long s20;
        long s21;
        long s22;
        long s23;
        long carry0;
        long carry1;
        long carry2;
        long carry3;
        long carry4;
        long carry5;
        long carry6;
        long carry7;
        long carry8;
        long carry9;
        long carry10;
        long carry11;
        long carry12;
        long carry13;
        long carry14;
        long carry15;
        long carry16;
        long carry17;
        long carry18;
        long carry19;
        long carry20;
        long carry21;
        long carry22;

        // @formatter:off
        s0  = a0*a0;
        s1  = a0*a1  + a1*a0;
        s2  = a0*a2  + a1*a1  + a2*a0;
        s3  = a0*a3  + a1*a2  + a2*a1  + a3*a0;
        s4  = a0*a4  + a1*a3  + a2*a2  + a3*a1  + a4*a0;
        s5  = a0*a5  + a1*a4  + a2*a3  + a3*a2  + a4*a1  + a5*a0;
        s6  = a0*a6  + a1*a5  + a2*a4  + a3*a3  + a4*a2  + a5*a1  + a6*a0;
        s7  = a0*a7  + a1*a6  + a2*a5  + a3*a4  + a4*a3  + a5*a2  + a6*a1  + a7*a0;
        s8  = a0*a8  + a1*a7  + a2*a6  + a3*a5  + a4*a4  + a5*a3  + a6*a2  + a7*a1  + a8*a0;
        s9  = a0*a9  + a1*a8  + a2*a7  + a3*a6  + a4*a5  + a5*a4  + a6*a3  + a7*a2  + a8*a1  + a9*a0;
        s10 = a0*a10 + a1*a9  + a2*a8  + a3*a7  + a4*a6  + a5*a5  + a6*a4  + a7*a3  + a8*a2  + a9*a1  + a10*a0;
        s11 = a0*a11 + a1*a10 + a2*a9  + a3*a8  + a4*a7  + a5*a6  + a6*a5  + a7*a4  + a8*a3  + a9*a2  + a10*a1  + a11*a0;
        s12 =          a1*a11 + a2*a10 + a3*a9  + a4*a8  + a5*a7  + a6*a6  + a7*a5  + a8*a4  + a9*a3  + a10*a2  + a11*a1;
        s13 =                   a2*a11 + a3*a10 + a4*a9  + a5*a8  + a6*a7  + a7*a6  + a8*a5  + a9*a4  + a10*a3  + a11*a2;
        s14 =                            a3*a11 + a4*a10 + a5*a9  + a6*a8  + a7*a7  + a8*a6  + a9*a5  + a10*a4  + a11*a3;
        s15 =                                     a4*a11 + a5*a10 + a6*a9  + a7*a8  + a8*a7  + a9*a6  + a10*a5  + a11*a4;
        s16 =                                              a5*a11 + a6*a10 + a7*a9  + a8*a8  + a9*a7  + a10*a6  + a11*a5;
        s17 =                                                       a6*a11 + a7*a10 + a8*a9  + a9*a8  + a10*a7  + a11*a6;
        s18 =                                                                a7*a11 + a8*a10 + a9*a9  + a10*a8  + a11*a7;
        s19 =                                                                         a8*a11 + a9*a10 + a10*a9  + a11*a8;
        s20 =                                                                                  a9*a11 + a10*a10 + a11*a9;
        s21 =                                                                                          a10*a11 + a11*a10;
        s22 =                                                                                                    a11*a11;

        carry0  = (s0  + (1<<20)) >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry2  = (s2  + (1<<20)) >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry4  = (s4  + (1<<20)) >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry12 = (s12 + (1<<20)) >> 21; s13 += carry12; s12 -= carry12 << 21;
        carry14 = (s14 + (1<<20)) >> 21; s15 += carry14; s14 -= carry14 << 21;
        carry16 = (s16 + (1<<20)) >> 21; s17 += carry16; s16 -= carry16 << 21;
        carry18 = (s18 + (1<<20)) >> 21; s19 += carry18; s18 -= carry18 << 21;
        carry20 = (s20 + (1<<20)) >> 21; s21 += carry20; s20 -= carry20 << 21;
        carry22 = (s22 + (1<<20)) >> 21; s23  = carry22; s22 -= carry22 << 21;

        carry1  = (s1  + (1<<20)) >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry3  = (s3  + (1<<20)) >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry5  = (s5  + (1<<20)) >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12 += carry11; s11 -= carry11 << 21;
        carry13 = (s13 + (1<<20)) >> 21; s14 += carry13; s13 -= carry13 << 21;
        carry15 = (s15 + (1<<20)) >> 21; s16 += carry15; s15 -= carry15 << 21;
        carry17 = (s17 + (1<<20)) >> 21; s18 += carry17; s17 -= carry17 << 21;
        carry19 = (s19 + (1<<20)) >> 21; s20 += carry19; s19 -= carry19 << 21;
        carry21 = (s21 + (1<<20)) >> 21; s22 += carry21; s21 -= carry21 << 21;
        // @formatter:on

        s11 += s23 * 666643;
        s12 += s23 * 470296;
        s13 += s23 * 654183;
        s14 -= s23 * 997805;
        s15 += s23 * 136657;
        s16 -= s23 * 683901;

        s10 += s22 * 666643;
        s11 += s22 * 470296;
        s12 += s22 * 654183;
        s13 -= s22 * 997805;
        s14 += s22 * 136657;
        s15 -= s22 * 683901;

        s9 += s21 * 666643;
        s10 += s21 * 470296;
        s11 += s21 * 654183;
        s12 -= s21 * 997805;
        s13 += s21 * 136657;
        s14 -= s21 * 683901;

        s8 += s20 * 666643;
        s9 += s20 * 470296;
        s10 += s20 * 654183;
        s11 -= s20 * 997805;
        s12 += s20 * 136657;
        s13 -= s20 * 683901;

        s7 += s19 * 666643;
        s8 += s19 * 470296;
        s9 += s19 * 654183;
        s10 -= s19 * 997805;
        s11 += s19 * 136657;
        s12 -= s19 * 683901;

        s6 += s18 * 666643;
        s7 += s18 * 470296;
        s8 += s18 * 654183;
        s9 -= s18 * 997805;
        s10 += s18 * 136657;
        s11 -= s18 * 683901;

        // @formatter:off
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry12 = (s12 + (1<<20)) >> 21; s13 += carry12; s12 -= carry12 << 21;
        carry14 = (s14 + (1<<20)) >> 21; s15 += carry14; s14 -= carry14 << 21;
        carry16 = (s16 + (1<<20)) >> 21; s17 += carry16; s16 -= carry16 << 21;

        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12 += carry11; s11 -= carry11 << 21;
        carry13 = (s13 + (1<<20)) >> 21; s14 += carry13; s13 -= carry13 << 21;
        carry15 = (s15 + (1<<20)) >> 21; s16 += carry15; s15 -= carry15 << 21;
        // @formatter:on

        s5 += s17 * 666643;
        s6 += s17 * 470296;
        s7 += s17 * 654183;
        s8 -= s17 * 997805;
        s9 += s17 * 136657;
        s10 -= s17 * 683901;

        s4 += s16 * 666643;
        s5 += s16 * 470296;
        s6 += s16 * 654183;
        s7 -= s16 * 997805;
        s8 += s16 * 136657;
        s9 -= s16 * 683901;

        s3 += s15 * 666643;
        s4 += s15 * 470296;
        s5 += s15 * 654183;
        s6 -= s15 * 997805;
        s7 += s15 * 136657;
        s8 -= s15 * 683901;

        s2 += s14 * 666643;
        s3 += s14 * 470296;
        s4 += s14 * 654183;
        s5 -= s14 * 997805;
        s6 += s14 * 136657;
        s7 -= s14 * 683901;

        s1 += s13 * 666643;
        s2 += s13 * 470296;
        s3 += s13 * 654183;
        s4 -= s13 * 997805;
        s5 += s13 * 136657;
        s6 -= s13 * 683901;

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = (s0  + (1<<20)) >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry2  = (s2  + (1<<20)) >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry4  = (s4  + (1<<20)) >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry6  = (s6  + (1<<20)) >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry8  = (s8  + (1<<20)) >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry10 = (s10 + (1<<20)) >> 21; s11 += carry10; s10 -= carry10 << 21;

        carry1  = (s1  + (1<<20)) >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry3  = (s3  + (1<<20)) >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry5  = (s5  + (1<<20)) >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry7  = (s7  + (1<<20)) >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry9  = (s9  + (1<<20)) >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry11 = (s11 + (1<<20)) >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        carry11 = s11 >> 21; s12  = carry11; s11 -= carry11 << 21;
        // @formatter:on

        s0 += s12 * 666643;
        s1 += s12 * 470296;
        s2 += s12 * 654183;
        s3 -= s12 * 997805;
        s4 += s12 * 136657;
        s5 -= s12 * 683901;

        // @formatter:off
        carry0  = s0  >> 21; s1  += carry0;  s0  -= carry0  << 21;
        carry1  = s1  >> 21; s2  += carry1;  s1  -= carry1  << 21;
        carry2  = s2  >> 21; s3  += carry2;  s2  -= carry2  << 21;
        carry3  = s3  >> 21; s4  += carry3;  s3  -= carry3  << 21;
        carry4  = s4  >> 21; s5  += carry4;  s4  -= carry4  << 21;
        carry5  = s5  >> 21; s6  += carry5;  s5  -= carry5  << 21;
        carry6  = s6  >> 21; s7  += carry6;  s6  -= carry6  << 21;
        carry7  = s7  >> 21; s8  += carry7;  s7  -= carry7  << 21;
        carry8  = s8  >> 21; s9  += carry8;  s8  -= carry8  << 21;
        carry9  = s9  >> 21; s10 += carry9;  s9  -= carry9  << 21;
        carry10 = s10 >> 21; s11 += carry10; s10 -= carry10 << 21;
        // @formatter:on

        byte[] result = new byte[32];
        result[0] = (byte) s0;
        result[1] = (byte) (s0 >> 8);
        result[2] = (byte) ((s0 >> 16) | (s1 << 5));
        result[3] = (byte) (s1 >> 3);
        result[4] = (byte) (s1 >> 11);
        result[5] = (byte) ((s1 >> 19) | (s2 << 2));
        result[6] = (byte) (s2 >> 6);
        result[7] = (byte) ((s2 >> 14) | (s3 << 7));
        result[8] = (byte) (s3 >> 1);
        result[9] = (byte) (s3 >> 9);
        result[10] = (byte) ((s3 >> 17) | (s4 << 4));
        result[11] = (byte) (s4 >> 4);
        result[12] = (byte) (s4 >> 12);
        result[13] = (byte) ((s4 >> 20) | (s5 << 1));
        result[14] = (byte) (s5 >> 7);
        result[15] = (byte) ((s5 >> 15) | (s6 << 6));
        result[16] = (byte) (s6 >> 2);
        result[17] = (byte) (s6 >> 10);
        result[18] = (byte) ((s6 >> 18) | (s7 << 3));
        result[19] = (byte) (s7 >> 5);
        result[20] = (byte) (s7 >> 13);
        result[21] = (byte) s8;
        result[22] = (byte) (s8 >> 8);
        result[23] = (byte) ((s8 >> 16) | (s9 << 5));
        result[24] = (byte) (s9 >> 3);
        result[25] = (byte) (s9 >> 11);
        result[26] = (byte) ((s9 >> 19) | (s10 << 2));
        result[27] = (byte) (s10 >> 6);
        result[28] = (byte) ((s10 >> 14) | (s11 << 7));
        result[29] = (byte) (s11 >> 1);
        result[30] = (byte) (s11 >> 9);
        result[31] = (byte) (s11 >> 17);
        return new Scalar(result);
    }

    public Scalar invert() {
        Scalar res = Scalar.ONE;
        for(int i = 255; i >= 0; i--){
            res = res.square();
            if (Constants.ModMinus2.testBit(i)) {
                res = this.multiplyAndAdd(res, Scalar.ZERO);
            }
        }
        return res;
    }

    public Scalar divide(Scalar a) {
        Scalar inv_a = a.invert();
        return this.multiplyAndAdd(inv_a, Scalar.ZERO);
    }

    /**
     * Writes this Scalar in radix 16, with coefficients in range $[-8, 8)$.
     *
     * @return 64 bytes, each between -8 and 7
     */
    byte[] toRadix16() {
        final byte[] e = new byte[64];
        int i;
        // Radix 16 notation
        for (i = 0; i < 32; i++) {
            e[2 * i + 0] = (byte) (s[i] & 15);
            e[2 * i + 1] = (byte) ((s[i] >> 4) & 15);
        }
        /* each e[i] is between 0 and 15 */
        /* e[63] is between 0 and 7 */
        int carry = 0;
        for (i = 0; i < 63; i++) {
            e[i] += carry;
            carry = e[i] + 8;
            carry >>= 4;
            e[i] -= carry << 4;
        }
        e[63] += carry;
        /* each e[i] is between -8 and 7 */
        return e;
    }

    public byte[] toRadix2w(int w) {
        long[] scalar64x4 = new long[4];
        ByteBuffer.wrap(s).order(ByteOrder.LITTLE_ENDIAN).asLongBuffer().get(scalar64x4);

        long radix = 1L << w;
        long windowMask = radix - 1;

        long carry = 0L;
        byte[] digits = new byte[43];
        int digitsCount = (256 + w - 1) / w;
        for (int i = 0; i < digitsCount; i++) {
            int bitOffset = i * w;
            int u64Idx = bitOffset / 64;
            int bitIdx = bitOffset % 64;

            long bitBuf;
            if (bitIdx < 64 - w || u64Idx == 3) {
                bitBuf = scalar64x4[u64Idx] >>> bitIdx;
            } else {
                bitBuf = (scalar64x4[u64Idx] >>> bitIdx) | (scalar64x4[1 + u64Idx] << (64 - bitIdx));
            }

            long coef = carry + (bitBuf & windowMask);

            carry = (coef + (radix / 2)) >>> w;
            digits[i] = (byte) ((coef - (carry << w)));
        }

        if (w == 8) {
            digits[digitsCount] += carry;
        } else {
            digits[digitsCount - 1] += (carry << w);
        }

        return digits;
    }

    /**
     * Compute a width-$w$ "Non-Adjacent Form" of this scalar.
     *
     * @return The byte array $naf$ in the above described form.
     */
    byte[] nonAdjacentForm() {
        byte[] naf = new byte[256];

        // Put each bit of this into a separate byte, 0 or 1
        for (int i = 0; i < 256; ++i) {
            naf[i] = (byte) (1 & (this.s[i >> 3] >> (i & 7)));
        }

        // Note: naf[i] will always be odd.
        for (int i = 0; i < 256; ++i) {
            if (naf[i] != 0) {
                for (int b = 1; b <= 6 && i + b < 256; ++b) {
                    // Accumulate bits if possible
                    if (naf[i + b] != 0) {
                        if (naf[i] + (naf[i + b] << b) <= 15) {
                            naf[i] += naf[i + b] << b;
                            naf[i + b] = 0;
                        } else if (naf[i] - (naf[i + b] << b) >= -15) {
                            naf[i] -= naf[i + b] << b;
                            for (int k = i + b; k < 256; ++k) {
                                if (naf[k] == 0) {
                                    naf[k] = 1;
                                    break;
                                }
                                naf[k] = 0;
                            }
                        } else
                            break;
                    }
                }
            }
        }

        return naf;
    }
}
