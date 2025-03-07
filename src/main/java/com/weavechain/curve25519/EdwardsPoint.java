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

/**
 * An EdwardsPoint represents a point on the Edwards form of Curve25519.
 */
public class EdwardsPoint implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final EdwardsPoint IDENTITY = new EdwardsPoint(FieldElement.ZERO, FieldElement.ONE, FieldElement.ONE,
            FieldElement.ZERO);

    transient FieldElement X;
    transient FieldElement Y;
    transient FieldElement Z;
    transient FieldElement T;

    /**
     * Only for internal use.
     */
    EdwardsPoint(FieldElement X, FieldElement Y, FieldElement Z, FieldElement T) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.T = T;
    }

    public EdwardsPoint copy() {
        return new EdwardsPoint(X, Y, Z, Y);
    }

    /**
     * Overrides class serialization to use the canonical encoded format.
     */
    private void writeObject(ObjectOutputStream out) throws IOException {
        out.write(this.compress().toByteArray());
    }

    /**
     * Overrides class serialization to use the canonical encoded format.
     */
    private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
        byte[] encoded = new byte[32];
        in.readFully(encoded);

        try {
            EdwardsPoint point = new CompressedEdwardsY(encoded).decompress();
            this.X = point.X;
            this.Y = point.Y;
            this.Z = point.Z;
            this.T = point.T;
        } catch (InvalidEncodingException iee) {
            throw new InvalidObjectException(iee.getMessage());
        }
    }

    @SuppressWarnings("unused")
    private void readObjectNoData() throws ObjectStreamException {
        throw new InvalidObjectException("Cannot deserialize EdwardsPoint from no data");
    }

    /**
     * Compress this point to CompressedEdwardsY format.
     *
     * @return the encoded point.
     */
    public CompressedEdwardsY compress() {
        FieldElement recip = this.Z.invert();
        FieldElement x = this.X.multiply(recip);
        FieldElement y = this.Y.multiply(recip);
        byte[] s = y.toByteArray();
        s[31] |= (x.isNegative() << 7);
        return new CompressedEdwardsY(s);
    }

    /**
     * Constant-time equality check.
     * <p>
     * Compares the encodings of the two EdwardsPoints.
     *
     * @return 1 if this and other are equal, 0 otherwise.
     */
    public int ctEquals(EdwardsPoint other) {
        return compress().ctEquals(other.compress());
    }

    /**
     * Constant-time selection between two EdwardsPoints.
     *
     * @param that the other point.
     * @param b    must be 0 or 1, otherwise results are undefined.
     * @return a copy of this if $b == 0$, or a copy of that if $b == 1$.
     */
    public EdwardsPoint ctSelect(EdwardsPoint that, int b) {
        return new EdwardsPoint(this.X.ctSelect(that.X, b), this.Y.ctSelect(that.Y, b), this.Z.ctSelect(that.Z, b),
                this.T.ctSelect(that.T, b));
    }

    /**
     * Equality check overridden to be constant-time.
     * <p>
     * Fails fast if the objects are of different types.
     *
     * @return true if this and other are equal, false otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof EdwardsPoint)) {
            return false;
        }

        EdwardsPoint other = (EdwardsPoint) obj;
        return ctEquals(other) == 1;
    }

    @Override
    public int hashCode() {
        // The general contract for the hashCode method states that equal objects must
        // have equal hash codes. Object equality is based on the encodings of the
        // points, not their internal representations (which may not be canonical).
        return compress().hashCode();
    }

    /**
     * Convert the representation of this point from extended coordinates to
     * projective coordinates.
     * <p>
     * Free.
     */
    ProjectivePoint toProjective() {
        return new ProjectivePoint(this.X, this.Y, this.Z);
    }

    /**
     * Convert to a ProjectiveNielsPoint.
     */
    ProjectiveNielsPoint toProjectiveNiels() {
        return new ProjectiveNielsPoint(this.Y.add(this.X), this.Y.subtract(this.X), this.Z,
                this.T.multiply(Constants.EDWARDS_2D));
    }

    /**
     * Dehomogenize to an AffineNielsPoint.
     */
    AffineNielsPoint toAffineNiels() {
        FieldElement recip = this.Z.invert();
        FieldElement x = this.X.multiply(recip);
        FieldElement y = this.Y.multiply(recip);
        FieldElement xy2D = x.multiply(y).multiply(Constants.EDWARDS_2D);
        return new AffineNielsPoint(y.add(x), y.subtract(x), xy2D);
    }

    /**
     * Point addition.
     *
     * @param Q the point to add to this one.
     * @return $P + Q$
     */
    public EdwardsPoint add(EdwardsPoint Q) {
        return this.add(Q.toProjectiveNiels()).toExtended();
    }

    /**
     * Point addition.
     *
     * @param Q the point to add to this one, in projective "Niels coordinates".
     * @return $P + Q$
     */
    CompletedPoint add(ProjectiveNielsPoint Q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PP = YPlusX.multiply(Q.YPlusX);
        FieldElement MM = YMinusX.multiply(Q.YMinusX);
        FieldElement TT2D = this.T.multiply(Q.T2D);
        FieldElement ZZ = this.Z.multiply(Q.Z);
        FieldElement ZZ2 = ZZ.add(ZZ);
        return new CompletedPoint(PP.subtract(MM), PP.add(MM), ZZ2.add(TT2D), ZZ2.subtract(TT2D));
    }

    EdwardsPoint addToExtended(ProjectiveNielsPoint Q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PP = YPlusX.multiply(Q.YPlusX);
        FieldElement MM = YMinusX.multiply(Q.YMinusX);
        FieldElement TT2D = this.T.multiply(Q.T2D);
        FieldElement ZZ = this.Z.multiply(Q.Z);
        FieldElement ZZ2 = ZZ.add(ZZ);

        FieldElement X = PP.subtract(MM);
        FieldElement Y = PP.add(MM);
        FieldElement Z = ZZ2.add(TT2D);
        FieldElement T = ZZ2.subtract(TT2D);

        return new EdwardsPoint(
                X.multiply(T),
                Y.multiply(Z),
                Z.multiply(T),
                X.multiply(Y)
        );
    }

    /**
     * Point addition.
     *
     * @param q the point to add to this one, in affine "Niels coordinates".
     * @return $P + q$
     */
    CompletedPoint add(AffineNielsPoint q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PP = YPlusX.multiply(q.yPlusx);
        FieldElement MM = YMinusX.multiply(q.yMinusx);
        FieldElement Txy2D = this.T.multiply(q.xy2D);
        FieldElement Z2 = this.Z.add(this.Z);
        return new CompletedPoint(PP.subtract(MM), PP.add(MM), Z2.add(Txy2D), Z2.subtract(Txy2D));
    }

    EdwardsPoint addToExtended(AffineNielsPoint q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PP = YPlusX.multiply(q.yPlusx);
        FieldElement MM = YMinusX.multiply(q.yMinusx);
        FieldElement Txy2D = this.T.multiply(q.xy2D);
        FieldElement Z2 = this.Z.add(this.Z);

        FieldElement X = PP.subtract(MM);
        FieldElement Y = PP.add(MM);
        FieldElement Z = Z2.add(Txy2D);
        FieldElement T = Z2.subtract(Txy2D);

        return new EdwardsPoint(
                X.multiply(T),
                Y.multiply(Z),
                Z.multiply(T),
                X.multiply(Y)
        );
    }

    /**
     * Point subtraction.
     *
     * @param Q the point to subtract from this one.
     * @return $P - Q$
     */
    public EdwardsPoint subtract(EdwardsPoint Q) {
        return this.subtract(Q.toProjectiveNiels()).toExtended();
    }

    /**
     * Point subtraction.
     *
     * @param Q the point to subtract from this one, in projective "Niels
     *          coordinates".
     * @return $P - Q$
     */
    CompletedPoint subtract(ProjectiveNielsPoint Q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PM = YPlusX.multiply(Q.YMinusX);
        FieldElement MP = YMinusX.multiply(Q.YPlusX);
        FieldElement TT2D = this.T.multiply(Q.T2D);
        FieldElement ZZ = Z.multiply(Q.Z);
        FieldElement ZZ2 = ZZ.add(ZZ);
        return new CompletedPoint(PM.subtract(MP), PM.add(MP), ZZ2.subtract(TT2D), ZZ2.add(TT2D));
    }

    /**
     * Point subtraction.
     *
     * @param q the point to subtract from this one, in affine "Niels coordinates".
     * @return $P - q$
     */
    CompletedPoint subtract(AffineNielsPoint q) {
        FieldElement YPlusX = this.Y.add(this.X);
        FieldElement YMinusX = this.Y.subtract(this.X);
        FieldElement PM = YPlusX.multiply(q.yMinusx);
        FieldElement MP = YMinusX.multiply(q.yPlusx);
        FieldElement Txy2D = this.T.multiply(q.xy2D);
        FieldElement Z2 = this.Z.add(this.Z);
        return new CompletedPoint(PM.subtract(MP), PM.add(MP), Z2.subtract(Txy2D), Z2.add(Txy2D));
    }

    /**
     * Point negation.
     *
     * @return $-P$
     */
    public EdwardsPoint negate() {
        return new EdwardsPoint(this.X.negate(), this.Y, this.Z, this.T.negate());
    }

    /**
     * Point doubling.
     *
     * @return $[2]P$
     */
    public EdwardsPoint dbl() {
        return this.toProjective().dbl().toExtended();
    }

    public EdwardsPoint dblInPlace() {
        FieldElement XX = this.X.square();
        FieldElement YY = this.Y.square();
        FieldElement ZZ2 = this.Z.squareAndDouble();
        FieldElement XPlusY = this.X.add(this.Y);
        FieldElement XPlusYSq = XPlusY.square();
        FieldElement YYPlusXX = YY.add(XX);
        FieldElement YYMinusXX = YY.subtract(XX);
        FieldElement X = XPlusYSq.subtract(YYPlusXX);
        FieldElement Y = YYPlusXX;
        FieldElement Z = YYMinusXX;
        FieldElement T = ZZ2.subtract(YYMinusXX);
        this.X = X.multiply(T);
        this.Y = Y.multiply(Z);
        this.Z = Z.multiply(T);
        this.T = X.multiply(Y);
        return this;
    }

    /**
     * Constant-time variable-base scalar multiplication.
     *
     * @param s the Scalar to multiply by.
     * @return $[s]P$
     */
    public EdwardsPoint multiply(final Scalar s) {
        // Construct a lookup table of [P,2P,3P,4P,5P,6P,7P,8P]
        final ProjectiveNielsPoint.LookupTable lookupTable = ProjectiveNielsPoint.buildLookupTable(this);

        // Compute
        //
        // s = s_0 + s_1*16^1 + ... + s_63*16^63,
        //
        // with -8 ≤ s_i < 8 for 0 ≤ i < 63 and -8 ≤ s_63 ≤ 8.
        final byte[] e = s.toRadix16();

        // Compute s*P as
        //
        // @formatter:off
        //    s*P = P*(s_0 +   s_1*16^1 +   s_2*16^2 + ... +   s_63*16^63)
        //    s*P =  P*s_0 + P*s_1*16^1 + P*s_2*16^2 + ... + P*s_63*16^63
        //    s*P = P*s_0 + 16*(P*s_1 + 16*(P*s_2 + 16*( ... + P*s_63)...))
        // @formatter:on
        //
        // We sum right-to-left.
        EdwardsPoint Q = EdwardsPoint.IDENTITY;
        for (int i = 63; i >= 0; i--) {
            Q = Q.multiplyByPow2InPlace(4);
            Q = Q.add(lookupTable.select(e[i])).toExtended();
        }
        return Q;
    }

    /**
     * Compute $r = [a]A + [b]B$ in variable time, where $B$ is the Ed25519
     * basepoint.
     *
     * @param a a Scalar.
     * @param A an EdwardsPoint.
     * @param b a Scalar.
     * @return $[a]A + [b]B$
     */
    public static EdwardsPoint vartimeDoubleScalarMultiplyBasepoint(final Scalar a, final EdwardsPoint A,
            final Scalar b) {
        final byte[] aNaf = a.nonAdjacentForm();
        final byte[] bNaf = b.nonAdjacentForm();

        ProjectiveNielsPoint.NafLookupTable tableA = ProjectiveNielsPoint.buildNafLookupTable(A);
        AffineNielsPoint.NafLookupTable tableB = Constants.AFFINE_ODD_MULTIPLES_OF_BASEPOINT;

        int i;
        for (i = 255; i >= 0; --i) {
            if (aNaf[i] != 0 || bNaf[i] != 0)
                break;
        }

        ProjectivePoint r = EdwardsPoint.IDENTITY.toProjective();
        for (; i >= 0; --i) {
            CompletedPoint t = r.dbl();

            if (aNaf[i] > 0) {
                t = t.toExtended().add(tableA.select(aNaf[i]));
            } else if (aNaf[i] < 0) {
                t = t.toExtended().subtract(tableA.select(-aNaf[i]));
            }

            if (bNaf[i] > 0) {
                t = t.toExtended().add(tableB.select(bNaf[i]));
            } else if (bNaf[i] < 0) {
                t = t.toExtended().subtract(tableB.select(-bNaf[i]));
            }

            r = t.toProjective();
        }

        return r.toExtended();
    }

    /**
     * Multiply by the cofactor.
     *
     * @return $[8]P$
     */
    public EdwardsPoint multiplyByCofactor() {
        return this.multiplyByPow2(3);
    }

    /**
     * Compute $[2^k]P$ by successive doublings.
     *
     * @param k the exponent of 2. Must be positive and non-zero.
     * @return $[2^k]P$
     */
    EdwardsPoint multiplyByPow2(int k) {
        if (!(k > 0)) {
            throw new IllegalArgumentException("Exponent must be positive and non-zero");
        }

        EdwardsPoint s = new EdwardsPoint(this.X, this.Y, this.Z, this.T);
        return s.multiplyByPow2InPlace(k);
    }

    EdwardsPoint multiplyByPow2InPlace(int k) {
        if (!(k > 0)) {
            throw new IllegalArgumentException("Exponent must be positive and non-zero");
        }

        if (k == 3) {
            FieldElement XX = this.X.square();
            FieldElement YY = this.Y.square();
            FieldElement ZZ2 = this.Z.squareAndDouble();
            FieldElement XPlusY = this.X.add(this.Y);
            FieldElement XPlusYSq = XPlusY.square();
            FieldElement YYPlusXX = YY.add(XX);
            FieldElement YYMinusXX = YY.subtract(XX);
            FieldElement X = XPlusYSq.subtract(YYPlusXX);
            FieldElement T = ZZ2.subtract(YYMinusXX);
            FieldElement X_2P = X.multiply(T);
            FieldElement Y_2P = YYPlusXX.multiply(YYMinusXX);
            FieldElement Z_2P = YYMinusXX.multiply(T);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ4 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ4.subtract(YYMinusXX2);
            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ8 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ8.subtract(YYMinusXX3);

            this.X = X3.multiply(T3);
            this.Y = YYPlusXX3.multiply(YYMinusXX3);
            this.Z = YYMinusXX3.multiply(T3);
            this.T = X3.multiply(YYPlusXX3);

            return this;
        } else if (k == 4) {
            FieldElement XX1 = this.X.square();
            FieldElement YY1 = this.Y.square();
            FieldElement ZZ2_1 = this.Z.squareAndDouble();
            FieldElement XPlusY1 = this.X.add(this.Y);
            FieldElement XPlusYSq1 = XPlusY1.square();
            FieldElement YYPlusXX1 = YY1.add(XX1);
            FieldElement YYMinusXX1 = YY1.subtract(XX1);
            FieldElement X1 = XPlusYSq1.subtract(YYPlusXX1);
            FieldElement T1 = ZZ2_1.subtract(YYMinusXX1);

            FieldElement X_2P = X1.multiply(T1);
            FieldElement Y_2P = YYPlusXX1.multiply(YYMinusXX1);
            FieldElement Z_2P = YYMinusXX1.multiply(T1);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ2_2 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ2_2.subtract(YYMinusXX2);

            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ2_3 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ2_3.subtract(YYMinusXX3);

            FieldElement X_8P = X3.multiply(T3);
            FieldElement Y_8P = YYPlusXX3.multiply(YYMinusXX3);
            FieldElement Z_8P = YYMinusXX3.multiply(T3);

            FieldElement XX4 = X_8P.square();
            FieldElement YY4 = Y_8P.square();
            FieldElement ZZ2_4 = Z_8P.squareAndDouble();
            FieldElement XPlusY4 = X_8P.add(Y_8P);
            FieldElement XPlusYSq4 = XPlusY4.square();
            FieldElement YYPlusXX4 = YY4.add(XX4);
            FieldElement YYMinusXX4 = YY4.subtract(XX4);
            FieldElement X4 = XPlusYSq4.subtract(YYPlusXX4);
            FieldElement T4 = ZZ2_4.subtract(YYMinusXX4);

            this.X = X4.multiply(T4);
            this.Y = YYPlusXX4.multiply(YYMinusXX4);
            this.Z = YYMinusXX4.multiply(T4);
            this.T = X4.multiply(YYPlusXX4);

            return this;
        } else if (k == 5) {
            FieldElement XX1 = this.X.square();
            FieldElement YY1 = this.Y.square();
            FieldElement ZZ2_1 = this.Z.squareAndDouble();
            FieldElement XPlusY1 = this.X.add(this.Y);
            FieldElement XPlusYSq1 = XPlusY1.square();
            FieldElement YYPlusXX1 = YY1.add(XX1);
            FieldElement YYMinusXX1 = YY1.subtract(XX1);
            FieldElement X1 = XPlusYSq1.subtract(YYPlusXX1);
            FieldElement T1 = ZZ2_1.subtract(YYMinusXX1);

            FieldElement X_2P = X1.multiply(T1);
            FieldElement Y_2P = YYPlusXX1.multiply(YYMinusXX1);
            FieldElement Z_2P = YYMinusXX1.multiply(T1);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ2_2 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ2_2.subtract(YYMinusXX2);

            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ2_3 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ2_3.subtract(YYMinusXX3);

            FieldElement X_8P = X3.multiply(T3);
            FieldElement Y_8P = YYPlusXX3.multiply(YYMinusXX3);
            FieldElement Z_8P = YYMinusXX3.multiply(T3);

            FieldElement XX4 = X_8P.square();
            FieldElement YY4 = Y_8P.square();
            FieldElement ZZ2_4 = Z_8P.squareAndDouble();
            FieldElement XPlusY4 = X_8P.add(Y_8P);
            FieldElement XPlusYSq4 = XPlusY4.square();
            FieldElement YYPlusXX4 = YY4.add(XX4);
            FieldElement YYMinusXX4 = YY4.subtract(XX4);
            FieldElement X4 = XPlusYSq4.subtract(YYPlusXX4);
            FieldElement T4 = ZZ2_4.subtract(YYMinusXX4);

            FieldElement X_16P = X4.multiply(T4);
            FieldElement Y_16P = YYPlusXX4.multiply(YYMinusXX4);
            FieldElement Z_16P = YYMinusXX4.multiply(T4);

            FieldElement XX5 = X_16P.square();
            FieldElement YY5 = Y_16P.square();
            FieldElement ZZ2_5 = Z_16P.squareAndDouble();
            FieldElement XPlusY5 = X_16P.add(Y_16P);
            FieldElement XPlusYSq5 = XPlusY5.square();
            FieldElement YYPlusXX5 = YY5.add(XX5);
            FieldElement YYMinusXX5 = YY5.subtract(XX5);
            FieldElement X5 = XPlusYSq5.subtract(YYPlusXX5);
            FieldElement T5 = ZZ2_5.subtract(YYMinusXX5);

            this.X = X5.multiply(T5);
            this.Y = YYPlusXX5.multiply(YYMinusXX5);
            this.Z = YYMinusXX5.multiply(T5);
            this.T = X5.multiply(YYPlusXX5);

            return this;
        } else if (k == 6) {
            FieldElement XX1 = this.X.square();
            FieldElement YY1 = this.Y.square();
            FieldElement ZZ2_1 = this.Z.squareAndDouble();
            FieldElement XPlusY1 = this.X.add(this.Y);
            FieldElement XPlusYSq1 = XPlusY1.square();
            FieldElement YYPlusXX1 = YY1.add(XX1);
            FieldElement YYMinusXX1 = YY1.subtract(XX1);
            FieldElement X1 = XPlusYSq1.subtract(YYPlusXX1);
            FieldElement T1 = ZZ2_1.subtract(YYMinusXX1);

            FieldElement X_2P = X1.multiply(T1);
            FieldElement Y_2P = YYPlusXX1.multiply(YYMinusXX1);
            FieldElement Z_2P = YYMinusXX1.multiply(T1);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ2_2 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ2_2.subtract(YYMinusXX2);

            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ2_3 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ2_3.subtract(YYMinusXX3);

            FieldElement X_8P = X3.multiply(T3);
            FieldElement Y_8P = YYPlusXX3.multiply(YYMinusXX3);
            FieldElement Z_8P = YYMinusXX3.multiply(T3);

            FieldElement XX4 = X_8P.square();
            FieldElement YY4 = Y_8P.square();
            FieldElement ZZ2_4 = Z_8P.squareAndDouble();
            FieldElement XPlusY4 = X_8P.add(Y_8P);
            FieldElement XPlusYSq4 = XPlusY4.square();
            FieldElement YYPlusXX4 = YY4.add(XX4);
            FieldElement YYMinusXX4 = YY4.subtract(XX4);
            FieldElement X4 = XPlusYSq4.subtract(YYPlusXX4);
            FieldElement T4 = ZZ2_4.subtract(YYMinusXX4);

            FieldElement X_16P = X4.multiply(T4);
            FieldElement Y_16P = YYPlusXX4.multiply(YYMinusXX4);
            FieldElement Z_16P = YYMinusXX4.multiply(T4);

            FieldElement XX5 = X_16P.square();
            FieldElement YY5 = Y_16P.square();
            FieldElement ZZ2_5 = Z_16P.squareAndDouble();
            FieldElement XPlusY5 = X_16P.add(Y_16P);
            FieldElement XPlusYSq5 = XPlusY5.square();
            FieldElement YYPlusXX5 = YY5.add(XX5);
            FieldElement YYMinusXX5 = YY5.subtract(XX5);
            FieldElement X5 = XPlusYSq5.subtract(YYPlusXX5);
            FieldElement T5 = ZZ2_5.subtract(YYMinusXX5);

            FieldElement X_32P = X5.multiply(T5);
            FieldElement Y_32P = YYPlusXX5.multiply(YYMinusXX5);
            FieldElement Z_32P = YYMinusXX5.multiply(T5);

            FieldElement XX6 = X_32P.square();
            FieldElement YY6 = Y_32P.square();
            FieldElement ZZ2_6 = Z_32P.squareAndDouble();
            FieldElement XPlusY6 = X_32P.add(Y_32P);
            FieldElement XPlusYSq6 = XPlusY6.square();
            FieldElement YYPlusXX6 = YY6.add(XX6);
            FieldElement YYMinusXX6 = YY6.subtract(XX6);
            FieldElement X6 = XPlusYSq6.subtract(YYPlusXX6);
            FieldElement T6 = ZZ2_6.subtract(YYMinusXX6);

            this.X = X6.multiply(T6);
            this.Y = YYPlusXX6.multiply(YYMinusXX6);
            this.Z = YYMinusXX6.multiply(T6);
            this.T = X6.multiply(YYPlusXX6);

            return this;
        } else if (k == 7) {
            FieldElement XX1 = this.X.square();
            FieldElement YY1 = this.Y.square();
            FieldElement ZZ2_1 = this.Z.squareAndDouble();
            FieldElement XPlusY1 = this.X.add(this.Y);
            FieldElement XPlusYSq1 = XPlusY1.square();
            FieldElement YYPlusXX1 = YY1.add(XX1);
            FieldElement YYMinusXX1 = YY1.subtract(XX1);
            FieldElement X1 = XPlusYSq1.subtract(YYPlusXX1);
            FieldElement T1 = ZZ2_1.subtract(YYMinusXX1);

            FieldElement X_2P = X1.multiply(T1);
            FieldElement Y_2P = YYPlusXX1.multiply(YYMinusXX1);
            FieldElement Z_2P = YYMinusXX1.multiply(T1);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ2_2 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ2_2.subtract(YYMinusXX2);

            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ2_3 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ2_3.subtract(YYMinusXX3);

            FieldElement X_8P = X3.multiply(T3);
            FieldElement Y_8P = YYPlusXX3.multiply(YYMinusXX3);
            FieldElement Z_8P = YYMinusXX3.multiply(T3);

            FieldElement XX4 = X_8P.square();
            FieldElement YY4 = Y_8P.square();
            FieldElement ZZ2_4 = Z_8P.squareAndDouble();
            FieldElement XPlusY4 = X_8P.add(Y_8P);
            FieldElement XPlusYSq4 = XPlusY4.square();
            FieldElement YYPlusXX4 = YY4.add(XX4);
            FieldElement YYMinusXX4 = YY4.subtract(XX4);
            FieldElement X4 = XPlusYSq4.subtract(YYPlusXX4);
            FieldElement T4 = ZZ2_4.subtract(YYMinusXX4);

            FieldElement X_16P = X4.multiply(T4);
            FieldElement Y_16P = YYPlusXX4.multiply(YYMinusXX4);
            FieldElement Z_16P = YYMinusXX4.multiply(T4);

            FieldElement XX5 = X_16P.square();
            FieldElement YY5 = Y_16P.square();
            FieldElement ZZ2_5 = Z_16P.squareAndDouble();
            FieldElement XPlusY5 = X_16P.add(Y_16P);
            FieldElement XPlusYSq5 = XPlusY5.square();
            FieldElement YYPlusXX5 = YY5.add(XX5);
            FieldElement YYMinusXX5 = YY5.subtract(XX5);
            FieldElement X5 = XPlusYSq5.subtract(YYPlusXX5);
            FieldElement T5 = ZZ2_5.subtract(YYMinusXX5);

            FieldElement X_32P = X5.multiply(T5);
            FieldElement Y_32P = YYPlusXX5.multiply(YYMinusXX5);
            FieldElement Z_32P = YYMinusXX5.multiply(T5);

            FieldElement XX6 = X_32P.square();
            FieldElement YY6 = Y_32P.square();
            FieldElement ZZ2_6 = Z_32P.squareAndDouble();
            FieldElement XPlusY6 = X_32P.add(Y_32P);
            FieldElement XPlusYSq6 = XPlusY6.square();
            FieldElement YYPlusXX6 = YY6.add(XX6);
            FieldElement YYMinusXX6 = YY6.subtract(XX6);
            FieldElement X6 = XPlusYSq6.subtract(YYPlusXX6);
            FieldElement T6 = ZZ2_6.subtract(YYMinusXX6);

            FieldElement X_64P = X6.multiply(T6);
            FieldElement Y_64P = YYPlusXX6.multiply(YYMinusXX6);
            FieldElement Z_64P = YYMinusXX6.multiply(T6);

            FieldElement XX7 = X_64P.square();
            FieldElement YY7 = Y_64P.square();
            FieldElement ZZ2_7 = Z_64P.squareAndDouble();
            FieldElement XPlusY7 = X_64P.add(Y_64P);
            FieldElement XPlusYSq7 = XPlusY7.square();
            FieldElement YYPlusXX7 = YY7.add(XX7);
            FieldElement YYMinusXX7 = YY7.subtract(XX7);
            FieldElement X7 = XPlusYSq7.subtract(YYPlusXX7);
            FieldElement T7 = ZZ2_7.subtract(YYMinusXX7);

            this.X = X7.multiply(T7);
            this.Y = YYPlusXX7.multiply(YYMinusXX7);
            this.Z = YYMinusXX7.multiply(T7);
            this.T = X7.multiply(YYPlusXX7);

            return this;
        } else if (k == 8) {
            FieldElement XX1 = this.X.square();
            FieldElement YY1 = this.Y.square();
            FieldElement ZZ2_1 = this.Z.squareAndDouble();
            FieldElement XPlusY1 = this.X.add(this.Y);
            FieldElement XPlusYSq1 = XPlusY1.square();
            FieldElement YYPlusXX1 = YY1.add(XX1);
            FieldElement YYMinusXX1 = YY1.subtract(XX1);
            FieldElement X1 = XPlusYSq1.subtract(YYPlusXX1);
            FieldElement T1 = ZZ2_1.subtract(YYMinusXX1);

            FieldElement X_2P = X1.multiply(T1);
            FieldElement Y_2P = YYPlusXX1.multiply(YYMinusXX1);
            FieldElement Z_2P = YYMinusXX1.multiply(T1);

            FieldElement XX2 = X_2P.square();
            FieldElement YY2 = Y_2P.square();
            FieldElement ZZ2_2 = Z_2P.squareAndDouble();
            FieldElement XPlusY2 = X_2P.add(Y_2P);
            FieldElement XPlusYSq2 = XPlusY2.square();
            FieldElement YYPlusXX2 = YY2.add(XX2);
            FieldElement YYMinusXX2 = YY2.subtract(XX2);
            FieldElement X2 = XPlusYSq2.subtract(YYPlusXX2);
            FieldElement T2 = ZZ2_2.subtract(YYMinusXX2);

            FieldElement X_4P = X2.multiply(T2);
            FieldElement Y_4P = YYPlusXX2.multiply(YYMinusXX2);
            FieldElement Z_4P = YYMinusXX2.multiply(T2);

            FieldElement XX3 = X_4P.square();
            FieldElement YY3 = Y_4P.square();
            FieldElement ZZ2_3 = Z_4P.squareAndDouble();
            FieldElement XPlusY3 = X_4P.add(Y_4P);
            FieldElement XPlusYSq3 = XPlusY3.square();
            FieldElement YYPlusXX3 = YY3.add(XX3);
            FieldElement YYMinusXX3 = YY3.subtract(XX3);
            FieldElement X3 = XPlusYSq3.subtract(YYPlusXX3);
            FieldElement T3 = ZZ2_3.subtract(YYMinusXX3);

            FieldElement X_8P = X3.multiply(T3);
            FieldElement Y_8P = YYPlusXX3.multiply(YYMinusXX3);
            FieldElement Z_8P = YYMinusXX3.multiply(T3);

            FieldElement XX4 = X_8P.square();
            FieldElement YY4 = Y_8P.square();
            FieldElement ZZ2_4 = Z_8P.squareAndDouble();
            FieldElement XPlusY4 = X_8P.add(Y_8P);
            FieldElement XPlusYSq4 = XPlusY4.square();
            FieldElement YYPlusXX4 = YY4.add(XX4);
            FieldElement YYMinusXX4 = YY4.subtract(XX4);
            FieldElement X4 = XPlusYSq4.subtract(YYPlusXX4);
            FieldElement T4 = ZZ2_4.subtract(YYMinusXX4);

            FieldElement X_16P = X4.multiply(T4);
            FieldElement Y_16P = YYPlusXX4.multiply(YYMinusXX4);
            FieldElement Z_16P = YYMinusXX4.multiply(T4);

            FieldElement XX5 = X_16P.square();
            FieldElement YY5 = Y_16P.square();
            FieldElement ZZ2_5 = Z_16P.squareAndDouble();
            FieldElement XPlusY5 = X_16P.add(Y_16P);
            FieldElement XPlusYSq5 = XPlusY5.square();
            FieldElement YYPlusXX5 = YY5.add(XX5);
            FieldElement YYMinusXX5 = YY5.subtract(XX5);
            FieldElement X5 = XPlusYSq5.subtract(YYPlusXX5);
            FieldElement T5 = ZZ2_5.subtract(YYMinusXX5);

            FieldElement X_32P = X5.multiply(T5);
            FieldElement Y_32P = YYPlusXX5.multiply(YYMinusXX5);
            FieldElement Z_32P = YYMinusXX5.multiply(T5);

            FieldElement XX6 = X_32P.square();
            FieldElement YY6 = Y_32P.square();
            FieldElement ZZ2_6 = Z_32P.squareAndDouble();
            FieldElement XPlusY6 = X_32P.add(Y_32P);
            FieldElement XPlusYSq6 = XPlusY6.square();
            FieldElement YYPlusXX6 = YY6.add(XX6);
            FieldElement YYMinusXX6 = YY6.subtract(XX6);
            FieldElement X6 = XPlusYSq6.subtract(YYPlusXX6);
            FieldElement T6 = ZZ2_6.subtract(YYMinusXX6);

            FieldElement X_64P = X6.multiply(T6);
            FieldElement Y_64P = YYPlusXX6.multiply(YYMinusXX6);
            FieldElement Z_64P = YYMinusXX6.multiply(T6);

            FieldElement XX7 = X_64P.square();
            FieldElement YY7 = Y_64P.square();
            FieldElement ZZ2_7 = Z_64P.squareAndDouble();
            FieldElement XPlusY7 = X_64P.add(Y_64P);
            FieldElement XPlusYSq7 = XPlusY7.square();
            FieldElement YYPlusXX7 = YY7.add(XX7);
            FieldElement YYMinusXX7 = YY7.subtract(XX7);
            FieldElement X7 = XPlusYSq7.subtract(YYPlusXX7);
            FieldElement T7 = ZZ2_7.subtract(YYMinusXX7);

            FieldElement X_128P = X7.multiply(T7);
            FieldElement Y_128P = YYPlusXX7.multiply(YYMinusXX7);
            FieldElement Z_128P = YYMinusXX7.multiply(T7);

            FieldElement XX8 = X_128P.square();
            FieldElement YY8 = Y_128P.square();
            FieldElement ZZ2_8 = Z_128P.squareAndDouble();
            FieldElement XPlusY8 = X_128P.add(Y_128P);
            FieldElement XPlusYSq8 = XPlusY8.square();
            FieldElement YYPlusXX8 = YY8.add(XX8);
            FieldElement YYMinusXX8 = YY8.subtract(XX8);
            FieldElement X8 = XPlusYSq8.subtract(YYPlusXX8);
            FieldElement T8 = ZZ2_8.subtract(YYMinusXX8);

            this.X = X8.multiply(T8);
            this.Y = YYPlusXX8.multiply(YYMinusXX8);
            this.Z = YYMinusXX8.multiply(T8);
            this.T = X8.multiply(YYPlusXX8);

            return this;
        } else {
            for (int i = 0; i < k; i++) {
                this.dblInPlace();
            }

            return this;
        }
    }

    EdwardsPoint multiplyByPow2loop(int k) {
        if (!(k > 0)) {
            throw new IllegalArgumentException("Exponent must be positive and non-zero");
        }
        EdwardsPoint s = new EdwardsPoint(this.X, this.Y, this.Z, this.T);

        for (int i = 0; i < k - 1; i++) {
            s = s.dblInPlace();
        }

        // Unroll last doubling so we can go directly to extended coordinates.
        return s.dblInPlace();
    }

    EdwardsPoint multiplyByPow2Ext(int k) {
        if (!(k > 0)) {
            throw new IllegalArgumentException("Exponent must be positive and non-zero");
        }
        ProjectivePoint s = this.toProjective();
        for (int i = 0; i < k - 1; i++) {
            s = s.dbl().toProjective();
        }
        // Unroll last doubling so we can go directly to extended coordinates.
        return s.dbl().toExtended();
    }

    /**
     * Determine if this point is the identity.
     *
     * @return true if this point is the identity, false otherwise.
     */
    public boolean isIdentity() {
        return this.ctEquals(EdwardsPoint.IDENTITY) == 1;
    }

    /**
     * Determine if this point is in the 8-torsion subgroup $(\mathcal E[8])$, and
     * therefore of small order.
     *
     * @return true if this point is of small order, false otherwise.
     */
    public boolean isSmallOrder() {
        return this.multiplyByCofactor().isIdentity();
    }

    /**
     * Determine if this point is contained in the prime-order subgroup $(\mathcal
     * E[\ell])$, and has no torsion component.
     *
     * @return true if this point has zero torsion component and is in the
     *         prime-order subgroup, false otherwise.
     */
    public boolean isTorsionFree() {
        return this.multiply(Constants.BASEPOINT_ORDER).isIdentity();
    }

    /**
     * For debugging.
     */
    String printInternalRepresentation() {
        String ir = "EdwardsPoint(\n";
        ir += "    X: " + this.X.printInternalRepresentation() + ",\n";
        ir += "    Y: " + this.Y.printInternalRepresentation() + ",\n";
        ir += "    Z: " + this.Z.printInternalRepresentation() + ",\n";
        ir += "    T: " + this.T.printInternalRepresentation() + ",\n";
        ir += ")";
        return ir;
    }
}
