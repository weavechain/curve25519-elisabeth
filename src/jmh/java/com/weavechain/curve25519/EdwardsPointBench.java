/*
 * This file is part of curve25519-elisabeth.
 * Copyright (c) 2019 Jack Grigg
 * See LICENSE for licensing information.
 */

package com.weavechain.curve25519;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.*;

@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 2, timeUnit = TimeUnit.SECONDS)
@Fork(1)
@State(Scope.Benchmark)
public class EdwardsPointBench {
    public EdwardsPoint P;
    public CompressedEdwardsY Penc;
    public EdwardsBasepointTable Pt;
    public EdwardsPoint Q;
    public Scalar a;
    public Scalar b;

    static Scalar randomScalar(Random r) {
        byte[] input = new byte[64];
        r.nextBytes(input);
        return Scalar.fromBytesModOrderWide(input);
    }

    @Setup
    public void prepare() {
        Random r = new Random();
        this.P = Constants.ED25519_BASEPOINT_TABLE.multiply(randomScalar(r));
        this.Penc = this.P.compress();
        this.Pt = new EdwardsBasepointTable(this.P);
        this.Q = Constants.ED25519_BASEPOINT_TABLE.multiply(randomScalar(r));
        this.a = randomScalar(r);
        this.b = randomScalar(r);
    }

    @Benchmark
    public EdwardsPoint decompress() throws InvalidEncodingException {
        return this.Penc.decompress();
    }

    @Benchmark
    public CompressedEdwardsY compress() {
        return this.P.compress();
    }

    @Benchmark
    public EdwardsPoint add() {
        return this.P.add(this.Q);
    }

    @Benchmark
    public EdwardsPoint dbl() {
        return this.P.dbl();
    }

    @Benchmark
    public EdwardsPoint variableBaseScalarMultiply() {
        return this.P.multiply(this.a);
    }

    @Benchmark
    public EdwardsPoint fixedBaseScalarMultiply() {
        return this.Pt.multiply(this.a);
    }

    @Benchmark
    public EdwardsPoint doubleScalarMultiplyBasepoint() {
        return EdwardsPoint.vartimeDoubleScalarMultiplyBasepoint(this.a, this.P, this.b);
    }
}
