/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package au.edu.adelaide.fxmr.om;

import java.util.Iterator;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Combinatorial utilities.
 *
 * @since 3.3
 */
public final class CombinatoricsUtils {

	/** All long-representable factorials */
	static final long[] FACTORIALS = new long[] { 1l, 1l, 2l, 6l, 24l, 120l, 720l, 5040l, 40320l, 362880l, 3628800l,
			39916800l, 479001600l, 6227020800l, 87178291200l, 1307674368000l, 20922789888000l, 355687428096000l,
			6402373705728000l, 121645100408832000l, 2432902008176640000l };

	/** Stirling numbers of the second kind. */
	static final AtomicReference<long[][]> STIRLING_S2 = new AtomicReference<long[][]>(null);

	/** Private constructor (class contains only static methods). */
	private CombinatoricsUtils() {
	}

	/**
	 * Returns an exact representation of the
	 * <a href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
	 * Coefficient</a>, "{@code n choose k}", the number of {@code k}-element
	 * subsets that can be selected from an {@code n}-element set.
	 * <p>
	 * <Strong>Preconditions</strong>:
	 * <ul>
	 * <li>{@code 0 <= k <= n } (otherwise {@code MathIllegalArgumentException}
	 * is thrown)</li>
	 * <li>The result is small enough to fit into a {@code long}. The largest
	 * value of {@code n} for which all coefficients are
	 * {@code  < Long.MAX_VALUE} is 66. If the computed value exceeds
	 * {@code Long.MAX_VALUE} a {@code MathArithMeticException} is thrown.</li>
	 * </ul>
	 * </p>
	 *
	 * @param n
	 *            the size of the set
	 * @param k
	 *            the size of the subsets to be counted
	 * @return {@code n choose k}
	 * @throws NotPositiveException
	 *             if {@code n < 0}.
	 * @throws NumberIsTooLargeException
	 *             if {@code k > n}.
	 * @throws MathArithmeticException
	 *             if the result is too large to be represented by a long
	 *             integer.
	 */
	public static long binomialCoefficient(final int n, final int k) {
		if (!checkBinomial(n, k))
			return -1;
		if ((n == k) || (k == 0)) {
			return 1;
		}
		if ((k == 1) || (k == n - 1)) {
			return n;
		}
		// Use symmetry for large k
		if (k > n / 2) {
			return binomialCoefficient(n, n - k);
		}

		// We use the formula
		// (n choose k) = n! / (n-k)! / k!
		// (n choose k) == ((n-k+1)*...*n) / (1*...*k)
		// which could be written
		// (n choose k) == (n-1 choose k-1) * n / k
		long result = 1;
		if (n <= 61) {
			// For n <= 61, the naive implementation cannot overflow.
			int i = n - k + 1;
			for (int j = 1; j <= k; j++) {
				result = result * i / j;
				i++;
			}
		} else {
			return -1;
		}
		return result;
	}


	/**
	 * Returns an iterator whose range is the k-element subsets of {0, ..., n -
	 * 1} represented as {@code int[]} arrays.
	 * <p>
	 * The arrays returned by the iterator are sorted in descending order and
	 * they are visited in lexicographic order with significance from right to
	 * left. For example, combinationsIterator(4, 2) returns an Iterator that
	 * will generate the following sequence of arrays on successive calls to
	 * {@code next()}:
	 * </p>
	 * <p>
	 * {@code [0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]}
	 * </p>
	 * <p>
	 * If {@code k == 0} an Iterator containing an empty array is returned and
	 * if {@code k == n} an Iterator containing [0, ..., n -1] is returned.
	 * </p>
	 *
	 * @param n
	 *            Size of the set from which subsets are selected.
	 * @param k
	 *            Size of the subsets to be enumerated.
	 * @return an {@link Iterator iterator} over the k-sets in n.
	 * @throws NotPositiveException
	 *             if {@code n < 0}.
	 * @throws NumberIsTooLargeException
	 *             if {@code k > n}.
	 */
	public static Iterator<int[]> combinationsIterator(int n, int k) {
		return new Combinations(n, k).iterator();
	}

	/**
	 * Check binomial preconditions.
	 *
	 * @param n
	 *            Size of the set.
	 * @param k
	 *            Size of the subsets to be counted.
	 * @throws NotPositiveException
	 *             if {@code n < 0}.
	 * @throws NumberIsTooLargeException
	 *             if {@code k > n}.
	 */
	public static boolean checkBinomial(final int n, final int k) {
		if (n < k) {
			return false;
		}
		if (n < 0) {
			return false;
		}
		return true;
	}
}
