/*
 * Copyright 2011-2014 JOptimizer
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */
package au.edu.adelaide.fxmr.joptimizer.util;

import cern.colt.function.IntIntDoubleFunction;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;

/**
 * Support class for recurrent algebra with Colt.
 * 
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 */
public class ColtUtils {

	/**
	 * Matrix-vector multiplication with diagonal matrix.
	 * 
	 * @param diagonalM
	 *            diagonal matrix M
	 * @param vector
	 * @return M.x
	 */

	/**
	 * Matrix-vector multiplication with diagonal matrix.
	 * 
	 * @param diagonalM
	 *            diagonal matrix M, in the form of a vector of its diagonal
	 *            elements
	 * @param vector
	 * @return M.x
	 */
	public static final DoubleMatrix1D diagonalMatrixMult(DoubleMatrix1D diagonalM, DoubleMatrix1D vector) {
		int n = diagonalM.size();
		DoubleMatrix1D ret = DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++) {
			ret.setQuick(i, diagonalM.getQuick(i) * vector.getQuick(i));
		}
		return ret;
	}

	/**
	 * Return diagonalU.A with diagonalU diagonal.
	 * 
	 * @param diagonal
	 *            matrix U, in the form of a vector of its diagonal elements
	 * @return U.A
	 */
	public static final DoubleMatrix2D diagonalMatrixMult(final DoubleMatrix1D diagonalU, DoubleMatrix2D A) {
		int r = diagonalU.size();
		int c = A.columns();
		final DoubleMatrix2D ret;
		if (A instanceof SparseDoubleMatrix2D) {
			ret = DoubleFactory2D.sparse.make(r, c);
			A.forEachNonZero(new IntIntDoubleFunction() {
				public double apply(int i, int j, double aij) {
					ret.setQuick(i, j, aij * diagonalU.getQuick(i));
					return aij;
				}
			});
		} else {
			ret = DoubleFactory2D.dense.make(r, c);
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c; j++) {
					ret.setQuick(i, j, A.getQuick(i, j) * diagonalU.getQuick(i));
				}
			}
		}

		return ret;
	}

	/**
	 * Return diagonalU.A.diagonalV with diagonalU and diagonalV diagonal.
	 * 
	 * @param diagonalU
	 *            diagonal matrix U, in the form of a vector of its diagonal
	 *            elements
	 * @param diagonalV
	 *            diagonal matrix V, in the form of a vector of its diagonal
	 *            elements
	 * @return U.A.V
	 */
	public static final DoubleMatrix2D diagonalMatrixMult(final DoubleMatrix1D diagonalU, DoubleMatrix2D A,
			final DoubleMatrix1D diagonalV) {
		int r = A.rows();
		int c = A.columns();
		final DoubleMatrix2D ret;
		if (A instanceof SparseDoubleMatrix2D) {
			ret = DoubleFactory2D.sparse.make(r, c);
			A.forEachNonZero(new IntIntDoubleFunction() {
				public double apply(int i, int j, double aij) {
					ret.setQuick(i, j, aij * diagonalU.getQuick(i) * diagonalV.getQuick(j));
					return aij;
				}
			});
		} else {
			ret = DoubleFactory2D.dense.make(r, c);
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c; j++) {
					ret.setQuick(i, j, A.getQuick(i, j) * diagonalU.getQuick(i) * diagonalV.getQuick(j));
				}
			}
		}

		return ret;
	}

	/**
	 * Return the sub-diagonal result of the multiplication. If A is sparse,
	 * returns a sparse matrix (even if, generally speaking, the multiplication
	 * of two sparse matrices is not sparse) because the result is at least 50%
	 * (aside the diagonal elements) sparse.
	 */
	public static void subdiagonalMultiply(final DoubleMatrix2D A, final DoubleMatrix2D B, DoubleMatrix2D ret) {
		final int r = A.rows();
		final int rc = A.columns();

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < i + 1; j++) {
				double s = 0;
				for (int k = 0; k < rc; k++) {
					s += A.getQuick(i, k) * B.getQuick(k, j);
				}
				ret.setQuick(i, j, s);
			}
		}
	}

	/**
	 * Returns v = A.a + beta*b. Useful in avoiding the need of the copy() in
	 * the colt api.
	 */
	public static final DoubleMatrix1D zMult(final DoubleMatrix2D A, final DoubleMatrix1D a, final DoubleMatrix1D b,
			final double beta) {
		int aRows = A.rows();
		final DoubleMatrix1D ret = DoubleFactory1D.dense.make(aRows);

		if (A instanceof SparseDoubleMatrix2D) {
			// if(1==2){
			// sparse matrix
			A.forEachNonZero(new IntIntDoubleFunction() {
				public double apply(int i, int j, double Aij) {
					ret.setQuick(i, ret.getQuick(i) + Aij * a.getQuick(j));
					return Aij;
				}
			});
			for (int i = 0; i < ret.size(); i++) {
				ret.setQuick(i, ret.getQuick(i) + beta * b.getQuick(i));
			}
		} else {
			// dense matrix
			int aCols = A.columns();
			for (int i = 0; i < aRows; i++) {
				// DoubleMatrix1D AI = A.viewRow(i);
				double vi = beta * b.getQuick(i);
				for (int j = 0; j < aCols; j++) {
					vi += A.getQuick(i, j) * a.getQuick(j);
				}
				ret.setQuick(i, vi);
			}
		}

		return ret;
	}

	public static final DoubleMatrix1D zMult1(final DoubleMatrix2D A, final DoubleMatrix1D a, final DoubleMatrix1D b) {

		int aRows = A.rows();
		final DoubleMatrix1D ret = DoubleFactory1D.dense.make(aRows);

		if (A instanceof SparseDoubleMatrix2D) {
			// if(1==2){
			// sparse matrix
			A.forEachNonZero(new IntIntDoubleFunction() {
				public double apply(int i, int j, double Aij) {
					ret.setQuick(i, ret.getQuick(i) + Aij * a.getQuick(j));
					return Aij;
				}
			});
			for (int i = 0; i < ret.size(); i++) {
				ret.setQuick(i, ret.getQuick(i) + b.getQuick(i));
			}
		} else {
			// dense matrix
			int aCols = A.columns();
			for (int i = 0; i < aRows; i++) {
				// DoubleMatrix1D AI = A.viewRow(i);
				double vi = b.getQuick(i);
				for (int j = 0; j < aCols; j++) {
					vi += A.getQuick(i, j) * a.getQuick(j);
				}
				ret.setQuick(i, vi);
			}
		}

		return ret;
	}

	/**
	 * Returns v = A[T].a + beta*b. Useful in avoiding the need of the copy() in
	 * the colt api.
	 */
	public static final DoubleMatrix1D zMultTranspose(final DoubleMatrix2D A, final DoubleMatrix1D a,
			final DoubleMatrix1D b, final double beta) {
		if (A.rows() != a.size() || A.columns() != b.size()) {
			throw new IllegalArgumentException("wrong matrices dimensions");
		}
		final DoubleMatrix1D ret = DoubleFactory1D.dense.make(A.columns());

		if (A instanceof SparseDoubleMatrix2D) {
			// if(1==2){
			A.forEachNonZero(new IntIntDoubleFunction() {
				public double apply(int i, int j, double Aij) {
					// log.debug(i + "," + j + ": " + Aij + ",
					// "+ret.getQuick(j)+", "+a.getQuick(i));
					ret.setQuick(j, ret.getQuick(j) + Aij * a.getQuick(i));
					return Aij;
				}
			});
			if (Double.compare(0., beta) != 0) {
				for (int i = 0; i < ret.size(); i++) {
					ret.setQuick(i, ret.getQuick(i) + beta * b.getQuick(i));
				}
			}
		} else {
			for (int i = 0; i < A.columns(); i++) {
				double vi = beta * b.getQuick(i);
				for (int j = 0; j < A.rows(); j++) {
					vi += A.getQuick(j, i) * a.getQuick(j);
				}
				ret.setQuick(i, vi);
			}
		}

		return ret;
	}

	/**
	 * Returns v = v1 + v2. Useful in avoiding the need of the copy() in the
	 * colt api.
	 */
	public static final DoubleMatrix1D add(DoubleMatrix1D v1, DoubleMatrix1D v2) {
		int n = v1.size();
		DoubleMatrix1D ret = DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++)
			ret.setQuick(i, v1.getQuick(i) + v2.getQuick(i));
		return ret;
	}

	/**
	 * Returns v = v1 + c*v2. Useful in avoiding the need of the copy() in the
	 * colt api.
	 */
	public static final DoubleMatrix1D add(DoubleMatrix1D v1, DoubleMatrix1D v2, double c) {
		if (v1.size() != v2.size()) {
			throw new IllegalArgumentException("wrong vectors dimensions");
		}
		DoubleMatrix1D ret = DoubleFactory1D.dense.make(v1.size());
		for (int i = 0; i < ret.size(); i++) {
			ret.setQuick(i, v1.getQuick(i) + c * v2.getQuick(i));
		}

		return ret;
	}

	/**
	 * Returns v = c * v1. Useful in avoiding the need of the copy() in the colt
	 * api.
	 */
	public static final DoubleMatrix1D scalarMult(DoubleMatrix1D v1, double c) {
		DoubleMatrix1D ret = DoubleFactory1D.dense.make(v1.size());
		for (int i = 0; i < ret.size(); i++) {
			ret.setQuick(i, c * v1.getQuick(i));
		}

		return ret;
	}
}
