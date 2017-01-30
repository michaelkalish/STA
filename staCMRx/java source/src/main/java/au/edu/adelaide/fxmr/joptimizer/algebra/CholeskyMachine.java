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
package au.edu.adelaide.fxmr.joptimizer.algebra;

import java.util.Arrays;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.joptimizer.util.Utils;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

/**
 * Cholesky L.L[T] factorization and inverse for symmetric and positive matrix:
 * 
 * Q = L.L[T], L lower-triangular
 * 
 * Just the subdiagonal elements of Q are used.
 * 
 * @author <a href="mailto:alberto.trivellato@gmail.com">alberto trivellato</a>
 */
public class CholeskyMachine {

	private int dim;
	private DoubleMatrix2D Q;
	private double[][] LDataQ;
	private double[][] LDataCur;

	private double[] y;
	private boolean failed;

	public boolean isFailed() {
		return failed;
	}

	/**
	 * 
	 * @param Q
	 *            the matrix to factorize
	 */
	public CholeskyMachine(DoubleMatrix2D Q) {
		this.dim = Q.rows();
		this.Q = Q;
		// this.rescaler = rescaler;

		y = new double[dim];
		LDataQ = new double[dim][];
		LDataCur = new double[dim][];
		for (int i = 0; i < dim; i++) {
			LDataQ[i] = new double[i + 1];
			LDataCur[i] = new double[i + 1];
		}

		factorizeQ();
	}

	public void rank1Update(DoubleMatrix1D x, double conSQRT) {
		// function [L] = cholupdate(L,x)
		int n = x.size();
		for (int i = 0; i < n; i++)
			x.setQuick(i, x.getQuick(i) * conSQRT);

		for (int k = 0; k < dim; k++) {
			double lkk = LDataCur[k][k];
			double xk = x.get(k);
			double r = Math.sqrt(lkk * lkk + xk * xk);
			double c = r / lkk;
			double s = xk / lkk;
			LDataCur[k][k] = r;

			for (int i = k + 1; i < dim; i++)
				LDataCur[i][k] = (LDataCur[i][k] + s * x.get(i)) / c;

			for (int i = k + 1; i < dim; i++)
				x.set(i, c * x.get(i) - s * LDataCur[i][k]);
		}
	}

	public void rank1Update(SimpleLinearConstraint fi, double conSQRT) {
		Arrays.fill(y, 0);
		y[fi.getNegIndex()] = -conSQRT;
		y[fi.getPosIndex()] = conSQRT;

		int start = Math.min(fi.getNegIndex(), fi.getPosIndex());
		for (int k = start; k < dim; k++) {
			double lkk = LDataCur[k][k];
			double xk = y[k];
			double r = Math.sqrt(lkk * lkk + xk * xk);
			double c = r / lkk;
			double s = xk / lkk;
			LDataCur[k][k] = r;

			for (int i = k + 1; i < dim; i++)
				LDataCur[i][k] = (LDataCur[i][k] + s * y[i]) / c;

			for (int i = k + 1; i < dim; i++)
				y[i] = c * y[i] - s * LDataCur[i][k];
		}
	}

	public void reset() {
		for (int i = 0; i < dim; i++)
			System.arraycopy(LDataQ[i], 0, LDataCur[i], 0, i + 1);
	}

	/**
	 * Cholesky factorization L of psd matrix, Q = L.LT. Construction of the
	 * matrix L.
	 */
	private void factorizeQ() {
		// double threshold = Math.pow(Utils.getDoubleMachineEpsilon(), 2);
		double threshold = Utils.RELATIVE_MACHINE_PRECISION;
		for (int i = 0; i < dim; i++) {
			double[] LDataI = LDataQ[i];
			// j < i
			for (int j = 0; j < i; j++) {
				double[] LDataJ = LDataQ[j];
				double sum = 0.0;
				for (int k = 0; k < j; k++) {
					sum += LDataI[k] * LDataJ[k];
				}
				LDataI[j] = 1.0 / LDataJ[j] * (Q.getQuick(i, j) - sum);
			}
			// j==i
			double sum = 0.0;
			for (int k = 0; k < i; k++) {
				sum += LDataI[k] * LDataI[k];
			}
			double d = Q.getQuick(i, i) - sum;
			if (!(d > threshold)) {
				failed = true;
			}
			LDataI[i] = Math.sqrt(d);
		}
	}

	public void solve(DoubleMatrix1D b, DoubleMatrix1D result) {
		// Solve L.y = b
		final double[] y = this.y;
		for (int i = 0; i < dim; i++) {
			double[] LI = LDataCur[i];
			double sum = 0;
			for (int j = 0; j < i; j++) {
				sum += LI[j] * y[j];
			}
			y[i] = (b.getQuick(i) - sum) / LI[i];
		}

		// Solve L[T].x = y
		for (int i = dim - 1; i > -1; i--) {
			double sum = 0;
			for (int j = dim - 1; j > i; j--) {
				sum += LDataCur[j][i] * result.getQuick(j);
			}
			result.setQuick(i, (y[i] - sum) / LDataCur[i][i]);
		}
	}

	public double[][] getLDataCur() {
		return LDataCur;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int r = 0; r < dim; r++) {
			for (int c = 0; c <= r; c++) {
				sb.append('\t');
				sb.append(LDataCur[r][c]);
			}
			for (int z = r + 1; z < dim; z++)
				sb.append("\t0");
			sb.append('\n');
		}

		return sb.toString();
	}

}
