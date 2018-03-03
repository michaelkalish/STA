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
package au.edu.adelaide.fxmr.joptimizer.optimizers;

import au.edu.adelaide.fxmr.joptimizer.algebra.CholeskyFactorization;
import au.edu.adelaide.fxmr.joptimizer.algebra.CholeskyMachine;
import au.edu.adelaide.fxmr.joptimizer.algebra.Matrix1NormRescaler;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.joptimizer.util.ColtUtils;
import au.edu.adelaide.fxmr.joptimizer.util.Utils;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.LUDecompositionQuick;
import cern.jet.math.Functions;
import cern.jet.math.Mult;

/**
 * Primal-dual interior-point method.
 * 
 * @see "S.Boyd and L.Vandenberghe, Convex Optimization, p. 609"
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 */
public class PrimalDualMethod extends OptimizationRequestHandler {
	private static final Mult NEGATE = Mult.mult(-1);
	private int iteration;

	public int optimize() throws Exception {
		int mIneq = getFi() == null ? 0 : getFi().length;
		int mEq = getFe() == null ? 0 : getFe().length;
		int dim = getDim();

		OptimizationResponse response = new OptimizationResponse();
		DoubleMatrix1D X0 = getInitialPoint();
		if (X0 == null || mIneq == 0) {
			return OptimizationResponse.FAILED;
		}

		// check X0 feasibility
		DoubleMatrix1D fiX0 = F1.make(mIneq);
		getFi(X0, fiX0);

		DoubleMatrix1D feX0 = F1.make(mEq);
		if (mEq > 0)
			rPri(X0, feX0);

		int maxIndex = Utils.getMaxIndex(fiX0);
		double maxValue = fiX0.getQuick(maxIndex);
		double rPriX0Norm = Math.sqrt(ALG.norm2(feX0));
		if (maxValue >= 0 || rPriX0Norm > getToleranceFeas()) {
			throw new Exception(
					"initial point must be strictly feasible, maxValue = " + maxValue + ", rPriX0Norm = " + rPriX0Norm);

		}

		DoubleMatrix1D L0 = getInitialLagrangian();
		if (L0 != null) {
			for (int j = 0; j < L0.size(); j++) {
				// must be >0
				if (L0.get(j) <= 0) {
					throw new IllegalArgumentException("initial lagrangian must be strictly > 0");
				}
			}
		} else {
			// must be >0
			L0 = F1.make(mIneq, Math.min(1, (double) dim / mIneq));
		}

		DoubleMatrix1D X = X0.copy();
		DoubleMatrix1D V = mEq == 0 ? null : feX0;
		DoubleMatrix1D L = L0.copy();
		DoubleMatrix1D rCentXLt = F1.make(L.size());
		// DoubleMatrix1D rCentX1L1t = F1.make(L.size());
		double t;
		iteration = 0;
		int endIter = getMaxIteration() + 1;

		SimpleLinearConstraint[] fis = getFi();
		SimpleLinearConstraint[] fes = getFe();

		// setup Cholesky or LUQuick?
		boolean useCM = mIneq <= 20 && mEq == 0;
		CholeskyMachine cm = null;
		LUDecompositionQuick lu = null;
		LUDecompositionQuick luZero = null;

		// ALG.rank( getHessF0())

		if (useCM) {
			cm = new CholeskyMachine(getHessF0());
			if (cm.isFailed()) {
				// Treading a fine line here - the hessian of the objective is
				// not positive definite!
				cm = null;
				useCM = false;
			}
		}
		if (!useCM) {
			lu = new LUDecompositionQuick();
			luZero = new LUDecompositionQuick(0);
		}

		// Save method calls
		double tolFeas = getToleranceFeas();
		double tol = getTolerance();

		// All of these used to be allocated every iteration!
		DoubleMatrix1D rPriX = F1.make(mEq);
		DoubleMatrix1D fiX = F1.make(mIneq);
		DoubleMatrix1D gradSum = F1.make(dim);
		DoubleMatrix1D rDualXLV = F1.make(dim);
		DoubleMatrix1D stepL = F1.make(mIneq);
		DoubleMatrix1D stepV = F1.make(mEq);
		DoubleMatrix1D X1 = F1.make(dim);
		DoubleMatrix1D L1 = F1.make(mIneq);
		DoubleMatrix1D V1 = F1.make(mEq);
		DoubleMatrix1D gradF0X = F1.make(dim);
		DoubleMatrix2D HessSum = F2.make(dim, dim);
		// DoubleMatrix2D HessSumOrig = F2.make(dim, dim);
		DoubleMatrix2D HInvAT = mEq == 0 ? null : F2.make(dim, mEq);
		DoubleMatrix1D HInvg = F1.make(dim);
		DoubleMatrix2D HessBase = getHessF0();

		DoubleMatrix2D matA = null;
		DoubleMatrix2D MenoSLower = null;
		if (mEq > 0) {
			matA = F2.make(mEq, dim);
			for (int i = 0; i < mEq; i++) {
				SimpleLinearConstraint eq = fes[i];
				matA.setQuick(i, eq.getPosIndex(), 1);
				matA.setQuick(i, eq.getNegIndex(), -1);
			}
			MenoSLower = F2.make(mEq, mEq);
		}

		while (true) {
			iteration++;

			// determine functions evaluations
			getGradF0(X, gradF0X);
			getFi(X, fiX);

			// determine t
			double surrDG = getSurrogateDualityGap(fiX, L);
			t = getMu() * mIneq / surrDG;

			// determine residuals
			if (mEq > 0)
				rPri(X, rPriX);
			rCent(fiX, L, t, rCentXLt);
			rDual(fis, gradF0X, fes, L, V, rDualXLV);
			// double rCentXLtNorm = Math.sqrt(ALG.norm2(rCentXLt));
			double rDualXLVNorm = Math.sqrt(ALG.norm2(rDualXLV));
			double rPriXNorm = mEq > 0 ? Math.sqrt(ALG.norm2(rPriX)) : 0;
			// double normRXLVt = Math.sqrt(rCentXLtNorm * rCentXLtNorm +
			// rDualXLVNorm * rDualXLVNorm);

			// exit condition
			if (rPriXNorm <= tolFeas && rDualXLVNorm <= tolFeas && surrDG <= tol) {
				response.setReturnCode(OptimizationResponse.SUCCESS);
				break;
			}

			// iteration limit condition
			if (iteration == endIter) {
				// System.err.println("Max iterations " + getMaxIteration() + "
				// limit reached ExitA:" + rDualXLVNorm + "<="
				// + getToleranceFeas() + ", ExitB: " + surrDG + "<=" +
				// getTolerance());
				response.setReturnCode(OptimizationResponse.WARN);
				response.setSolution(X.toArray());
				setOptimizationResponse(response);
				return response.getReturnCode();
				// throw new Exception("Max iterations limit reached");
			}

			// compute primal-dual search direction
			// a) prepare 11.55 system
			// DoubleMatrix2D HessSum = getHessF0().copy();
			gradSum.assign(gradF0X);
			for (int j = 0; j < mIneq; j++) {
				SimpleLinearConstraint fi = fis[j];
				int pos = fi.getPosIndex();
				int neg = fi.getNegIndex();

				double adder = 1. / (-t * fiX.getQuick(j));
				gradSum.setQuick(pos, gradSum.getQuick(pos) + adder);
				gradSum.setQuick(neg, gradSum.getQuick(neg) - adder);
			}
			if (mEq > 0) {
				for (int j = 0; j < mEq; j++) {
					SimpleLinearConstraint fi = fes[j];
					int pos = fi.getPosIndex();
					int neg = fi.getNegIndex();

					double v = V.getQuick(j);
					gradSum.setQuick(pos, gradSum.getQuick(pos) + v);
					gradSum.setQuick(neg, gradSum.getQuick(neg) - v);
				}
			}

			if (useCM) {
				cm.reset();
				for (int j = 0; j < mIneq; j++) {
					final double c = -L.getQuick(j) / fiX.getQuick(j);
					SimpleLinearConstraint fi = fis[j];
					// tmpCon.assign(fi.gradient());

					cm.rank1Update(fi, Math.sqrt(c));
				}
				cm.solve(gradSum, HInvg);
			} else {
				HessSum.assign(HessBase);
				for (int j = 0; j < mIneq; j++) {
					final double c = -L.getQuick(j) / fiX.getQuick(j);

					SimpleLinearConstraint fi = fis[j];
					int pos = fi.getPosIndex();
					int neg = fi.getNegIndex();

					HessSum.setQuick(pos, pos, HessSum.getQuick(pos, pos) + c);
					HessSum.setQuick(neg, neg, HessSum.getQuick(neg, neg) + c);
					HessSum.setQuick(pos, neg, HessSum.getQuick(pos, neg) - c);
					HessSum.setQuick(neg, pos, HessSum.getQuick(neg, pos) - c);
				}

				// b) solving 11.55 system
				lu.decompose(HessSum); // same as Hpd
				HInvg.assign(gradSum);

				if (lu.det() == 0) {
					// Singular matrix to tolerance...Make the tolerance ZERO
					// with a new LU decomposition
					luZero.decompose(HessSum);
					if (luZero.isNonsingular()) {
						luZero.solve(HInvg);
					} else {
						// This will not be the best inverse going around...
						try {
							DoubleMatrix2D inverse = Algebra.ZERO.inverse(HessSum);
							DoubleMatrix1D tmp = F1.make(dim);
							inverse.zMult(HInvg, tmp);
							HInvg.assign(tmp);
						} catch (IllegalArgumentException e) {
							// Completely singular - oh what to do?
							System.out.println(rPriXNorm + " - " + rDualXLVNorm + " - " + surrDG);
							return OptimizationResponse.FAILED;
						}

					}
				} else {
					lu.solve(HInvg);
				}
			}

			DoubleMatrix1D stepX;
			if (mEq > 0) {
				HInvAT.assign(matA.viewDice());
				// stepV also needs to be calculated - uses rPriX as RHS
				lu.solve(HInvAT);
				ColtUtils.subdiagonalMultiply(matA, HInvAT, MenoSLower);
				DoubleMatrix1D AHInvg = ALG.mult(matA, HInvg);

				CholeskyFactorization MSFact = new CholeskyFactorization(MenoSLower, new Matrix1NormRescaler());
				MSFact.factorize();

				rPriX.assign(AHInvg, Functions.minus);
				stepV = MSFact.solve(rPriX);
				
				stepX = HInvg.assign(ALG.mult(HInvAT, stepV), Functions.plus).assign(NEGATE);
			} else {
				// Solving KKT system via elimination
				stepX = HInvg.assign(NEGATE);
			}

			// c) solving for L
			DoubleMatrix1D a2 = rCentXLt.assign(fiX, Functions.div);
			makeStepL(fis, stepX, L, stepL, fiX, a2);

			// line search and update
			// a) sMax computation
			double sMax = Double.MAX_VALUE;
			for (int j = 0; j < mIneq; j++) {
				if (stepL.get(j) < 0) {
					sMax = Math.min(-L.getQuick(j) / stepL.getQuick(j), sMax);
				}
			}
			sMax = Math.min(1, sMax);
			double s = 0.99 * sMax;

			// b) backtracking with f
			int cnt = 0;
			boolean areAllNegative = true;
			while (cnt < 500) {
				cnt++;
				// X1 = X + s*stepX
				add(X, stepX, s, X1);
				areAllNegative = true;
				for (int i = 0; areAllNegative && i < mIneq; i++)
					areAllNegative = fis[i].value(X1) < 0;

				if (areAllNegative)
					break;
				s = getBeta() * s;
			}

			if (!areAllNegative) {
				// exited from the feasible region
				System.err.println("Optimization failed: impossible to remain within the faesible region");
				throw new Exception("Optimization failed: impossible to remain within the faesible region");
			}

			add(L, stepL, s, L1);
			if (mEq != 0)
				add(V, stepV, s, V1);

			// c) backtracking with norm - This was removed because it never
			// happens?!
			// double previousNormRX1L1V1t = Double.NaN;
			// cnt = 0;
			// while (cnt < 500) {
			// cnt++;
			// add(X, stepX, s, X1);
			// add(L, stepL, s, L1);
			//
			// if (isInDomainF0(X1)) {
			// getFi(X1, fiX1);
			// getGradF0(X1, gradF0X1);
			//
			// rCent(fiX1, L1, t, rCentX1L1t);
			// rDual(fis, gradF0X1, L1, rDualX1L1V1);
			// double normRX1L1V1t = Math.sqrt(ALG.norm2(rCentX1L1t) +
			// ALG.norm2(rDualX1L1V1));
			// if (normRX1L1V1t <= (1 - getAlpha() * s) * normRXLVt) {
			// break;
			// }
			//
			// if (!Double.isNaN(previousNormRX1L1V1t)) {
			// if (previousNormRX1L1V1t <= normRX1L1V1t) {
			// break;
			// }
			// }
			// previousNormRX1L1V1t = normRX1L1V1t;
			// }
			//
			// s = getBeta() * s;
			// }

			// update
			X.assign(X1);
			L.assign(L1);
			if (mEq != 0)
				V.assign(V1);
		}

		response.setSolution(X.toArray());

		setOptimizationResponse(response);
		return response.getReturnCode();
	}

	/**
	 * Returns result = v1 + c*v2.
	 */
	private void add(DoubleMatrix1D v1, DoubleMatrix1D v2, double c, DoubleMatrix1D result) {
		int n = result.size();
		for (int i = 0; i < n; i++)
			result.setQuick(i, v1.getQuick(i) + c * v2.getQuick(i));

		// result.assign(v2);
		// //Save object creation
		// mult.multiplicator = c;
		// result.assign(mult);
		// result.assign(v1, Functions.plus);
	}

	private void makeStepL(SimpleLinearConstraint[] fis, DoubleMatrix1D stepX, DoubleMatrix1D L, DoubleMatrix1D result,
			DoubleMatrix1D fiX, DoubleMatrix1D a2) {
		// a2 - gradfi * stepX .* L ./ fiX
		int i = 0;
		for (SimpleLinearConstraint fi : fis) {
			double v = a2.getQuick(i)
					- (stepX.getQuick(fi.getPosIndex()) - stepX.getQuick(fi.getNegIndex())) * L.get(i) / fiX.get(i);
			result.setQuick(i++, v);
		}
	}

	/**
	 * Surrogate duality gap.
	 * 
	 * @see "Convex Optimization, 11.59"
	 */
	private double getSurrogateDualityGap(DoubleMatrix1D fiX, DoubleMatrix1D L) {
		return -ALG.mult(fiX, L);
	}

	/**
	 * ColtUtils.zMultTranspose(getA(), V, ColtUtils.zMultTranspose(GradFiX, L,
	 * gradF0X, 1.), 1.)
	 * 
	 * @see "Convex Optimization, p. 610"
	 */
	private void rDual(SimpleLinearConstraint[] fis, DoubleMatrix1D gradF0X, SimpleLinearConstraint[] fes,
			DoubleMatrix1D L, DoubleMatrix1D V, DoubleMatrix1D result) {

		result.assign(gradF0X);
		int j = 0;
		for (SimpleLinearConstraint fi : fis) {
			double curL = L.getQuick(j++);
			int i = fi.getPosIndex();
			result.setQuick(i, result.getQuick(i) + curL);
			i = fi.getNegIndex();
			result.setQuick(i, result.getQuick(i) - curL);
		}

		if (V != null) {
			j = 0;
			for (SimpleLinearConstraint fe : fes) {
				double curV = V.getQuick(j++);
				int i = fe.getPosIndex();
				result.setQuick(i, result.getQuick(i) + curV);
				i = fe.getNegIndex();
				result.setQuick(i, result.getQuick(i) - curV);
			}
		}
	}

	/**
	 * @param rCentXLt
	 * @see "Convex Optimization, p. 610"
	 */
	private void rCent(DoubleMatrix1D fiX, DoubleMatrix1D L, double t, DoubleMatrix1D out) {
		int n = out.size();
		double oneOnT = -1.0 / t;
		for (int i = 0; i < n; i++)
			out.setQuick(i, -L.getQuick(i) * fiX.getQuick(i) + oneOnT);
	}

	/**
	 * rPri := Ax where A is the equality matrix
	 */
	protected void rPri(DoubleMatrix1D X, DoubleMatrix1D out) {
		int i = 0;
		for (SimpleLinearConstraint c : request.getA())
			out.setQuick(i++, c.value(X));
	}

	public int getIteration() {
		return iteration;
	}
}
