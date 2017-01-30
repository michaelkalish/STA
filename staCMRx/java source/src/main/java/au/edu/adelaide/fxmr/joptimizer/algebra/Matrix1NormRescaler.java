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

import au.edu.adelaide.fxmr.joptimizer.util.ColtUtils;
import cern.colt.function.IntIntDoubleFunction;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * Calculates the matrix rescaling factors so that the 1-norm of each row and
 * each column of the scaled matrix asymptotically converges to one.
 * 
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 * @see Daniel Ruiz,
 *      "A scaling algorithm to equilibrate both rows and columns norms in matrices"
 * @see Philip A. Knight, Daniel Ruiz, Bora Ucar
 *      "A Symmetry Preserving Algorithm for Matrix Scaling"
 */
public class Matrix1NormRescaler implements MatrixRescaler {

	private double eps = 1.e-3;

	public Matrix1NormRescaler() {
	}

	/**
	 * Scaling factors for symmetric (not singular) matrices. Just the
	 * subdiagonal elements of the matrix are required.
	 * 
	 * @see Daniel Ruiz,
	 *      "A scaling algorithm to equilibrate both rows and columns norms in matrices"
	 * @see Philip A. Knight, Daniel Ruiz, Bora Ucar
	 *      "A Symmetry Preserving Algorithm for Matrix Scaling"
	 */
	public DoubleMatrix1D getMatrixScalingFactorsSymm(DoubleMatrix2D A) {
		DoubleFactory1D F1 = DoubleFactory1D.dense;
		DoubleFactory2D F2 = DoubleFactory2D.sparse;
		int dim = A.columns();
		DoubleMatrix1D D1 = F1.make(dim, 1);
		DoubleMatrix2D AK = A.copy();
		DoubleMatrix2D DR = F2.identity(dim);
		DoubleMatrix1D DRInv = F1.make(dim);
		// log.debug("eps : " + eps);
		int maxIteration = 50;
		for (int k = 0; k <= maxIteration; k++) {
			double normR = -Double.MAX_VALUE;
			for (int i = 0; i < dim; i++) {
				// double dri = ALG.normInfinity(AK.viewRow(i));
				double dri = this.getRowInfinityNorm(AK, i);
				DR.setQuick(i, i, Math.sqrt(dri));
				DRInv.setQuick(i, 1. / Math.sqrt(dri));
				normR = Math.max(normR, Math.abs(1 - dri));
				if (Double.isNaN(normR)) {
					throw new IllegalArgumentException("matrix is singular");
				}
			}

			// log.debug("normR: " + normR);
			if (normR < eps) {
				break;
			}

			for (int i = 0; i < dim; i++) {
				double prevD1I = D1.getQuick(i);
				double newD1I = prevD1I * DRInv.getQuick(i);
				D1.setQuick(i, newD1I);
			}
			// logger.debug("D1: " + ArrayUtils.toString(D1.toArray()));

			AK = ColtUtils.diagonalMatrixMult(DRInv, AK, DRInv);
		}

		return D1;
	}

	/**
	 * 
	 * @param ASymm
	 *            symm matrix filled in its subdiagonal elements
	 * @param r
	 *            the index of the row
	 * @return
	 */
	private double getRowInfinityNorm(final DoubleMatrix2D ASymm, final int r) {

		final double[] maxValueHolder = new double[] { -Double.MAX_VALUE };

		IntIntDoubleFunction myFunct = new IntIntDoubleFunction() {
			public double apply(int i, int j, double pij) {
				// logger.warn("(" + i + "," + j + ")=" + pij);
				maxValueHolder[0] = Math.max(maxValueHolder[0], Math.abs(pij));
				return pij;
			}
		};

		// view A row from starting element to diagonal
		DoubleMatrix2D AR = ASymm.viewPart(r, 0, 1, r + 1);
		AR.forEachNonZero(myFunct);
		// view A col from diagonal to final element
		DoubleMatrix2D AC = ASymm.viewPart(r, r, ASymm.rows() - r, 1);
		AC.forEachNonZero(myFunct);

		return maxValueHolder[0];
	}

}
