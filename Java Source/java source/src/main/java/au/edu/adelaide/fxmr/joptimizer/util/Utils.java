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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

/**
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 */
public class Utils {

	public static final Double RELATIVE_MACHINE_PRECISION;

	static {
		double eps = 1.;
		do {
			eps /= 2.;
		} while ((double) (1. + (eps / 2.)) != 1.);

		RELATIVE_MACHINE_PRECISION = eps;
	}

	/**
	 * Calculate the scaled residual <br>
	 * ||Ax-b||_oo/( ||A||_oo . ||x||_oo + ||b||_oo ), with <br>
	 * ||x||_oo = max(||x[i]||)
	 */
	public static double calculateScaledResidual(DoubleMatrix2D A, DoubleMatrix1D x, DoubleMatrix1D b) {
		double residual = -Double.MAX_VALUE;
		double nix = Algebra.DEFAULT.normInfinity(x);
		double nib = Algebra.DEFAULT.normInfinity(b);
		if (Double.compare(nix, 0.) == 0 && Double.compare(nib, 0.) == 0) {
			return 0;
		} else {
			double num = Algebra.DEFAULT.normInfinity(ColtUtils.zMult(A, x, b, -1));
			double den = Algebra.DEFAULT.normInfinity(A) * nix + nib;
			residual = num / den;
			// log.debug("scaled residual: " + residual);

			return residual;
		}
	}

	/**
	 * Get the index of the maximum entry.
	 */
	public static int getMaxIndex(DoubleMatrix1D v) {
		int maxIndex = -1;
		double maxValue = -Double.MAX_VALUE;
		for (int i = 0; i < v.size(); i++) {
			if (v.getQuick(i) > maxValue) {
				maxIndex = i;
				maxValue = v.getQuick(i);
			}
		}
				
		return maxIndex;
	}

}
