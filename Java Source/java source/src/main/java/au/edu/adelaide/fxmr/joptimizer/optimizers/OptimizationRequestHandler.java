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

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

abstract class OptimizationRequestHandler {
	protected OptimizationRequest request;
	private OptimizationResponse response;
	private int dim = -1;
	protected Algebra ALG = Algebra.DEFAULT;
	protected DoubleFactory1D F1 = DoubleFactory1D.dense;
	protected DoubleFactory2D F2 = DoubleFactory2D.dense;

	public void setOptimizationRequest(OptimizationRequest request) {
		this.request = request;
	}

	protected void setOptimizationResponse(OptimizationResponse response) {
		this.response = response;
	}

	public OptimizationResponse getOptimizationResponse() {
		return this.response;
	}

	/**
	 * Number of variables.
	 */
	protected final int getDim() {
		if (dim < 0) {
			dim = this.request.getF0().getDim();
		}
		return dim;
	}

	protected DoubleMatrix1D getInitialPoint() {
		return request.getInitialPoint();
	}

	protected DoubleMatrix1D getInitialLagrangian() {
		return request.getInitialLagrangian();
	}

	protected final int getMaxIteration() {
		return request.getMaxIteration();
	}

	protected final double getTolerance() {
		return request.getTolerance();
	}

	protected final double getToleranceFeas() {
		return request.getToleranceFeas();
	}

	protected final double getAlpha() {
		return request.getAlpha();
	}

	protected final double getBeta() {
		return request.getBeta();
	}

	protected final double getMu() {
		return request.getMu();
	}

	/**
	 * Objective function domain.
	 */
	protected boolean isInDomainF0(DoubleMatrix1D X) {
		double F0X = request.getF0().value(X);
		return !Double.isInfinite(F0X) && !Double.isNaN(F0X);
	}

	/**
	 * Objective function gradient at X.
	 * 
	 * @param gradF0X
	 */
	protected void getGradF0(DoubleMatrix1D X, DoubleMatrix1D result) {
		request.getF0().gradient(X, result);
	}

	/**
	 * Objective function hessian at X.
	 */
	protected DoubleMatrix2D getHessF0() {
		return request.getF0().hessian();
	}

	/**
	 * Inequality functions.
	 */
	protected SimpleLinearConstraint[] getFi() {
		return request.getFi();
	}
	
	/**
	 * Equality functions.
	 */
	protected SimpleLinearConstraint[] getFe() {
		return request.getA();
	}

	/**
	 * Inequality functions values at X.
	 */
	protected void getFi(DoubleMatrix1D x, DoubleMatrix1D result) {
		SimpleLinearConstraint[] fis = request.getFi();
		int n = fis.length;
		for (int i = 0; i < n; i++)
			result.setQuick(i, fis[i].value(x));
	}
}
