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
package au.edu.adelaide.fxmr.joptimizer.functions;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 */
public interface  TwiceDifferentiableMultivariateRealFunction {

	/**
	 * Evaluation of the function at point X.
	 */
	//public double value(double[] X);
	double value(DoubleMatrix1D x);
	
	/**
	 * Function gradient at point X.
	 */
	DoubleMatrix1D gradient(DoubleMatrix1D x);
	
	public void gradient(DoubleMatrix1D x, DoubleMatrix1D result);
	/**
	 * Function hessian at point X.
	 */
	public DoubleMatrix2D hessian();
	
	/**
	 * Dimension of the function argument.
	 */
	public int getDim();

	
}
