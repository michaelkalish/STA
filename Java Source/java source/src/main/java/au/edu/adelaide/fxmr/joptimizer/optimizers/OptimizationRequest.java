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

import au.edu.adelaide.fxmr.joptimizer.functions.QuadraticMultivariateRealFunction;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;

/**
 * Optimization problem. Setting the field's values you define an optimization
 * problem.
 * 
 * @see "S.Boyd and L.Vandenberghe, Convex Optimization"
 * @author alberto trivellato (alberto.trivellato@gmail.com)
 */
public class OptimizationRequest {

	/**
	 * Maximum number of iteration in the search algorithm. Not mandatory,
	 * default is provided.
	 */
	private int maxIteration = JOptimizer.DEFAULT_MAX_ITERATION;

	/**
	 * Tolerance for the minimum value. Not mandatory, default is provided.
	 * NOTE: as a golden rule, do not ask for more accuracy than you really
	 * need.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double tolerance = JOptimizer.DEFAULT_TOLERANCE;

	/**
	 * Tolerance for the constraints satisfaction. Not mandatory, default is
	 * provided. NOTE: as a golden rule, do not ask for more accuracy than you
	 * really need.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double toleranceFeas = JOptimizer.DEFAULT_FEASIBILITY_TOLERANCE;

	/**
	 * Tolerance for inner iterations in the barrier-method. NB: it makes sense
	 * only for barrier method. Not mandatory, default is provided. NOTE: as a
	 * golden rule, do not ask for more accuracy than you really need.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double toleranceInnerStep = JOptimizer.DEFAULT_TOLERANCE_INNER_STEP;

	/**
	 * Calibration parameter for line search. Not mandatory, default is
	 * provided.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double alpha = JOptimizer.DEFAULT_ALPHA;

	/**
	 * Calibration parameter for line search. Not mandatory, default is
	 * provided.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double beta = JOptimizer.DEFAULT_BETA;

	/**
	 * Calibration parameter for line search. Not mandatory, default is
	 * provided.
	 * 
	 * @see "Convex Optimization, p. 11.7.3"
	 */
	private double mu = JOptimizer.DEFAULT_MU;

	/**
	 * Acceptable tolerance for KKT system resolution. Not mandatory, default is
	 * provided.
	 */
	private double toleranceKKT = JOptimizer.DEFAULT_KKT_TOLERANCE;

	/**
	 * Should matrix rescaling be disabled? Rescaling is involved in LP
	 * presolving and in the solution of the KKT systems associated with the
	 * problem. It is an heuristic process, in some situations it could be
	 * useful to turn off this feature.
	 */
	private boolean rescalingDisabled = false;

	/**
	 * The objective function to minimize. Mandatory.
	 */
	private QuadraticMultivariateRealFunction f0;

	/**
	 * Feasible starting point for the minimum search. It must be feasible. Not
	 * mandatory.
	 */
	private DoubleMatrix1D initialPoint = null;

	/**
	 * Not-feasible starting point for the minimum search. It does not have to
	 * be feasible. This provide the possibility to give the algorithm a
	 * starting point even if it does not satisfies inequality constraints. The
	 * algorithm will search a feasible point starting from here. Not mandatory.
	 */
	private DoubleMatrix1D notFeasibleInitialPoint = null;

	/**
	 * Starting point for the Lagrangian multipliers. Must have the same
	 * dimension of the inequalities constraints array. Not mandatory, but very
	 * useful in some case.
	 */
	private DoubleMatrix1D initialLagrangian = null;

	/**
	 * Inequalities constraints array. Not mandatory.
	 * 
	 * @see "Convex Optimization, 11.1"
	 */
	private SimpleLinearConstraint[] fi;

	/**
	 * Equalities constraints matrix, in the simple form of xp - xn = 0.
	 * 
	 * The vector b = 0 in our case.
	 * 
	 * @see "Convex Optimization, 11.1"
	 */
	private SimpleLinearConstraint[] A;

	public int getMaxIteration() {
		return maxIteration;
	}

	public void setMaxIteration(int maxIteration) {
		this.maxIteration = maxIteration;
	}

	double getTolerance() {
		return tolerance;
	}

	public void setTolerance(double tolerance) {
		this.tolerance = tolerance;
	}

	public double getToleranceFeas() {
		return toleranceFeas;
	}

	public void setToleranceFeas(double toleranceFeas) {
		this.toleranceFeas = toleranceFeas;
	}

	public double getToleranceInnerStep() {
		return toleranceInnerStep;
	}

	public void setToleranceInnerStep(double toleranceInnerStep) {
		this.toleranceInnerStep = toleranceInnerStep;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public double getMu() {
		return mu;
	}

	public void setMu(double mu) {
		this.mu = mu;
	}


	public double getToleranceKKT() {
		return toleranceKKT;
	}

	public void setToleranceKKT(double toleranceKKT) {
		this.toleranceKKT = toleranceKKT;
	}

	public QuadraticMultivariateRealFunction getF0() {
		return f0;
	}

	public void setF0(QuadraticMultivariateRealFunction f0) {
		this.f0 = f0;
	}

	public DoubleMatrix1D getInitialPoint() {
		return initialPoint;
	}

	public void setInitialPoint(double[] initialPoint) {
		this.initialPoint = DoubleFactory1D.dense.make(initialPoint);
	}

	public DoubleMatrix1D getNotFeasibleInitialPoint() {
		return notFeasibleInitialPoint;
	}

	public void setNotFeasibleInitialPoint(double[] notFeasibleInitialPoint) {
		this.notFeasibleInitialPoint = DoubleFactory1D.dense.make(notFeasibleInitialPoint);
	}

	public DoubleMatrix1D getInitialLagrangian() {
		return initialLagrangian;
	}

	public void setInitialLagrangian(double[] initialLagrangian) {
		this.initialLagrangian = DoubleFactory1D.dense.make(initialLagrangian);
	}

	public SimpleLinearConstraint[] getFi() {
		return fi;
	}

	public void setFi(SimpleLinearConstraint[] fi) {
		this.fi = fi;
	}

	public void setRescalingDisabled(boolean rescalingDisabled) {
		this.rescalingDisabled = rescalingDisabled;
	}

	public boolean isRescalingDisabled() {
		return rescalingDisabled;
	}

	public SimpleLinearConstraint[] getA() {
		return A;
	}

	public void setA(SimpleLinearConstraint[] A) {
		this.A = A;
	}
}
