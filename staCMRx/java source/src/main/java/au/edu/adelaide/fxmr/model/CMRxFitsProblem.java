package au.edu.adelaide.fxmr.model;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * A CMRxProblem with the extra information required to fit the solution
 * (covariance matrices and size of the sample)
 * 
 * 
 */
public class CMRxFitsProblem extends CMRxProblem {
	private int[] n;
	private DoubleMatrix2D[] cov;

	public CMRxFitsProblem(double[][] means, DoubleMatrix2D[] weights, HashSet<SimpleLinearConstraint>[] adj,
			DoubleMatrix2D model, int[] n, DoubleMatrix2D[] cov) {
		super(means, weights, adj, model);
		this.n = n;
		this.cov = cov;
	}

	public int[] getN() {
		return n;
	}

	public DoubleMatrix2D[] getCov() {
		return cov;
	}
}
