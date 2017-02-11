package au.edu.adelaide.fxmr.model.mr;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MRProblemMaker {
	private HashSet<SimpleLinearConstraint> curMat = new HashSet<>();
	private double[][] weights;
	private double[] y;

	public MRProblemMaker(double[] y) {
		this.y = y;
	}

	public void setWeightArray(double[][] weights) {
		this.weights = weights;
	}

	public void addConstraint(int pos, int neg) {
		// Note: -1s are because java is zero indexed, MATLAB/R isnt
		curMat.add(new SimpleLinearConstraint(pos - 1, neg - 1));
	}

	public void addConstraints(int[] constraints) {
		// Note: -1s are because java is zero indexed, MATLAB/R isnt
		for (int i = 1; i < constraints.length; i++) {
			curMat.add(new SimpleLinearConstraint(constraints[i - 1] - 1, constraints[i] - 1));
		}
	}

	public MRProblem getProblem() {
		return new MRProblem(y, new DenseDoubleMatrix2D(weights), curMat);
	}
}
