package au.edu.adelaide.fxmr.model;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class CMRProblem {
	private int[][] rangeSet;
	private double[][] means;
	private int nVar;
	private DoubleMatrix2D[] weights;

	public CMRProblem(double[][] means, DoubleMatrix2D[] weights, int[][] rangeSet) {
		this.means = means;
		this.rangeSet = rangeSet;
		this.weights = weights;
		this.nVar = means.length;
		forceSymetry();
	}

	private void forceSymetry() {
		for (DoubleMatrix2D w : weights) {
			DoubleMatrix2D check = w.copy();
			check.assign(w.viewDice(), Functions.minus);
			double normCheck = Algebra.ZERO.normInfinity(check);
			if (normCheck > 1e-14) {
				// Force symmetry
				check.assign(w);
				check.assign(w.viewDice(), Functions.plus);
				check.assign(Functions.div(2));
				w.assign(check);
			}
		}
	}

	public CMRProblem(CombinedStatsSTA[] stats, int[][] rangeSet) {
		nVar = stats.length;
		means = new double[nVar][];
		weights = new DoubleMatrix2D[nVar];
		for (int i = 0; i < nVar; i++) {
			means[i] = stats[i].getMeans().toArray();
			weights[i] = stats[i].getWeights();
		}
		this.rangeSet = rangeSet;
	}

	public int[][] getRangeSet() {
		return rangeSet;
	}

	public double[][] getMeans() {
		return means;
	}

	public int getNVar() {
		return nVar;
	}

	public DoubleMatrix2D[] getWeights() {
		return weights;
	}
}
