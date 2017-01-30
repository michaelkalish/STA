package au.edu.adelaide.fxmr.model;

import cern.colt.matrix.DoubleMatrix2D;

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
