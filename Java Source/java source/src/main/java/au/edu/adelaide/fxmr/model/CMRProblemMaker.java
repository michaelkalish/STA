package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class CMRProblemMaker {
	private ArrayList<Object> means = new ArrayList<>();
	private ArrayList<Object> weights = new ArrayList<>();
	private ArrayList<Object> rangeSet = new ArrayList<>();

	public CMRProblemMaker() {
	}

	public void addMeanArray(double[] means) {
		this.means.add(means);
	}

	public void addWeightArray(double[][] weight) {
		this.weights.add(weight);
	}

	public void addRangeSet(int[] rangeSet) {
		this.rangeSet.add(rangeSet);
	}

	public CMRProblem getProblem() {
		if (means.size() != weights.size())
			return null;

		int[][] dRangeSet = new int[rangeSet.size()][];
		DoubleMatrix2D[] dWeights = new DoubleMatrix2D[weights.size()];
		double[][] dMeans = new double[means.size()][];

		for (int i = 0; i < weights.size(); i++) {
			dWeights[i] = new DenseDoubleMatrix2D((double[][]) weights.get(i));
			dMeans[i] = (double[]) means.get(i);
		}

		for (int i = 0; i < rangeSet.size(); i++) {
			dRangeSet[i] = (int[]) rangeSet.get(i);
		}

		return new CMRProblem(dMeans, dWeights, dRangeSet);
	}
}
