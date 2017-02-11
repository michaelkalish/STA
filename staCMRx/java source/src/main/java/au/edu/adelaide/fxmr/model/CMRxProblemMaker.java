package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;
import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Allows MATLAB to create CMRxProblems
 * 
 * 
 */
public class CMRxProblemMaker {
	protected ArrayList<Object> means = new ArrayList<>();
	protected ArrayList<Object> weights = new ArrayList<>();
	protected ArrayList<HashSet<SimpleLinearConstraint>> matAs = new ArrayList<>();
	protected double[][] model;

	public CMRxProblemMaker() {
	}

	public void setModel(double[][] model) {
		this.model = model;
	}

	/**
	 * Specific method for R that allows a vector to be set (interpreted as a
	 * matrix)
	 * 
	 * @param model
	 */
	public void setModel(double[] model) {
		int n = model.length;
		this.model = new double[n][];
		for (int i = 0; i < n; i++)
			this.model[i] = new double[] { model[i] };
	}

	public void addMeanArray(double[] means) {
		this.means.add(means);
	}

	public void addWeightArray(double[][] weight) {
		this.weights.add(weight);
	}

	/**
	 * Add adjacency matrix as a full matrix
	 * 
	 * @param adj
	 */
	public void addAdj(double[][] adj) {
		HashSet<SimpleLinearConstraint> curMat = new HashSet<>();

		for (int i = 0; i < adj.length; i++)
			// assume square!
			for (int j = 0; j < adj.length; j++)
				if (adj[i][j] > 0)
					curMat.add(new SimpleLinearConstraint(i, j));

		matAs.add(curMat);
	}

	public void dupeAdj(int nvar) {
		if (matAs.isEmpty())
			throw new IllegalArgumentException("No adjacency matrix exists to copy");
		while (matAs.size() < nvar)
			matAs.add(matAs.get(0));
	}

	public int initAdj() {
		HashSet<SimpleLinearConstraint> curMat = new HashSet<>();
		matAs.add(curMat);
		return matAs.size() - 1;
	}

	/**
	 * Add adjacencies as a matlab cell of arrays. This generally requires
	 * initAdj to have been called beforehand
	 * 
	 * @return
	 */
	public void addAdj(int nCond, int index, int[] e) {
		HashSet<SimpleLinearConstraint> curMat = matAs.get(index);

		int rl = e.length;
		if (rl >= 2) {
			// New code, fewer constraints
			int rlm1 = rl - 1;
			for (int i = 0; i < rlm1; i++)
				// Take away 1 to convert from matlab 1-index to java 0-index
				curMat.add(new SimpleLinearConstraint(e[i] - 1, e[i + 1] - 1));
		}
	}

	@SuppressWarnings("unchecked")
	public CMRxProblem getProblem() {
		if (means.size() != weights.size())
			return null;

		int nvar = weights.size();
		DoubleMatrix2D[] dWeights = new DoubleMatrix2D[nvar];
		double[][] dMeans = new double[means.size()][];

		for (int i = 0; i < nvar; i++) {
			dWeights[i] = new DenseDoubleMatrix2D((double[][]) weights.get(i));
			dMeans[i] = (double[]) means.get(i);
		}

		HashSet<SimpleLinearConstraint>[] dMatAs = null;
		if (matAs.size() > 0) {
			// @SuppressWarnings("unchecked")
			dMatAs = new HashSet[nvar];
			matAs.toArray(dMatAs);
		}

		return new CMRxProblem(dMeans, dWeights, dMatAs, new DenseDoubleMatrix2D(model));
	}
}
