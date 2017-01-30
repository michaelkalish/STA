package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;
import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Allows MATLAB to create CMRxFitsProblems
 * 
 * 
 */
public class CMRxFitsProblemMaker extends CMRxProblemMaker {
	private int[] n;
	private ArrayList<DoubleMatrix2D> covs = new ArrayList<>();

	public void setN(int[] n) {
		this.n = n;
	}

	public void addCov(double[][] cov) {
		covs.add(new DenseDoubleMatrix2D(cov));
	}

	@SuppressWarnings("unchecked")
	public CMRxFitsProblem getProblem() {
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

		DoubleMatrix2D[] cov = new DoubleMatrix2D[covs.size()];
		covs.toArray(cov);

		if (model == null) {
			// This implies the CMRxFitsProblem will be used in non-coupled
			// mode.
			model = new double[][]{{1}};
		}

		return new CMRxFitsProblem(dMeans, dWeights, dMatAs, new DenseDoubleMatrix2D(model), n, cov);
	}
}
