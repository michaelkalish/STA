package au.edu.adelaide.fxmr.model.util;

import java.util.Random;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.jet.math.Functions;

/**
 * Based on simple wikipedia algorithm
 * 
 * <a href=
 * "https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution">
 * here</a>
 * 
 * 
 */
public class MultivariateNormalDistribution {
	private DoubleMatrix2D L;
	private DoubleMatrix1D tmpRands;
	private DoubleMatrix1D tmpResult;
	private Random rand;
	private double[] means;
	private DoubleMatrix1D meansVec;

	public MultivariateNormalDistribution(double[] means, DoubleMatrix2D cov) {
		this(means, cov, System.currentTimeMillis());
	}

	public MultivariateNormalDistribution(double[] means, DoubleMatrix2D cov, long seed) {
		int n = means.length;
		meansVec = new DenseDoubleMatrix1D(means);
		CholeskyDecomposition chol = new CholeskyDecomposition(cov);

		this.rand = new Random(seed);
		this.L = chol.getL();
		this.tmpRands = new DenseDoubleMatrix1D(n);
		this.tmpResult = new DenseDoubleMatrix1D(n);
		this.means = means;
	}

	public double[] sample(double[] ans) {
		int n = means.length;
		for (int i = 0; i < n; i++)
			tmpRands.setQuick(i, rand.nextGaussian());
		L.zMult(tmpRands, tmpResult);

		for (int i = 0; i < n; i++)
			ans[i] = tmpResult.getQuick(i) + means[i];

		return ans;
	}

	public void fill(DoubleMatrix2D ans) {
		int n = means.length;
		int rows = ans.rows();
		for (int r = 0; r < rows; r++) {
			DoubleMatrix1D row = ans.viewRow(r);

			for (int i = 0; i < n; i++)
				tmpRands.setQuick(i, rand.nextGaussian());
			L.zMult(tmpRands, row);
			row.assign(meansVec, Functions.plus);
		}
	}
}
