package au.edu.adelaide.fxmr.model;

import java.util.HashSet;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.data.Subject;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Allows MATLAB to create CMRxProblems
 * 
 * 
 */
public class CMRxGMProblemMaker extends CMRxProblemMaker {
	private GeneralModel gm;
	private double shrink;
	private int subjNum = 1;
	private double[][] shrinkages;

	public CMRxGMProblemMaker() {
	}

	public void setShrink(double shrink) {
		this.shrink = shrink;
	}

	public GeneralModel getGm() {
		return gm;
	}

	public void setGM(double[][] data) {
		gm = new GeneralModel();
		for (double[] row : data)
			gm.addData(new Subject(row));
	}

	public void addCell(int cond, int var, double[][] data) {
		if (gm == null)
			gm = new GeneralModel();

		for (double[] row : data)
			gm.addData(new Subject(subjNum++, cond, var, row));
	}

	@SuppressWarnings("unchecked")
	public CMRxProblem getProblem() {
		CombinedStatsSTA[] stats = shrink < 0 ? gm.calcStats() : gm.calcStats(shrink);

		int nvar = stats.length;
		shrinkages = new double[nvar][];

		DoubleMatrix2D[] dWeights = new DoubleMatrix2D[nvar];
		double[][] dMeans = new double[nvar][];

		for (int i = 0; i < nvar; i++) {
			shrinkages[i] = stats[i].getShrinkages();
			dWeights[i] = stats[i].getWeights();
			dMeans[i] = stats[i].getMeans().toArray();
		}

		HashSet<SimpleLinearConstraint>[] dMatAs = null;
		if (matAs.size() > 0) {
			// @SuppressWarnings("unchecked")
			dMatAs = new HashSet[nvar];
			matAs.toArray(dMatAs);
		}

		return new CMRxProblem(dMeans, dWeights, dMatAs, new DenseDoubleMatrix2D(model));
	}

	public double[][] getShrinkages() {
		return shrinkages;
	}
}
