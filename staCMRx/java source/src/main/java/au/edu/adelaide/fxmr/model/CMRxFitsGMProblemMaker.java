package au.edu.adelaide.fxmr.model;

import java.util.HashSet;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.data.Subject;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Allows MATLAB to create CMRxFitsGMProblems
 * 
 * 
 */
public class CMRxFitsGMProblemMaker extends CMRxProblemMaker {
	private GeneralModel gm;
	private double shrink;

	public void setShrink(double shrink) {
		this.shrink = shrink;
	}

	public GeneralModel getGm() {
		return gm;
	}

	private int subjNum = 1;

	public void setGM(double[][] data) {
		gm = new GeneralModel();
		for (double[] row : data)
			gm.addData(new Subject(row));
	}

	public void addCell(int group, int var, double[][] data) {
		if (gm == null)
			gm = new GeneralModel();

		for (double[] row : data)
			gm.addData(new Subject(subjNum++, group, var, row));
	}

	public Fits solve(int nSample, int proc) {
		return solve(nSample, proc, false);
	}

	@SuppressWarnings("unchecked")
	public Fits solve(int nSample, int proc, boolean cheapP) {
		HashSet<SimpleLinearConstraint>[] dMatAs = null;
		if (matAs.size() > 0) {
			dMatAs = new HashSet[matAs.size()];
			matAs.toArray(dMatAs);
		}

		return new CMRxGMFits(nSample, gm, shrink, new DenseDoubleMatrix2D(model), dMatAs, proc, cheapP, false);
	}

	@SuppressWarnings("unchecked")
	public Fits solve(int nSample, int proc, boolean cheapP, boolean onlySTAMR) {
		HashSet<SimpleLinearConstraint>[] dMatAs = null;
		if (matAs.size() > 0) {
			dMatAs = new HashSet[matAs.size()];
			matAs.toArray(dMatAs);
		}

		if (model == null) {
			// This implies the CMRxFitsProblem will be used in non-coupled
			// mode.
			model = new double[][]{{1}};
		}

		return new CMRxGMFits(nSample, gm, shrink, new DenseDoubleMatrix2D(model), dMatAs, proc, cheapP, onlySTAMR);
	}
}
