package au.edu.adelaide.fxmr.model.bin;

import java.util.ArrayList;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class BinCMRProblemMaker {
	private ArrayList<Object> rangeSet = new ArrayList<>();
	private BinModel model;
	/**
	 * Optional model
	 */
	protected double[][] cmrModel;

	public BinCMRProblemMaker(int nSubj, int nVar) {
		model = new BinModel(nSubj, nVar);
	}

	public void setElement(int subj, int var, int[][] values) {
		model.setElement(new BinElement(values), subj, var);
	}

	public void setModel(double[][] model) {
		this.cmrModel = model;
	}

	/**
	 * Specific method for R that allows a vector to be set (interpreted as a
	 * matrix)
	 * 
	 * @param model
	 */
	public void setModel(double[] model) {
		int n = model.length;
		this.cmrModel = new double[n][];
		for (int i = 0; i < n; i++)
			this.cmrModel[i] = new double[] { model[i] };
	}

	/**
	 * Will convert the single array into a 2D array!
	 */
	public void setElement(int subj, int var, int[] values) {
		if (values.length % 2 != 0)
			throw new IllegalArgumentException("Values must have an even number of elements");
		int n = values.length / 2;
		int[][] gv = new int[2][n];
		for (int i = 0; i < n; i++) {
			gv[0][i] = values[i];
			gv[1][i] = values[i + n];
		}

		model.setElement(new BinElement(gv), subj, var);
	}

	public void addRangeSet(int[] r) {
		for (int i = 0; i < r.length; i++)
			r[i]--;

		this.rangeSet.add(r);
	}

	public BinProblem[] getProblems() {
		int nSubj = model.getnSubj();
		BinProblem[] problems = new BinProblem[nSubj];

		DoubleMatrix2D cmrModelMatrix = null;
		if (cmrModel != null)
			cmrModelMatrix = new DenseDoubleMatrix2D(cmrModel);

		for (int i = 0; i < nSubj; i++) {
			problems[i] = new BinProblem(model, i, makeRangeSet(), cmrModelMatrix);
		}

		return problems;
	}

	private int[][] makeRangeSet() {
		if (rangeSet.isEmpty())
			return null;

		int[][] dRangeSet = new int[rangeSet.size()][];
		for (int i = 0; i < rangeSet.size(); i++)
			dRangeSet[i] = (int[]) rangeSet.get(i);

		return dRangeSet;
	}

	public BinBaseProblem getBaseProblem() {
		if (cmrModel == null)
			return new BinBaseProblem(model, makeRangeSet());
		return new BinBaseProblem(model, makeRangeSet(), new DenseDoubleMatrix2D(cmrModel));
	}
}
