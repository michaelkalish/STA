package au.edu.adelaide.fxmr.model.bin;

import java.util.ArrayList;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;

public class BinCMRProblemMaker {
	private ArrayList<Object> rangeSet = new ArrayList<>();
	private BinModel model;

	public BinCMRProblemMaker(int nSubj, int nVar) {
		model = new BinModel(nSubj, nVar);
	}

	public void setElement(int subj, int var, int[][] values) {
		model.setElement(new BinElement(values), subj, var);
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

		for (int i = 0; i < nSubj; i++)
			problems[i] = new BinProblem(model, i, makeRangeSet());

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
		return new BinBaseProblem(model, makeRangeSet());
	}
}
