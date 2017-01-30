package au.edu.adelaide.fxmr.model.bin;

import au.edu.adelaide.fxmr.data.BinModel;

public class BinProblem extends BinBaseProblem {
	private int subjectIndex;

	public BinProblem(BinModel model, int subjectIndex, int[][] rangeSet) {
		super(model, rangeSet);
		this.subjectIndex = subjectIndex;
	}

	public int getSubjectIndex() {
		return subjectIndex;
	}
}
