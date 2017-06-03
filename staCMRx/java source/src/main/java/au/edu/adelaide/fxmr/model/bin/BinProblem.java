package au.edu.adelaide.fxmr.model.bin;

import au.edu.adelaide.fxmr.data.BinModel;
import cern.colt.matrix.DoubleMatrix2D;

public class BinProblem extends BinBaseProblem {
	private int subjectIndex;

	public BinProblem(BinModel model, int subjectIndex, int[][] rangeSet, DoubleMatrix2D cmrModel) {
		super(model, rangeSet, cmrModel);
		this.subjectIndex = subjectIndex;
	}

	public BinProblem(BinModel model, int subjectIndex, int[][] rangeSet) {
		super(model, rangeSet);
		this.subjectIndex = subjectIndex;
	}

	public int getSubjectIndex() {
		return subjectIndex;
	}
}
