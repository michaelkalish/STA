package au.edu.adelaide.fxmr.model.bin;

import au.edu.adelaide.fxmr.data.BinModel;

public class BinBaseProblem {
	private int[][] rangeSet;
	private BinModel model;

	public BinBaseProblem(BinModel model, int[][] rangeSet) {
		this.model = model;
		this.rangeSet = rangeSet;
	}

	public int[][] getRangeSet() {
		return rangeSet;
	}

	public BinModel getModel() {
		return model;
	}
}
