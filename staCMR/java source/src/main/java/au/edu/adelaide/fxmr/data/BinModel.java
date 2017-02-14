package au.edu.adelaide.fxmr.data;

import au.edu.adelaide.fxmr.model.bin.BinSolution;

/**
 * Stores a binomial model
 */
public class BinModel {
	/**
	 * nsubj x nvar array of elements each consisting of counts (hits, misses)
	 */
	private BinElement[][] data;
	/**
	 * Number of Variables in data
	 */
	private int nVar;
	/**
	 * Number of Subjects in data
	 */
	private int nSubj;

	public BinModel(int nSubj, int nVar) {
		this.data = new BinElement[nSubj][nVar];
		this.nSubj = nSubj;
		this.nVar = nVar;
	}

	public BinElement setElement(BinElement element, int subj, int var) {
		return data[subj][var] = element;
	}

	public BinElement[][] getData() {
		return data;
	}

	public int getnVar() {
		return nVar;
	}

	public int getnSubj() {
		return nSubj;
	}

	public BinModel resample() {
		BinModel newModel = new BinModel(nSubj, nVar);
		for (int s = 0; s < nSubj; s++)
			for (int v = 0; v < nVar; v++)
				newModel.data[s][v] = data[s][v].resample();
		return newModel;
	}

	public BinModel resample(BinSolution[] meanXBases) {
		BinModel newModel = new BinModel(nSubj, nVar);
		for (int s = 0; s < nSubj; s++)
			for (int v = 0; v < nVar; v++)
				newModel.data[s][v] = data[s][v].resample(meanXBases[s].getXStar(v));
		return newModel;
	}
}
