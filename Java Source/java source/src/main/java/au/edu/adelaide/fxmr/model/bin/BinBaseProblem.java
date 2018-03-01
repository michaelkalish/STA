package au.edu.adelaide.fxmr.model.bin;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.om.OMUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.set.hash.TIntHashSet;

public class BinBaseProblem {
	private int[][] rangeSet;
	private BinModel model;
	private int[][] cv;
	private DoubleMatrix2D cmrModel;
	private TIntHashSet infeasZones;

	/**
	 * No model defined - this is generally used for default model version.
	 * 
	 * @param model
	 * @param rangeSet
	 */
	public BinBaseProblem(BinModel model, int[][] rangeSet) {
		this.model = model;
		this.rangeSet = rangeSet;
	}

	public BinBaseProblem(BinModel model, int[][] rangeSet, DoubleMatrix2D cmrModel) {
		this.model = model;
		this.rangeSet = rangeSet;
		this.cmrModel = cmrModel;
		setup();
	}

	private void setup() {
		if (cmrModel == null) {
			// Default model
			cmrModel = new DenseDoubleMatrix2D(model.getnVar(), 1);
			cmrModel.assign(1);
		}

		DoubleMatrix2D checked = OMUtil.checkRank(cmrModel);
		int[][] vectors = OMUtil.vectors(checked);

		infeasZones = new TIntHashSet();
		if (vectors != null)
			infeasZones.addAll(OMUtil.zEncode(OMUtil.signAugment(vectors)));

		int[][] cvBase = OMUtil.signClosure(OMUtil.covectors(checked));
		int[] d = OMUtil.dim(cvBase);
		int dn = d.length;
		int max = -Integer.MAX_VALUE;
		int n = 0;
		for (int dv : d)
			if (dv > max) {
				max = dv;
				n = 1;
			} else if (dv == max) {
				n++;
			}

		cv = new int[n * 2][];
		int len = cvBase[0].length;
		int ci = 0;
		for (int i = 0; i < dn; i++) {
			if (d[i] == max) {
				cv[ci] = cvBase[i];
				cv[ci + n] = new int[len];
				for (int j = 0; j < len; j++)
					cv[ci + n][j] = -cvBase[i][j];
				ci++;
			}
		}
	}

	public int[][] getRangeSet() {
		return rangeSet;
	}

	public BinModel getModel() {
		return model;
	}

	/**
	 * Get the covectors, calculated from the CMRModel
	 * 
	 * @return
	 */
	public int[][] getCv() {
		if (cv == null)
			setup();

		return cv;
	}

	public TIntHashSet getInfeasZones() {
		if (infeasZones == null)
			setup();
		return infeasZones;
	}

	public DoubleMatrix2D getCmrModel() {
		return cmrModel;
	}

	public void setCmrModel(DoubleMatrix2D cmrModel) {
		this.cmrModel = cmrModel;
		this.cv = null;
		this.infeasZones = null;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();

		
		for (BinElement[] d : model.getData()) {
			sb.append("{");
			for (BinElement d2 : d) {
				sb.append("{");
				int[] h = d2.getHits();
				int[] m = d2.getMisses();
				int n = h.length;
				for (int i = 0; i < n; i++) {
					sb.append("{");
					sb.append(h[i]);
					sb.append(",");
					sb.append(m[i]);
					sb.append("},");
				}
				sb.append("},");
			}
			
		}
		sb.append("}");
		return sb.toString();
	}
	
	
	public String toMatlabString(){
		StringBuilder sb = new StringBuilder();

		
		for (BinElement[] d : model.getData()) {
			sb.append("{");
			for (BinElement d2 : d) {
				sb.append("[");
				int[] h = d2.getHits();
				int[] m = d2.getMisses();
				int n = h.length;
				for (int i = 0; i < n; i++) {
					sb.append("[");
					sb.append(h[i]);
					sb.append(",");
					sb.append(m[i]);
					sb.append("];");
				}
				sb.append("],");
			}
			
		}
		sb.append("}");
		return sb.toString();
	}
}