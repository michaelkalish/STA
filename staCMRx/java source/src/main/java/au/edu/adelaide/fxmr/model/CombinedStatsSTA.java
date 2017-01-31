package au.edu.adelaide.fxmr.model;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

/**
 * Used by the GeneralModel, cumulate statistics for all data from multiple
 * StatsSTA objects (one for each between subjects condition)
 * 
 * 
 * 
 */
public class CombinedStatsSTA {
	private DoubleMatrix1D means;
	private int[] n;
	private double[] shrinkages;
	private double[] lmValues;
	private DoubleMatrix2D covariances;
	private DoubleMatrix2D regCovariances;
	private DoubleMatrix2D weights;
	private StatsSTA[] sstaSource;
	private int variable;

	public CombinedStatsSTA(StatsSTA[] ssta, int variable) {
		this.variable = variable;
		sstaSource = ssta;

		int nSSTA = ssta.length;

		DoubleMatrix1D[] means = new DoubleMatrix1D[nSSTA];
		DoubleMatrix2D[][] covariances = new DoubleMatrix2D[nSSTA][nSSTA];
		DoubleMatrix2D[][] regCovariances = new DoubleMatrix2D[nSSTA][nSSTA];
		DoubleMatrix2D[][] weights = new DoubleMatrix2D[nSSTA][nSSTA];

		shrinkages = new double[nSSTA];
		lmValues = new double[nSSTA];
		n = new int[nSSTA];

		for (int i = 0; i < nSSTA; i++) {
			StatsSTA cur = ssta[i];
			means[i] = cur.getMeans();
			covariances[i][i] = cur.getCovariance();
			regCovariances[i][i] = cur.getRegCovariance();
			weights[i][i] = cur.getWeights();
			shrinkages[i] = cur.getShrinkage();
			lmValues[i] = cur.getLMValue();
			n[i] = cur.getN();
		}
		this.means = DoubleFactory1D.dense.make(means);
		this.covariances = DoubleFactory2D.dense.compose(covariances);
		this.regCovariances = DoubleFactory2D.dense.compose(regCovariances);
		this.weights = DoubleFactory2D.dense.compose(weights);
	}

	public StatsSTA[] getSstaSource() {
		return sstaSource;
	}

	public DoubleMatrix1D getMeans() {
		return means;
	}

	public DoubleMatrix2D getCovariances() {
		return covariances;
	}

	public int[] getN() {
		return n;
	}

	public double[] getShrinkages() {
		return shrinkages;
	}

	public DoubleMatrix2D getRegCovariances() {
		return regCovariances;
	}

	public double[] getLmValues() {
		return lmValues;
	}

	public DoubleMatrix2D getWeights() {
		return weights;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("N\t");
		for (int i = 0; i < n.length; i++) {
			sb.append('\t');
			sb.append(n[i]);
		}
		sb.append('\n');
		sb.append("Shrinkage");
		for (int i = 0; i < shrinkages.length; i++) {
			sb.append('\t');
			sb.append(shrinkages[i]);
		}
		sb.append('\n');
		sb.append("Loftus & Masson");
		for (int i = 0; i < lmValues.length; i++) {
			sb.append('\t');
			sb.append(lmValues[i]);
		}
		return sb.toString();
	}

	public int getVariable() {
		return variable;
	}

	public Badness getWorst() {
		Badness worst = Badness.NONE;
		for (StatsSTA s : sstaSource) {
			worst = Badness.worst(s.getBad(), worst);
		}
		return worst;
	}
}
