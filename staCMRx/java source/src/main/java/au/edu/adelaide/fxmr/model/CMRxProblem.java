package au.edu.adelaide.fxmr.model;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.om.OMUtil;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import gnu.trove.set.hash.TIntHashSet;

public class CMRxProblem {
	private double[][] means;
	private double[][] meansOriginal;
	private int nVar;
	private int nVarOriginal;
	private int nCond;
	private DoubleMatrix2D model;
	private DoubleMatrix2D modelOriginal;
	/**
	 * Optional value that defines the sample size for each variable
	 */
	private DoubleMatrix2D[] weights;
	private DoubleMatrix2D[] weightsOriginal;
	/**
	 * Permitted
	 */
	private int[][] cv;

	/**
	 * Forbidden
	 */
	private TIntHashSet infeasZones;
	private HashSet<SimpleLinearConstraint>[] adj;
	/**
	 * Which rows from the original problem are we solving?
	 * 
	 * If null, the original model had no zero rows.
	 */
	private int[] notNullRows;
	private int[] nullRows;
	private double fAddWeightedMeans;
	private double[] weightedMeans;

	public CMRxProblem(double[][] means, DoubleMatrix2D[] weights, HashSet<SimpleLinearConstraint>[] adj,
			DoubleMatrix2D model) {
		this.means = meansOriginal = means;
		this.adj = adj;
		this.weights = weightsOriginal = weights;
		this.nVar = means.length;
		this.nVarOriginal = nVar;
		this.nCond = means[0].length;
		this.model = modelOriginal = model;
		setup();
	}

	public double[][] getMeansOriginal() {
		return meansOriginal;
	}

	public DoubleMatrix2D[] getWeightsOriginal() {
		return weightsOriginal;
	}

	public int getnVarOriginal() {
		return nVarOriginal;
	}

	@SuppressWarnings("unchecked")
	private void setup() {
		if (nVar != weights.length)
			throw new IllegalArgumentException("Weights has " + weights.length + " dimension(s) but nVar = " + nVar);

		if (adj != null && nVar != adj.length) {
			// If there are no adjs, we can infer what to do
			boolean empty = true;
			for (HashSet<SimpleLinearConstraint> a : adj)
				empty = empty && a.isEmpty();

			if (!empty)
				throw new IllegalArgumentException("Adjacency constraints has " + adj.length + " dimension(s) but nVar = " + nVar);

			// Recreate empty with correct size
			adj = new HashSet[nVar];
			for (int i = 0; i < nVar; i++)
				adj[i] = new HashSet<>();
		}

		if (nVar != means.length)
			throw new IllegalArgumentException("Means has " + means.length + " dimension(s) but nVar = " + nVar);

		if (model == null) {
			// If this is the case, there wont be null rows
			model = new DenseDoubleMatrix2D(nVar, 1).assign(1);
		} else {
			// Check for null rows in model
			int rows = model.rows();
			int cols = model.columns();
			notNullRows = new int[rows];
			nullRows = new int[rows];
			int nNotNull = 0;
			int nNull = 0;

			for (int r = 0; r < rows; r++) {
				boolean notNull = false;
				for (int c = 0; c < cols; c++) {
					if (model.getQuick(r, c) != 0) {
						notNull = true;
						break;
					}
				}
				if (notNull)
					notNullRows[nNotNull++] = r;
				else
					nullRows[nNull++] = r;
			}
			if (nNull > 0) {
				// Resize problem and calculate weighted means to be inserted
				// back into solution
				int[] tmp = new int[nNotNull];
				System.arraycopy(notNullRows, 0, tmp, 0, nNotNull);
				notNullRows = tmp;

				tmp = new int[nNull];
				System.arraycopy(nullRows, 0, tmp, 0, nNull);
				nullRows = tmp;

				double[][] subMeans = new double[nNotNull][];
				DoubleMatrix2D[] subWeights = new DoubleMatrix2D[nNotNull];
				HashSet<SimpleLinearConstraint>[] subAdj = null;
				if (adj != null)
					subAdj = new HashSet[nNotNull];
				// Create sub-problem
				try {
					for (int r = 0; r < nNotNull; r++) {
						subMeans[r] = means[notNullRows[r]];
						subWeights[r] = weights[notNullRows[r]];
						if (adj != null)
							subAdj[r] = adj[notNullRows[r]];
					}
				} catch (ArrayIndexOutOfBoundsException e) {
					throw new IllegalArgumentException("Model dimension does not match the input data");
				}

				// calculate weighted means
				weightedMeans = new double[nNull];
				fAddWeightedMeans = 0;
				for (int r = 0; r < nNull; r++) {
					int row = nullRows[r];
					DoubleMatrix1D w = DoubleFactory2D.dense.diagonal(weights[row]);
					double[] m = means[row];
					double z = 0;
					for (int i = 0; i < nCond; i++)
						z += m[i] * w.getQuick(i);
					z /= w.zSum();
					weightedMeans[r] = z;

					for (int i = 0; i < nCond; i++) {
						double d = m[i] - z;
						fAddWeightedMeans += d * w.getQuick(i) * d;
					}
				}

				// As far as CMRxSolver is concerned, nVar is smaller
				nVar = nNotNull;
				model = model.viewSelection(notNullRows, null);
				means = subMeans;
				weights = subWeights;
				adj = subAdj;
			} else {
				nullRows = null;
			}
		}

		DoubleMatrix2D checked = OMUtil.checkRank(model);
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

		// Ensure adj exists
		if (adj == null || adj.length != nVar) {
			adj = new HashSet[nVar];
			for (int i = 0; i < nVar; i++)
				adj[i] = new HashSet<>();
		}
	}

	// public CMRxProblem(CombinedStatsSTA[] stats, int[][] rangeSet,
	// DoubleMatrix2D model) {
	// nVar = stats.length;
	// nCond = stats[0].getMeans().size();
	// means = new double[nVar][];
	// weights = new DoubleMatrix2D[nVar];
	// for (int i = 0; i < nVar; i++) {
	// means[i] = stats[i].getMeans().toArray();
	// weights[i] = stats[i].getWeights();
	// }
	// this.adj = adj;
	// this.model = model;
	// setup();
	// }

	public double[][] getMeans() {
		return means;
	}

	public int getNVar() {
		return nVar;
	}

	public int getNCond() {
		return nCond;
	}

	public DoubleMatrix2D[] getWeights() {
		return weights;
	}

	public TIntHashSet getInfeasZones() {
		return infeasZones;
	}

	public HashSet<SimpleLinearConstraint>[] getAdj() {
		return adj;
	}

	public String toString() {
		return "nVar  = " + nVar + ", nCond = " + nCond;
	}

	public String toMatlab() {
		StringBuilder sb = new StringBuilder();
		sb.append("ysJ={");
		for (int i = 0; i < nVar; i++) {
			sb.append("struct(");
			sb.append("'means',");
			sb.append(Arrays.toString(means[i]));
			sb.append(",'weights',[");
			for (int row = 0; row < weights[i].rows(); row++) {
				sb.append(weights[i].get(row, 0));
				for (int col = 1; col < weights[i].columns(); col++) {
					sb.append(',');
					sb.append(weights[i].get(row, col));
				}
				if (row != weights[i].rows() - 1)
					sb.append(';');
			}
			sb.append("]) ");
		}
		sb.append("};");

		sb.append("\nE={");
		if (adj != null) {
			for (HashSet<SimpleLinearConstraint> a : adj) {
				for (SimpleLinearConstraint slc : a) {
					sb.append("[" + (slc.getPosIndex() + 1) + "," + (slc.getNegIndex() + 1) + "]");
					sb.append(",");
				}
			}
		}
		sb.append("};");
		return sb.toString();
	}

	public DoubleMatrix2D getModel() {
		return model;
	}

	/**
	 * Get the covectors that are permitted, based on the model
	 * 
	 * @return
	 */

	public int[][] getCv() {
		return cv;
	}

	/**
	 * If this problem has null rows, insert the weighted means into the
	 * solution and return.
	 * 
	 * Otherwise, create a new solution and return it.
	 * 
	 * @param fBar
	 * @param xBar
	 * @param iter
	 * @param seconds
	 * @param fBarReductions
	 * @param free
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public CMRSolution createSolution(double fBar, double[][] xBar, List<CMRIter> iter, double seconds,
			HashSet<SimpleLinearConstraint>[] adjBar, int mrCalls, int fBarReductions) {
		if (nullRows == null)
			return new CMRSolution(fBar, xBar, iter, seconds, adjBar, mrCalls, fBarReductions);
		double[][] finalXBar = new double[nullRows.length + xBar.length][];
		for (int r = 0; r < notNullRows.length; r++) {
			finalXBar[notNullRows[r]] = xBar[r];
		}

		for (int r = 0; r < nullRows.length; r++) {
			double[] tmp = new double[nCond];
			Arrays.fill(tmp, weightedMeans[r]);
			finalXBar[nullRows[r]] = tmp;
		}

		HashSet<SimpleLinearConstraint>[] finalAdjBar = null;
		if (adjBar != null) {
			finalAdjBar = new HashSet[nullRows.length + xBar.length];
			for (int r = 0; r < notNullRows.length; r++)
				finalAdjBar[notNullRows[r]] = adjBar[r];
		}

		return new CMRSolution(fBar + fAddWeightedMeans, finalXBar, iter, seconds, finalAdjBar, mrCalls, fBarReductions);
	}

	public boolean hasBadCondition() {
		for (DoubleMatrix2D w : weights) {
			if (new SingularValueDecomposition(w).cond() > 1e25)
				return true;
		}
		return false;
	}

	public DoubleMatrix2D getModelOriginal() {
		return modelOriginal;
	}

	public double getfAddWeightedMeans() {
		return fAddWeightedMeans;
	}
}
