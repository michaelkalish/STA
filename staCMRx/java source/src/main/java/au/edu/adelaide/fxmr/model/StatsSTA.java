package au.edu.adelaide.fxmr.model;

import javax.sound.midi.SysexMessage;

import au.edu.adelaide.fxmr.model.mr.MRUtil;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;

/**
 * This class calculates a series of stats for a given Nsubjects x Nconditions
 * matrix
 * 
 * Once the contructor is called, the resulting object has means, covariance
 * (std and shrunk), weights, shrinkage and the Loftus & Masson within subjects
 * error based on Mean Square residual.
 * 
 * 
 */
public class StatsSTA {
	private DoubleMatrix1D means;
	private int n;
	private int[][] nMat;
	private int nWithinCond;
	private double shrinkage;
	private DoubleMatrix2D covariance;
	private DoubleMatrix2D regCovariance;
	private DoubleMatrix2D weights;
	private double lmValue;
	private boolean autoShrink = false;
	private Badness bad = Badness.NONE;

	public StatsSTA(DoubleMatrix2D data) {
		this.autoShrink = true;
		calculate(data);
	}

	public StatsSTA(DoubleMatrix2D data, double shrinkage) {
		this.shrinkage = shrinkage;
		this.autoShrink = false;
		calculate(data);
	}

	private void calculate(DoubleMatrix2D data) {
		int nCol = nWithinCond = data.columns();
		means = MRUtil.colMeans(data);
		n = data.rows();

		nMat = new int[nCol][nCol];
		for (int c1 = 0; c1 < nCol; c1++) {
			DoubleMatrix1D col1 = data.viewColumn(c1);
			for (int c2 = c1; c2 < nCol; c2++) {
				DoubleMatrix1D col2 = data.viewColumn(c2);
				int count = 0;
				for (int r = 0; r < n; r++)
					if (!Double.isNaN(col1.get(r)) && !Double.isNaN(col2.get(r)))
						count++;
				nMat[c1][c2] = nMat[c2][c1] = count;
			}
		}

		if (n > 1) {
			covariance = MRUtil.covariance(data);
			ShrinkDiagonal sd = new ShrinkDiagonal(data);

			EigenvalueDecomposition eig = new EigenvalueDecomposition(covariance);
			double cond = MRUtil.cond(eig);

			boolean posDef = true;
			for (double s : eig.getRealEigenvalues().toArray())
				if (s <= 0) {
					posDef = false;
					break;
				}

			if (cond < 1e6 && posDef) {
				// Normal behaviour
				if (autoShrink)
					sd.shrink();
				else
					sd.shrink(shrinkage);
			} else {
				sd.shrink(1);
				bad = Badness.DIAGONALISED;
			}
			shrinkage = sd.getShrinkage();
			regCovariance = sd.getResult();

			try {
				eig = new EigenvalueDecomposition(regCovariance);
			} catch (ArrayIndexOutOfBoundsException e) {
				System.out.println(regCovariance);
				System.out.println(data);
				e.printStackTrace();
			}
			cond = MRUtil.cond(eig);
			posDef = true;
			for (double s : eig.getRealEigenvalues().toArray())
				if (s <= 0) {
					posDef = false;
					break;
				}

			if (cond > 1e6 || !posDef) {
				MRUtil.forceDiagonal(regCovariance);
				bad = Badness.DIAGONALS_HAD_ZEROS;
			}

			weights = Algebra.ZERO.inverse(regCovariance);
			for (int c = 0; c < nCol; c++) {
				for (int r = 0; r < nCol; r++) {
					weights.setQuick(r, c, weights.getQuick(r, c) * nMat[r][c]);
				}
			}

			DoubleMatrix2D check = weights.copy();
			check.assign(weights.viewDice(), Functions.minus);
			double normCheck = Algebra.ZERO.normInfinity(check);
			if (normCheck > 1e-6) {
				// This means the regCovariance matrix isn't full rank!
				throw new IllegalArgumentException("Resulting weights matrix is non-symmetric");
			}
		} else {
			regCovariance = covariance = new DenseDoubleMatrix2D(nCol, nCol);
			weights = DoubleFactory2D.dense.identity(nCol);
			shrinkage = -1;
		}

		DoubleMatrix1D rowMeans = MRUtil.rowMeans(data);
		double mean = 0;
		int nM = 0;
		for (int row = n; --row >= 0;)
			for (int column = nCol; --column >= 0;) {
				double d = data.getQuick(row, column);
				if (!Double.isNaN(d)) {
					mean += d;
					nM++;
				}
			}
		mean /= nM;

		DoubleMatrix2D ya = new DenseDoubleMatrix2D(n, nCol);

		for (int row = n; --row >= 0;)
			for (int column = nCol; --column >= 0;) {
				double d = data.getQuick(row, column);
				if (!Double.isNaN(d))
					ya.setQuick(row, column,
							d - means.getQuick(column) - rowMeans.getQuick(row) + mean);
			}
		ya.assign(Functions.square);
		double ss = ya.zSum();
		double df = (n - 1) * (nCol - 1);
		lmValue = ss / df;
	}

	public int getNWithinCond() {
		return nWithinCond;
	}

	public DoubleMatrix2D getRegCovariance() {
		return regCovariance;
	}

	public DoubleMatrix1D getMeans() {
		return means;
	}

	public int getN() {
		return n;
	}

	public double getShrinkage() {
		return shrinkage;
	}

	public DoubleMatrix2D getCovariance() {
		return covariance;
	}

	public double getLMValue() {
		return lmValue;
	}

	public boolean isAutoShrink() {
		return autoShrink;
	}

	public int[][] getnMat() {
		return nMat;
	}

	public Badness getBad() {
		return bad;
	}

	public DoubleMatrix2D getWeights() {
		return weights;
	}

}
