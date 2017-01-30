package au.edu.adelaide.fxmr.model;

import au.edu.adelaide.fxmr.model.mr.MRUtil;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.Property;
import cern.jet.math.Functions;

/**
 * This class shrinks matrices toward a diagonal.
 * 
 * 
 * x (t*n): t iid observations on n random variables sigma (n*n): invertible
 * covariance matrix estimator
 * 
 * Shrinks towards diagonal matrix if shrink is specified, then this constant is
 * used for shrinkage
 * 
 * 
 *
 */
public class ShrinkDiagonal {
	/**
	 * Matrix to shrink
	 */
	private DoubleMatrix2D data;
	private DoubleMatrix2D result;
	private double shrinkage;

	public ShrinkDiagonal(DoubleMatrix2D data) {
		this.data = data;
	}

	/**
	 * Shrink matrix data by the optimal amount
	 * 
	 * @param shrinkage
	 */
	public void shrink() {
		int nRow = data.rows();
		int nCol = data.columns();

		// sample = dM' * dM / nRow
		DoubleMatrix2D sample = MRUtil.covariance(data);

		if (Property.ZERO.isDiagonal(sample)){
			shrinkage = 0;
			result = finalShrink(sample, null, nCol);
			return;
		}
		
		DoubleMatrix2D prior = new DenseDoubleMatrix2D(nCol, nCol);
		for (int i = nCol; --i >= 0;)
			prior.setQuick(i, i, sample.getQuick(i, i));

		DoubleMatrix1D means = MRUtil.colMeans(data);
		DoubleMatrix2D phiMat = data.copy();
		for (int c = 0; c < nCol; c++) {
			double m = means.get(c);
			for (int r = 0; r < nRow; r++) {
				if (Double.isNaN(phiMat.getQuick(r, c))) {
					phiMat.setQuick(r, c, 0);
				} else {
					double v = phiMat.getQuick(r, c) - m;
					phiMat.setQuick(r, c, v * v);
				}
			}
		}

		phiMat = phiMat.zMult(phiMat, null, 1.0 / nRow, 0, true, false);
		DoubleMatrix2D s2 = sample.copy();
		s2.assign(Functions.square);
		double phi = phiMat.assign(s2, Functions.minus).zSum();

		double rho = 0;
		for (int i = nCol; --i >= 0;)
			rho += phiMat.getQuick(i, i);

		double gamma = Algebra.ZERO.normF(sample.copy().assign(prior, Functions.minus));
		gamma *= gamma;

		double kappa = (phi - rho) / gamma;
		shrinkage = Math.max(0, Math.min(1, kappa / nRow));
		result = finalShrink(sample, prior, nCol);
	}

	private DoubleMatrix2D finalShrink(DoubleMatrix2D sample, DoubleMatrix2D prior, int t) {
		if (shrinkage == 0)
			return sample;
		if (shrinkage == 1)
			return prior;

		double oms = 1.0 - shrinkage;
		for (int row = t; --row >= 0;)
			for (int column = t; --column >= 0;)
				prior.setQuick(row, column,
						shrinkage * prior.getQuick(row, column) + oms * sample.getQuick(row, column));
		return prior;
	}

	/**
	 * Shrink matrix data by the specified amount
	 * 
	 * @param shrinkage
	 */
	public void shrink(double shrinkage) {
		this.shrinkage = shrinkage;
		int nCol = data.columns();
		DoubleMatrix2D sample = MRUtil.covariance(data);

		DoubleMatrix2D prior = new DenseDoubleMatrix2D(nCol, nCol);
		for (int i = nCol; --i >= 0;)
			prior.setQuick(i, i, sample.getQuick(i, i));

		result = finalShrink(sample, prior, nCol);
	}

	public double getShrinkage() {
		return shrinkage;
	}

	public DoubleMatrix2D getResult() {
		return result;
	}
}
