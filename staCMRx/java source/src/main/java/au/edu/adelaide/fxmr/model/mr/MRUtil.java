package au.edu.adelaide.fxmr.model.mr;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Statistic;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

/**
 * Various utility methods
 * 
 * 
 */
public class MRUtil {
	/**
	 * Calculate the means of each column
	 * 
	 * @param matrix
	 * @return
	 */
	public static DoubleMatrix1D colMeans(DoubleMatrix2D matrix) {
		int nCol = matrix.columns();
		int nRow = matrix.rows();
		final double[] ret = new double[nCol];

		for (int column = nCol; --column >= 0;) {
			int count = 0;
			double sum = 0;
			for (int r = 0; r < nRow; r++) {
				double v = matrix.get(r, column);
				if (!Double.isNaN(v)) {
					count++;
					sum += v;
				}
			}
			ret[column] = sum / count;
		}

		return new DenseDoubleMatrix1D(ret);
	}

	/**
	 * Calculate the means of each row
	 * 
	 * @param matrix
	 * @return
	 */
	public static DoubleMatrix1D rowMeans(DoubleMatrix2D matrix) {
		int nCol = matrix.columns();
		int nRow = matrix.rows();
		final double[] ret = new double[nRow];

		for (int r = nRow; --r >= 0;) {
			int count = 0;
			double sum = 0;
			for (int c = 0; c < nCol; c++) {
				double v = matrix.get(r, c);
				if (!Double.isNaN(v)) {
					count++;
					sum += v;
				}
			}
			ret[r] = sum / count;
		}

		return new DenseDoubleMatrix1D(ret);
	}

	// Note: this is slightly innacurate in some extreme cases...
	// public static DoubleMatrix2D covarianceOld(DoubleMatrix2D matrix) {
	// int rows = matrix.rows();
	// int columns = matrix.columns();
	// DoubleMatrix2D covariance = new DenseDoubleMatrix2D(columns, columns);
	//
	// DoubleMatrix1D[] cols = new DoubleMatrix1D[columns];
	// for (int i = columns; --i >= 0;)
	// cols[i] = matrix.viewColumn(i);
	//
	// for (int i = columns; --i >= 0;) {
	// DoubleMatrix1D colI = cols[i];
	// for (int j = i + 1; --j >= 0;) {
	// DoubleMatrix1D colJ = cols[j];
	// int count = 0;
	// double sumI = 0;
	// double sumJ = 0;
	// double sumOfProducts = 0;
	// for (int row = 0; row < rows; row++) {
	// double vi = colI.get(row);
	// double vj = colJ.get(row);
	// if (!Double.isNaN(vi) && !Double.isNaN(vj)) {
	// sumI += vi;
	// sumJ += vj;
	// count++;
	// sumOfProducts += vi * vj;
	// }
	// }
	//
	// double cov = count == 0 ? 0 : (sumOfProducts - sumI * sumJ / count) /
	// count;
	// covariance.setQuick(i, j, cov);
	// covariance.setQuick(j, i, cov); // symmetric
	// }
	// }
	// return covariance;
	// }

	public static DoubleMatrix2D covariance(DoubleMatrix2D matrix) {
		int rows = matrix.rows();
		int columns = matrix.columns();
		DoubleMatrix2D covariance = new DenseDoubleMatrix2D(columns, columns);

		DoubleMatrix1D[] cols = new DoubleMatrix1D[columns];
		for (int i = columns; --i >= 0;)
			cols[i] = matrix.viewColumn(i);

		for (int i = columns; --i >= 0;) {
			DoubleMatrix1D colI = cols[i];
			for (int j = i + 1; --j >= 0;) {
				DoubleMatrix1D colJ = cols[j];
				double meanJ = 0;
				double meanI = 0;
				int count = 0;

				for (int row = 0; row < rows; row++) {
					double vi = colI.get(row);
					double vj = colJ.get(row);
					if (!Double.isNaN(vi) && !Double.isNaN(vj)) {
						count++;
						meanI += vi;
						meanJ += vj;
					}
				}
				meanI /= count;
				meanJ /= count;

				double sumOfProducts = 0;
				for (int row = 0; row < rows; row++) {
					double vi = colI.get(row);
					double vj = colJ.get(row);
					if (!Double.isNaN(vi) && !Double.isNaN(vj))
						sumOfProducts += (vi - meanI) * (vj - meanJ);
				}

				double cov = count == 0 ? 0 : sumOfProducts / count;
				covariance.setQuick(i, j, cov);
				covariance.setQuick(j, i, cov); // symmetric
			}
		}
		return covariance;
	}

	/**
	 * Average diagonal * identity
	 * 
	 * Assumes input is square
	 * 
	 * @param matrix
	 */
	public static void forceDiagonal(DoubleMatrix2D matrix) {
		int rows = matrix.rows();

		double sum = 0;
		for (int i = rows; --i >= 0;)
			sum += matrix.getQuick(i, i);

		sum /= rows;
		matrix.assign(0);
		for (int i = rows; --i >= 0;)
			matrix.setQuick(i, i, sum);
	}

	/**
	 * Calculate the condition number from the eigenvalue decomposition. We
	 * assume the imaginary parts of the numbers are all zero (safe assumption
	 * since we're dealing with covariance matrices).
	 * 
	 * This is useful if we need the eigenvalues anyway.
	 * 
	 * @param eig
	 * @return
	 */
	public static double cond(EigenvalueDecomposition eig) {
		double max = 0;
		double min = Double.MAX_VALUE;

		DoubleMatrix1D v = eig.getRealEigenvalues();
		int n = v.size();
		for (int i = 0; i < n; i++) {
			double val = Math.abs(v.get(i));
			if (val > max)
				max = val;
			if (val < min)
				min = val;
		}

		return max / min;
	}
}
