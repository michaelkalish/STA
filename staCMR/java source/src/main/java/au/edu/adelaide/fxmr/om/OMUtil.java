package au.edu.adelaide.fxmr.om;

import java.util.Arrays;
import java.util.Iterator;

import au.edu.adelaide.fxmr.model.CMRUtil;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.LUDecomposition;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;

public class OMUtil {

	/**
	 * returns reduced form of a if not full rank
	 * 
	 * assume a.rows() >= a.columns()
	 * 
	 * @param a
	 * @return
	 */
	public static DoubleMatrix2D checkRank(DoubleMatrix2D a) {
		int rank = Algebra.ZERO.rank(a);
		if (rank < a.columns()) {
			LUDecomposition lu = new LUDecomposition(a.viewDice());
			DoubleMatrix2D b = lu.getU().viewDice();

			int[] cols = new int[rank];
			for (int i = 0; i < rank; i++)
				cols[i] = i;

			return b.viewSelection(null, cols);
		} else {
			return a;
		}
	}

	public static DoubleMatrix1D minors(DoubleMatrix2D x) {
		if (x.columns() > x.rows())
			x = x.viewDice();

		int nChoosekLen = (int) CombinatoricsUtils.binomialCoefficient(x.rows(), x.columns());
		Iterator<int[]> iter = CombinatoricsUtils.combinationsIterator(x.rows(), x.columns());

		DoubleMatrix1D m = new DenseDoubleMatrix1D(nChoosekLen);
		int i = 0;
		while (iter.hasNext()) {
			DoubleMatrix2D subMatrix = x.viewSelection(iter.next(), null);
			double det = Algebra.ZERO.det(subMatrix);
			if (Math.abs(det) > 1e-10)
				m.setQuick(i, det);
			i++;
		}

		return m;
	}

	/**
	 * returns orthogonal complement of n x m matrix x
	 * 
	 * @param x
	 * @return
	 */
	public static DoubleMatrix2D orthocom(DoubleMatrix2D x) {
		// TODO: probably dont need checkrank...
		x = checkRank(x);
		if (x.rows() > x.columns()) {
			// if m = n-1 then orthog complement is normal
			if (x.rows() == x.columns() + 1) {
				DoubleMatrix1D subMat = genCrossProd(x);
				DoubleMatrix2D ret = new DenseDoubleMatrix2D(1, subMat.size());
				ret.viewRow(0).assign(subMat);
				return ret;
			} else {
				// otherwise find a good basis for x
				// int nChoosekLen = (int)
				// CombinatoricsUtils.binomialCoefficient(x.rows(),
				// x.columns());
				Iterator<int[]> iter = CombinatoricsUtils.combinationsIterator(x.rows(), x.columns());

				double dmax = 0;
				int[] rBasis = null;
				while (iter.hasNext()) {
					int[] cur = iter.next();
					DoubleMatrix2D subMatrix = x.viewSelection(cur, null);
					double det = Math.abs(Algebra.ZERO.det(subMatrix));
					if (det > dmax) {
						dmax = det;
						rBasis = cur;
					}
				}
				int[] rr = setDiff(x.rows(), rBasis);

				DoubleMatrix2D xbinv = Algebra.ZERO.inverse(x.viewSelection(rBasis, null));
				DoubleMatrix2D b = x.zMult(xbinv, x.copy());
				DoubleMatrix2D xperp = new DenseDoubleMatrix2D(rr.length, x.rows());
				xperp.viewSelection(null, rBasis).assign(b.viewSelection(rr, null));
				int n = rr.length;
				for (int i = 0; i < n; i++)
					xperp.set(i, rr[i], -1);
				return xperp;
			}
		}
		//
		return null;
	}

	private static int[] setDiff(int max, int[] maxInds) {
		int n = max - maxInds.length;
		int[] ret = new int[n];
		int maxIndsInd = 0;
		int j = 0;
		for (int i = 0; i < max; i++) {
			if (maxIndsInd < maxInds.length && i == maxInds[maxIndsInd]) {
				// Not in
				maxIndsInd++;
			} else {
				ret[j++] = i;
			}
		}
		return ret;
	}

	/**
	 * returns a row vector equaling the generalized cross product of matrix x;
	 * x must be n by (n-1) matrix
	 * 
	 * gencrossprod(x)*x = null vector;
	 * 
	 * @param x
	 * @return
	 */
	private static DoubleMatrix1D genCrossProd(DoubleMatrix2D x) {
		DoubleMatrix1D m = minors(x).viewFlip();
		int n = m.size();
		for (int i = 0; i < n; i++)
			m.setQuick(i, m.getQuick(i) * (((i + 1) % 2) * 2 - 1));

		return m;
	}

	/**
	 * returns list of circuits of oriented matroid specified by matrix a
	 * ordered by dimension
	 * 
	 * @param a
	 *            n * k matrix n >= k
	 * @return
	 */
	public static int[][] circuits(DoubleMatrix2D a) {
		// DoubleMatrix2D x = checkRank(a);
		int[][] c = null;
		int rows = a.rows();
		int cols = a.columns();
		if (a.size() > 0) {
			if (cols < rows) {
				long nChoosekLen = CombinatoricsUtils.binomialCoefficient(rows, cols + 1);
				Iterator<int[]> iter = CombinatoricsUtils.combinationsIterator(rows, cols + 1);

				c = new int[(int) nChoosekLen][rows];
				int i = 0;
				while (iter.hasNext()) {
					int[] curIter = iter.next();
					DoubleMatrix2D xsub = a.viewSelection(curIter, null);
					DoubleMatrix1D cp = genCrossProd(xsub);

					for (int j = 0; j < cols + 1; j++)
						c[i][curIter[j]] = (int) Math.signum(cp.get(j));

					i++;
				}
			}
		}

		if (c != null) {
			rows = c.length;
			int[] zAbss = zEncode(c);
			TIntHashSet unique = new TIntHashSet();
			for (int i = 0; i < rows; i++)
				if (zAbss[i] != 0)
					unique.add(Math.abs(zAbss[i]));
			zAbss = unique.toArray();
			c = zDecode(zAbss, c[0].length);
			Arrays.sort(c, new NegAbsIntArrayComparitor());
		}
		return c;
	}

	/**
	 * returns sign vectors of length n from vector of zone numbers z;
	 * 
	 * @param d
	 * @return
	 */
	public static int[][] zDecode(int[] z, int n) {
		int rows = z.length;
		int[][] d = new int[rows][n];

		long[] powers = CMRUtil.REVERSE_POW_3[n];
		long sumP = CMRUtil.SUM_REVERSE_POW_3[n];

		for (int k = 0; k < rows; k++) {
			long r = z[k] + sumP;
			long u = floorDiv(r, powers[0]) - 1;
			d[k][0] = (int) u;
			for (int i = 1; i < n; i++) {
				r -= (1 + u) * powers[i - 1];
				u = floorDiv(r, powers[i]) - 1;
				d[k][i] = (int) u;
			}
		}
		return d;
	}

	/**
	 * returns ONE sign vector of length n from zone number z;
	 * 
	 * @param d
	 * @return
	 */
	public static int[] zDecode(int z, int n) {
		int[] d = new int[n];
		long[] powers = CMRUtil.REVERSE_POW_3[n];
		long sumP = CMRUtil.SUM_REVERSE_POW_3[n];

		long r = z + sumP;
		long u = floorDiv(r, powers[0]) - 1;
		d[0] = (int) u;
		for (int i = 1; i < n; i++) {
			r -= (1 + u) * powers[i - 1];
			u = floorDiv(r, powers[i]) - 1;
			d[i] = (int) u;
		}

		return d;
	}

	/**
	 * returns zone number of matrix of difference vectors, d
	 * 
	 * @param d
	 * @return
	 */
	public static int[] zEncode(int[][] d) {
		int n = d[0].length;
		int rows = d.length;
		int[] s = new int[rows];

		long[] powers = CMRUtil.REVERSE_POW_3[n];

		for (int i = 0; i < rows; i++) {
			int sum = 0;
			for (int j = 0; j < n; j++) {
				double v = d[i][j];
				if (v < 0)
					sum -= powers[j];
				else if (v > 0)
					sum += powers[j];
			}

			s[i] = sum;
		}
		return s;
	}

	/**
	 * returns zone number for one difference vector, d
	 * 
	 * @param d
	 * @return
	 */
	public static int zEncode(int[] d) {
		int n = d.length;

		long[] powers = CMRUtil.REVERSE_POW_3[n];
		int sum = 0;
		for (int j = 0; j < n; j++) {
			double v = d[j];
			if (v < 0)
				sum -= powers[j];
			else if (v > 0)
				sum += powers[j];
		}

		return sum;
	}

	/**
	 * generates all sign vectors of length n. excludes null vector
	 * 
	 * @param n
	 * @return
	 */
	public static int[][] allSignVectors(int n) {
		int nn = (int) ((CMRUtil.POW_3[n] - 1) / 2);

		int[] z = new int[nn];
		for (int i = 0; i < nn; i++)
			z[i] = i + 1;

		return zDecode(z, n);
	}

	/**
	 * checks for sign concordance between each row of matrix x and each row of
	 * matrix z. x is a n * k matrix and z is a m * k matrix. status is a n * m
	 * logical matrix where 1=sign concordant, 0=not. let w = u o v be the
	 * Hadamard (or componentwise product) of vectors u and v, then u and v are
	 * sign concordant iff: (1) sign(w_i) >= 0 for all i (2) sign(w_i) > 0 for
	 * at least one i
	 * 
	 * if u and v are not sign concordant then they are orthogonal (in oriented
	 * matroid sense)
	 * 
	 * @param x
	 * @param z
	 * @return
	 */
	public static boolean[][] signCon(int[][] x, int[][] z) {
		if (x.length == 0 || z.length == 0 || x[0].length != z[0].length)
			return null;
		int n = x[0].length;
		int xRows = x.length;
		int zRows = z.length;

		boolean[][] status = new boolean[x.length][z.length];

		for (int i = 0; i < zRows; i++) {
			int[] curZRow = z[i];
			for (int j = 0; j < xRows; j++) {
				int[] curXRow = x[j];
				boolean negOk = true;
				boolean posOk = false;
				for (int c = 0; c < n && negOk; c++) {
					int cur = curZRow[c] * curXRow[c];
					negOk = cur >= 0;
					posOk = posOk || cur > 0;
				}
				status[j][i] = negOk && posOk;
			}
		}
		return status;
	}

	/**
	 * returns list of covectors of oriented matroid specified by matrix a
	 * ordered by dimension
	 * 
	 * @param a
	 * @return
	 */
	public static int[][] covectors(DoubleMatrix2D a) {
		// DoubleMatrix2D x = checkRank(a);
		// if (x.size() == 0)
		// return null;

		int[][] c = circuits(a);
		int[][] d = allSignVectors(a.rows());
		if (c == null)
			return d;
		boolean[][] pos = signCon(d, c);
		if (pos == null)
			throw new IllegalArgumentException("Unable to create positive covectors");
		for (int[] cc : c) {
			int n = cc.length;
			for (int i = 0; i < n; i++)
				cc[i] *= -1;
		}
		boolean[][] neg = signCon(d, c);
		if (neg == null)
			throw new IllegalArgumentException("Unable to create negative covectors");

		int n = pos.length;
		int nCol = neg[0].length;
		int count = 0;
		boolean[] falseN = new boolean[n];
		for (int i = 0; i < n; i++) {
			boolean[] posI = pos[i];
			boolean[] negI = neg[i];
			boolean allFalse = true;

			for (int j = 0; j < nCol && allFalse; j++)
				allFalse = !posI[j] && !negI[j];

			falseN[i] = allFalse;
			if (allFalse)
				count++;
		}

		int[][] cv = new int[count][];
		int j = 0;
		for (int i = 0; i < n; i++)
			if (falseN[i])
				cv[j++] = d[i];
		return cv;
	}

	/**
	 * returns the closure of a set of a sign vectors for two sign vectors s and
	 * t, if t_i = s_i or t_i = 0, then t is an element of the sign closure of s
	 * 
	 * @param s
	 * @return
	 */
	public static int[][] signClosure(int[][] s) {
		int n = s[0].length;

		TIntHashSet unique = new TIntHashSet();
		TIntArrayList nonZeroIndices = new TIntArrayList(n);
		int[] tmp = new int[n];
		for (int[] cur : s) {
			// Add the vector itself
			unique.add(Math.abs(zEncode(cur)));

			// For every non-zero, make it zero then add it.
			nonZeroIndices.clear();
			for (int i = 0; i < n; i++)
				if (cur[i] != 0)
					nonZeroIndices.add(i);
			int count = nonZeroIndices.size();
			for (int k = 1; k <= count; k++) {
				Iterator<int[]> iter = CombinatoricsUtils.combinationsIterator(count, k);
				while (iter.hasNext()) {
					System.arraycopy(cur, 0, tmp, 0, n);
					int[] comb = iter.next();
					for (int i = 0; i < k; i++)
						tmp[nonZeroIndices.get(comb[i])] = 0;
					unique.add(Math.abs(zEncode(tmp)));
				}
			}
		}
		unique.remove(0);
		int[] uniqueZs = unique.toArray();
		Arrays.sort(uniqueZs);
		int[][] zs = zDecode(uniqueZs, n);

		return zs;
	}

	/**
	 * returns the "augment" of a set of a sign vector for two sign vectors s
	 * and t, if t_i = s_i or s_i = 0, then t is an element of the sign augment
	 * of s
	 * 
	 * @param vectors
	 * @return
	 */
	public static int[][] signAugment(int[][] s) {
		int n = s[0].length;
		TIntHashSet unique = new TIntHashSet();

		int[] tmp = new int[n];
		int[] tmpNeg = new int[n];
		boolean[] zeroIndices = new boolean[n];
		for (int[] cur : s) {
			unique.add(Math.abs(zEncode(cur)));

			int count = 0;
			for (int i = 0; i < n; i++) {
				zeroIndices[i] = cur[i] == 0;
				if (zeroIndices[i])
					count++;
			}

			// not called often, no cache needed
			int[][] u = allSignVectors(count);
			for (int[] curU : u) {
				int iz = 0;
				for (int i = 0; i < n; i++) {
					if (zeroIndices[i]) {
						tmp[i] = curU[iz];
						tmpNeg[i] = -curU[iz++];
					} else {
						tmp[i] = cur[i];
						tmpNeg[i] = cur[i];
					}
				}
				unique.add(Math.abs(zEncode(tmp)));
				unique.add(Math.abs(zEncode(tmpNeg)));
			}
		}
		unique.remove(0);

		int[] uniqueZs = unique.toArray();
		Arrays.sort(uniqueZs);
		int[][] zs = zDecode(uniqueZs, n);

		return zs;
	}

	/**
	 * returns list of vectors of matrix a
	 * 
	 * tolerance is 1e-10
	 * 
	 * @param a
	 * @return
	 */
	public static int[][] vectors(DoubleMatrix2D a) {
		double tol = 1e-10;
		// TODO: I probably dont need check rank!
		DoubleMatrix2D x = checkRank(a);

		if (x == null || x.size() == 0)
			return null;
		x = orthocom(x);
		if (x == null || x.size() == 0)
			return null;

		int rows = x.rows();
		int columns = x.columns();
		for (int row = rows; --row >= 0;)
			for (int column = columns; --column >= 0;)
				if (Math.abs(x.getQuick(row, column)) < tol)
					x.setQuick(row, column, 0);

		return covectors(x.viewDice());
	}

	/**
	 * returns dimension of sign vector z
	 * 
	 * @param allSignVectors
	 * @return
	 */
	public static int[] dim(int[][] z) {
		int l = z.length;
		int n = z[0].length;
		int[] diag = new int[l];
		for (int i = 0; i < l; i++) {
			int sum = 0;
			for (int j = 0; j < n; j++) {
				sum += z[i][j] * z[i][j];
			}
			diag[i] = sum;
		}

		return diag;
	}

	/**
	 * Copied from jdk 1.8 (can be removed later...)
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static int floorDiv(int x, int y) {
		int r = x / y;
		// if the signs are different and modulo not zero, round down
		if ((x ^ y) < 0 && (r * y != x)) {
			r--;
		}
		return r;
	}

	public static long floorDiv(long x, long y) {
		long r = x / y;
		// if the signs are different and modulo not zero, round down
		if ((x ^ y) < 0 && (r * y != x)) {
			r--;
		}
		return r;
	}
}
