package au.edu.adelaide.fxmr.eqoptimiser;

import au.edu.adelaide.fxmr.joptimizer.functions.QuadraticMultivariateRealFunction;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.LUDecompositionQuick;
/**
 * Simple system with equalities. Use of this is very rare
 * @author luke
 *
 */
public class EQSolver {
	private SimpleLinearConstraint[] eq;
	private QuadraticMultivariateRealFunction obj;

	public EQSolver(QuadraticMultivariateRealFunction objectiveFunction, SimpleLinearConstraint[] eqMat) {
		this.obj = objectiveFunction;
		this.eq = eqMat;
	}

	public int getIterations() {
		return 1;
	}

	public double[] solve() {
		DoubleMatrix2D p = obj.getP();
		DoubleMatrix1D c = obj.getQ();

		int nq = obj.getP().columns();
		int n = nq + eq.length;
		DoubleMatrix2D lh = new DenseDoubleMatrix2D(n, n);
		DoubleMatrix1D rhs = new DenseDoubleMatrix1D(n);

		for (int r = 0; r < nq; r++) {
			rhs.setQuick(r, -c.getQuick(r));
			for (int col = 0; col < nq; col++)
				lh.setQuick(r, col, p.getQuick(r, col));
		}

		int i = nq;
		for (SimpleLinearConstraint con : eq) {
			// rhs.setQuick(n+i, 0); - nope, already 0
			lh.setQuick(i, con.getPosIndex(), 1);
			lh.setQuick(i, con.getNegIndex(), -1);

			lh.setQuick(con.getPosIndex(), i, 1);
			lh.setQuick(con.getNegIndex(), i, -1);
			i++;
		}

		LUDecompositionQuick lu = new LUDecompositionQuick();
		lu.decompose(lh);
		lu.solve(rhs);

		double[] solLag = new double[nq];
		for (i = 0; i < nq; i++) {
			solLag[i] = rhs.getQuick(i);
		}

		return solLag;
	}

}
