package au.edu.adelaide.fxmr.model.mr;

/**
 * MRSolver based on oj! algo's QP capability
 * 
 * Of note: as of V39, oj! algo doesn't work. This may change with V40 but JOptimizer works.
 * 
 * This class is disabled - it didn't work and I changed the underlying data
 * structure. It should be an easy fix once oj! releases a new version for
 * testing.
 * 
 * 
 */
public class MRSolverOJAlgo implements MRSolver {
	public MRSolverOJAlgo() {
	}

	public MRSolution solve(MRProblem p) {
		return null;
	}
	// final ExpressionsBasedModel tmpModel = buildModel(p);
	// System.out.println(tmpModel);
	// final long tmpBefore = System.nanoTime();
	//
	// final Result tmpResult = tmpModel.minimise();
	// final long tmpAfter = System.nanoTime();
	// System.out.println(tmpResult);
	// final double timeS = (tmpAfter - tmpBefore) / 1_000_000_000.0;
	//
	// int n = p.getN();
	// double[] x = new double[n];
	// for (int i = 0; i < n; i++) {
	// x[i] = tmpResult.doubleValue(i);
	// }
	//
	// return new MRSolution(p, x, timeS, 0);
	// }
	//
	// private ExpressionsBasedModel buildModel(MRProblem p) {
	// int n = p.getN();
	//
	// DoubleMatrix1D q = new DenseDoubleMatrix1D(n);
	// p.getWeights().zMult(new DenseDoubleMatrix1D(p.getY()), q, -1, 0, true);
	//
	// final Variable[] tmpVariables = new Variable[n];
	// for (int i = 0; i < n; i++) {
	// tmpVariables[i] = Variable.make("V" + Integer.toString(i)).lower(ZERO);
	// tmpVariables[i].weight(q.get(i));
	// }
	//
	// final ExpressionsBasedModel retVal = new
	// ExpressionsBasedModel(tmpVariables);
	// p.getAdj().forEachNonZero(new IntIntDoubleFunction() {
	// @Override
	// public double apply(int row, int col, double val) {
	// Expression e = retVal.addExpression(row + "," + col);
	// e.upper(ZERO);
	// e.set(tmpVariables[row], 1);
	// e.set(tmpVariables[col], -1);
	// return val;
	// }
	// });
	//
	// final Expression quads = retVal.addExpression("Quads");
	// p.getWeights().forEachNonZero(new IntIntDoubleFunction() {
	// @Override
	// public double apply(int row, int col, double val) {
	// quads.set(tmpVariables[row], tmpVariables[col], val);
	// return val;
	// }
	// });
	// quads.weight(0.5);
	//
	// return retVal;
	// }

	@Override
	public int getCalls() {
		return 0;
	}


}
