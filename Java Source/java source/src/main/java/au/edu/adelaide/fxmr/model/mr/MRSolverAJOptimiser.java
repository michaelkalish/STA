package au.edu.adelaide.fxmr.model.mr;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.QuadraticMultivariateRealFunction;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.joptimizer.optimizers.OptimizationRequest;
import au.edu.adelaide.fxmr.joptimizer.optimizers.OptimizationResponse;
import au.edu.adelaide.fxmr.joptimizer.optimizers.PrimalDualMethod;
import cern.colt.Arrays;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

/**
 * Solver based on a heavily modified JOptimiser (which can only solve problems
 * of the MR type)
 */
public class MRSolverAJOptimiser extends MRSolver {
	private static final double FEAS_DIST = 0.01;
	private static final double FEAS_DIST_2 = 1;
	// default is 0.55 - for simple problems, 0.7 generally goes faster
	private static final double BETA = 0.7;
	// If it fails the first time, do as many iterations as it takes!
	private static final int ITER_2 = 100000;
	// This beta should be used when the solution wasn't found (more iterations
	// will be required as well)
	private static final double BETA_2 = 0.1;

	private int free = 0;

	private boolean allowCyclicProblems = true;
	
	//private static double maxTime = 0;
	

	public MRSolution solve(MRProblem p) {// throws CatastrophicMRFailure {
		if (p.constraintsAlreadySatisfied()) {
			// Trivial case
			return new MRSolution(p, p.getY(), 0, 0);
		}

		int n = p.getN();
		HashSet<SimpleLinearConstraint> ineq = new HashSet<>();
		HashSet<SimpleLinearConstraint> eq = new HashSet<>();
		double[] initial = p.findFeasibleStart(FEAS_DIST, ineq, eq);

		if (!eq.isEmpty() && !allowCyclicProblems)
			return null;

		// Only non-trivial problems count as a "call"
		calls++;

		// We assume weights is positive definite and symmetric. If not, this
		// will probably fail!
		DoubleMatrix2D weights = p.getWeights();

		// inequalities
		int nIneq = ineq.size();
		SimpleLinearConstraint[] inequalities = new SimpleLinearConstraint[nIneq];
		ineq.toArray(inequalities);

		// Equalities
		SimpleLinearConstraint[] eqMat = null;
		if (!eq.isEmpty()) {
			eqMat = new SimpleLinearConstraint[eq.size()];
			eq.toArray(eqMat);
		}

		// Create the objective function
		DoubleMatrix1D q = new DenseDoubleMatrix1D(n);
		weights.zMult(new DenseDoubleMatrix1D(p.getY()), q, -1, 0, true);

		QuadraticMultivariateRealFunction objectiveFunction = new QuadraticMultivariateRealFunction(
				weights.toArray(), q.toArray(), 0);

		// optimization problem
		OptimizationRequest or = new OptimizationRequest();
		or.setBeta(BETA);
		or.setF0(objectiveFunction);
		or.setInitialPoint(initial);
		or.setFi(inequalities);
		or.setA(eqMat);
		or.setToleranceFeas(tolInitFeas);
		or.setTolerance(tolInit);

		// optimization
		PrimalDualMethod pdm = new PrimalDualMethod();
		pdm.setOptimizationRequest(or);
		int returnCode;

		MRSolution s = null;
		long start = System.nanoTime();
		try {
			returnCode = pdm.optimize();

			if (returnCode == OptimizationResponse.SUCCESS) {
				double[] soln = pdm.getOptimizationResponse().getSolution();
				s = new MRSolution(p, soln, (double) (System.nanoTime() - start) / 1_000_000_000.0, pdm.getIteration());
				
				// if (s.getTimeS() > maxTime){
				// maxTime = s.getTimeS();
				// System.out.println(maxTime);
				// System.out.println(objectiveFunction.getP());
				// System.out.println(objectiveFunction.getQ());
				// System.out.println(objectiveFunction.getR());
				// System.out.println(eq);
				// System.out.println(ineq);
				// System.out.println(Arrays.toString(initial));
				// System.out.println(pdm.getIteration());
				// }
				
			}
		} catch (Exception e) {
			// Not too worried at this point - s will be null
			// e.printStackTrace();
			// System.out.println(p.toMatlabString());
		}

		boolean failedOnce = s == null;

		if (failedOnce) {
			// Something went wrong, try more conservative values
			ineq.clear();
			eq.clear();
			initial = p.findFeasibleStart(FEAS_DIST_2, ineq, eq);
			or.setToleranceFeas(tol2Feas);
			or.setTolerance(tol2);
			or.setBeta(BETA_2);
			or.setMaxIteration(ITER_2);
			or.setInitialPoint(initial);

			try {
				returnCode = pdm.optimize();
				if (returnCode == OptimizationResponse.SUCCESS) {
					double[] soln = pdm.getOptimizationResponse().getSolution();
					s = new MRSolution(p, soln, (double) (System.nanoTime() - start) / 1_000_000_000.0,
							pdm.getIteration());
				}
			} catch (Exception e) {
				// NOTE: return value will be null!

				// if (!failQuietly) {
				// System.err.println(
				// "Exception in optimisation with maxIter=" + ITER_2 + ",
				// infNormCheck=" + normCheck + ", n=" + p.getN()
				// + ", nIneq="
				// + inequalities.length + (eqMat == null ? ", no equalities" :
				// ", nEq=" + eqMat.length));
				// e.printStackTrace();
				// }
			}
		}

		// if (s == null && !failQuietly) {
		// // Something went wrong AGAIN - no idea what to do? :(
		// System.err.println("Failed optimisation with maxIter=" + ITER_2 + ",
		// infNormCheck=" + normCheck + ", infNormCheck2="
		// + normCheck2 + ", n=" + p.getN() + ", nIneq="
		// + inequalities.length + (eqMat == null ? ", no equalities" : ", nEq="
		// + eqMat.length));
		// }

		// if (s == null)
		// when s == null, MR failed catastrophically
		// throw new CatastrophicMRFailure();

		return s;
	}

	public int getFree() {
		return free;
	}

	public int getCalls() {
		return calls;
	}

	public boolean isAllowCyclicProblems() {
		return allowCyclicProblems;
	}

	public void setAllowCyclicProblems(boolean allowCyclicProblems) {
		this.allowCyclicProblems = allowCyclicProblems;
	}
}
