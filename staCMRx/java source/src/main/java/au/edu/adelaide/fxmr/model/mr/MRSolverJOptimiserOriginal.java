package au.edu.adelaide.fxmr.model.mr;

/**
 * Solver based on JOptimiser directly (uncomment to test)
 * 
 * 
 */
public class MRSolverJOptimiserOriginal extends MRSolver {
	public MRSolution solve(MRProblem p) {
		return null;
/*
		if (p.getMatIneq().isEmpty()) {
			// Trivial case
			return new MRSolution(p, p.getY(), 0, 0);
		}
		calls++;

		DoubleMatrix2D weights = p.getWeights();
		DoubleMatrix2D check = p.getWeights().copy();
		check.assign(p.getWeights().viewDice(), Functions.minus);
		double normCheck = Algebra.ZERO.normInfinity(check);
		if (normCheck > 1e-15) {
			// Force symmetry
			weights = check;
			weights.assign(p.getWeights());
			weights.assign(p.getWeights().viewDice(), Functions.plus);
			weights.assign(Functions.div(2));
			if (normCheck > 1e-4) {
				// Tell the user how bad things are
				System.err.println("Weights matrix in MRSolverAJOptimiser is not symmetric! InfNorm(W-W^T)=" + normCheck + " forcing");
			}
		}
		
		// Find initial feasible point - also find cycles
		HashSet<SimpleLinearConstraint> ineq = new HashSet<>();
		HashSet<SimpleLinearConstraint> eq = new HashSet<>();
		double[] initial = p.findFeasibleStart(MRSolverAJOptimiser.FEAS_DIST, ineq, eq);
		if (initial == null)
			return null;

		// Create the objective function
		int n = p.getN();
		DoubleMatrix1D q = new DenseDoubleMatrix1D(n);
		weights.zMult(new DenseDoubleMatrix1D(p.getY()), q, -1, 0, true);

		ConvexMultivariateRealFunction objectiveFunction = new PDQuadraticMultivariateRealFunction(
				weights.toArray(), q.toArray(), 0);

		// inequalities
		int nIneq = ineq.size();
		LinearMultivariateRealFunction[] inequalities = new LinearMultivariateRealFunction[nIneq];
		int i = 0;
		for (SimpleLinearConstraint c : ineq) {
			double[] qVector = new double[c.getDim()];
			qVector[c.getPosIndex()] = 1;
			qVector[c.getNegIndex()] = -1;

			LinearMultivariateRealFunction f = new LinearMultivariateRealFunction(qVector, 0);
			inequalities[i++] = f;
		}

		// equalities
		int nEq = 0;
		DoubleMatrix2D eqA = null;
		if (!eq.isEmpty()) {
			nEq = eq.size();
			eqA = new DenseDoubleMatrix2D(nEq, n);
			i = 0;
			for (SimpleLinearConstraint c : eq) {
				eqA.set(i, c.getPosIndex(), 1);
				eqA.set(i, c.getNegIndex(), -1);
				i++;
			}
		}

		// optimization problem
		OptimizationRequest or = new OptimizationRequest();
		or.setBeta(MRSolverAJOptimiser.BETA);
		or.setF0(objectiveFunction);

		or.setInitialPoint(initial);
		or.setFi(inequalities);
		if (!eq.isEmpty()) {
			or.setA(eqA);
			or.setB(new DenseDoubleMatrix1D(nEq));
		}
		or.setToleranceFeas(MRSolverAJOptimiser.TOL_INIT_FEAS);
		or.setTolerance(MRSolverAJOptimiser.TOL_INIT);

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
				s = new MRSolution(p, soln, (double) (System.nanoTime() - start) / 1_000_000_000.0, 0);
			}
		} catch (Exception e) {
			e.printStackTrace();
			// Silent fail
		}

		if (s == null) {
			// Something went wrong, try more conservative values
			ineq.clear();
			eq.clear();
			initial = p.findFeasibleStart(MRSolverAJOptimiser.FEAS_DIST_2, ineq, eq);
			or.setToleranceFeas(MRSolverAJOptimiser.TOL_2_FEAS);
			or.setTolerance(MRSolverAJOptimiser.TOL_2);
			or.setBeta(MRSolverAJOptimiser.BETA_2);
			or.setMaxIteration(MRSolverAJOptimiser.ITER_2);
			or.setInitialPoint(initial);

			try {
				returnCode = pdm.optimize();
				if (returnCode == OptimizationResponse.SUCCESS) {
					double[] soln = pdm.getOptimizationResponse().getSolution();
					s = new MRSolution(p, soln, (double) (System.nanoTime() - start) / 1_000_000_000.0, 0);
				} else {
					System.err.println("Failed optimisation with " + MRSolverAJOptimiser.ITER_2 + " iterations!");
				}
			} catch (Exception e) {
				e.printStackTrace();
				// return null

				// System.err.println("Optimisation warning...");
				// System.out.println("\n\n");
				// System.out.println(p.toMatlabString());
				//
				// System.out.println("errF = abs(fit-" + s.getfVal() +
				// ")");
				// System.out.println("errX = abs(x-" +
				// Arrays.toString(soln) + "')");
			}
		}

		// previousSolnMap.put(p, s);
		return s;
*/
	}
}
