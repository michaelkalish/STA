package au.edu.adelaide.fxmr.model.bin;

import java.util.HashSet;
import java.util.TreeSet;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.CMRUtil;
import au.edu.adelaide.fxmr.model.CMRxSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.set.hash.TIntHashSet;

public class BinCMRSolver {
	public BinCMRSolver() {
	}

	public BinSolution[] solve(BinBaseProblem problem) {
		int n = problem.getModel().getnSubj();
		BinSolution[] solutions = new BinSolution[n];
		for (int i = 0; i < n; i++)
			solutions[i] = solve(new BinProblem(problem.getModel(), i, problem.getRangeSet()));

		return solutions;
	}

	public BinSolution[] solve(BinProblem[] problems) {
		int n = problems.length;
		BinSolution[] solutions = new BinSolution[n];
		for (int i = 0; i < n; i++)
			solutions[i] = solve(problems[i]);

		return solutions;
	}

	public BinSolution solve(BinProblem problem) {
		long start = System.nanoTime();

		BinModel model = problem.getModel();
		int nvar = model.getnVar();
		double gFloor = Double.NEGATIVE_INFINITY;
		double gBar = Double.POSITIVE_INFINITY;
		double fBar = Double.POSITIVE_INFINITY;
		double[][] xBar = null;

		BinElement[] data = model.getData()[problem.getSubjectIndex()];
		int ncond = data[0].getMeans().length;
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);

		HashSet<SimpleLinearConstraint>[] adjBar = null;
		TIntHashSet infeas = CMRUtil.calcInfeasibleZoneNumbers(nvar);
		int iter = 0;

		TreeSet<BinTrial> remaining = new TreeSet<>();
		MRSolver mrSolver = new MRSolverAJOptimiser();

		// Add first CMRTrial
		remaining.add(new BinTrial(problem.getRangeSet(), nvar, data));

		while (!remaining.isEmpty() && gFloor < gBar) {
			BinTrial current = remaining.pollFirst();
			gFloor = current.getG2();
			iter++;

			if (gFloor < gBar) {
				double[][] xPrimes = current.solveMRs(mrSolver);
				if (xPrimes == null) {
					// failed optimisation
					return null;
				}

				double gFit = current.getG2();
				int[] feas = CMRxSolver.isFeasible3n(xPrimes, infeas, tmpVolumes, tmpZoneNumbers);

				if (gFit < gBar) {
					// Solution is better than current best
					if (feas == null) {
						// Solution is feasible
						gBar = gFit;
						fBar = current.getF();

						xBar = xPrimes;
						adjBar = current.getAdjs();
					} else {
						// Solution not feasible - make it so!
						int negIndex = feas[0];
						int posIndex = feas[1];
						remaining.add(current.split(posIndex, negIndex, 0));
						remaining.add(current.split(negIndex, posIndex, 1));
					}
				}
			}
		}

		return new BinSolution(fBar, gBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls());
	}
}
