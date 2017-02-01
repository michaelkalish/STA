package au.edu.adelaide.fxmr.model.bin;

import java.util.HashSet;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;

/**
 * Same as BinCMRSolver but without the coupling (mostly for use with
 * binCMRFits)
 *
 */
public class BinMRSolver extends BinCMRSolver {
	@Override
	public BinSolution solve(BinProblem problem) {
		long start = System.nanoTime();

		BinModel model = problem.getModel();
		int nvar = model.getnVar();

		BinElement[] data = model.getData()[problem.getSubjectIndex()];

		MRSolver mrSolver = new MRSolverAJOptimiser();
		BinTrial base = new BinTrial(problem.getRangeSet(), nvar, data);

		double[][] xBar = base.solveMRs(mrSolver);
		double fBar = base.getF();
		double gBar = base.getG2();
		HashSet<SimpleLinearConstraint>[] adjBar = base.getAdjs();

		return new BinSolution(fBar, gBar, xBar, 1, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls());
	}
}
