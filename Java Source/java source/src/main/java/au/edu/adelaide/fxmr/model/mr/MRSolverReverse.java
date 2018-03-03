package au.edu.adelaide.fxmr.model.mr;

public class MRSolverReverse extends MRSolver {
	private MRSolverAJOptimiser internalSolver = new MRSolverAJOptimiser();

	@Override
	public MRSolution solve(MRProblem p) {
		// Copy tolerance over
		internalSolver.setTolerance(tolInit, tol2);

		MRSolution solBase = internalSolver.solve(p);
		if (solBase == null)
			return null;

		MRProblem inverted = p.invertConstraints();

		MRSolution solInverted = internalSolver.solve(inverted);
		if (solInverted == null)
			return null;

		return new MRSolution(solBase, solInverted);
	}

	public void setAllowCyclicProblems(boolean allowCyclicProblems) {
		if (internalSolver == null)
			internalSolver = new MRSolverAJOptimiser();
		internalSolver.setAllowCyclicProblems(allowCyclicProblems);
	}
}
