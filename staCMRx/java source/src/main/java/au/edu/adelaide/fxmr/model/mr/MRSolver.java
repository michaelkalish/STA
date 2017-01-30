package au.edu.adelaide.fxmr.model.mr;

public interface MRSolver {
	public MRSolution solve(MRProblem p);// throws CatastrophicMRFailure;

	public int getCalls();
}
