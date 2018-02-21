package au.edu.adelaide.fxmr.model;

import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import gnu.trove.map.hash.TIntObjectHashMap;

public class TrialCache {
	private TIntObjectHashMap<MRSolution> map = new TIntObjectHashMap<>();
	private int nVar;

	public TrialCache(int nVar) {
		this.nVar = nVar;
	}

	public void clear() {
		map.clear();
	}

	public MRSolution solve(int i, MRSolver solver, MRProblem mrProblem) {
		if (mrProblem.getNewestSingleConstraint() == null)
			return solver.solve(mrProblem);

		int hash = (i + 1) * (mrProblem.getNewestSingleConstraint().getPosIndex() + nVar + 1);

		if (map.containsKey(hash))
			return map.get(hash);

		MRSolution s = solver.solve(mrProblem);
		map.put(hash, s);
		return s;
	}
}
