package au.edu.adelaide.fxmr.model.bin;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class ParBinCMRxSolver extends BinCMRxSolver {
	public ParBinCMRxSolver() {
	}

	public BinSolution[] solve(BinBaseProblem problem, int proc) {
		int n = problem.getModel().getnSubj();
		BinProblem[] problems = new BinProblem[n];
		for (int i = 0; i < n; i++)
			problems[i] = new BinProblem(problem.getModel(), i, problem.getRangeSet());
		return solve(problems, proc);
	}

	public BinSolution[] solve(final BinProblem[] problems, int proc) {
		int n = problems.length;
		final BinSolution[] solutions = new BinSolution[n];

		int procBase = proc < 1 ? Runtime.getRuntime().availableProcessors() : proc;
		if (n >= procBase || onlyFeas) {
			// More (or equal) problems than processors. Easy split. onlyFeas isn't done in parallel anyway.
			ExecutorService pool = Executors.newFixedThreadPool(procBase);
			for (int i = 0; i < n; i++) {
				final int fi = i;
				pool.execute(new Runnable() {
					@Override
					public void run() {
						solutions[fi] = solve(problems[fi], 1);
					}
				});
			}
			pool.shutdown();
			try {
				pool.awaitTermination(10000, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		} else {
			// Complicated split, run sequentially with L-parallelisation
			// (procBase > 1)
			for (int i = 0; i < n; i++)
				solutions[i] = solve(problems[i], procBase);
		}

		return solutions;
	}

	public BinSolution[] solve(BinProblem[] problems) {
		return solve(problems, 0);
	}

	public BinSolution solve(BinProblem problem) {
		return solve(problem, 0);
	}

	public BinSolution solve(BinProblem problem, int proc) {
		if (proc == 1) {
			return solve(problem, null);
		} else {
			ParBinCMRxSolverImpl parSolver = new ParBinCMRxSolverImpl();
			return parSolver.solve(problem, proc);
		}
	}

}
