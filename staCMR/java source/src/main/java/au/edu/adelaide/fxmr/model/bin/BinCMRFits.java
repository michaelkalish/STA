package au.edu.adelaide.fxmr.model.bin;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;

public class BinCMRFits {
	private BinBaseProblem problem;
	private BinModel model;
	private double[] baseFitDiff;
	private double[][] fits;
	private double[] p;
	private int nSubj;
	private boolean notCoupled;

	/**
	 * 
	 * @param nSample
	 * @param problem
	 * @param proc
	 */

	public BinCMRFits(int nSample, BinBaseProblem problem) {
		this(nSample, problem, Runtime.getRuntime().availableProcessors(), false);
	}

	public BinCMRFits(int nSample, BinBaseProblem problem, int proc, boolean notCoupled) {
		this.problem = problem;
		this.model = problem.getModel();
		this.nSubj = model.getnSubj();
		this.notCoupled = notCoupled;
		this.baseFitDiff = calcFitDiff(problem, notCoupled ? new BinMRSolver() : new BinCMRSolver());
		

		fits = new double[nSample][];

		ExecutorService pool = Executors.newFixedThreadPool(proc < 1 ? Runtime.getRuntime().availableProcessors() : proc);

		for (int i = 0; i < nSample; i++)
			pool.execute(new FitJob(i));
		pool.shutdown();
		try {
			pool.awaitTermination(10000, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		p = new double[nSubj];
		for (int s = 0; s < nSubj; s++) {
			int count = 0;
			double baseFit = baseFitDiff[s];
			for (int i = 0; i < nSample; i++)
				if (fits[i][s] >= baseFit)
					count++;

			p[s] = (double) count / nSample;
		}
	}

	private double[] calcFitDiff(BinBaseProblem problemLoc, BinCMRSolver solver) {
		BinModel modelLoc = problemLoc.getModel();
		int nVar = modelLoc.getnVar();
		int nSubj = modelLoc.getnSubj();

		double[] fitDiff = new double[nSubj];

		if (problemLoc.getRangeSet() != null && problemLoc.getRangeSet().length != 0  && !notCoupled) {
			MRSolver mrSolver = new MRSolverAJOptimiser();
			for (int s = 0; s < nSubj; s++) {
				double g2Val = 0;
				for (int v = 0; v < nVar; v++) {
					BinElement data = modelLoc.getData()[s][v];
					MRProblem mrp = new MRProblem(data.getMeans(), data.getN(), problemLoc.getRangeSet());

					double[] x = mrSolver.solve(mrp).getxVector();
					BinTrial.clampZeroOne(x);
					g2Val += BinTrial.mleBN(data, x);
				}

				fitDiff[s] -= g2Val;
			}
		}

		BinSolution[] solns = solver.solve(problemLoc);
		for (int s = 0; s < nSubj; s++)
			fitDiff[s] += solns[s].getG2Star();

		return fitDiff;
	}

	public class FitJob implements Runnable {

		private int index;

		public FitJob(int i) {
			this.index = i;
		}

		@Override
		public void run() {
			BinCMRSolver solver = notCoupled ? new BinMRSolver() : new BinCMRSolver();
			BinModel r1 = model.resample();

			BinBaseProblem pr1 = new BinBaseProblem(r1, problem.getRangeSet());
			BinSolution[] solnR1 = solver.solve(pr1);

			BinModel r2 = model.resample(solnR1);
			BinBaseProblem pr2 = new BinBaseProblem(r2, problem.getRangeSet());

			fits[index] = calcFitDiff(pr2, solver);
		}
	}

	public double[] getP() {
		return p;
	}

	public double[][] getFits() {
		return fits;
	}

	public double[] getBaseFitDiff() {
		return baseFitDiff;
	}
}
