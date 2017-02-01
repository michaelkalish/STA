package au.edu.adelaide.fxmr.model;

import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.util.MultivariateNormalDistribution;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class CMRxFits implements Fits {
	private double[] fits;
	private double[] times;
	private double p;
	private double dataFit;
	private MultivariateNormalDistribution[] dists;
	private int nVar;
	private CMRxSolver solver;
	private CMRxFitsProblem problem;
	private double shrink;
	private Badness[] badnesses;
	private boolean cheapP;
	private HashSet<SimpleLinearConstraint>[] adj;

	public double getDataFit() {
		return dataFit;
	}

	public double[] getFits() {
		return fits;
	}

	public double getP() {
		return p;
	}

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc) {
		this(nSample, problem, shrink, proc, false, false, 0, 0);
	}

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP) {
		this(nSample, problem, shrink, proc, cheapP, false, 0, 0);
	}

	/**
	 * Do a parametric fit (original data not available)
	 * 
	 * @param nSample
	 * @param problem
	 * @param proc
	 */
	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1,
			double mrTol2) {
		this.problem = problem;
		this.shrink = shrink;
		this.cheapP = cheapP;
		this.adj = problem.getAdj();
		nVar = problem.getnVarOriginal();
		if (onlySTAMR)
			solver = new STASolver();
		else
			solver = new CMRxSolver();
		solver.setMrTolerance1(mrTol1);
		solver.setMrTolerance2(mrTol2);

		CMRSolution initSoln = solver.solve(problem);
		dataFit = initSoln.getFStar();
		initSoln = null;

		HashSet<SimpleLinearConstraint>[] adj = problem.getAdj();
		if (adj != null && adj.length > 0 && !onlySTAMR) {
			MRSolverAJOptimiser mrSolver = new MRSolverAJOptimiser();
			mrSolver.setTolerance(mrTol1, mrTol2);
			// Take away MR from dataFit
			for (int v = 0; v < problem.getNVar(); v++) {
				if (adj[v].size() > 0) {
					MRProblem mrp = new MRProblem(problem.getMeans()[v], problem.getWeights()[v], adj[v]);
					dataFit -= mrSolver.solve(mrp).getfVal();
				}
			}
		}

		fits = new double[nSample];
		times = new double[nSample];
		badnesses = new Badness[nSample];
		dists = new MultivariateNormalDistribution[nVar];
		for (int i = 0; i < nVar; i++) {
			dists[i] = new MultivariateNormalDistribution(problem.getMeansOriginal()[i], problem.getCov()[i],
					System.currentTimeMillis() - i * 99999);
		}

		ExecutorService pool = Executors.newFixedThreadPool(proc < 1 ? Runtime.getRuntime().availableProcessors() : proc);

		for (int i = 0; i < nSample; i++)
			pool.execute(new FitJob(i));
		pool.shutdown();

		try {
			pool.awaitTermination(10000, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		int numGood = 0;
		for (int s = 0; s < nSample; s++)
			if (fits[s] >= dataFit)
				numGood++;

		p = (double) numGood / nSample;
	}

	/**
	 * Parametric bootstrab using recycled MultivariateNormalDistribution
	 * objects
	 * 
	 * @param dists
	 * @param ybs
	 */
	private void parametricBootstrap(MultivariateNormalDistribution[] dists, DoubleMatrix2D[] ybs) {
		int n = dists.length;
		for (int i = 0; i < n; i++)
			dists[i].fill(ybs[i]);
	}

	public class FitJob implements Runnable {
		private int index;

		public FitJob(int i) {
			this.index = i;
		}

		@Override
		public void run() {
			long start = System.nanoTime();
			Badness worst = Badness.NONE;
			double[][] tmpMeans = new double[nVar][];
			DoubleMatrix2D[] tmpWeights = new DoubleMatrix2D[nVar];
			MultivariateNormalDistribution[] tmpDists = new MultivariateNormalDistribution[nVar];

			int nCond = problem.getNCond();
			DenseDoubleMatrix2D[] ybs = new DenseDoubleMatrix2D[nVar];
			for (int i = 0; i < nVar; i++) {
				ybs[i] = new DenseDoubleMatrix2D(problem.getN()[i], nCond);
			}

			parametricBootstrap(dists, ybs);
			// Fit to 1D model??
			for (int i = 0; i < nVar; i++) {
				StatsSTA sta;
				if (shrink < 0)
					sta = new StatsSTA(ybs[i]);
				else
					sta = new StatsSTA(ybs[i], shrink);

				worst = Badness.worst(worst, sta.getBad());
				tmpMeans[i] = sta.getMeans().toArray();
				tmpWeights[i] = sta.getWeights();
			}
			CMRxProblem tmpProblem = new CMRxProblem(tmpMeans, tmpWeights, adj, problem.getModelOriginal());
			CMRSolution tmpSolution = solver.solve(tmpProblem);

			// bootstrap sample from 1D model data
			for (int i = 0; i < nVar; i++) {
				tmpDists[i] = new MultivariateNormalDistribution(tmpSolution.getXStar()[i], problem.getCov()[i],
						System.currentTimeMillis() - i * 99999);
			}
			parametricBootstrap(tmpDists, ybs);
			// calculate actual fit
			for (int i = 0; i < nVar; i++) {
				StatsSTA sta;
				if (shrink < 0)
					sta = new StatsSTA(ybs[i]);
				else
					sta = new StatsSTA(ybs[i], shrink);

				worst = Badness.worst(worst, sta.getBad());
				tmpMeans[i] = sta.getMeans().toArray();
				tmpWeights[i] = sta.getWeights();
			}
			tmpProblem = new CMRxProblem(tmpMeans, tmpWeights, adj, problem.getModelOriginal());
			tmpSolution = solver.solve(tmpProblem, null, cheapP, dataFit);

			fits[index] = tmpSolution.getFStar();
			times[index] = (System.nanoTime() - start) / 1000000000.0;

			badnesses[index] = worst;
		}
	}

	@Override
	public Badness[] getBadnesses() {
		return badnesses;
	}

	@Override
	public double[] getTimes() {
		return times;
	}
}
