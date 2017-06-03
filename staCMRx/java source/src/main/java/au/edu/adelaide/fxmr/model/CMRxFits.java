package au.edu.adelaide.fxmr.model;

import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.mr.MRSolverReverse;
import au.edu.adelaide.fxmr.model.ui.StatusFrame;
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
	private MRSolver mrSolver;
	private boolean onlySTAMR;
	private long seed;
	private StatusFrame sf;
	private AtomicBoolean running = new AtomicBoolean(true);
	private boolean[] complete;
	private long nextUpdate;
	private int nSample;

	public double getDataFit() {
		return dataFit;
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

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR) {
		this(nSample, problem, shrink, proc, cheapP, onlySTAMR, 0, 0);
	}

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1,
			double mrTol2) {
		this(nSample, problem, shrink, proc, cheapP, onlySTAMR, mrTol1, mrTol2, false, false);
	}

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1,
			double mrTol2, boolean approximate) {
		this(nSample, problem, shrink, proc, cheapP, onlySTAMR, mrTol1, mrTol2, approximate, false);
	}

	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1,
			double mrTol2, boolean approximate, boolean reverse) {
		this(nSample, problem, shrink, proc, cheapP, onlySTAMR, mrTol1, mrTol2, approximate, false, -1, false);
	}

	/**
	 * Do a parametric fit (original data not available)
	 * 
	 * @param nSample
	 * @param problem
	 * @param proc
	 */
	public CMRxFits(int nSample, CMRxFitsProblem problem, double shrink, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1,
			double mrTol2, boolean approximate, boolean reverse, long seed, boolean showStatus) {
		this.problem = problem;
		this.shrink = shrink;
		this.cheapP = cheapP;
		this.adj = problem.getAdj();
		this.onlySTAMR = onlySTAMR;
		this.nSample = nSample;
		nVar = problem.getnVarOriginal();

		if (showStatus) {
			sf = new StatusFrame("CMRxFits", running);
			sf.updateStatus("Init...");
		}

		this.seed = seed == -1 ? System.currentTimeMillis() : seed;

		dists = new MultivariateNormalDistribution[nVar];
		for (int i = 0; i < nVar; i++) {
			dists[i] = new MultivariateNormalDistribution(problem.getMeansOriginal()[i], problem.getCov()[i]);
			if (!dists[i].isSymPosDef()) {
				// This will probably be called from MATLAB so report it
				// one-indexed.
				System.err.println("Input covariance matrix " + (i + 1) + " is not positive definite?! (1 indexed, nvar=\" + nVar + \")");
				return;
			}
		}

		if (onlySTAMR)
			solver = new STASolver(reverse);
		else {
			solver = new CMRxSolver();
			solver.setOnlyFeas(approximate);
		}
		solver.setMrTolerance1(mrTol1);
		solver.setMrTolerance2(mrTol2);

		CMRSolution initSoln = solver.solve(problem);

		dataFit = initSoln.getFStar();
		initSoln = null;

		HashSet<SimpleLinearConstraint>[] adj = problem.getAdj();
		if (adj != null && adj.length > 0 && !onlySTAMR) {
			mrSolver = reverse ? new MRSolverReverse() : new MRSolverAJOptimiser();
			mrSolver.setTolerance(mrTol1, mrTol2);
			// Take away MR from dataFit
			for (int v = 0; v < problem.getNVar(); v++) {
				if (adj[v].size() > 0) {
					MRProblem mrp = new MRProblem(problem.getMeans()[v], problem.getWeights()[v], adj[v]);
					dataFit -= mrSolver.solve(mrp).getfVal();
				}
			}
		}

		proc = proc < 1 ? Runtime.getRuntime().availableProcessors() : proc;
		if (showStatus) {
			sf.updateStatus("Ready, Starting " + proc + " thread" + (proc == 1 ? "" : "s") + ".");
		}

		fits = new double[nSample];
		if (showStatus)
			complete = new boolean[nSample];

		times = new double[nSample];
		badnesses = new Badness[nSample];
		ExecutorService pool = Executors.newFixedThreadPool(proc);

		if (showStatus)
			sf.setCancelText("Accept Current Fit Data");
		for (int i = 0; i < nSample; i++)
			pool.execute(new FitJob(i));
		pool.shutdown();

		try {
			pool.awaitTermination(10000, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		int countDone = 0;
		if (running.get())
			countDone = nSample;
		else
			for (int i = 0; i < nSample; i++)
				if (complete[i])
					countDone++;

		int count = 0;
		for (int i = 0; i < nSample; i++)
			if ((!showStatus || complete[i]) && fits[i] >= dataFit)
				count++;

		p = (double) count / countDone;

		if (showStatus)
			sf.setFinished();

	}

	/**
	 * Parametric bootstrab using recycled MultivariateNormalDistribution
	 * objects
	 * 
	 * @param dists
	 * @param ybs
	 * @param rand
	 */
	private void parametricBootstrap(MultivariateNormalDistribution[] dists, DoubleMatrix2D[] ybs, Random rand) {
		int n = dists.length;
		for (int i = 0; i < n; i++)
			dists[i].fill(ybs[i], rand);
	}

	public class FitJob implements Runnable {
		private int index;

		public FitJob(int i) {
			this.index = i;
		}

		@Override
		public void run() {
			if (!running.get())
				return;

			long start = System.nanoTime();
			Random rand = new Random(seed + index);

			Badness worst = Badness.NONE;
			double[][] tmpMeans = new double[nVar][];
			DoubleMatrix2D[] tmpWeights = new DoubleMatrix2D[nVar];
			MultivariateNormalDistribution[] tmpDists = new MultivariateNormalDistribution[nVar];

			int nCond = problem.getNCond();
			DenseDoubleMatrix2D[] ybs = new DenseDoubleMatrix2D[nVar];
			for (int i = 0; i < nVar; i++) {
				ybs[i] = new DenseDoubleMatrix2D(problem.getN()[i], nCond);
			}

			parametricBootstrap(dists, ybs, rand);
			// Fit to 1D model
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
				tmpDists[i] = new MultivariateNormalDistribution(tmpSolution.getXStar()[i], problem.getCov()[i]);
			}
			parametricBootstrap(tmpDists, ybs, rand);
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

			if (adj != null && adj.length > 0 && !onlySTAMR) {
				// Take away MR from dataFit
				for (int v = 0; v < problem.getNVar(); v++) {
					if (adj[v].size() > 0) {
						MRProblem mrp = new MRProblem(tmpMeans[v], tmpWeights[v], adj[v]);
						fits[index] -= mrSolver.solve(mrp).getfVal();
					}
				}
			}

			times[index] = (System.nanoTime() - start) / 1000000000.0;
			badnesses[index] = worst;

			if (sf != null) {
				complete[index] = true;
				updateSF();
			}
		}
	}

	private void updateSF() {
		if (System.currentTimeMillis() < nextUpdate) {
			return;
		}
		nextUpdate = System.currentTimeMillis() + 1000;
		StringBuilder sb = new StringBuilder();
		int countDone = 0;

		for (int i = 0; i < nSample; i++)
			if (complete[i])
				countDone++;

		double min = Double.MAX_VALUE;
		double max = 0;

		int count = 0;
		for (int i = 0; i < nSample; i++)
			if (complete[i]) {
				double fis = fits[i];

				if (fis >= dataFit)
					count++;
				if (fis < min)
					min = fis;
				if (fis > max)
					max = fis;
			}
		double curP = (double) count / countDone;

		sb.append("Completed ");
		sb.append(countDone);
		sb.append(" / ");
		sb.append(nSample);
		sb.append("\n\ndataFit\t");
		sb.append(dataFit);

		sb.append("\n\nmin\t");
		sb.append(min);

		sb.append("\nmax\t");
		sb.append(max);

		sb.append("\n\nP\t");
		sb.append(curP);
		sf.updateStatus(sb.toString());
	}

	@Override
	public Badness[] getBadnesses() {
		int countDone = 0;
		if (running.get())
			countDone = nSample;
		else
			for (int i = 0; i < nSample; i++)
				if (complete[i])
					countDone++;
		if (countDone == nSample)
			return badnesses;

		Badness[] subBadnesses = new Badness[countDone];
		int j = 0;
		for (int i = 0; i < nSample; i++)
			if (complete[i])
				subBadnesses[j++] = badnesses[i];
		return subBadnesses;
	}

	@Override
	public double[] getTimes() {
		return subGet(times);
	}

	private double[] subGet(double[] vals) {
		int countDone = 0;
		if (running.get())
			countDone = nSample;
		else
			for (int i = 0; i < nSample; i++)
				if (complete[i])
					countDone++;
		if (countDone == nSample)
			return vals;

		double[] subVals = new double[countDone];
		int j = 0;
		for (int i = 0; i < nSample; i++)
			if (complete[i])
				subVals[j++] = vals[i];
		return subVals;
	}

	public double[] getFits() {
		return subGet(fits);
	}
}
