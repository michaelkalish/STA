package au.edu.adelaide.fxmr.model;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class CMRxGMFits implements Fits {
	private double[] fits;
	private Badness[] badnesses;
	private double p;
	private double dataFit;
	private CMRxSolver solver;
	private GeneralModel gm;
	private int nVar;
	private int nCond;
	private Random rand = new Random();
	private DoubleMatrix2D model;
	private boolean autoShrink;
	private double shrinkage;
	private HashSet<SimpleLinearConstraint>[] adj;
	private double[][] csSTAmeans;
	private double[] times;
	private boolean cheapP;
	private FitListener listener;
	private boolean running = true;
	private long nextUpdate;

	@Override
	public double getDataFit() {
		return dataFit;
	}

	@Override
	public double[] getFits() {
		return fits;
	}

	public Badness[] getBadnesses() {
		return badnesses;
	}

	@Override
	public double getP() {
		return p;
	}

	public CMRxGMFits(int nSample, GeneralModel gm, DoubleMatrix2D model, HashSet<SimpleLinearConstraint>[] adj,
			int proc) {
		this(nSample, gm, 0, true, model, adj, proc, false, null, false, MRSolver.TOL1, MRSolver.TOL2);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrinkage, DoubleMatrix2D model,
			HashSet<SimpleLinearConstraint>[] adj, int proc, boolean cheapP, boolean onlySTAMR) {
		this(nSample, gm, shrinkage, shrinkage < 0, model, adj, proc, cheapP, null, onlySTAMR, MRSolver.TOL1, MRSolver.TOL2);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrinkage, DoubleMatrix2D model,
			HashSet<SimpleLinearConstraint>[] adj, int proc, boolean cheapP, FitListener listner, double mrTol1, double mrTol2) {
		this(nSample, gm, shrinkage, shrinkage < 0, model, adj, proc, cheapP, listner, false, mrTol1, mrTol2);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrink, DenseDoubleMatrix2D denseDoubleMatrix2D,
			HashSet<SimpleLinearConstraint>[] dMatAs, int proc, boolean cheapP2, boolean onlySTAMR, double mrTol1, double mrTol2) {
		this(nSample, gm, shrink, onlySTAMR, denseDoubleMatrix2D, dMatAs, proc, cheapP2, null, onlySTAMR, mrTol1, mrTol2);
	}

	/**
	 * Do a parametric fit (original data not available)
	 * 
	 * @param nSample
	 * @param gm
	 * @param adj
	 * @param cheapP
	 */
	private CMRxGMFits(int nSample, GeneralModel gm, double shrinkage, boolean autoShrink, DoubleMatrix2D model,
			HashSet<SimpleLinearConstraint>[] adj, int proc, boolean cheapP, FitListener listener, boolean onlySTAMR, double mrTolerance1,
			double mrTolerance2) {
		this.gm = gm;
		this.model = model;
		this.shrinkage = shrinkage;
		this.autoShrink = autoShrink;
		this.adj = adj;
		this.cheapP = cheapP;
		this.listener = listener;
		int[] depVar = gm.getNDepVar();
		int[] conds = gm.getNBSConds();
		nVar = depVar.length;
		nCond = conds.length;
		if (onlySTAMR)
			solver = new STASolver();
		else
			solver = new CMRxSolver();
		solver.setEasyFail(true);
		solver.setMrTolerance1(mrTolerance1);
		solver.setMrTolerance2(mrTolerance2);

		double[][][][] ybSource = new double[nVar][nCond][][];
		int[] nVarSamples = new int[nVar];
		for (int v = 0; v < nVar; v++) {
			for (int c = 0; c < nCond; c++) {
				ybSource[v][c] = gm.getBlock(depVar[v], conds[c]);
				nVarSamples[v] += ybSource[v][c].length;
			}
		}

		CombinedStatsSTA[] csSTA = autoShrink ? gm.calcStats() : gm.calcStats(shrinkage);

		csSTAmeans = new double[nVar][];
		DoubleMatrix2D[] weights = new DoubleMatrix2D[nVar];
		for (int v = 0; v < nVar; v++) {
			csSTAmeans[v] = csSTA[v].getMeans().toArray();
			weights[v] = csSTA[v].getWeights();
		}

		CMRxProblem initProblem = new CMRxProblem(csSTAmeans, weights, adj, model);
		// Use parallel solver for first solution
		CMRSolution initSoln;
		if (onlySTAMR) {
			initSoln = solver.solve(initProblem);
		} else {
			ParCMRxSolver parSolver = new ParCMRxSolver();
			parSolver.setMrTolerance1(mrTolerance1);
			parSolver.setMrTolerance2(mrTolerance2);
			initSoln = parSolver.solve(initProblem);
		}

		dataFit = initSoln.getFStar();
		initSoln = null;

		if (adj != null && adj.length > 0 && !onlySTAMR) {
			MRSolverAJOptimiser mrSolver = new MRSolverAJOptimiser();
			mrSolver.setTolerance(mrTolerance1, mrTolerance2);
			// Take away MR from dataFit
			for (int v = 0; v < nVar; v++) {
				MRProblem problem = new MRProblem(csSTAmeans[v], weights[v], adj[v]);
				problem.forceSymetry();
				dataFit -= mrSolver.solve(problem).getfVal();
			}
		}

		fits = new double[nSample];
		Arrays.fill(fits, -1);
		times = new double[nSample];
		badnesses = new Badness[nSample];
		ExecutorService pool = Executors.newFixedThreadPool(proc < 1 ? Runtime.getRuntime().availableProcessors() : proc);

		nextUpdate = System.currentTimeMillis() + 3000;
		for (int i = 0; i < nSample; i++)
			pool.execute(new FitJob(i));
		pool.shutdown();

		try {
			pool.awaitTermination(10000, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			// This can happen when the user cancels
		}

		int numGood = 0;
		int finalNSample = 0;
		for (int s = 0; s < nSample; s++) {
			if (fits[s] != -1)
				finalNSample++;
			if (fits[s] >= dataFit)
				numGood++;
		}

		p = (double) numGood / finalNSample;

		if (listener != null) {
			listener.updateStatus(fits, dataFit);
			listener.setFinished();
		}
	}

	public class FitJob implements Runnable {
		private int index;
		private Badness worst = Badness.NONE;

		public FitJob(int i) {
			this.index = i;
		}

		@Override
		public void run() {
			long start = System.nanoTime();
			boolean needSoln = true;
			while (needSoln && running) {
				// Bootstrap from ybSource
				GeneralModel tmpGeneralModel = gm.bootStrap(rand);

				if (!tmpGeneralModel.sanityCheck()) {
					// Resample!
					continue;
				}

				// Solve bootstrapped model
				CMRxProblem tmpProblem = makeCMRxProblem(tmpGeneralModel);
				if (tmpProblem == null || !running)
					continue;

				CMRSolution tmpSolution = solver.solve(tmpProblem);
				if (tmpSolution == null || !running)
					continue;

				// "movedata" = original data - means + tmpSolutionX
				tmpGeneralModel = gm.moveFrom(csSTAmeans, tmpSolution.getXStar()).bootStrap(rand);
				if (!tmpGeneralModel.sanityCheck()) {
					// Resample
					continue;
				}

				tmpProblem = makeCMRxProblem(tmpGeneralModel);
				if (tmpProblem == null || !running)
					continue;

				tmpSolution = solver.solve(tmpProblem, null, cheapP, dataFit);
				if (tmpSolution == null || !running) {
					continue;
				}

				fits[index] = tmpSolution.getFStar();
				needSoln = false;
			}

			times[index] = (System.nanoTime() - start) / 1000000000.0;

			System.out.println(times[index]);

			badnesses[index] = worst;

			if (listener != null && System.currentTimeMillis() > nextUpdate) {
				running = listener.updateStatus(fits, dataFit);
				nextUpdate = System.currentTimeMillis() + 1000;
			}
		}

		private CMRxProblem makeCMRxProblem(GeneralModel gm) {
			try {
				CombinedStatsSTA[] csSTA = autoShrink ? gm.calcStats() : gm.calcStats(shrinkage);
				int nVar = csSTA.length;

				double[][] means = new double[nVar][];
				DoubleMatrix2D[] weights = new DoubleMatrix2D[nVar];
				for (int v = 0; v < nVar; v++) {
					means[v] = csSTA[v].getMeans().toArray();
					weights[v] = csSTA[v].getWeights();
					worst = Badness.worst(csSTA[v].getWorst(), worst);
				}

				return new CMRxProblem(means, weights, adj, model);
			} catch (IllegalArgumentException e) {
				// This occurs when the problem is singular - just make another
				// problem
				return null;
			}
		}
	}

	@Override
	public double[] getTimes() {
		return times;
	}
}
