package au.edu.adelaide.fxmr.model;

import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.mr.MRSolverReverse;
import au.edu.adelaide.fxmr.model.ui.StatusFrame;
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
	private DoubleMatrix2D model;
	private boolean autoShrink;
	private double shrinkage;
	private HashSet<SimpleLinearConstraint>[] adj;
	private double[][] csSTAmeans;
	private double[] times;
	private boolean cheapP;
	private AtomicBoolean running = new AtomicBoolean(true);
	private long nextUpdate;
	private MRSolver mrSolver;
	private boolean onlySTAMR;
	private long seed;
	private StatusFrame sf;
	private boolean[] complete;
	private int nSample;

	@Override
	public double getDataFit() {
		return dataFit;
	}

	@Override
	public double getP() {
		return p;
	}

	// TODO: this overloading is a mess - clean it up one day!
	public CMRxGMFits(int nSample, GeneralModel gm, DoubleMatrix2D model, HashSet<SimpleLinearConstraint>[] adj,
			int proc) {
		this(nSample, gm, 0, true, model, adj, proc, false, false, MRSolver.TOL1, MRSolver.TOL2, false, false, -1, false);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrinkage, DoubleMatrix2D model,
			HashSet<SimpleLinearConstraint>[] adj, int proc, boolean cheapP, boolean onlySTAMR) {
		this(nSample, gm, shrinkage, shrinkage < 0, model, adj, proc, cheapP, onlySTAMR, MRSolver.TOL1, MRSolver.TOL2, false, false, -1, false);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrink, DoubleMatrix2D denseDoubleMatrix2D,
			HashSet<SimpleLinearConstraint>[] dMatAs, int proc, boolean cheapP, boolean onlySTAMR, double mrTol1, double mrTol2) {
		this(nSample, gm, shrink, denseDoubleMatrix2D, dMatAs, proc, cheapP, onlySTAMR, mrTol1, mrTol2, false);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrink, DoubleMatrix2D denseDoubleMatrix2D,
			HashSet<SimpleLinearConstraint>[] dMatAs, int proc, boolean cheapP2, boolean onlySTAMR, double mrTol1, double mrTol2,
			boolean approximate) {
		this(nSample, gm, shrink, shrink < 0, denseDoubleMatrix2D, dMatAs, proc, cheapP2, onlySTAMR, mrTol1, mrTol2, approximate,
				false, -1, false);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrink, DoubleMatrix2D denseDoubleMatrix2D,
			HashSet<SimpleLinearConstraint>[] dMatAs, int proc, boolean cheapP2, boolean onlySTAMR, double mrTol1, double mrTol2,
			boolean approximate, boolean reverse) {
		this(nSample, gm, shrink, shrink < 0, denseDoubleMatrix2D, dMatAs, proc, cheapP2, onlySTAMR, mrTol1, mrTol2, approximate,
				reverse, -1, false);
	}

	public CMRxGMFits(int nSample, GeneralModel gm, double shrink, DenseDoubleMatrix2D denseDoubleMatrix2D, HashSet<SimpleLinearConstraint>[] dMatAs, int proc,
			boolean cheapP, boolean onlySTAMR, double mrTol1, double mrTol2, boolean approximate, boolean reverse, long seed, boolean showStatus) {
		this(nSample, gm, shrink, shrink < 0, denseDoubleMatrix2D, dMatAs, proc, cheapP, onlySTAMR, mrTol1, mrTol2, approximate,
				reverse, seed, showStatus);
	}

	/**
	 * Do a parametric fit (original data available)
	 * 
	 * @param nSample
	 * @param gm
	 * @param adj
	 * @param cheapP
	 */
	private CMRxGMFits(int nSample, GeneralModel gm, double shrinkage, boolean autoShrink, DoubleMatrix2D model,
			HashSet<SimpleLinearConstraint>[] adj, int proc, boolean cheapP, boolean onlySTAMR, double mrTolerance1,
			double mrTolerance2, boolean approximate, boolean reverse, long seed, boolean showStatus) {
		this.gm = gm;
		this.model = model;
		this.shrinkage = shrinkage;
		this.autoShrink = autoShrink;
		this.adj = adj;
		this.cheapP = cheapP;
		this.onlySTAMR = onlySTAMR;
		int[] depVar = gm.getNDepVar();
		int[] conds = gm.getNBSConds();
		this.seed = seed == -1 ? System.currentTimeMillis() : seed;
		this.nSample = nSample;

		if (showStatus) {
			sf = new StatusFrame("CMRxGMFits", running);
			sf.updateStatus("Init...");
		}

		nVar = depVar.length;
		nCond = conds.length;
		if (onlySTAMR)
			solver = new STASolver(reverse);
		else {
			solver = new CMRxSolver();
			solver.setOnlyFeas(approximate);
		}
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
		if (onlySTAMR || approximate) {
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
			mrSolver = reverse ? new MRSolverReverse() : new MRSolverAJOptimiser();
			mrSolver.setTolerance(mrTolerance1, mrTolerance2);
			// Take away MR from dataFit
			for (int v = 0; v < nVar; v++) {
				MRProblem problem = new MRProblem(csSTAmeans[v], weights[v], adj[v]);
				problem.forceSymetry();
				dataFit -= mrSolver.solve(problem).getfVal();
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
			// This can happen when the user cancels
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
			Random rand = new Random(seed + index);
			while (needSoln && running.get()) {
				// Bootstrap from ybSource
				GeneralModel tmpGeneralModel = gm.bootStrap(rand);

				if (!tmpGeneralModel.sanityCheck()) {
					// Resample!
					continue;
				}

				// Solve bootstrapped model
				CMRxProblem tmpProblem = makeCMRxProblem(tmpGeneralModel);
				if (tmpProblem == null || !running.get())
					continue;

				CMRSolution tmpSolution = solver.solve(tmpProblem);
				if (tmpSolution == null || !running.get())
					continue;

				// "movedata" = original data - means + tmpSolutionX
				tmpGeneralModel = gm.moveFrom(csSTAmeans, tmpSolution.getXStar()).bootStrap(rand);
				if (!tmpGeneralModel.sanityCheck()) {
					// Resample
					continue;
				}

				tmpProblem = makeCMRxProblem(tmpGeneralModel);
				if (tmpProblem == null || !running.get())
					continue;

				tmpSolution = solver.solve(tmpProblem, null, cheapP, dataFit);
				if (tmpSolution == null || !running.get()) {
					continue;
				}

				fits[index] = tmpSolution.getFStar();

				if (adj != null && adj.length > 0 && !onlySTAMR) {
					// Take away MR from fit
					for (int v = 0; v < nVar; v++) {
						MRProblem problem = new MRProblem(tmpProblem.getMeans()[v], tmpProblem.getWeights()[v], adj[v]);
						problem.forceSymetry();
						fits[index] -= mrSolver.solve(problem).getfVal();
					}
				}

				needSoln = false;
			}

			times[index] = (System.nanoTime() - start) / 1000000000.0;
			badnesses[index] = worst;
			
			if (sf != null) {
				complete[index] = true;
				updateSF();
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
