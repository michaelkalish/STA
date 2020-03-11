package au.edu.adelaide.fxmr.model.bin;

import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.ui.StatusFrame;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

public class BinCMRxFits {
	private BinBaseProblem problem;
	private BinModel model;
	private double[] baseFitDiff;
	private double[][] fits;
	private double[][][][] xStars;
	private boolean[] complete;
	private AtomicBoolean running = new AtomicBoolean(true);
	private double[] p;
	private int nSubj;
	private StatusFrame sf;
	private long nextUpdate = 0;
	private int nSample;
	private BinCMRxSolver solver;
	private static final String FORMAT_SF = "%3d %15.8f %15.8f %15.8f %15.8f\n";

	/**
	 * 
	 * @param nSample
	 * @param problem
	 * @param proc
	 */

	public BinCMRxFits(int nSample, BinBaseProblem problem) {
		this(nSample, problem, 0, false, false);
	}

	public BinCMRxFits(int nSample, BinBaseProblem problem, int proc, boolean approximate, boolean showStatus) {
		this(nSample, problem, proc, approximate, showStatus, -1);
	}

	public BinCMRxFits(int nSample, BinBaseProblem problem, int proc, boolean approximate, boolean showStatus,
			long seed) {
		if (showStatus) {
			sf = new StatusFrame("BinCMRxFits", running);
			sf.updateStatus("Init...");
		}

		solver = new BinCMRxSolver();
		solver.setOnlyFeas(approximate);

		this.problem = problem;
		this.model = problem.getModel();
		this.nSubj = model.getnSubj();
		this.baseFitDiff = calcFitDiff(problem, -1);
		this.nSample = nSample;

		proc = proc < 1 ? Runtime.getRuntime().availableProcessors() : proc;
		if (showStatus) {
			sf.updateStatus("Ready, Starting " + proc + " thread" + (proc == 1 ? "" : "s") + ".");
		}

		fits = new double[nSample][];
		xStars = new double[nSample][][][];
		if (showStatus)
			complete = new boolean[nSample];

		ExecutorService pool = Executors.newFixedThreadPool(proc);

		Random r = new Random(seed == -1 ? System.currentTimeMillis() : seed);

		for (int i = 0; i < nSample; i++)
			pool.execute(new FitJob(i, r.nextInt()));
		pool.shutdown();

		if (showStatus)
			sf.setCancelText("Accept Current Fit Data");
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

		p = new double[nSubj];
		for (int s = 0; s < nSubj; s++) {
			int count = 0;
			double baseFit = baseFitDiff[s];
			for (int i = 0; i < nSample; i++)
				if ((!showStatus || complete[i]) && fits[i][s] >= baseFit)
					count++;

			p[s] = (double) count / countDone;
		}

		if (showStatus)
			sf.setFinished();
	}

	private double[] calcFitDiff(BinBaseProblem problemLoc, int index) {
		BinModel modelLoc = problemLoc.getModel();
		int nVar = modelLoc.getnVar();
		int nSubj = modelLoc.getnSubj();

		double[] fitDiff = new double[nSubj];
		double[][][] curXStars = new double[nSubj][][];

		if (problemLoc.getRangeSet() != null && problemLoc.getRangeSet().length != 0) {
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

		BinSolution[] solns = solver.solve(problemLoc, running);
		for (int s = 0; s < nSubj; s++) {
			fitDiff[s] += solns[s].getG2Star();
			curXStars[s] = solns[s].getXStar();
		}

		if (index != -1) {
			fits[index] = fitDiff;
			xStars[index] = curXStars;
		}

		return fitDiff;
	}

	public class FitJob implements Runnable {
		private int index;
		private RandomEngine binRand;

		public FitJob(int i, int seed) {
			this.index = i;
			this.binRand = new MersenneTwister(seed);
		}

		@Override
		public void run() {
			if (!running.get())
				return;

			BinModel r1 = model.resample(binRand);

			if (!running.get())
				return;

			BinBaseProblem pr1 = new BinBaseProblem(r1, problem.getRangeSet(), problem.getCmrModel());
			// OLD code: BinBaseProblem pr1 = new BinBaseProblem(r1,
			// problem.getRangeSet());

			// System.out.println("int[][][] pr1 = " + pr1.toString());

			BinSolution[] solnR1 = solver.solve(pr1, running);

			if (!running.get())
				return;

			BinModel r2 = model.resample(solnR1, binRand);
			BinBaseProblem pr2 = new BinBaseProblem(r2, problem.getRangeSet(), problem.getCmrModel());
			// OLD code: BinBaseProblem pr2 = new BinBaseProblem(r2,
			// problem.getRangeSet());

			// System.out.println("int[][][] pr2 = " + pr2.toString());
			if (!running.get())
				return;

			calcFitDiff(pr2, index);

			if (sf != null) {
				complete[index] = true;
				updateSF();
			}
		}
	}

	public double[] getP() {
		return p;
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

		double[] curP = new double[nSubj];
		double[] min = new double[nSubj];
		double[] max = new double[nSubj];

		for (int s = 0; s < nSubj; s++) {
			int count = 0;
			min[s] = Double.MAX_VALUE;
			double baseFit = baseFitDiff[s];
			for (int i = 0; i < nSample; i++)
				if (complete[i]) {
					double fis = fits[i][s];

					if (fis >= baseFit)
						count++;
					if (fis < min[s])
						min[s] = fis;
					if (fis > max[s])
						max[s] = fis;
				}
			curP[s] = (double) count / countDone;
		}

		sb.append("Completed ");
		sb.append(countDone);
		sb.append(" / ");
		sb.append(nSample);
		sb.append("\n\nN        baseFit         min             max             P\n");
		for (int s = 0; s < nSubj; s++) {
			sb.append(String.format(FORMAT_SF, s + 1, baseFitDiff[s], min[s], max[s], curP[s]));
		}
		sf.updateStatus(sb.toString());
	}

	public double[][] getFits() {
		int countDone = 0;
		if (running.get())
			countDone = nSample;
		else
			for (int i = 0; i < nSample; i++)
				if (complete[i])
					countDone++;
		if (countDone == nSample)
			return fits;

		double[][] subFits = new double[countDone][];
		int j = 0;
		for (int i = 0; i < nSample; i++)
			if (complete[i])
				subFits[j++] = fits[i];
		return subFits;
	}

	public double[] getBaseFitDiff() {
		return baseFitDiff;
	}

	public double[][][][] getXStars() {
		return xStars;
	}
}
