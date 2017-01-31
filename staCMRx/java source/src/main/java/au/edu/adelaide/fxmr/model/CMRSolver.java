package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.set.hash.TIntHashSet;

public class CMRSolver {
	private static final double ZERO_TOL = 1e-6;

	private CMRListener listener;

	private double mrTolerance1;
	private double mrTolerance2;

	public double getMrTolerance1() {
		return mrTolerance1;
	}

	public void setMrTolerance1(double mrTolerance1) {
		this.mrTolerance1 = mrTolerance1;
	}

	public double getMrTolerance2() {
		return mrTolerance2;
	}

	public void setMrTolerance2(double mrTolerance2) {
		this.mrTolerance2 = mrTolerance2;
	}

	public CMRSolver() {
	}

	public CMRSolver(CMRListener listener) {
		this.listener = listener;
	}

	public CMRSolution solve(CMRProblem problem) {
		long start = System.nanoTime();

		int nvar = problem.getNVar();
		int fBarReductions = 0;
		double[][] means = problem.getMeans();
		DoubleMatrix2D[] weights = problem.getWeights();
		double fFloor = Double.NEGATIVE_INFINITY;
		double fBar = Double.POSITIVE_INFINITY;
		double[][] xBar = null;
		HashSet<SimpleLinearConstraint>[] adjBar = null;
		TIntHashSet infeas = CMRUtil.calcInfeasibleZoneNumbers(nvar);
		ArrayList<CMRIter> iter = new ArrayList<>();

		TreeSet<CMRTrial> remaining = new TreeSet<>();

		MRSolver mrSolver = new MRSolverAJOptimiser();
		if (mrTolerance1 != 0)
			mrSolver.setTolerance(mrTolerance1, mrTolerance2);

		// Add first CMRTrial
		remaining.add(new CMRTrial(problem.getRangeSet(), nvar, weights, means));

		while (!remaining.isEmpty() && fFloor < fBar) {
			CMRTrial current = remaining.pollFirst();
			fFloor = current.getF();
			iter.add(new CMRIter(fFloor, fBar, remaining.size()));

			if (fFloor < fBar) {
				double[][] xPrimes = current.solveMRs(mrSolver);
				if (xPrimes == null) {
					// failed optimisation
					return null;
				}
				double fit = current.getF();
				int[] feas = isFeasible4(xPrimes, infeas);

				// for (int i = 0; i < xPrimes[0].length; i++)
				// System.out.println(xPrimes[0][i] + "\t" + xPrimes[1][i]);

				// if (feas == null)
				// System.out.println("NO (in)FEAS" + " - " + fit);
				// else
				// System.out.println(Arrays.toString(feas) + " - " + fit);

				if (fit < fBar)
					// Solution is better than current best
					if (feas == null) {
					// Solution is feasible
					fBar = fit;
					fBarReductions++;
					xBar = xPrimes;
					adjBar = current.getAdjs();
					} else {
					// Solution not feasible - make it so!
					int negIndex = feas[0];
					int posIndex = feas[1];
					remaining.add(current.split(posIndex, negIndex, 0));
					remaining.add(current.split(negIndex, posIndex, 1));
					}
			}
		}

		if (listener != null) {
			// Tell GUI what's going on
			listener.step(xBar, fBar);
		}

		return new CMRSolution(fBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls(), fBarReductions);
	}

	/**
	 * modified version of isFeasible3n - given arbitrary ordered y vectors,
	 * will determine the largest non coupled monotonic pair.
	 * 
	 * @param y
	 *            2d array of means from StatsSTA
	 * @param infeas
	 *            zone numbers corresponding to infeasible points. For example,
	 *            [+1 -1] corresponds to zone number 2
	 * @return null if feasible or the largest inversion if not (ZERO INDEXED!)
	 */
	int[] isFeasible4(double[][] y, TIntHashSet infeas) {
		int nvar = y.length;
		int ny = y[0].length;

		int[][] zoneNumbers = new int[ny][ny];
		DoubleMatrix2D volumes = new DenseDoubleMatrix2D(ny, ny);
		volumes.assign(1);

		for (int i = 0; i < nvar; i++) {
			double[] curY = y[i];

			for (int row = ny; --row >= 0;)
				for (int column = ny; --column >= 0;) {
					double diff = curY[row] - curY[column];
					if (Math.abs(diff) > ZERO_TOL) {
						volumes.setQuick(row, column, volumes.getQuick(row, column) * diff);
						zoneNumbers[row][column] += (int) (Math.signum(diff)) * CMRUtil.POW_3[nvar - i - 1];
					}
				}
		}

		double maxVol = 0;
		int maxRow = -1;
		int maxColumn = -1;

		for (int row = ny; --row >= 0;) {
			for (int column = ny; --column >= 0;) {
				int curZone = Math.abs(zoneNumbers[row][column]);
				if (infeas.contains(curZone)) {
					double curVol = Math.abs(volumes.getQuick(row, column));
					if (curVol > maxVol) {
						maxVol = curVol;
						maxRow = row;
						maxColumn = column;
					}
				}
			}
		}
		if (maxRow != -1)
			return new int[] { maxColumn, maxRow };
		return null;
	}
}
