package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

import au.edu.adelaide.fxmr.joptimizer.algebra.CholeskyMachine;
import au.edu.adelaide.fxmr.model.mr.CatastrophicMRFailure;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * This class contains disabled tests for comparing times of JOptimiser. They
 * are only really used to measure performance differences as opposed to actual
 * unit tests.
 * 
 * Just add an @ Test annotation to time JOptimiser.
 * 
 * 
 */
public class JoptimiserTests {

	@Test
	public void makeSureAnnotationIsImportedTest() {
	}

	//@Test
	public void averagePerformanceTest() throws CatastrophicMRFailure {
		int n = 630;
		// Avg= 0.10855242468000005
		// Avg= 0.11103568089999998
		// Avg= 0.10755750362
		// Avg= 0.10110787809800005 (500)

		double[] y = new double[n];

		int nRange = n / 2;

		int[][] ranges = new int[2][nRange];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
		}

		for (int i = 1; i <= nRange; i++) {
			int im1 = i - 1;
			ranges[0][im1] = im1;
			ranges[1][im1] = im1 + n / 2;
		}

		MRProblem p = new MRProblem(y, ranges);

		MRSolverAJOptimiser solver = new MRSolverAJOptimiser();
		MRSolution solution = solver.solve(p);

		int nI = 5000;
		for (int i = 0; i < nI; i++) {
			solution = solver.solve(p);
			System.out.println(solution.getTimeS() + ",\t" + solution.getIterations());
		}

		// System.out.println("Avg=\t" + sum / nI);
	}

	private MRProblem problemN(int n, int c) {
		double[] y = new double[n];

		int[][] ranges = new int[1][c];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
		}

		for (int i = 1; i <= c; i++) {
			int im1 = i - 1;
			ranges[0][im1] = im1;
		}

		return new MRProblem(y, ranges);
	}

	@Test
	public void choleskyMachineTest() {
		double[][] exp = { { 1.41421356237310 }, { 0.707106781186548, 1.22474487139159 },
				{ 0.707106781186548, 2.04124145231932, 1.15470053837925 },
				{ 0, 3.26598632371090, 2.88675134594813, 1.41421356237309 } };

		double[][] matA = { { 1, 1, 1, 1 }, { 1, 2, 3, 4 }, { 1, 3, 6, 10 }, { 1, 4, 10, 20 } };

		CholeskyMachine machine = new CholeskyMachine(new DenseDoubleMatrix2D(matA));
		double[][] L = machine.getLDataCur();
		machine.reset();

		DenseDoubleMatrix1D x = new DenseDoubleMatrix1D(4);
		x.set(0, 1);
		x.set(3, -1);
		machine.rank1Update(x, 1);

		for (int i = 0; i < exp.length; i++) {
			assertArrayEquals(exp[i], L[i], 1e-14);
		}

		exp = new double[][] { { 1.24498995979887 }, { 0.803219328902499, 1.16397539049476 },
				{ 0.803219328902499, 2.02310008347898, 1.12334534400814 },
				{ 0.361448698006124, 3.18707547397374, 2.90374173828519, 1.35284466190516 } };
		machine = new CholeskyMachine(new DenseDoubleMatrix2D(matA));
		L = machine.getLDataCur();
		machine.reset();

		x = new DenseDoubleMatrix1D(4);
		x.set(0, 1);
		x.set(3, -1);
		machine.rank1Update(x, Math.sqrt(.55));

		for (int i = 0; i < exp.length; i++) {
			assertArrayEquals(exp[i], L[i], 1e-14);
		}
	}
}
