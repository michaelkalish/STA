package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import org.junit.Test;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.CatastrophicMRFailure;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.mr.MRUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MRTests {

	@Test
	public void solnTest() throws CatastrophicMRFailure {
		int n = 50;
		double[] y = new double[n];
		int[][] ranges = new int[1][n];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
			ranges[0][im1] = im1;
		}

		MRProblem p = new MRProblem(y, ranges);

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution solution = solver.solve(p);

		double[] xv = solution.getxVector();
		for (int i = 0; i < n - 1; i++) {
			// Ensure feasibility
			assertTrue(xv[i] <= xv[i + 1]);
		}
		// System.out.println(" E=263.4844337812948");
		// System.out.println("JOptimiser=" + solution.getfVal());
		// System.out.println(solution.getTimeS());
		assertEquals(263.4844337812948, solution.getfVal(), 1e-5);

		// Solve with oj! (oj! is broken)
		// solver = new MRSolverOJAlgo(1e-10);
		// solution = solver.solve(p);
		//
		// xv = solution.getxVector();
		// for (int i = 0; i < n - 1; i++) {
		// // Ensure feasibility
		// assertTrue(xv[i] <= xv[i + 1]);
		// }
		// System.out.println(" E=263.4844337812948");
		// System.out.println("oj!=" + solution.getfVal());
		// System.out.println(Arrays.toString(y));
		// System.out.println(Arrays.toString(xv));
		// // assertEquals(263.4844337812948, solution.getfVal(), 1e-5);
	}

	@Test
	public void badStartTest() {
		double[] y = { 0.9568965517241379, 0.853448275862069, 0.9921875, 0.890625, 0.9814814814814815, 0.9259259259259259,
				0.9714285714285714, 0.8357142857142857, 0.925, 0.7916666666666666, 0.8972222222222223, 0.7972222222222223,
				0.9396551724137931, 0.8189655172413793, 0.9615384615384616, 0.8269230769230769, 0.9435483870967742, 0.7983870967741935,
				0.9583333333333334, 0.8541666666666666 };
		int[][] rangeSet = { { 15, 17 }, { 1, 3 }, { 3, 17 }, { 17, 15 }, { 4, 0 }, { 3, 13 }, { 5, 13 }, { 13, 1 } };
		MRProblem p = new MRProblem(y, rangeSet);

		HashSet<SimpleLinearConstraint> eq = new HashSet<>();
		HashSet<SimpleLinearConstraint> ineq = new HashSet<>();

		double[] start = p.findFeasibleStart(1, ineq, eq);
		assertTrue(p.testFeas(start));
	}

	@Test
	public void doublePointDAGTest() {
		double[] y = { 0.330833333333333, 0.455, 0.534583333333333, 0.549166666666667, 0.28359375, 0.303125, 0.31796875,
				0.309765625 };
		int[][] rangeSet = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 1, 7 } };
		MRProblem p = new MRProblem(y, rangeSet);
		HashSet<SimpleLinearConstraint> ineqCon = new HashSet<>();
		HashSet<SimpleLinearConstraint> eqCon = new HashSet<>();
		assertTrue(p.testFeas(p.findFeasibleStart(0.01, ineqCon, eqCon)));
	}

	@Test
	public void solnRandomTest() throws CatastrophicMRFailure {
		Random rand = new Random(1337);
		double secs = 0;

		for (int r = 0; r < 10; r++) {
			int n = 10 + r / 20;

			double[] y = new double[n];
			int[][] ranges = new int[1][n];
			for (int i = 1; i <= n; i++) {
				int im1 = i - 1;
				y[im1] = rand.nextGaussian();
				ranges[0][im1] = im1;
			}

			MRProblem p = new MRProblem(y, ranges);

			HashSet<SimpleLinearConstraint> ineqCon = new HashSet<>();
			HashSet<SimpleLinearConstraint> eqCon = new HashSet<>();
			double[] initial = p.findFeasibleStart(0.01, ineqCon, eqCon);
			assertTrue(p.testFeas(initial));

			MRSolver solver = new MRSolverAJOptimiser();
			MRSolution solution = solver.solve(p);

			if (solution == null) {
				System.out.println(Arrays.toString(y));
			} else {
				secs += solution.getTimeS();
				double[] xv = solution.getxVector();
				for (int i = 0; i < n - 1; i++) {
					// Ensure feasibility
					assertTrue(xv[i] <= xv[i + 1]);
				}
			}
		}
		// System.out.println(secs);
	}

	@Test
	public void badProblemTest() throws CatastrophicMRFailure {
		double[][] weights = {
				{ 173.307369, 79.14989, 26.897914, 117.854693, 55.739509, 43.013188, 36.16501, 38.71727, 36.150232,
						55.906534, 67.877555, 4.227007 },
				{ 79.14989, 177.65603, 12.19409, 105.591342, 117.749225, 28.133631, 0.247612, 4.073618, 19.379166,
						56.511404, 33.051736, 7.850132 },
				{ 26.897914, 12.19409, 215.055741, 84.277943, 115.508946, 104.33712, 118.805283, 110.773987, 117.127714,
						96.671801, 93.955705, 36.752727 },
				{ 117.854693, 105.591342, 84.277943, 364.820278, 147.488265, 87.927012, 67.833489, 74.229795,
						106.003233, 137.817021, 174.081407, 13.274677 },
				{ 55.739509, 117.749225, 115.508946, 147.488265, 527.354692, 122.481887, -6.620381, 57.110282,
						125.051188, 161.387466, 97.088822, 56.715019 },
				{ 43.013188, 28.133631, 104.33712, 87.927012, 122.481887, 170.137358, 115.665601, 88.949684, 112.797962,
						82.449866, 83.281253, 44.616607 },
				{ 36.16501, 0.247612, 118.805283, 67.833489, -6.620381, 115.665601, 470.626178, 142.644778, 117.814058,
						28.224717, 89.057712, 52.137383 },
				{ 38.71727, 4.073618, 110.773987, 74.229795, 57.110282, 88.949684, 142.644778, 177.642192, 94.223102,
						70.790498, 103.529847, 33.436765 },
				{ 36.150232, 19.379166, 117.127714, 106.003233, 125.051188, 112.797962, 117.814058, 94.223102,
						201.610931, 94.487843, 102.150004, 40.397261 },
				{ 55.906534, 56.511404, 96.671801, 137.817021, 161.387466, 82.449866, 28.224717, 70.790498, 94.487843,
						196.265109, 127.273198, 24.683237 },
				{ 67.877555, 33.051736, 93.955705, 174.081407, 97.088822, 83.281253, 89.057712, 103.529847, 102.150004,
						127.273198, 363.761642, 18.753787 },
				{ 4.227007, 7.850132, 36.752727, 13.274677, 56.715019, 44.616607, 52.137383, 33.436765, 40.397261,
						24.683237, 18.753787, 44.273899 } };

		double[] means = { 78.72100621871694, 75.25525264520714, 6.6402333440731764, 90.0355933088869,
				40.37485195861801, 24.309873738591715, 26.66177277906867, 34.18957087473541, 6.19308278575369,
				92.2218631244646, 80.53179366777084, 54.89136683289008 };

		HashSet<SimpleLinearConstraint> cons = new HashSet<>();
		cons.add(new SimpleLinearConstraint(means.length, 11, 8));
		cons.add(new SimpleLinearConstraint(means.length, 0, 5));
		cons.add(new SimpleLinearConstraint(means.length, 7, 2));
		cons.add(new SimpleLinearConstraint(means.length, 1, 8));

		MRProblem p = new MRProblem(means, new DenseDoubleMatrix2D(weights), cons);
		// System.out.println(p.toMatlabString());
		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution s = solver.solve(p);
		assertEquals(623312.662257973, s.getfVal(), 1e-6);

		// Test the tolerance parameter
		solver.setTolerance(1.0);
		s = solver.solve(p);
		assertEquals(623312.662257973, s.getfVal(), 1e-2);
	}

	@Test
	public void equalityConstraintsOneTest() throws CatastrophicMRFailure {
		int n = 100;
		double[] y = new double[n];
		int[][] ranges = new int[2][n];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
			ranges[0][im1] = im1;
		}
		ranges[1] = new int[] { 79, 19 };

		MRProblem p = new MRProblem(y, ranges);

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution solution = solver.solve(p);

		assertEquals(22474.8106931395, solution.getfVal(), 1e-5);
	}

	@Test
	public void equalityConstraintsTinyATest() throws CatastrophicMRFailure {
		int n = 10;
		double[] y = new double[n];
		int[][] ranges = new int[2][n];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
			ranges[0][im1] = im1;
		}
		ranges[1] = new int[] { 5, 1 };

		MRProblem p = new MRProblem(y, ranges);

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution solution = solver.solve(p);

		assertEquals(5.2836542015136, solution.getfVal(), 1e-8);
	}

	/**
	 * Test using the Original JOptimizer (normally not applicabble)
	 * 
	 * @throws CatastrophicMRFailure
	 */
	// @Test
	public void equalityConstraintsTinyOriginalTest() throws CatastrophicMRFailure {
		int n = 10;
		double[] y = new double[n];
		int[][] ranges = new int[2][n];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
			ranges[0][im1] = im1;
		}
		ranges[1] = new int[] { 5, 1 };

		MRProblem p = new MRProblem(y, ranges);

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution solution = solver.solve(p);

		assertEquals(5.2836542015136, solution.getfVal(), 1e-11);
	}

	/**
	 * Tests overlapping equality cycles
	 * 
	 * @throws CatastrophicMRFailure
	 */
	@Test
	public void realDataEqualityTest() throws CatastrophicMRFailure {
		int n = 16;
		double[] wBase = { 0.046913, 0.027412, 0.053279, 0.066098, 0.73046, 81.632653, 0.138408, 0.170874, 0.114387,
				0.11562, 0.052893, 0.065564, 2.163332, 1.189061, 0.145159, 0.066098 };
		double[] y = { 65.6, 48.6, 71.9, 62.4, 88.9, 99.0, 78.0, 89.8, 80.4, 87.7, 80.6, 74.9, 95.7, 96.4, 82.7, 76.6 };

		// cons = [x11 - x7, x6 - x15, x7 - x14, x4 - x10, x7 - x10, x8 - x11,
		// x9 - x11, x7 - x8, x10 - x11]

		HashSet<SimpleLinearConstraint> matA = new HashSet<>();
		matA.add(new SimpleLinearConstraint(n, 11, 7));
		matA.add(new SimpleLinearConstraint(n, 6, 15));
		matA.add(new SimpleLinearConstraint(n, 7, 14));
		matA.add(new SimpleLinearConstraint(n, 4, 10));
		matA.add(new SimpleLinearConstraint(n, 7, 10));
		matA.add(new SimpleLinearConstraint(n, 8, 11));
		matA.add(new SimpleLinearConstraint(n, 9, 11));
		matA.add(new SimpleLinearConstraint(n, 7, 8));
		matA.add(new SimpleLinearConstraint(n, 10, 11));

		DoubleMatrix2D weights = new DenseDoubleMatrix2D(n, n);
		for (int i = 0; i < n; i++)
			weights.set(i, i, wBase[i]);

		MRProblem p = new MRProblem(y, weights, matA);

		double[] init = p.findFeasibleStart(10, new HashSet<SimpleLinearConstraint>(),
				new HashSet<SimpleLinearConstraint>());
		assertTrue(p.testFeas(init));

		MRSolverAJOptimiser solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);
		assertEquals(23.3255481711711, sol.getfVal(), 1e-5);
	}

	/**
	 * This test currently fails (due to no inequality constraints). This case
	 * should never arise in the context of CMR(x) and is therefore not handled
	 * 
	 * @throws CatastrophicMRFailure
	 */
	// @Test
	public void onlyEqualityTest() throws CatastrophicMRFailure {
		int n = 16;
		double[] wBase = { 0.046913, 0.027412, 0.053279, 0.066098, 0.73046, 81.632653, 0.138408, 0.170874, 0.114387,
				0.11562, 0.052893, 0.065564, 2.163332, 1.189061, 0.145159, 0.066098 };
		double[] y = { 65.6, 48.6, 71.9, 62.4, 88.9, 99.0, 78.0, 89.8, 80.4, 87.7, 80.6, 74.9, 95.7, 96.4, 82.7, 76.6 };

		// cons = [x11 - x7, x6 - x15, x7 - x14, x4 - x10, x7 - x10, x8 - x11,
		// x9 - x11, x7 - x8, x10 - x11]

		HashSet<SimpleLinearConstraint> matA = new HashSet<>();
		matA.add(new SimpleLinearConstraint(n, 0, 1));
		matA.add(new SimpleLinearConstraint(n, 1, 2));
		matA.add(new SimpleLinearConstraint(n, 2, 0));
		matA.add(new SimpleLinearConstraint(n, 1, 5));
		matA.add(new SimpleLinearConstraint(n, 5, 0));

		DoubleMatrix2D weights = new DenseDoubleMatrix2D(n, n);
		for (int i = 0; i < n; i++)
			weights.set(i, i, wBase[i]);

		MRProblem p = new MRProblem(y, weights, matA);

		double[] init = p.findFeasibleStart(10, new HashSet<SimpleLinearConstraint>(),
				new HashSet<SimpleLinearConstraint>());

		assertTrue(p.testFeas(init));

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);

		double[] expected = { 98.9462780558617,
				98.9462780558691,
				98.9462780558618,
				62.3999999995485,
				88.8999999993536,
				98.9462780558666,
				77.9999999994338,
				89.799999999347,
				80.3999999994162,
				87.6999999993625,
				80.5999999994147,
				74.8999999994566,
				95.6999999993036,
				96.3999999992985,
				82.6999999993992,
				76.5999999994441 };

		assertEquals(160.857798628857, sol.getfVal(), 1e-10);
		assertArrayEquals(expected, sol.getxVector(), 1e-10);
	}

	@Test
	public void brickOfEqualityTest() throws CatastrophicMRFailure {
		int n = 16;
		double[] wBase = { 0.046913, 0.027412, 0.053279, 0.066098, 0.73046, 81.632653, 0.138408, 0.170874, 0.114387,
				0.11562, 0.052893, 0.065564, 2.163332, 1.189061, 0.145159, 0.066098 };
		double[] y = { 65.6, 48.6, 71.9, 62.4, 88.9, 99.0, 78.0, 89.8, 80.4, 87.7, 80.6, 74.9, 95.7, 96.4, 82.7, 76.6 };

		HashSet<SimpleLinearConstraint> matA = new HashSet<>();
		int nMaxC = 14;
		for (int i = 0; i < nMaxC; i++)
			for (int j = i + 1; j < nMaxC; j++)
				matA.add(new SimpleLinearConstraint(n, j, i));

		matA.add(new SimpleLinearConstraint(n, 0, nMaxC - 1));
		matA.add(new SimpleLinearConstraint(n, nMaxC, nMaxC - 1));

		DoubleMatrix2D weights = new DenseDoubleMatrix2D(n, n);
		for (int i = 0; i < n; i++)
			weights.set(i, i, wBase[i]);

		MRProblem p = new MRProblem(y, weights, matA);

		double[] init = p.findFeasibleStart(10, new HashSet<SimpleLinearConstraint>(),
				new HashSet<SimpleLinearConstraint>());

		assertTrue(p.testFeas(init));

		MRSolver solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);

		assertEquals(527.5135336748779, sol.getfVal(), 1e-10);
	}

	/**
	 * Tests only a triple cycle
	 * 
	 * @throws CatastrophicMRFailure
	 */
	@Test
	public void tripleEqualityTest() throws CatastrophicMRFailure {
		int n = 16;
		double[] wBase = { 0.046913, 0.027412, 0.053279, 0.066098, 0.73046, 81.632653, 0.138408, 0.170874, 0.114387,
				0.11562, 0.052893, 0.065564, 2.163332, 1.189061, 0.145159, 0.066098 };
		double[] y = { 65.6, 48.6, 71.9, 62.4, 88.9, 99.0, 78.0, 89.8, 80.4, 87.7, 80.6, 74.9, 95.7, 96.4, 82.7, 76.6 };

		// cons = [x11 - x7, x6 - x15, x7 - x14, x4 - x10, x7 - x10, x8 - x11,
		// x9 - x11, x7 - x8, x10 - x11]

		HashSet<SimpleLinearConstraint> matA = new HashSet<>();
		matA.add(new SimpleLinearConstraint(n, 0, 1));
		matA.add(new SimpleLinearConstraint(n, 1, 2));
		matA.add(new SimpleLinearConstraint(n, 2, 0));

		matA.add(new SimpleLinearConstraint(n, 1, 5));
		matA.add(new SimpleLinearConstraint(n, 5, 0));

		matA.add(new SimpleLinearConstraint(n, 1, 9));
		matA.add(new SimpleLinearConstraint(n, 9, 5));

		matA.add(new SimpleLinearConstraint(n, 15, 14));

		DoubleMatrix2D weights = new DenseDoubleMatrix2D(n, n);
		for (int i = 0; i < n; i++)
			weights.set(i, i, wBase[i]);

		MRProblem p = new MRProblem(y, weights, matA);

		double[] init = p.findFeasibleStart(10, new HashSet<SimpleLinearConstraint>(),
				new HashSet<SimpleLinearConstraint>());

		assertTrue(p.testFeas(init));

		MRSolverAJOptimiser solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);
		assertEquals(175.460623672191, sol.getfVal(), 1e-10);
	}

	/**
	 * Tests only a triple cycle
	 * 
	 * @throws CatastrophicMRFailure
	 */
	@Test
	public void quadrupleEqualityTest() throws CatastrophicMRFailure {
		int n = 30;
		double[] y = new double[n];
		int[][] ranges = new int[4][n];
		for (int i = 1; i <= n; i++) {
			int im1 = i - 1;
			y[im1] = (double) i + Math.sin(i) * n / 10;
			ranges[0][im1] = im1;
		}
		ranges[1] = new int[] { 5, 1 };
		ranges[2] = new int[] { 7, 3 };
		ranges[3] = new int[] { 19, 7 };

		MRProblem p = new MRProblem(y, ranges);

		double[] initial = p.findFeasibleStart(1, new HashSet<SimpleLinearConstraint>(),
				new HashSet<SimpleLinearConstraint>());

		assertTrue(p.testFeas(initial));
		MRSolverAJOptimiser solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);

		// Strictly speaking, MRSolverAJOptimiser will only give 6 decimal
		// places...
		assertEquals(667.766045406726, sol.getfVal(), 1e-3);
	}

	public void evilDataTest() throws CatastrophicMRFailure {
		double[] y = { 0.4893617021276593, 0.2978723404255321, 0.9361702127659575, 0.6879432624113477, 0.9290780141843976,
				0.5319148936170214, 0.8581560283687948, 0.5602836879432626, 0.9503546099290782, 0.7375886524822698, 0.9432624113475184,
				0.8226950354609927 };
		DoubleMatrix2D w = new DenseDoubleMatrix2D(new double[][] {
				{ 499.5046983020609, -242.8198555351167, -504.9319114464205, 56.97854804305997, -80.81179861180473, 65.72218377195566,
						110.8904445383659, -17.32660652645919, -93.31749940819267, -77.91018388464946, 306.1559519544281,
						76.67323780554035 },
				{ -242.8198555351167, 587.6777274189657, 253.1854534226979, -126.6738683948532, 65.52694526170335, -158.7317593366188,
						-218.7938455655018, 26.18889702669812, 59.34396173943481, 121.3305783697826, -163.2750439178449,
						-27.27510303019637 },
				{ -504.9319114464206, 253.1854534226978, 3365.64424007323, -45.56771176954812, -298.7988470458105, -170.946725049071,
						-8.384323834324263, -63.74159750176837, -756.3463641903754, 78.7436666124878, -478.3057368844885, 219.50683684999 },
				{ 56.97854804305997, -126.6738683948531, -45.5677117695481, 519.0954066922626, -67.66240662856458, -154.5381932081344,
						249.5899196389197, -75.41323727572726, -217.5692844735258, -4.914356536381657, -154.0349063307196,
						-46.20186664509577 },
				{ -80.81179861180472, 65.52694526170332, -298.7988470458104, -67.66240662856454, 2973.655452383485, -116.8208483852999,
						244.028207909786, 222.7685060730922, 95.22551521825147, -193.3124776645275, -1642.55071987272, 18.5966076725447 },
				{ 65.72218377195567, -158.7317593366188, -170.9467250490711, -154.5381932081344, -116.8208483853, 509.3494057079729,
						-339.741807304412, -38.40840458414962, 341.654722334674, -25.83465972398282, 384.4469093242443,
						-25.82499144561633 },
				{ 110.8904445383659, -218.7938455655018, -8.384323834324134, 249.5899196389196, 244.028207909786, -339.741807304412,
						2412.309072937001, -167.4932448943234, -1096.396244778376, -500.1939062366727, -1033.094524908262,
						608.8267394785238 },
				{ -17.32660652645918, 26.18889702669808, -63.74159750176836, -75.41323727572718, 222.7685060730922, -38.40840458414963,
						-167.4932448943234, 628.0642170525244, -274.7502881563543, -241.7647997745418, -177.3126473188264,
						-231.9087344690153 },
				{ -93.31749940819265, 59.3439617394348, -756.3463641903754, -217.5692844735257, 95.22551521825142, 341.654722334674,
						-1096.396244778376, -274.7502881563544, 4143.278525846607, 244.519748104209, 376.2186585955296, -211.277083050609 },
				{ -77.91018388464946, 121.3305783697826, 78.74366661248767, -4.91435653638166, -193.3124776645275, -25.83465972398284,
						-500.1939062366726, -241.7647997745417, 244.519748104209, 761.478870057487, 206.4775845970088, -334.3803354432869 },
				{ 306.1559519544282, -163.2750439178449, -478.3057368844885, -154.0349063307196, -1642.55071987272, 384.4469093242443,
						-1033.094524908262, -177.3126473188264, 376.2186585955299, 206.4775845970088, 7.781822460176559e+31,
						6.712185478291982 },
				{ 76.67323780554031, -27.27510303019634, 219.5068368499901, -46.2018666450958, 18.59660767254466, -25.82499144561631,
						608.8267394785237, -231.9087344690153, -211.277083050609, -334.3803354432868, 6.712185478291982, 913.7217568147602 }
		});

		MRProblem p = new MRProblem(y, w, new int[][] { { 10, 5 } });

		MRSolverAJOptimiser solver = new MRSolverAJOptimiser();
		MRSolution sol = solver.solve(p);
		assertEquals(1533.30319481554, sol.getfVal(), 1e-10);
	}

	@Test
	public void covTest() {
		DoubleMatrix2D m = new DenseDoubleMatrix2D(new double[][] { { 1, Double.NaN, 2 },
				{ Double.NaN, 33, 2 },
				{ 4, 3, 2 },
				{ 4, 3, 88 }
		});

		DoubleMatrix2D c = MRUtil.covariance(m);
		DenseDoubleMatrix2D expected = new DenseDoubleMatrix2D(new double[][] { { 2, 0, 28.66666666666667 }, { 0, 200, -286.6666666666667 },
				{ 28.66666666666667, -286.6666666666667, 1386.75 }
		});

		assertTrue(c.equals(expected));
	}
}
