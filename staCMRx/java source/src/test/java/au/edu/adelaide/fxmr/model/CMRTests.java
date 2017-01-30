package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Arrays;

import org.junit.Test;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import gnu.trove.set.hash.TIntHashSet;

public class CMRTests {

	private static final double[][][] WEIGHT_BASE = {
			{ { 4909.21477287896, -2498.46049296539, 287.220544418912, -74.8884624739565, 0, 0, 0, 0 },
					{ -2498.46049296539, 3232.67047265783, -1204.08834504430, -379.675841004848, 0, 0, 0, 0 },
					{ 287.220544418912, -1204.08834504430, 2727.57207512732, -1349.75597766002, 0, 0, 0, 0 },
					{ -74.8884624739565, -379.675841004848, -1349.75597766002, 1780.38690501671, 0, 0, 0, 0 },
					{ 0, 0, 0, 0, 3192.72426569220, -1210.68772553161, 436.261656263469, -356.435549731305 },
					{ 0, 0, 0, 0, -1210.68772553161, 3460.62316379029, -1937.45853848169, -669.580010152516 },
					{ 0, 0, 0, 0, 436.261656263469, -1937.45853848169, 4287.83548539312, -2319.49072628248 },
					{ 0, 0, 0, 0, -356.435549731305, -669.580010152516, -2319.49072628248, 3033.35421927641 } },
			{ { 5917.79594564311, -2643.97092422918, -307.525233751586, -204.057587400608, 0, 0, 0, 0 },
					{ -2643.97092422918, 3345.70933946505, -1229.73039263390, -12.6525122869357, 0, 0, 0, 0 },
					{ -307.525233751586, -1229.73039263390, 4018.69584375822, -2498.85535725609, 0, 0, 0, 0 },
					{ -204.057587400608, -12.6525122869357, -2498.85535725609, 2771.05064789733, 0, 0, 0, 0 },
					{ 0, 0, 0, 0, 6023.17355125993, -954.924303309050, -777.728064087585, 226.297642437044 },
					{ 0, 0, 0, 0, -954.924303309050, 4904.97576716846, -1456.50320237336, -1186.22427010323 },
					{ 0, 0, 0, 0, -777.728064087585, -1456.50320237336, 4180.29811680008, -1616.19733790714 },
					{ 0, 0, 0, 0, 226.297642437044, -1186.22427010323, -1616.19733790714, 3192.47268741585 } } };

	private static final double[][] MEANS = {
			{ 0.367647058823529, 0.467647058823529, 0.575735294117647, 0.611764705882353, 0.344485294117647,
					0.443382352941177, 0.508088235294118, 0.516911764705882 },
			{ 0.330833333333333, 0.455000000000000, 0.534583333333333, 0.549166666666667, 0.283593750000000,
					0.303125000000000, 0.317968750000000, 0.309765625000000 } };

	@Test
	public void isFeasible3nTest() {
		double[][] y = {
				{ 0.367647058823716, 0.467647058823891, 0.575735294117454, 0.611764705886237, 0.344485294143725,
						0.443382352976257, 0.508088235295267, 0.516911764839725 },
				{ 0.330833333333334, 0.455000000000000, 0.534583333333334, 0.549166666666667, 0.283044120061215,
						0.303360912877124, 0.314923129852751, 0.314923129852751 } };
		TIntHashSet infeas = CMRUtil.calcInfeasibleZoneNumbers(2);

		CMRSolver cmrSolver = new CMRSolver();

		int[] idx = cmrSolver.isFeasible4(y, infeas);
		assertArrayEquals(new int[] { 1, 7 }, idx);

		idx = cmrSolver.isFeasible4(new double[][] { { 1, 2, 3, 4, 5 }, { 1, 2, 3, 4, 5 } }, infeas);
		assertNull(idx);
		idx = cmrSolver.isFeasible4(new double[][] { { 1, 2, 3, 4, 5 }, { 1, 2, 4, 3, 5 } }, infeas);
		assertArrayEquals(new int[] { 2, 3 }, idx);

		infeas = CMRUtil.calcInfeasibleZoneNumbers(4);
		idx = cmrSolver.isFeasible4(new double[][] { { 1, 2, 3, 4, 5 }, { 1, 2, 3, 4, 5 }, { 1, 2, 3, 4, 5 } }, infeas);
		assertNull(idx);

		idx = cmrSolver.isFeasible4(new double[][] { { 1, 2, 3, 4, 5 }, { 1, 2, 4, 3, 5 }, { 1, 2, 4, 3, 5 } }, infeas);
		assertArrayEquals(new int[] { 2, 3 }, idx);
	}

	@Test
	public void getInfeasibleTest() {
		int[][] expected = { { 2 }, { 2, 5, 6, 7, 8, 11 },
				{ 2, 5, 6, 7, 8, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 29, 32, 33, 34, 35, 38 } };

		for (int i = 0; i < expected.length; i++) {
			int[] result = CMRUtil.calcInfeasibleZoneNumbers(i + 2).toArray();
			Arrays.sort(result);
			assertArrayEquals(expected[i], result);
		}
	}

	@Test
	public void cmrWithInitialOrderingTest() {

		int[][] rangeSet = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };

		DoubleMatrix2D[] weights = new DoubleMatrix2D[WEIGHT_BASE.length];
		for (int i = 0; i < WEIGHT_BASE.length; i++)
			weights[i] = new SparseDoubleMatrix2D(WEIGHT_BASE[i]);

		CMRProblem problem = new CMRProblem(MEANS, weights, rangeSet);
		CMRSolver solver = new CMRSolver();
		CMRSolution solution = solver.solve(problem);

		assertEquals(1.749306505018210, solution.getFStar(), 1e-8);

		double[][] expectedXStar = {
				{ 0.375896260823335, 0.485035432099519, 0.589841619860198, 0.626514203714855, 0.335286414933073,
						0.41859656830624, 0.48058130590643, 0.485035432095858 },
				{ 0.314984500052892, 0.435834962941795, 0.51669581542477, 0.531781607126321, 0.285597306945978,
						0.314984500052892, 0.322682298846644, 0.322682298846644 } };

		for (int i = 0; i < expectedXStar.length; i++) {
			assertArrayEquals(expectedXStar[i], solution.getXStar()[i], 1e-8);
		}
	}

	@Test
	public void cmrWithNoInitialOrderTest() {
		int[][] rangeSet = null;

		DoubleMatrix2D[] weights = new DoubleMatrix2D[WEIGHT_BASE.length];
		for (int i = 0; i < WEIGHT_BASE.length; i++)
			weights[i] = new SparseDoubleMatrix2D(WEIGHT_BASE[i]);

		CMRProblem problem = new CMRProblem(MEANS, weights, rangeSet);
		CMRSolver solver = new CMRSolver();
		CMRSolution solution = solver.solve(problem);

		assertEquals(1.65296751647382, solution.getFStar(), 1e-5);

		double[][] expectedXStar = {
				{ 0.375138676211893, 0.483438530482393, 0.588546132695485, 0.625159649348022, 0.335195015927247,
						0.419855256393548, 0.48343853050931, 0.483438530482387 },
				{ 0.314854809066627, 0.435678135463463, 0.516549441884311, 0.531639345202168, 0.286141954638518,
						0.314854809066627, 0.325648033366418, 0.317831087527433 } };

		for (int i = 0; i < expectedXStar.length; i++) {
			assertArrayEquals(expectedXStar[i], solution.getXStar()[i], 1e-6);
		}
	}
}
