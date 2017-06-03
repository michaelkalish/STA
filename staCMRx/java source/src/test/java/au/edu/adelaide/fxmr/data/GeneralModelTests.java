package au.edu.adelaide.fxmr.data;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.HashSet;

import org.junit.Test;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.CMRSolution;
import au.edu.adelaide.fxmr.model.CMRxGMFits;
import au.edu.adelaide.fxmr.model.CMRxGMProblemMaker;
import au.edu.adelaide.fxmr.model.CMRxProblem;
import au.edu.adelaide.fxmr.model.CMRxProblemMaker;
import au.edu.adelaide.fxmr.model.CMRxSolver;
import au.edu.adelaide.fxmr.model.CombinedStatsSTA;
import au.edu.adelaide.fxmr.model.ParCMRxSolver;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolverReverse;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class GeneralModelTests {
	private static final double TOL = 1e-15;

	public static int[][] delayInitData = { { 1, 2, 1 }, { 2, 1, 2 }, { 3, 1, 1 }, { 4, 2, 2 }, { 5, 2, 1 }, { 6, 1, 2 }, { 7, 1, 1 },
			{ 8, 2, 2 }, { 9, 2, 1 },
			{ 10, 2, 1 }, { 11, 1, 1 }, { 12, 2, 2 }, { 13, 2, 1 }, { 14, 1, 2 }, { 15, 1, 1 }, { 16, 2, 2 }, { 17, 2, 1 }, { 18, 1, 2 },
			{ 19, 1, 1 },
			{ 20, 2, 2 }, { 21, 2, 1 }, { 22, 2, 1 }, { 23, 1, 1 }, { 24, 2, 2 }, { 25, 2, 1 }, { 26, 1, 2 }, { 27, 1, 1 }, { 28, 2, 2 },
			{ 29, 2, 1 },
			{ 30, 1, 2 }, { 31, 1, 1 }, { 32, 2, 2 }, { 33, 2, 1 }, { 34, 1, 2 }, { 35, 1, 1 }, { 36, 2, 2 }, { 37, 2, 1 }, { 38, 1, 2 },
			{ 39, 1, 1 },
			{ 40, 2, 2 }, { 41, 2, 1 }, { 42, 1, 2 }, { 43, 1, 1 }, { 44, 2, 2 }, { 45, 2, 1 }, { 46, 1, 2 }, { 47, 1, 1 }, { 48, 2, 2 },
			{ 49, 2, 1 },
			{ 50, 1, 2 }, { 51, 1, 1 }, { 52, 2, 2 }, { 53, 2, 1 }, { 54, 1, 2 }, { 55, 1, 1 }, { 56, 2, 2 }, { 57, 2, 1 }, { 58, 1, 2 },
			{ 59, 1, 1 },
			{ 60, 2, 2 }, { 61, 2, 1 }, { 62, 1, 2 }, { 63, 1, 1 }, { 64, 2, 2 }, { 66, 1, 2 }, { 67, 1, 1 }, { 68, 2, 2 }, { 69, 2, 1 },
			{ 71, 1, 2 },
			{ 72, 1, 1 }, { 73, 2, 1 }, { 74, 1, 2 }, { 75, 1, 1 }, { 76, 2, 2 }, { 77, 2, 1 }, { 78, 1, 2 }, { 79, 1, 1 }, { 80, 2, 2 },
			{ 81, 2, 1 },
			{ 82, 1, 2 }, { 83, 1, 1 }, { 84, 2, 2 }, { 85, 2, 1 }, { 86, 1, 2 }, { 87, 1, 1 }, { 91, 1, 1 }, { 92, 2, 2 }, { 94, 1, 2 },
			{ 95, 1, 1 },
			{ 96, 2, 2 }, { 97, 2, 1 }, { 98, 1, 2 }, { 99, 1, 1 }, { 100, 2, 2 }, { 101, 2, 1 }, { 102, 1, 2 }, { 103, 1, 1 },
			{ 105, 2, 1 }, { 106, 1, 2 },
			{ 107, 1, 1 }, { 108, 2, 2 }, { 109, 2, 1 }, { 110, 1, 2 }, { 111, 1, 1 }, { 112, 2, 2 }, { 113, 2, 1 }, { 115, 1, 1 },
			{ 116, 2, 2 }, { 117, 2, 1 },
			{ 118, 1, 2 }, { 120, 2, 2 }, { 121, 2, 1 }, { 123, 1, 1 }, { 124, 2, 2 }, { 125, 2, 1 }, { 126, 1, 2 }, { 127, 1, 1 },
			{ 128, 2, 2 }, { 129, 2, 1 },
			{ 130, 1, 2 }, { 131, 1, 1 }, { 132, 2, 2 }, { 133, 2, 1 }, { 134, 1, 2 }, { 135, 1, 1 }, { 136, 2, 2 }, { 137, 2, 1 },
			{ 138, 1, 2 }, { 139, 1, 1 },
			{ 140, 2, 2 } };

	public static double[][] delayData = { { 0.3375, 0.8875, 0.9375, 1 }, { 0.2875, 0.475, 0.7, 0.7625 }, { 0.525, 0.6, 0.9, 0.925 },
			{ 0.35, 0.3, 0.575, 0.45 },
			{ 0.2375, 0.2625, 0.225, 0.2875 }, { 0.2125, 0.2875, 0.6875, 0.7625 }, { 0.5125, 0.85, 0.9625, 0.8875 },
			{ 0.4375, 0.55, 0.4625, 0.475 },
			{ 0.2, 0.225, 0.275, 0.25 }, { 0.275, 0.1375, 0.25, 0.2625 }, { 0.3125, 0.2, 0.3375, 0.25 }, { 0.225, 0.2375, 0.2375, 0.3125 },
			{ 0.2625, 0.4375, 0.425, 0.4875 }, { 0.25, 0.3625, 0.275, 0.25 }, { 0.3, 0.3125, 0.5, 0.3625 }, { 0.3, 0.25, 0.2375, 0.2375 },
			{ 0.55, 0.8125, 0.8625, 0.9375 }, { 0.225, 0.325, 0.2375, 0.3125 }, { 0.3, 0.4125, 0.3875, 0.575 },
			{ 0.2375, 0.225, 0.275, 0.225 },
			{ 0.3875, 0.7375, 0.85, 0.8 }, { 0.375, 0.4625, 0.45, 0.4125 }, { 0.375, 0.7, 1, 1 }, { 0.2625, 0.2875, 0.25, 0.325 },
			{ 0.275, 0.3125, 0.5625, 0.7 },
			{ 0.375, 0.5375, 0.6375, 0.6875 }, { 0.175, 0.35, 0.5375, 0.35 }, { 0.2875, 0.4, 0.3875, 0.3875 },
			{ 0.2125, 0.2375, 0.275, 0.275 },
			{ 0.175, 0.2875, 0.4, 0.55 }, { 0.5625, 0.4875, 0.4, 0.5375 }, { 0.225, 0.25, 0.25, 0.2875 }, { 0.2625, 0.8625, 0.975, 0.925 },
			{ 0.2875, 0.2125, 0.2, 0.275 }, { 0.325, 0.475, 0.525, 0.55 }, { 0.3375, 0.175, 0.25, 0.225 }, { 0.275, 0.2125, 0.225, 0.2375 },
			{ 0.175, 0.3, 0.575, 0.6875 }, { 0.35, 0.5625, 0.5125, 0.9 }, { 0.2375, 0.225, 0.225, 0.225 }, { 0.4125, 0.65, 0.8625, 0.9125 },
			{ 0.55, 0.6875, 0.6875, 0.775 }, { 0.275, 0.3, 0.3125, 0.2875 }, { 0.225, 0.3375, 0.375, 0.3375 }, { 0.25, 0.25, 0.25, 0.25 },
			{ 0.65, 0.8875, 0.925, 0.9125 }, { 0.3125, 0.525, 0.55, 0.6125 }, { 0.25, 0.25, 0.2375, 0.25 }, { 0.2, 0.2375, 0.3625, 0.25 },
			{ 0.225, 0.25, 0.2375, 0.2125 }, { 0.225, 0.175, 0.2625, 0.275 }, { 0.2, 0.6125, 0.75, 0.8375 }, { 0.2, 0.2625, 0.25, 0.55 },
			{ 0.45, 0.4625, 0.6, 0.525 }, { 0.25, 0.3375, 0.8125, 0.9625 }, { 0.275, 0.25, 0.175, 0.325 }, { 0.4375, 0.3, 0.225, 0.2625 },
			{ 0.3125, 0.425, 0.4875, 0.3375 }, { 0.35, 0.425, 0.8, 0.6625 }, { 0.225, 0.2, 0.25, 0.25 }, { 0.3875, 0.4375, 0.725, 0.85 },
			{ 0.35, 0.6375, 0.7, 0.85 }, { 0.7, 0.85, 0.8625, 0.925 }, { 0.2875, 0.275, 0.3375, 0.4375 }, { 0.275, 0.2375, 0.225, 0.3375 },
			{ 0.3625, 0.5875, 0.7875, 0.9375 }, { 0.425, 0.325, 0.375, 0.2875 }, { 0.475, 0.875, 0.9375, 0.925 },
			{ 0.2125, 0.325, 0.2625, 0.25 },
			{ 0.5, 0.8625, 0.8875, 0.925 }, { 0.3875, 0.3375, 0.35, 0.2625 }, { 0.2625, 0.5, 0.6625, 0.5875 },
			{ 0.225, 0.375, 0.375, 0.3875 },
			{ 0.2875, 0.2125, 0.3, 0.35 }, { 0.3, 0.3375, 0.3875, 0.4125 }, { 0.475, 0.7625, 0.8375, 0.7625 }, { 0.3375, 0.4, 0.4625, 0.4 },
			{ 0.2, 0.25, 0.2125, 0.1875 }, { 0.4125, 0.25, 0.3625, 0.3125 }, { 0.425, 0.5625, 0.6375, 0.6375 },
			{ 0.25, 0.4125, 0.5375, 0.5 },
			{ 0.5375, 0.5375, 0.475, 0.45 }, { 0.375, 0.6375, 0.9, 0.85 }, { 0.4625, 0.675, 0.675, 0.7375 },
			{ 0.2875, 0.1625, 0.3875, 0.425 },
			{ 0.35, 0.3125, 0.4125, 0.4 }, { 0.1875, 0.3375, 0.2, 0.1125 }, { 0.4375, 0.6625, 0.65, 0.7625 }, { 0.3375, 0.25, 0.4125, 0.4 },
			{ 0.3, 0.2375, 0.3625, 0.2375 }, { 0.3, 0.3125, 0.375, 0.3375 }, { 0.2375, 0.3375, 0.3, 0.2125 },
			{ 0.5125, 0.5625, 0.4625, 0.5125 },
			{ 0.2625, 0.3625, 0.3875, 0.2625 }, { 0.275, 0.1875, 0.275, 0.2875 }, { 0.4875, 0.575, 0.825, 0.775 },
			{ 0.5125, 0.9, 0.95, 0.95 },
			{ 0.35, 0.1875, 0.3625, 0.325 }, { 0.225, 0.6, 0.6, 0.4875 }, { 0.7375, 0.975, 0.9875, 0.975 },
			{ 0.375, 0.425, 0.4125, 0.4375 },
			{ 0.6125, 0.8625, 0.975, 0.9625 }, { 0.3125, 0.35, 0.525, 0.525 }, { 0.375, 0.425, 0.65, 0.775 }, { 0.2, 0.275, 0.175, 0.1875 },
			{ 0.625, 0.95, 0.925, 0.9875 }, { 0.3125, 0.375, 0.375, 0.525 }, { 0.3375, 0.3125, 0.3375, 0.1875 }, { 0.3, 0.175, 0.3, 0.35 },
			{ 0.4375, 0.6875, 0.6125, 0.6 }, { 0.3, 0.175, 0.2125, 0.2375 }, { 0.25, 0.325, 0.3125, 0.25 },
			{ 0.275, 0.2875, 0.4625, 0.7125 },
			{ 0.2625, 0.2625, 0.3, 0.1875 }, { 0.2, 0.3125, 0.275, 0.2375 }, { 0.375, 0.4, 0.6625, 0.8 }, { 0.3, 0.3375, 0.5875, 0.8375 },
			{ 0.2375, 0.25, 0.275, 0.3375 }, { 0.225, 0.25, 0.4125, 0.1875 }, { 0.2375, 0.275, 0.25, 0.3125 },
			{ 0.3, 0.3125, 0.4125, 0.2625 },
			{ 0.2625, 0.4, 0.3375, 0.425 }, { 0.825, 0.8375, 0.7, 0.825 }, { 0.3, 0.3625, 0.5125, 0.4125 },
			{ 0.4625, 0.5125, 0.4125, 0.5375 },
			{ 0.2, 0.2375, 0.2375, 0.2 }, { 0.2625, 0.5125, 0.4375, 0.4125 }, { 0.2375, 0.2, 0.45, 0.375 }, { 0.2125, 0.2875, 0.35, 0.275 },
			{ 0.3375, 0.275, 0.3, 0.2375 } };

	@Test
	public void GeneralModelTest() {
		double[] expectedMeans = { 0.367647058823529, 0.467647058823529, 0.575735294117647, 0.611764705882353, 0.344485294117647,
				0.443382352941177,
				0.508088235294118, 0.516911764705882 };
		double[][] expectedCov = { { 0.017197772491349, 0.022730860726644, 0.017445663927336, 0.018836505190311, 0, 0, 0, 0 },
				{ 0.022730860726644, 0.045423875432526, 0.038870296280277, 0.040642571366782, 0, 0, 0, 0 },
				{ 0.017445663927336, 0.038870296280277, 0.049971885813149, 0.048984915657439, 0, 0, 0, 0 },
				{ 0.018836505190311, 0.040642571366782, 0.048984915657439, 0.063271518166090, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0.019294036548443, 0.024292549740484, 0.022369971885813, 0.024808336937716 },
				{ 0, 0, 0, 0, 0.024292549740484, 0.065397383217993, 0.064194150086505, 0.066844452854671 },
				{ 0, 0, 0, 0, 0.022369971885813, 0.064194150086505, 0.072691933391003, 0.074183066608997 },
				{ 0, 0, 0, 0, 0.024808336937716, 0.066844452854671, 0.074183066608997, 0.083270977508651 } };
		double[][] expectedWeights1 = { { 4909.21477287896, -2498.46049296539, 287.220544418912, -74.8884624739565, 0, 0, 0, 0 },
				{ -2498.46049296539, 3232.67047265783, -1204.0883450443, -379.675841004848, 0, 0, 0, 0 },
				{ 287.220544418912, -1204.0883450443, 2727.57207512732, -1349.75597766002, 0, 0, 0, 0 },
				{ -74.8884624739565, -379.675841004848, -1349.75597766002, 1780.38690501671, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 3192.7242656922, -1210.68772553161, 436.261656263469, -356.435549731305 },
				{ 0, 0, 0, 0, -1210.68772553161, 3460.62316379029, -1937.45853848169, -669.580010152516 },
				{ 0, 0, 0, 0, 436.261656263469, -1937.45853848169, 4287.83548539312, -2319.49072628248 },
				{ 0, 0, 0, 0, -356.435549731305, -669.580010152516, -2319.49072628248, 3033.35421927641 } };
		double[][] expectedWeights2 = { { 5917.79594564311, -2643.97092422918, -307.525233751586, -204.057587400608, 0, 0, 0, 0 },
				{ -2643.97092422918, 3345.70933946505, -1229.73039263390, -12.6525122869357, 0, 0, 0, 0 },
				{ -307.525233751586, -1229.73039263390, 4018.69584375822, -2498.85535725609, 0, 0, 0, 0 },
				{ -204.057587400608, -12.6525122869357, -2498.85535725609, 2771.05064789733, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 6023.17355125993, -954.924303309050, -777.728064087585, 226.297642437044 },
				{ 0, 0, 0, 0, -954.924303309050, 4904.97576716846, -1456.50320237336, -1186.22427010323 },
				{ 0, 0, 0, 0, -777.728064087585, -1456.50320237336, 4180.29811680008, -1616.19733790714 },
				{ 0, 0, 0, 0, 226.297642437044, -1186.22427010323, -1616.19733790714, 3192.47268741585 } };

		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++) {
			if (i < 10)
				gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);
			else
				gm.addData(new Subject(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]));
		}

		CombinedStatsSTA stats[] = gm.calcStats();
		assertEquals(2, stats.length);
		CombinedStatsSTA s1 = stats[0];

		assertEquals(8, s1.getMeans().size());
		assertArrayEquals(expectedMeans, s1.getMeans().toArray(), TOL);

		DoubleMatrix2D covs = s1.getCovariances();
		for (int i = 0; i < covs.rows(); i++) {
			assertArrayEquals(expectedCov[i], covs.viewRow(i).toArray(), TOL);
		}

		assertEquals(covs.columns(), Algebra.ZERO.rank(covs));

		DoubleMatrix2D weights1 = s1.getWeights();
		DoubleMatrix2D weights2 = stats[1].getWeights();
		for (int i = 0; i < covs.rows(); i++) {
			assertArrayEquals(expectedWeights1[i], weights1.viewRow(i).toArray(), 1e-10);
			assertArrayEquals(expectedWeights2[i], weights2.viewRow(i).toArray(), 1e-10);
		}
	}

	@Test
	public void gmMakerTest() {
		int n = delayData.length;
		double[][] rawDelayData = new double[n][];
		for (int i = 0; i < n; i++) {
			int[] in = delayInitData[i];
			double[] dn = delayData[i];
			rawDelayData[i] = new double[] { in[0], in[1], in[2], dn[0], dn[1], dn[2], dn[3] };
		}

		CMRxGMProblemMaker maker = new CMRxGMProblemMaker();
		maker.setModel(new double[] { 1, 1 });
		maker.setShrink(-1);
		maker.setGM(rawDelayData);

		CMRxProblem problem = maker.getProblem();

		ParCMRxSolver solver = new ParCMRxSolver();

		CMRSolution sol = solver.solve(problem);
		assertEquals(1.65296751579206, sol.getFStar(), 1e-5);
	}

	@Test
	public void staMRTest() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		double[] x1Expected = { 0.307307683091782, 0.336186028243926, 0.411231543807238, 0.411231543807375, 0.411231543807457,
				0.524785168266143,
				0.583048612432285, 0.600042950685741 };
		double[] x2Expected = { 0.243969803017848, 0.304387501017471, 0.313993948960525, 0.313993948960525, 0.313993948960525,
				0.315929574797554,
				0.327067491570140, 0.327067491570140 };

		MRSolution[] soln = gm.calcStaMR(null, null);
		double sum = 0;
		for (int i = 0; i < soln.length; i++)
			sum += soln[i].getfVal();

		assertEquals(73.6403692693032, sum, 1e-8);
		assertArrayEquals(x1Expected, soln[0].getxVector(), 1e-10);
		assertArrayEquals(x2Expected, soln[1].getxVector(), 1e-10);
	}

	@Test
	public void staReverseMRTest() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		CombinedStatsSTA[] s = gm.calcStats();

		MRProblem p = new MRProblem(s[0].getMeans().toArray(), s[0].getWeights(), new int[][] { { 1, 2, 3 }, { 7, 6, 5 } });

		MRSolverReverse revSolver = new MRSolverReverse();
		MRSolution sol = revSolver.solve(p);

		assertEquals(12.93488115156335, sol.getfVal(), 1e-8);
	}

	@Test
	public void datafitTest() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		DoubleMatrix2D model = new DenseDoubleMatrix2D(4, 1);
		model.assign(1);

		int nvar = gm.getNDepVar().length;

		@SuppressWarnings("unchecked")
		HashSet<SimpleLinearConstraint>[] adj = new HashSet[nvar];

		HashSet<SimpleLinearConstraint> a = new HashSet<>();
		a.add(new SimpleLinearConstraint(0, 1));
		a.add(new SimpleLinearConstraint(1, 2));
		a.add(new SimpleLinearConstraint(2, 3));

		a.add(new SimpleLinearConstraint(4, 5));
		a.add(new SimpleLinearConstraint(5, 6));
		a.add(new SimpleLinearConstraint(6, 7));

		a.add(new SimpleLinearConstraint(4, 0));

		a.add(new SimpleLinearConstraint(5, 1));

		a.add(new SimpleLinearConstraint(6, 2));

		a.add(new SimpleLinearConstraint(7, 3));

		for (int i = 0; i < nvar; i++) {
			adj[i] = a;
		}

		CMRxGMFits fits = new CMRxGMFits(5, gm, -1, model, adj, -1, false, false, 0, 0, false);
		
		assertEquals(1.5771779458404, fits.getDataFit(), 1e-9);
	}

	@Test
	public void seedRepeatTest() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		DenseDoubleMatrix2D model = new DenseDoubleMatrix2D(4, 1);
		model.assign(1);

		int nvar = gm.getNDepVar().length;

		@SuppressWarnings("unchecked")
		HashSet<SimpleLinearConstraint>[] adj = new HashSet[nvar];

		HashSet<SimpleLinearConstraint> a = new HashSet<>();
		a.add(new SimpleLinearConstraint(0, 1));
		a.add(new SimpleLinearConstraint(1, 2));
		a.add(new SimpleLinearConstraint(2, 3));
		a.add(new SimpleLinearConstraint(4, 5));
		a.add(new SimpleLinearConstraint(5, 6));
		a.add(new SimpleLinearConstraint(6, 7));
		a.add(new SimpleLinearConstraint(4, 0));
		a.add(new SimpleLinearConstraint(5, 1));
		a.add(new SimpleLinearConstraint(6, 2));
		a.add(new SimpleLinearConstraint(7, 3));

		for (int i = 0; i < nvar; i++) {
			adj[i] = a;
		}

		CMRxGMFits fits = new CMRxGMFits(1000, gm, -1, model, adj, -1, false, false, 0, 0, false, false, 242343l, false);
		assertEquals(1.5771779458404, fits.getDataFit(), 1e-9);
		assertEquals(0.177, fits.getP(), 1e-9);
		assertEquals(0.23202811288387803, fits.getFits()[1], 1e-9);
		assertEquals(0.633664696462513, fits.getFits()[10], 1e-9);
	}

	@Test
	public void datafit2Test() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		int nvar = gm.getNDepVar().length;
		DoubleMatrix2D model = new DenseDoubleMatrix2D(nvar, 1);
		model.assign(1);

		CMRxGMFits fits = new CMRxGMFits(100, gm, -1, model, null, -1, false, false, 0, 0, true);

		assertEquals(1.6529675157507, fits.getDataFit(), 1e-9);
	}

	@Test
	public void approximateSpeedTest() {
		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		int nvar = gm.getNDepVar().length;
		DoubleMatrix2D model = new DenseDoubleMatrix2D(nvar, 1);
		model.assign(1);

		long start1 = System.nanoTime();
		new CMRxGMFits(20, gm, -1, model, null, -1, false, false, 0, 0, false);
		long end1 = System.nanoTime();

		long start2 = System.nanoTime();
		new CMRxGMFits(20, gm, -1, model, null, -1, false, false, 0, 0, true);
		long end2 = System.nanoTime();

		assertTrue(end1 - start1 > end2 - start2);
	}

	@Test
	public void strangeModelTest() {
		CMRxProblemMaker maker = new CMRxProblemMaker();

		GeneralModel gm = new GeneralModel();
		for (int i = 0; i < delayData.length; i++)
			gm.addData(delayInitData[i][0], delayInitData[i][1], delayInitData[i][2], delayData[i]);

		CombinedStatsSTA[] stats = gm.calcStats();

		maker.addWeightArray(stats[0].getWeights().toArray());
		maker.addWeightArray(stats[1].getWeights().toArray());

		maker.addMeanArray(stats[0].getMeans().toArray());
		maker.addMeanArray(stats[1].getMeans().toArray());

		maker.setModel(new double[][] { { 1 }, { 1 } });

		int index = maker.initAdj();
		maker.addAdj(1, index, new int[] { 1, 2, 3, 4 });
		maker.addAdj(1, index, new int[] { 5, 6, 7, 8 });
		maker.addAdj(1, index, new int[] { 5, 1 });
		maker.addAdj(1, index, new int[] { 6, 2 });
		maker.addAdj(1, index, new int[] { 7, 3 });
		maker.addAdj(1, index, new int[] { 8, 4 });
		maker.dupeAdj(2);

		CMRxSolver solver = new CMRxSolver();

		CMRSolution sol = solver.solve(maker.getProblem());

		assertEquals(1.74930650492367, sol.getFStar(), 1e-6);
	}
}
