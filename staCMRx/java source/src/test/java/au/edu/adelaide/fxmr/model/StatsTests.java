package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Tests various stats methods
 * 
 */
public class StatsTests {

	private static final double TOL = 1e-15;
	private double[][] zeroCov = { { 0.017197772491349, 0.022730860726644, 0.017445663927336, 0.018836505190311 },
			{ 0.022730860726644, 0.045423875432526, 0.038870296280277, 0.040642571366782 },
			{ 0.017445663927336, 0.038870296280277, 0.049971885813149, 0.048984915657439 },
			{ 0.018836505190311, 0.040642571366782, 0.048984915657439, 0.063271518166090 } };
	private double[][] optimalCov = { { 0.017197772491349, 0.021549498513419, 0.016538982548385, 0.017857539381289 },
			{ 0.021549498513419, 0.045423875432526, 0.036850139639726, 0.038530306519515 },
			{ 0.016538982548385, 0.036850139639726, 0.049971885813149, 0.046439084724259 },
			{ 0.017857539381289, 0.038530306519515, 0.046439084724259, 0.063271518166090 } };

	/**
	 * Tests the ShrinkDiag class's with complete shrinkage
	 */
	@Test
	public void shrinkFullTest() {
		double[] expected = { 0.017197772491349, 0.045423875432526, 0.049971885813149, 0.063271518166090 };

		ShrinkDiagonal sd = new ShrinkDiagonal(makeTestMatrix());
		sd.shrink(1);

		assertEquals(1, sd.getShrinkage(), TOL);

		double[][] actuals = sd.getResult().toArray();
		for (int i = 0; i < expected.length; i++)
			assertEquals(expected[i], actuals[i][i], TOL);
	}

	/**
	 * Tests the ShrinkDiag class's with complete shrinkage
	 */
	@Test
	public void shrinkZeroTest() {
		double[][] expected = zeroCov;

		ShrinkDiagonal sd = new ShrinkDiagonal(makeTestMatrix());
		sd.shrink(0);

		assertEquals(0, sd.getShrinkage(), TOL);

		double[][] actuals = sd.getResult().toArray();
		for (int i = 0; i < expected.length; i++)
			assertArrayEquals(expected[i], actuals[i], TOL);
	}

	/**
	 * Tests the ShrinkDiag class's optimal shrinkage method
	 */
	@Test
	public void shrinkAutoTest() {
		double[][] expected = optimalCov;

		ShrinkDiagonal sd = new ShrinkDiagonal(makeTestMatrix());
		sd.shrink();

		assertEquals(0.051971732501974, sd.getShrinkage(), TOL);

		double[][] actuals = sd.getResult().toArray();
		for (int i = 0; i < expected.length; i++)
			assertArrayEquals(expected[i], actuals[i], TOL);
	}

	@Test
	public void statSTAMiniTest() {
		DoubleMatrix2D d = makeTestMatrix().viewPart(0, 0, 1, 4);
		StatsSTA s = new StatsSTA(d);
		assertArrayEquals(new double[] { 0.525000000000000, 0.600000000000000, 0.900000000000000, 0.925000000000000 },
				s.getMeans().toArray(), 1e-15);
		assertEquals(-1, s.getShrinkage(), 1e-15);
		assertEquals(1, s.getN());

		double[][] expectedCov = new double[4][4];

		double[][] actuals = s.getCovariance().toArray();
		for (int i = 0; i < expectedCov.length; i++)
			assertArrayEquals(expectedCov[i], actuals[i], TOL);

		assertTrue(Double.isNaN(s.getLMValue()));
	}

	@Test
	public void statSTATest() {
		DoubleMatrix2D d = makeTestMatrix();
		StatsSTA s = new StatsSTA(d);
		assertArrayEquals(new double[] { 0.367647058823529, 0.467647058823529, 0.575735294117647, 0.611764705882353 },
				s.getMeans().toArray(), 1e-15);
		assertEquals(0.051971732501974, s.getShrinkage(), 1e-15);
		assertEquals(34, s.getN());

		double[][] expectedCov = zeroCov;
		double[][] actuals = s.getCovariance().toArray();
		for (int i = 0; i < expectedCov.length; i++)
			assertArrayEquals(expectedCov[i], actuals[i], TOL);

		expectedCov = optimalCov;
		actuals = s.getRegCovariance().toArray();
		for (int i = 0; i < expectedCov.length; i++)
			assertArrayEquals(expectedCov[i], actuals[i], TOL);

		double[][] expectedW = { { 4909.21477287896, -2498.46049296539, 287.220544418912, -74.8884624739565 },
				{ -2498.46049296539, 3232.67047265783, -1204.0883450443, -379.675841004848 },
				{ 287.220544418912, -1204.0883450443, 2727.57207512732, -1349.75597766002 },
				{ -74.8884624739565, -379.675841004848, -1349.75597766002, 1780.38690501671 } };
		actuals = s.getWeights().toArray();
		for (int i = 0; i < expectedW.length; i++)
			for (int j = 0; j < expectedW[i].length; j++) {
				double tol = 1e-10;
				assertEquals(expectedW[i][j], actuals[i][j], tol);
			}

		assertEquals(0.0130997474747475, s.getLMValue(), TOL);
	}

	private DoubleMatrix2D makeTestMatrix() {
		double[][] data = { { 0.525, 0.6, 0.9, 0.925 }, { 0.5125, 0.85, 0.9625, 0.8875 }, { 0.3125, 0.2, 0.3375, 0.25 },
				{ 0.3, 0.3125, 0.5, 0.3625 }, { 0.3, 0.4125, 0.3875, 0.575 }, { 0.375, 0.7, 1, 1 },
				{ 0.175, 0.35, 0.5375, 0.35 }, { 0.5625, 0.4875, 0.4, 0.5375 }, { 0.325, 0.475, 0.525, 0.55 },
				{ 0.35, 0.5625, 0.5125, 0.9 }, { 0.275, 0.3, 0.3125, 0.2875 }, { 0.3125, 0.525, 0.55, 0.6125 },
				{ 0.225, 0.175, 0.2625, 0.275 }, { 0.25, 0.3375, 0.8125, 0.9625 }, { 0.35, 0.425, 0.8, 0.6625 },
				{ 0.7, 0.85, 0.8625, 0.925 }, { 0.3625, 0.5875, 0.7875, 0.9375 }, { 0.5, 0.8625, 0.8875, 0.925 },
				{ 0.225, 0.375, 0.375, 0.3875 }, { 0.3375, 0.4, 0.4625, 0.4 }, { 0.25, 0.4125, 0.5375, 0.5 },
				{ 0.2875, 0.1625, 0.3875, 0.425 }, { 0.35, 0.3125, 0.4125, 0.4 }, { 0.3375, 0.25, 0.4125, 0.4 },
				{ 0.5125, 0.5625, 0.4625, 0.5125 }, { 0.5125, 0.9, 0.95, 0.95 }, { 0.7375, 0.975, 0.9875, 0.975 },
				{ 0.375, 0.425, 0.65, 0.775 }, { 0.3125, 0.375, 0.375, 0.525 }, { 0.275, 0.2875, 0.4625, 0.7125 },
				{ 0.3, 0.3375, 0.5875, 0.8375 }, { 0.3, 0.3125, 0.4125, 0.2625 }, { 0.4625, 0.5125, 0.4125, 0.5375 },
				{ 0.2125, 0.2875, 0.35, 0.275 } };

		return new DenseDoubleMatrix2D(data);
	}

}
