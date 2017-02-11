package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class RealDataTests {
	private DoubleMatrix2D model1 = new DenseDoubleMatrix2D(new double[][] { { 1 }, { 0 }, { 1 }, { 0 } });
	private DoubleMatrix2D model3 = new DenseDoubleMatrix2D(new double[][] { { 1, 0 }, { 0, 0 }, { 0, 1 }, { 0, 0 } });

	private DoubleMatrix2D model6 = new DenseDoubleMatrix2D(
			new double[][] { { 0, 1 }, { 0, -1 }, { 1, 1 }, { 1, -1 } });
	private DoubleMatrix2D model4 = new DenseDoubleMatrix2D(
			new double[][] { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } });
	private DoubleMatrix2D model9 = new DenseDoubleMatrix2D(new double[][] { { 1, 1 }, { 1, 0 }, { 1, 1 }, { 1, 0 } });
	private DoubleMatrix2D model8 = new DenseDoubleMatrix2D(
			new double[][] { { 0, 1, 0 }, { 0, -1, 0 }, { 1, 0, 1 }, { 1, 0, -1 } });

	private double[][] getMeansSK1() {
		return new double[][] {
				{ 99.2, 99.8, 99.4, 99.6, 92.7, 99.2, 97.3, 91.4, 76.9, 79.2, 70.5, 53.2, 74.2, 57.6, 83.6, 66.4 },
				{ 58.7, 62.7, 52, 49.8, 74.1, 64.6, 64.7, 65.3, 77, 69.4, 79.5, 83.1, 99.6, 84.3, 89.6, 75.1 },
				{ 94.2, 99.3, 95.5, 92.9, 61.7, 49.8, 86.4, 55.6, 93.9, 99.4, 93.8, 83.1, 67.3, 80.2, 82.3, 63.7 },
				{ 65.6, 48.6, 71.9, 62.4, 88.9, 99, 78, 89.8, 80.4, 87.7, 80.6, 74.9, 95.7, 96.4, 82.7, 76.6 } };
	}

	private DoubleMatrix2D[] getWeightsSK1() {
		DoubleMatrix2D[] weights = new DoubleMatrix2D[4];
		weights[0] = DoubleFactory2D.dense.diagonal(
				new DenseDoubleMatrix1D(new double[] { 250, 250, 160, 160, 0.0990074503106359, 250, 0.657462195923734,
						0.0881659282770173, 0.0237953599048186, 0.0334124093688396, 0.0243865264441396,
						0.015439186972414, 0.0221453287197232, 0.016, 0.0362897372623022, 0.0196655866981972 }));
		weights[1] = DoubleFactory2D.dense.diagonal(new DenseDoubleMatrix1D(new double[] { 0.0167965600644988,
				0.017146114904689, 0.0239118608807934, 0.0160641925132831, 0.0221453287197232, 0.0180309320639557,
				0.0179546107440391, 0.0198411714227608, 0.0240292195309496, 0.0174336757597814, 0.0330294622803541,
				0.0322830578512397, 160, 0.0354308390022676, 0.0403124212648022, 0.025 }));
		weights[2] = DoubleFactory2D.dense.diagonal(new DenseDoubleMatrix1D(new double[] { 0.425124880433627, 62.5,
				2.26757369614512, 0.324648973297622, 0.0424407686023194, 0.0277008310249307, 0.182615047479912,
				0.0418931515170557, 1.07497984412792, 160, 0.324648973297622, 0.0485619589894256, 0.0356426821118289,
				0.0811622433244055, 0.177777777777778, 0.0776261910768693 }));
		weights[3] = DoubleFactory2D.dense.diagonal(new DenseDoubleMatrix1D(new double[] { 0.0469131169074873,
				0.027411529289219, 0.0532793435984869, 0.0660982219578293, 0.730460189919649, 81.6326530612245,
				0.13840830449827, 0.170874449997864, 0.114387028510967, 0.115620302925194, 0.0528925619834711,
				0.0655640971004278, 2.16333153055706, 1.18906064209275, 0.145158949049209, 0.0660982219578293 }));
		return weights;
	}

	@Test
	public void model9Test() {
		CMRxProblem problem = new CMRxProblem(getMeansSK1(), getWeightsSK1(), null, model9);
		// System.out.println(problem.toMatlab());
		ParCMRxSolver solver = new ParCMRxSolver();

		CMRSolution sol = solver.solve(problem, false);
		// With Fast implementation
		assertEquals(136.3192329017816, sol.getFStar(), 1e-5);
		// With JOptimizer original
		// assertEquals(136.31923290183096, sol.getFStar(), 1e-15);
	}

	@Test
	public void model6Test() {
		CMRxProblem problem = new CMRxProblem(getMeansSK1(), getWeightsSK1(), null, model6);
		// System.out.println(problem.toMatlab());
		ParCMRxSolver solver = new ParCMRxSolver();
		CMRSolution sol = solver.solve(problem, false);
		assertEquals(102.03152394, sol.getFStar(), 1e-5);
	}

	@Test
	public void tolerance6Test() {
		double tol = 0.5;
		
		CMRxProblem problem = new CMRxProblem(getMeansSK1(), getWeightsSK1(), null, model6);
		// System.out.println(problem.toMatlab());
		ParCMRxSolver solver = new ParCMRxSolver();
		solver.setTolerance(tol);
		CMRSolution sol = solver.solve(problem, false);
		// assertEquals(107.6, sol.getFStar(), 1e-1); // for 0.2

		CMRxSolver s2 = new CMRxSolver();
		s2.setTolerance(tol);
		CMRSolution sol2 = s2.solve(problem, false);
		assertEquals(sol.getFStar(), sol2.getFStar(), 1e-15);
	}

	@Test
	public void model1Test() {
		CMRxProblem problem = new CMRxProblem(getMeansSK1(), getWeightsSK1(), null, model1);
		// System.out.println(problem.toMatlab());
		ParCMRxSolver solver = new ParCMRxSolver();

		CMRSolution sol = solver.solve(problem, false);
		assertEquals(1000.90753441429, sol.getFStar(), 1e-5);
		assertEquals(4, sol.getXStar().length);
		assertEquals(4, sol.getAdjs().length);
	}
}
