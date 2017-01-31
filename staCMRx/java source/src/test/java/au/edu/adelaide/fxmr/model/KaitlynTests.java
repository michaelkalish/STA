package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.*;

import java.util.DoubleSummaryStatistics;

import org.junit.Test;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.model.mr.MRUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Statistic;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class KaitlynTests {

	double[][][][] data = { {
			{ { 1, 0.25 }, { 1, 1 }, { 0.75, 1 }, { 0.75, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 },
					{ 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.75 }, { 1, 1 }, { 1, 0.5 },
					{ 1, 0.75 }, { 0.75, 0 }, { 0.5, 0.25 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.25 }
			},
			{ { 1, 0 }, { 1, 1 }, { 1, 1 }, { 0.25, 0.5 }, { 0, 0.5 }, { 0.75, 1 }, { 0, 0 }, { 1, 1 }, { 0, 0 }, { 1, 1 }, { 1, 1 },
					{ 0.75, 1 }, { 0.5, 0.5 }, { 0, 0.5 }, { 1, 0.75 }, { 0, 0 }, { 1, 1 }, { 0, 0 }, { 0, 0 }, { 0, 0.75 }, { 0.25, 0 },
					{ 1, 0 }, { 0, 0.25 }, { 0, 0 }, { 0, 0 }, { 0.5, 0.25 }, { 0, 0 }, { 0, 0 }, { 0.75, 0.25 }
			},
			{ { 1, 1 }, { 1, 0 }, { 0.75, 0.5 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.25 }, { 1, 1 }, { 1, 0 },
					{ 1, 1 }, { 1, 1 }, { 1, 0.25 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.75 }, { 1, 1 }, { 1, 1 },
					{ 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.75 }
			},
			{ { 1, 1 }, { 1, 0.5 }, { 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 0.75, 1 }, { 1, 0.75 }, { 0.75, 1 }, { 1, 0 },
					{ 1, 0.75 }, { 0, 0 }, { 1, 0.25 }, { 0, 0 }, { 0, 0 }, { 0.5, 0.75 }, { 0, 0 }, { 0.75, 0.5 }, { 0, 0 }, { 0, 0 },
					{ 0, 0 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 0.75, 1 }, { 0, 0 }
			} }, {
					{ { 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.25 }, { 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 },
							{ 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 0.75, 0.25 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 },
							{ 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0 }, { 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 },
							{ 1, 1 }
					},
					{ { 0.75, 0.5 }, { 0.75, 0.75 }, { 1, 1 }, { 0.5, 1 }, { 0.75, 0.5 }, { 0.5, 0.75 }, { 1, 1 }, { 0, 0 }, { 0, 0 },
							{ 0.75, 0.75 }, { 0, 0.25 }, { 0.75, 0.5 }, { 0, 0 }, { 1, 1 }, { 0.75, 0.5 }, { 0, 0.25 }, { 1, 1 }, { 0, 0 },
							{ 0.75, 0.75 }, { 0.75, 1 }, { 0, 0 }, { 0, 0.25 }, { 0, 0 }, { 0.5, 0.75 }, { 0, 0 }, { 0, 0 }, { 0, 0 },
							{ 1, 1 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0.5, 1 }
					},
					{ { 0.75, 1 }, { 0.75, 0.5 }, { 0.75, 1 }, { 1, 0 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0.25 }, { 1, 1 }, { 1, 1 },
							{ 1, 0.75 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 }, { 1, 1 },
							{ 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 0 }, { 1, 1 }, { 1, 1 }, { 1, 1 }, { 1, 1 }
					},
					{ { 1, 0.25 }, { 0.75, 0.75 }, { 1, 1 }, { 0.75, 0 }, { 1, 1 }, { 1, 0.75 }, { 1, 1 }, { 0.5, 0 }, { 1, 1 },
							{ 1, 0.75 }, { 0, 0 }, { 1, 0.75 }, { 0, 0 }, { 0, 0 }, { 1, 1 }, { 1, 0.75 }, { 1, 0 }, { 1, 0 }, { 1, 0.75 },
							{ 0.25, 0 }, { 0.5, 0.75 }, { 0, 0 }, { 1, 0.75 }, { 0, 0 }, { 1, 0 }, { 0.25, 0 }, { 0, 0 }, { 1, 1 }, { 1, 1 }
					}
			}
	};

	double[][][] allModels = { { { 1, 0, 1 }, { 1, 0, 0 }, { 0, 1, 1 }, { 0, 1, 0 }
	},
			{ { 1, 0, 1 }, { 1, 0, -1 }, { 0, 1, 1 }, { 0, 1, -1 }
			},
			{ { 1, 0, 1, 0 }, { 1, 0, 0, 0 }, { 0, 1, 0, 1 }, { 0, 1, 0, 0 }
			},
			{ { 1, 0, 1, 0 }, { 1, 0, -1, 0 }, { 0, 1, 0, 1 }, { 0, 1, 0, -1 }
			},
			{ { 1, 1 }, { 1, 0 }, { 1, 1 }, { 1, 0 }
			},
			{ { 1, 1 }, { 1, -1 }, { 1, 1 }, { 1, -1 }
			},
			{ { 1, 1, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 0, 0 }
			},
			{ { 1, 1, 0 }, { 1, -1, 0 }, { 1, 0, 1 }, { 1, 0, -1 }
			},
			{ { 0, 1 }, { 0, 0 }, { 1, 1 }, { 1, 0 }
			},
			{ { 0, 1 }, { 0, -1 }, { 1, 1 }, { 1, -1 }
			},
			{ { 0, 1, 0 }, { 0, 0, 0 }, { 1, 0, 1 }, { 1, 0, 0 }
			},
			{ { 0, 1, 0 }, { 0, -1, 0 }, { 1, 0, 1 }, { 1, 0, -1 }
			},
			{ { 1 }, { 0 }, { 1 }, { 0 }
			},
			{ { 1 }, { -1 }, { 1 }, { -1 }
			},
			{ { 1, 0 }, { 0, 0 }, { 0, 1 }, { 0, 0 }
			},
			{ { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }
			}
	};

	@Test
	public void cmrxTest() {
		// CMRxGMProblemMaker maker = new CMRxGMProblemMaker();
		CMRxSolver solver = new CMRxSolver();
		CMRxGMProblemMaker maker = new CMRxGMProblemMaker();
		for (int var = 0; var < 4; var++)
			for (int cond = 0; cond < 2; cond++)
				maker.addCell(cond, var, data[cond][var]);

		maker.setModel(new double[] { 1, 1, 1, 1 });
		maker.setShrink(-1);

		CMRxProblem problem = maker.getProblem();

		DoubleMatrix2D expectedW1 = new DenseDoubleMatrix2D(
				new double[][] { { 2295.093787019131, -192.6346977006495, 0, 0 }, { -192.6346977006495, 360.8893057326901, 0, 0 },
						{ 0, 0, 16941.05779131645, -121.2245706372899 }, { 0, 0, -121.2245706372899, 514.8752858145195 }
				});
		assertEquals(expectedW1, problem.getWeights()[0]);
		CMRSolution sol = solver.solve(problem);
		assertEquals(3.22374073887459, sol.getFStar(), 1e-5);
	}

	@Test
	public void fitTest() {
		CMRxFitsGMProblemMaker maker = new CMRxFitsGMProblemMaker();
		for (int var = 0; var < 4; var++)
			for (int cond = 0; cond < 2; cond++)
				maker.addCell(cond, var, data[cond][var]);

		maker.setShrink(-1);
		for (double[][] m : allModels) {
			long start = System.nanoTime();
			maker.setModel(m);
			Fits f = maker.solve(2, -1, false);
			// System.out.print(f.getP() + "\t" + (System.nanoTime() - start) /
			// 1000000000.0);

			start = System.nanoTime();
			f = maker.solve(2, -1, true);
			// System.out.println("\t" + f.getP() + "\t" + (System.nanoTime() -
			// start) / 1000000000.0);
		}
	}

	public void zeroModelTest() {
		// TODO: a model of all zeroes?
		CMRxFitsGMProblemMaker maker = new CMRxFitsGMProblemMaker();
		for (int var = 0; var < 4; var++)
			for (int cond = 0; cond < 2; cond++)
				maker.addCell(cond, var, data[cond][var]);

		maker.setModel(new double[] { 1, 1, 1, 1 });
		maker.setShrink(-1);

		// maker.solve(10000000, -1);
	}

	@Test
	public void evilSampleTest() {
		DenseDoubleMatrix2D m = new DenseDoubleMatrix2D(27, 2);
		for (int i = 0; i < 27; i++) {
			m.set(i, 0, 0.9970173329376881);
			m.set(i, 1, 1.0273944040692427);
		}

		DoubleMatrix2D c1 = MRUtil.covariance(m);
		double[][] expectedCov = { { 1.10933564796705e-31, 1.4791141972894e-31 },
				{ 1.4791141972894e-31, 1.97215226305253e-31 } };
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				assertEquals(expectedCov[i][j], c1.get(i, j), 1e-45);

		try {
			new StatsSTA(m);
			assertTrue(true);
		} catch (Exception e) {
			assertTrue(false);
		}
	}
}
