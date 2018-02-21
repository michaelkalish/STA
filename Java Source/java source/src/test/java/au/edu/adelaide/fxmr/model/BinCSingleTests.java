package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.model.bin.BinBaseProblem;
import au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker;
import au.edu.adelaide.fxmr.model.bin.BinCMRxFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRxSolver;
import au.edu.adelaide.fxmr.model.bin.BinSolution;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class BinCSingleTests {

	private int[][][] data = {
			{ { 101, 10 }, { 24, 87 }, { 38, 73 }, { 28, 83 }, { 37, 74 }, { 41, 70 }, { 41, 70 }, { 53, 58 }, { 59, 52 }, { 86, 25 }, { 26, 85 }, { 38, 73 },
					{ 25, 86 }, { 33, 78 }, { 38, 73 }, { 39, 72 }, { 52, 59 }, { 53, 58 }, { 35, 76 }, { 29, 82 }, { 36, 75 }, { 90, 21 }, { 18, 93 },
					{ 37, 74 }, { 20, 91 }, { 28, 83 }, { 36, 75 }, { 40, 71 }, { 50, 61 }, { 57, 54 }
			},
			{ { 89, 22 }, { 52, 59 }, { 58, 53 }, { 29, 82 }, { 79, 32 }, { 51, 60 }, { 77, 34 }, { 107, 4 }, { 78, 33 }, { 87, 24 }, { 54, 57 }, { 64, 47 },
					{ 29, 82 }, { 68, 43 }, { 54, 57 }, { 61, 50 }, { 63, 48 }, { 56, 55 }, { 73, 38 }, { 69, 42 }, { 44, 67 }, { 89, 22 }, { 53, 58 },
					{ 54, 57 }, { 22, 89 }, { 62, 49 }, { 37, 74 }, { 69, 42 }, { 78, 33 }, { 68, 43 }
			},
			{ { 57, 54 }, { 10, 101 }, { 26, 85 }, { 21, 90 }, { 27, 84 }, { 33, 78 }, { 30, 81 }, { 45, 66 }, { 48, 63 }, { 38, 73 }, { 19, 92 }, { 27, 84 },
					{ 20, 91 }, { 31, 80 }, { 34, 77 }, { 32, 79 }, { 38, 73 }, { 43, 68 }, { 18, 93 }, { 22, 89 }, { 36, 75 }, { 49, 62 }, { 12, 99 },
					{ 27, 84 }, { 18, 93 }, { 23, 88 }, { 31, 80 }, { 30, 81 }, { 37, 74 }, { 46, 65 }
			},
			{ { 90, 21 }, { 68, 43 }, { 79, 32 }, { 52, 59 }, { 104, 7 }, { 67, 44 }, { 77, 34 }, { 109, 2 }, { 71, 40 }, { 91, 20 }, { 76, 35 }, { 84, 27 },
					{ 61, 50 }, { 82, 29 }, { 62, 49 }, { 72, 39 }, { 74, 37 }, { 82, 29 }, { 75, 36 }, { 84, 27 }, { 63, 48 }, { 88, 23 }, { 68, 43 },
					{ 75, 36 }, { 17, 94 }, { 77, 34 }, { 66, 45 }, { 84, 27 }, { 84, 27 }, { 66, 45 }
			}
	};

	double[][] m22 = { { 1, 0 }, { 0, 1 }, { -0.8170432915914517, 0.7011001349274827 }, { 2.28878612327359, 1.792193175830165 }
	};

	int[][][] evilData = { { { 101, 10 }, { 17, 94 }, { 39, 72 }, { 25, 86 }, { 46, 65 }, { 43, 68 }, { 44, 67 }, { 62, 49 }, { 47, 64 }, { 80, 31 },
			{ 19, 92 }, { 39, 72 }, { 20, 91 }, { 26, 85 }, { 34, 77 }, { 40, 71 }, { 62, 49 }, { 52, 59 }, { 28, 83 }, { 29, 82 }, { 33, 78 }, { 82, 29 },
			{ 25, 86 }, { 41, 70 }, { 23, 88 }, { 31, 80 }, { 37, 74 }, { 36, 75 }, { 55, 56 }, { 52, 59 }
	}, { { 97, 14 }, { 47, 64 }, { 56, 55 }, { 32, 79 }, { 80, 31 }, { 49, 62 }, { 69, 42 }, { 108, 3 }, { 78, 33 }, { 87, 24 }, { 58, 53 }, { 67, 44 },
			{ 33, 78 }, { 62, 49 }, { 54, 57 }, { 61, 50 }, { 63, 48 }, { 56, 55 }, { 69, 42 }, { 63, 48 }, { 40, 71 }, { 92, 19 }, { 50, 61 }, { 52, 59 },
			{ 31, 80 }, { 60, 51 }, { 36, 75 }, { 61, 50 }, { 77, 34 }, { 69, 42 }
	}, { { 59, 52 }, { 10, 101 }, { 26, 85 }, { 18, 93 }, { 21, 90 }, { 32, 79 }, { 19, 92 }, { 41, 70 }, { 60, 51 }, { 40, 71 }, { 23, 88 }, { 23, 88 },
			{ 19, 92 }, { 29, 82 }, { 31, 80 }, { 27, 84 }, { 30, 81 }, { 41, 70 }, { 19, 92 }, { 19, 92 }, { 33, 78 }, { 49, 62 }, { 17, 94 }, { 21, 90 },
			{ 16, 95 }, { 22, 89 }, { 31, 80 }, { 25, 86 }, { 34, 77 }, { 43, 68 }
	}, { { 91, 20 }, { 62, 49 }, { 78, 33 }, { 50, 61 }, { 107, 4 }, { 77, 34 }, { 85, 26 }, { 111, 0 }, { 67, 44 }, { 94, 17 }, { 68, 43 }, { 80, 31 },
			{ 59, 52 }, { 91, 20 }, { 59, 52 }, { 73, 38 }, { 78, 33 }, { 94, 17 }, { 78, 33 }, { 82, 29 }, { 61, 50 }, { 91, 20 }, { 70, 41 }, { 81, 30 },
			{ 14, 97 }, { 81, 30 }, { 62, 49 }, { 78, 33 }, { 85, 26 }, { 68, 43 }

			}
	};

	@Test
	public void fitsxFINDEvilDataTest() {
		BinModel binModel = makeBinModel(data);
		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(m22);
		BinBaseProblem p = new BinBaseProblem(binModel, null, cmrModel);
		int nSample = 5000;
		new BinCMRxFits(nSample, p, -1, true, true);
	}

	@Test
	public void evilDataTest1() {
		BinCMRProblemMaker maker = new BinCMRProblemMaker(1, 4);
		maker.setModel(new double[][] { { 1 }, { 1 }, { 1 }, { 1 } });
		for (int i = 0; i < 4; i++) {
			maker.setElement(0, i, evilData[i]);
		}

		BinCMRxSolver solver = new BinCMRxSolver();
		solver.setOnlyFeas(true);
		BinSolution[] sol = solver.solve(maker.getProblems());
		assertTrue(sol[0] != null);
	}

	// @Test
	public void evilDataTest() {
		BinModel binModel = makeBinModel(evilData);
		BinBaseProblem p = new BinBaseProblem(binModel, null);

		BinCMRxSolver solver = new BinCMRxSolver();
		solver.setOnlyFeas(true);
		BinSolution[] solnR1 = solver.solve(p);
	}

	private BinModel makeBinModel(int[][][] data) {
		BinModel model = new BinModel(1, data.length);
		int i = 0;
		for (int[][] d : data)
			model.setElement(new BinElement(d), 0, i++);
		return model;
	}

}
