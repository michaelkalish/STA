package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.model.bin.BinBaseProblem;
import au.edu.adelaide.fxmr.model.bin.BinCMRxFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRxSolver;
import au.edu.adelaide.fxmr.model.bin.BinSolution;
import cern.colt.Arrays;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class BinSingleTests {

	private int[][][] data = { { { 88, 33 }, { 20, 101 }, { 32, 89 }, { 25, 96 }, { 34, 87 }, { 34, 87 }, { 36, 85 }, { 44, 77 }, { 50, 71 }, { 75, 46 },
			{ 22, 99 }, { 33, 88 }, { 22, 99 }, { 30, 91 }, { 31, 90 }, { 32, 89 }, { 43, 78 }, { 45, 76 }, { 33, 88 }, { 26, 95 }, { 29, 92 }, { 78, 43 },
			{ 16, 105 }, { 31, 90 }, { 16, 105 }, { 27, 94 }, { 32, 89 }, { 33, 88 }, { 43, 78 }, { 49, 72 }
	}, { { 75, 46 }, { 45, 76 }, { 51, 70 }, { 26, 95 }, { 70, 51 }, { 43, 78 }, { 65, 56 }, { 93, 28 }, { 68, 53 }, { 74, 47 }, { 46, 75 }, { 55, 66 },
			{ 27, 94 }, { 62, 59 }, { 44, 77 }, { 51, 70 }, { 51, 70 }, { 45, 76 }, { 65, 56 }, { 59, 62 }, { 38, 83 }, { 76, 45 }, { 48, 73 }, { 48, 73 },
			{ 19, 102 }, { 53, 68 }, { 33, 88 }, { 58, 63 }, { 68, 53 }, { 58, 63 }
	}, { { 51, 70 }, { 8, 113 }, { 24, 97 }, { 20, 101 }, { 25, 96 }, { 29, 92 }, { 27, 94 }, { 38, 83 }, { 42, 79 }, { 32, 89 }, { 15, 106 }, { 26, 95 },
			{ 18, 103 }, { 29, 92 }, { 28, 93 }, { 27, 94 }, { 30, 91 }, { 35, 86 }, { 15, 106 }, { 19, 102 }, { 31, 90 }, { 39, 82 }, { 10, 111 }, { 24, 97 },
			{ 18, 103 }, { 21, 100 }, { 27, 94 }, { 25, 96 }, { 32, 89 }, { 36, 85 }
	}, { { 76, 45 }, { 59, 62 }, { 67, 54 }, { 48, 73 }, { 90, 31 }, { 58, 63 }, { 66, 55 }, { 94, 27 }, { 60, 61 }, { 78, 43 }, { 66, 55 }, { 72, 49 },
			{ 52, 69 }, { 76, 45 }, { 54, 67 }, { 63, 58 }, { 63, 58 }, { 69, 52 }, { 62, 59 }, { 70, 51 }, { 54, 67 }, { 75, 46 }, { 57, 64 }, { 69, 52 },
			{ 14, 107 }, { 71, 50 }, { 56, 65 }, { 74, 47 }, { 74, 47 }, { 57, 64 }
	} };
	private int[][] hitsExample = { { 88, 18, 31, 23, 32, 35, 40, 36, 53, 75, 25, 33, 19, 33, 36, 33, 42, 37, 38, 23, 27, 78, 20, 28, 13, 19, 31, 42, 44, 55 },
			{ 80, 50, 53, 29, 67, 42, 69, 95, 66, 79, 48, 55, 26, 59, 45, 52, 49, 45, 70, 64, 40, 72, 54, 52, 15, 52, 41, 67, 65, 57 },
			{ 49, 10, 29, 20, 34, 22, 31, 38, 36, 26, 14, 27, 17, 34, 23, 31, 33, 34, 9, 14, 40, 45, 9, 18, 13, 30, 34, 25, 28, 41 },
			{ 72, 64, 73, 46, 89, 53, 59, 97, 52, 69, 70, 77, 46, 77, 56, 63, 61, 84, 66, 69, 46, 71, 56, 68, 15, 76, 49, 75, 72, 53 }
	};
	private int[][] missesExample = {
			{ 33, 103, 90, 98, 89, 86, 81, 85, 68, 46, 96, 88, 102, 88, 85, 88, 79, 84, 83, 98, 94, 43, 101, 93, 108, 102, 90, 79, 77, 66 },
			{ 41, 71, 68, 92, 54, 79, 52, 26, 55, 42, 73, 66, 95, 62, 76, 69, 72, 76, 51, 57, 81, 49, 67, 69, 106, 69, 80, 54, 56, 64 },
			{ 72, 111, 92, 101, 87, 99, 90, 83, 85, 95, 107, 94, 104, 87, 98, 90, 88, 87, 112, 107, 81, 76, 112, 103, 108, 91, 87, 96, 93, 80 },
			{ 49, 57, 48, 75, 32, 68, 62, 24, 69, 52, 51, 44, 75, 44, 65, 58, 60, 37, 55, 52, 75, 50, 65, 53, 106, 45, 72, 46, 49, 68 }
	};

	@Test
	public void singleBinApproxTest() {
		BinModel binModel = makeSingleBinModel();

		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { 1 }, { 1 }, { 1 } });
		BinBaseProblem p = new BinBaseProblem(binModel, null, cmrModel);

		BinCMRxSolver solver = new BinCMRxSolver();
		solver.setOnlyFeas(true);
		BinSolution[] sol = solver.solve(p);

		assertEquals(198.35249070689395, sol[0].getG2Star(), 1e-16);
	}

	@Test
	public void fitsxApproxTest() {
		BinModel binModel = makeBinModel();

		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { 1 }, { 1 }, { 1 } });
		BinBaseProblem p = new BinBaseProblem(binModel, null, cmrModel);

		int nSample = 10;

		BinCMRxFits fits = new BinCMRxFits(nSample, p, -1, true, true);
		// Assuming we dont cancel...
		assertTrue(fits.getFits().length == nSample);
		assertTrue(fits.getXStars().length == nSample);
	}

	/**
	 * Bug that john found where the cmrModel wasn't being passed through!
	 */
	@Test
	public void fitsxCMRModel() {
		BinModel binModel = makeBinModel();
		int nSample = 1;

		DoubleMatrix2D cmrModel2 = new DenseDoubleMatrix2D(new double[][] { { 1, 0 }, { 0, 1 }, { 1, 0 }, { 0, 1 } });
		BinBaseProblem p2 = new BinBaseProblem(binModel, null, cmrModel2);
		BinCMRxFits fits2 = new BinCMRxFits(nSample, p2, -1, true, false);

		DoubleMatrix2D cmrModel1 = new DenseDoubleMatrix2D(new double[][] { { 1 }, { 1 }, { 1 }, { 1 } });
		BinBaseProblem p1 = new BinBaseProblem(binModel, null, cmrModel1);
		BinCMRxFits fits1 = new BinCMRxFits(nSample, p1, -1, true, false);

		assertTrue(Math.abs(fits1.getBaseFitDiff()[0] - fits2.getBaseFitDiff()[0]) > 1);
	}

	private BinModel makeBinModel() {
		BinModel model = new BinModel(1, data.length);
		int i = 0;
		for (int[][] d : data) {

			model.setElement(new BinElement(d), 0, i++);
		}
		return model;
	}

	private BinModel makeSingleBinModel() {
		BinModel model = new BinModel(1, 4);

		for (int i = 0; i < 4; i++)
			model.setElement(new BinElement(hitsExample[i], missesExample[i]), 0, i);
		return model;
	}
}
