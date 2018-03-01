package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.model.bin.BinBaseProblem;
import au.edu.adelaide.fxmr.model.bin.BinCMRFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRSolver;
import au.edu.adelaide.fxmr.model.bin.BinCMRxFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRxSolver;
import au.edu.adelaide.fxmr.model.bin.BinProblem;
import au.edu.adelaide.fxmr.model.bin.BinSolution;
import au.edu.adelaide.fxmr.model.bin.ParBinCMRxSolver;
import cern.colt.Arrays;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class BinTests {
	private double[][][] xExpected = {
			{ { 0.717948717948718, 0.6923076923076923, 0.8717948717948718, 0.4487179487179487, 0.6410256410256411, 0.717948717948718 },
					{ 0.6794871794871795, 0.6730769230769569, 0.7692307692307693, 0.6153846153846154, 0.6730769230768893,
							0.7564102564102564 },
			},
			{ { 0.7307692307692307, 0.6025641025641025, 0.7307692307692307, 0.5897435897435898, 0.6410256410256411,
					0.6538461538461539 },
					{ 0.7008547070830424, 0.6538461538461539, 0.7307692307692307, 0.6025641025641025, 0.7008546884038074,
							0.7008547070772528 },
			},
			{ { 0.6923076923076923, 0.6923076923076923, 0.9102564102564102, 0.6025641025641025, 0.7051282051282053,
					0.6923076923076923 },
					{ 0.6025641025641025, 0.6666666666666666, 0.8076923076923077, 0.5769230769230769, 0.8076923076923077,
							0.7948717948717948 },
			},
			{ { 0.5854700854700782, 0.7179487179441275, 0.8076923076969649, 0.5854700854700806, 0.5854700854700778,
					0.6153846153845684 },
					{ 0.6025641025641025, 0.6698717996648093, 0.6698717996779072, 0.669871780477728, 0.5512820512794809,
							0.6698717996693057 },
			},
			{ { 0.6410256410256411, 0.6410256410256411, 0.6410256410256411, 0.5641025639645181, 0.56410256424061, 0.6666666666666666 },
					{ 0.7307692307692307, 0.7115384615600647, 0.7564102564102564, 0.5512820512768887, 0.7115384615220209,
							0.8205128205128205 },
			},
			{ { 0.6474358974120321, 0.7820512820512907, 0.7820512820512821, 0.6282051282051282, 0.6923076923076871,
					0.6474358974597595 },
					{ 0.5897435897435898, 0.6752136756468904, 0.7564102564102564, 0.5512820512820513, 0.675213674347686,
							0.6752136756464493 },
			},
			{ { 0.576923076923077, 0.6794871794871795, 0.717948717948718, 0.5256410256410257, 0.576923076923077, 0.6153846153846154 },
					{ 0.666666671305421, 0.7435897436571827, 0.8333333333333333, 0.666666660226635, 0.7307692325705082,
							0.7435897435223046 },
			},
			{ { 0.5384615379277006, 0.5384615379277006, 0.8205128207595457, 0.6538461525394378, 0.6538461562205457,
					0.6923076920609672 },
					{ 0.6623931623970206, 0.6623931623854452, 0.8076923076923133, 0.6623931623970214, 0.7307692307692307,
							0.8076923076923021 },
			},
			{ { 0.5961538461538462, 0.6346153846153845, 0.6346153846153847, 0.5961538461538462, 0.6923076923076923,
					0.6923076923076923 },
					{ 0.6410256428554484, 0.6538461533302703, 0.7435897424825051, 0.5512820494522439, 0.7435897450756739,
							0.8589743591115507 },
			},
			{ { 0.5384615384616484, 0.6794871794873101, 0.632478632417269, 0.632478632417269, 0.4615384615383516, 0.6324786326012289 },
					{ 0.5064102564102571, 0.6282051282051282, 0.5512820512820513, 0.5512820512820513, 0.5064102564102558,
							0.6282051282051282 },
			},
			{ { 0.6410256410256412, 0.6794871794871795, 0.8717948717948719, 0.4230769230769231, 0.6794871794871795,
					0.7435897435897435 },
					{ 0.6538461538461539, 0.7179487179113061, 0.858974358974359, 0.5512820512820513, 0.8205128204858109,
							0.8205128205772421 },
			},
			{ { 0.5128205128205128, 0.8205128209482262, 0.7435897435997461, 0.5384615386153406, 0.538461538267354, 0.6666666662616407 },
					{ 0.5256410256410257, 0.7147435935752421, 0.7147435935498484, 0.7147435783133656, 0.5897435897314114,
							0.7147435935480813 },
			},
			{ { 0.7521367521367521, 0.7521367521367521, 0.8717948717948717, 0.5256410256410257, 0.717948717948718, 0.7521367521367521 },
					{ 0.7307692307692307, 0.7692307692307693, 0.858974358974359, 0.7051282051282051, 0.7051282051282051,
							0.8461538461538463 },
			},
			{ { 0.5769230769230769, 0.717948717948718, 0.7307692307692307, 0.5769230769230769, 0.6666666666666666, 0.576923076923077 },
					{ 0.6346153803858957, 0.6346153888796556, 0.7435897435897436, 0.4999999999898454, 0.634615388820246,
							0.6346153803858957 },
			},
			{ { 0.8504273504211976, 0.9487179487179493, 0.9743589743589743, 0.6794871794871795, 0.8504273504212087,
					0.8504273504396447 },
					{ 0.6794871794869027, 0.7905982918448843, 0.8076923076923077, 0.6538461538461539, 0.7905982918445264,
							0.790598288105738 },
			},
			{ { 0.5641025641025597, 0.6794871794871795, 0.6987179487179487, 0.6538461538461583, 0.6538461538461539,
					0.6987179487179488 },
					{ 0.717948717948718, 0.7435897435897435, 0.8076923076923077, 0.717948717948718, 0.7435897435897435,
							0.9102564102564102 },
			},
			{ { 0.634615384498341, 0.6025641025641025, 0.6923076923076923, 0.5769230769230769, 0.6346153847323316, 0.7307692307693274 },
					{ 0.6153846153813314, 0.5897435897435898, 0.7222222204350568, 0.5641025641025641, 0.7222222231151475,
							0.7222222231197463 },
			},
			{ { 0.6196581196581187, 0.717948717948718, 0.7820512820512821, 0.4999999999999999, 0.6196581196581216, 0.6196581196581187 },
					{ 0.647435897435903, 0.713675213570946, 0.7136752138837794, 0.6474358974358817, 0.6538461538461539,
							0.7136752135709259 },
			}, };

	private double[] fExpected = { 0.00641025641032409,
			0.0598290722973402,
			0.205128205128205,
			0.710470114363538,
			0.0833333339996483,
			0.681623935140027,
			0.358974381671449,
			0.521367528760678,
			0.346153851340185,
			0.194444444935008,
			0.435897436080298,
			1.08653850133574,
			0.64957264957265,
			0.448717965883704,
			0.170940173446912,
			0.160256410256411,
			0.245726499564599,
			0.79273504315248 };

	private double[] gExpected = { 0.029133053262594,
			0.284404641248016,
			1.143832639091350,
			3.095654947167363,
			0.385529842613039,
			3.082865966584057,
			1.600538974642498,
			2.782001501548158,
			1.531509070781395,
			0.816885267161333,
			2.064096253252863,
			5.312586340042671,
			3.309665018119363,
			1.867676699107720,
			1.044493600569392,
			0.781225269023814,
			1.212706913473313,
			3.709622484334310 };

	@Test
	public void binCMRBigTest() {
		BinModel binModel = makeBinModel();

		BinCMRSolver solver = new BinCMRSolver();

		BinProblem[] problems = new BinProblem[binModel.getnSubj()];
		for (int i = 0; i < binModel.getnSubj(); i++) {
			problems[i] = new BinProblem(binModel, i, null);
		}

		BinSolution[] solns = solver.solve(problems);

		for (int i = 0; i < binModel.getnSubj(); i++) {
			assertEquals(gExpected[i], solns[i].getG2Star(), 1e-6);
			assertEquals(fExpected[i], solns[i].getFStar(), 1e-6);

			assertArrayEquals(xExpected[i][0], solns[i].getXStar(0), 1e-7);
			assertArrayEquals(xExpected[i][1], solns[i].getXStar(1), 1e-7);
		}
	}

	@Test
	public void binCMRBigTestDefModelCMRx() {
		BinModel binModel = makeBinModel();

		BinCMRxSolver solver = new BinCMRxSolver();

		BinProblem[] problems = new BinProblem[binModel.getnSubj()];
		for (int i = 0; i < binModel.getnSubj(); i++) {
			problems[i] = new BinProblem(binModel, i, null);
		}

		BinSolution[] solns = solver.solve(problems);

		for (int i = 0; i < binModel.getnSubj(); i++) {
			assertEquals(gExpected[i], solns[i].getG2Star(), 1e-6);
			assertEquals(fExpected[i], solns[i].getFStar(), 1e-6);

			assertArrayEquals(xExpected[i][0], solns[i].getXStar(0), 1e-7);
			assertArrayEquals(xExpected[i][1], solns[i].getXStar(1), 1e-7);
		}
	}

	@Test
	public void binCMRBigTestDefModelParCMRx() {
		BinModel binModel = makeBinModel();

		ParBinCMRxSolver solver = new ParBinCMRxSolver();

		BinProblem[] problems = new BinProblem[binModel.getnSubj()];
		for (int i = 0; i < binModel.getnSubj(); i++) {
			problems[i] = new BinProblem(binModel, i, null);
		}

		BinSolution[] solns = solver.solve(problems);

		for (int i = 0; i < binModel.getnSubj(); i++) {
			assertEquals(gExpected[i], solns[i].getG2Star(), 1e-6);
			assertEquals(fExpected[i], solns[i].getFStar(), 1e-6);

			assertArrayEquals(xExpected[i][0], solns[i].getXStar(0), 1e-7);
			assertArrayEquals(xExpected[i][1], solns[i].getXStar(1), 1e-7);
		}
	}

	@Test
	public void binCMRBigTestDefModelSimpleCMRx() {
		BinModel binModel = makeBinModel();

		BinCMRxSolver solver = new BinCMRxSolver();
		solver.setOnlyFeas(true);
		BinCMRxSolver solverFull = new BinCMRxSolver();

		BinProblem[] problems = new BinProblem[binModel.getnSubj()];
		for (int i = 0; i < binModel.getnSubj(); i++) {
			problems[i] = new BinProblem(binModel, i, null);
		}

		BinSolution[] solns = solver.solve(problems);
		BinSolution[] solnFull = solverFull.solve(problems);

		for (int i = 0; i < binModel.getnSubj(); i++) {
			assertTrue(solnFull[i].getG2Star() <= solns[i].getG2Star());
		}
	}

	@Test
	public void bigTestSingle() {
		BinModel model = new BinModel(1, 4);
		model.setElement(
				new BinElement(
						new int[] { 48, 18, 4, 9, 6, 39, 41, 47, 50, 64, 29, 10, 16, 28, 35, 41 },
						new int[] { 57, 87, 101, 96, 99, 66, 64, 58, 55, 41, 76, 95, 89, 77, 70, 64 }),
				0, 0);
		model.setElement(
				new BinElement(
						new int[] { 53, 42, 19, 14, 11, 12, 9, 14, 12, 61, 39, 26, 19, 19, 16, 6 },
						new int[] { 52, 63, 86, 91, 94, 93, 96, 91, 93, 44, 66, 79, 86, 86, 89, 99 }),
				0, 1);
		model.setElement(
				new BinElement(
						new int[] { 33, 23, 30, 24, 24, 42, 40, 41, 41, 51, 34, 27, 40, 36, 41, 34 },
						new int[] { 72, 82, 75, 81, 81, 63, 65, 64, 64, 54, 71, 78, 65, 69, 64, 71 }),
				0, 2);
		model.setElement(
				new BinElement(
						new int[] { 55, 39, 28, 25, 16, 17, 14, 20, 19, 69, 44, 30, 26, 24, 14, 11 },
						new int[] { 50, 66, 77, 80, 89, 88, 91, 85, 86, 36, 61, 75, 79, 81, 91, 94 }),
				0, 3);
		BinCMRxSolver solver = new BinCMRxSolver();

		BinBaseProblem bbp = new BinBaseProblem(model, null);

		BinSolution[] solns = solver.solve(bbp);

		assertEquals(96.16838217424265, solns[0].getG2Star(), 1e-8);
		assertEquals(16.91508967382458, solns[0].getFStar(), 1e-8);
	}

	@Test
	public void bigTestSinglePar() {
		BinModel model = new BinModel(1, 4);
		model.setElement(
				new BinElement(
						new int[] { 48, 18, 4, 9, 6, 39, 41, 47, 50, 64, 29, 10, 16, 28, 35, 41 },
						new int[] { 57, 87, 101, 96, 99, 66, 64, 58, 55, 41, 76, 95, 89, 77, 70, 64 }),
				0, 0);
		model.setElement(
				new BinElement(
						new int[] { 53, 42, 19, 14, 11, 12, 9, 14, 12, 61, 39, 26, 19, 19, 16, 6 },
						new int[] { 52, 63, 86, 91, 94, 93, 96, 91, 93, 44, 66, 79, 86, 86, 89, 99 }),
				0, 1);
		model.setElement(
				new BinElement(
						new int[] { 33, 23, 30, 24, 24, 42, 40, 41, 41, 51, 34, 27, 40, 36, 41, 34 },
						new int[] { 72, 82, 75, 81, 81, 63, 65, 64, 64, 54, 71, 78, 65, 69, 64, 71 }),
				0, 2);
		model.setElement(
				new BinElement(
						new int[] { 55, 39, 28, 25, 16, 17, 14, 20, 19, 69, 44, 30, 26, 24, 14, 11 },
						new int[] { 50, 66, 77, 80, 89, 88, 91, 85, 86, 36, 61, 75, 79, 81, 91, 94 }),
				0, 3);
		ParBinCMRxSolver solver = new ParBinCMRxSolver();

		BinBaseProblem bbp = new BinBaseProblem(model, null);

		BinSolution[] solns = solver.solve(bbp);
		
		assertEquals(96.16838217424265, solns[0].getG2Star(), 1e-8);
		assertEquals(16.91508967382458, solns[0].getFStar(), 1e-8);
	}

	@Test
	public void binCMRBigTestAltModelCMRx() {
		BinModel binModel = makeBinModel();
		BinCMRxSolver solver = new BinCMRxSolver();
		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { -1 } });
		BinProblem[] problems = new BinProblem[binModel.getnSubj()];
		for (int i = 0; i < binModel.getnSubj(); i++) {
			problems[i] = new BinProblem(binModel, i, null, cmrModel);
		}

		BinSolution[] solns = solver.solve(problems);
		double[] gExpectedModel1m1 = { 6.28673559071436,
				3.9410308627660497,
				10.880375799829133,
				3.7746413483525596,
				3.301589152119811,
				5.670013487918262,
				7.293199453433908,
				10.088324830794589,
				4.592401406191948,
				5.5051253496781705,
				28.266117830617855,
				9.884639265850378,
				9.866141389995416,
				4.925393061276791,
				10.060606362669596,
				4.1299219906412015,
				4.920868907974132,
				3.0364746800668994 };

		for (int i = 0; i < binModel.getnSubj(); i++) {
			assertEquals(gExpectedModel1m1[i], solns[i].getG2Star(), 1e-6);
		}
	}

	@Test
	public void fitsTest() {
		BinModel binModel = makeBinModel();
		BinBaseProblem p = new BinBaseProblem(binModel, new int[][] { { 3, 4 } });
		BinCMRFits fits = new BinCMRFits(10, p);

		assertEquals(18, fits.getP().length);
		// Difficult to test random output!
		// System.out.println(Arrays.toString(fits.getP()));
	}

	@Test
	public void fitsxTest() {
		BinModel binModel = makeBinModel();

		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { 1 } });
		BinBaseProblem p = new BinBaseProblem(binModel, null, cmrModel);

		int nSample = 5;

		BinCMRxFits fits = new BinCMRxFits(nSample, p, -1, false, false);
		assertEquals(18, fits.getP().length);

		// Assuming we cancel...
		assertTrue(fits.getFits().length == nSample);

		// Difficult to test random output!
		// System.out.println(Arrays.toString(fits.getP()));
	}

	//@Test
	public void fitsxCancelTest() {
		BinModel binModel = makeBinModel();

		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { -1 } });
		BinBaseProblem p = new BinBaseProblem(binModel, new int[][] { { 3, 4 } }, cmrModel);

		int nSample = 100000;

		BinCMRxFits fits = new BinCMRxFits(nSample, p, -1, false, true);
		assertEquals(18, fits.getP().length);

		// Assuming we cancel...
		assertTrue(fits.getFits().length < nSample);

		// Difficult to test random output!
		// System.out.println(Arrays.toString(fits.getP()));
	}

	@Test
	public void fitsxApproxTest() {
		BinModel binModel = makeBinModel();

		DoubleMatrix2D cmrModel = new DenseDoubleMatrix2D(new double[][] { { 1 }, { -1 } });
		BinBaseProblem p = new BinBaseProblem(binModel, new int[][] { { 3, 4 } }, cmrModel);

		int nSample = 100;

		BinCMRxFits fits = new BinCMRxFits(nSample, p, -1, true, false);
		assertEquals(18, fits.getP().length);

		// Assuming we dont cancel...
		assertTrue(fits.getFits().length == nSample);

		// Difficult to test random output!
		// System.out.println(Arrays.toString(fits.getP()));
	}

	@Test
	public void mrNoCTest() {
		BinModel binModel = makeBinModel();
		BinBaseProblem p = new BinBaseProblem(binModel, new int[][] { { 3, 4 } });
		BinCMRFits f = new BinCMRFits(10, p, -1, true);
		// System.out.println(Arrays.toString(f.getBaseFitDiff()));
		// Difficult to test random output!

		// System.out.println(Arrays.toString(fits.getP()));
	}

	@Test
	public void fitsETest() {
		double[] expectedDatafit = { 0.00975472359921814, 0.156111096940363, 1.02048305451518, 1.71678713485153, 0.199714048709331,
				3.08286603183291, 1.60053887224234, 2.75224597024015, 1.28261832234251, 0.685697670388934, 2.06409625478214,
				5.3125863017906, 3.30966501986503, 1.62227018798227, 1.04449368628261, 0.781225173717857, 1.18506096852275,
				3.70962252422463 };

		BinModel binModel = makeBinModel();
		BinBaseProblem p = new BinBaseProblem(binModel, new int[][] { { 0, 1, 2 } });
		BinCMRFits f = new BinCMRFits(10, p);

		assertArrayEquals(expectedDatafit, f.getBaseFitDiff(), 1e-6);
	}

	@Test
	public void valueCloseToZeroTest() {
		BinModel model = new BinModel(1, 2);
		model.setElement(new BinElement(new int[] { 2, 54, 68, 35, 50, 56 }, new int[] { 1, 24, 10, 43, 28, 22 }), 0, 0);
		model.setElement(new BinElement(new int[] { 0, 52, 60, 48, 53, 59 }, new int[] { 25, 26, 18, 30, 25, 19 }), 0, 1);

		BinBaseProblem p = new BinBaseProblem(model, new int[][] { { 0, 1, 2 } });
		BinCMRFits f = new BinCMRFits(10000, p);

		// System.out.println(Arrays.toString(f.getP()));

		// assertArrayEquals(expectedDatafit, f.getBaseFitDiff(), 1e-6);
	}

	@Test
	public void zeroMissTest() {
		BinModel model = new BinModel(1, 2);
		model.setElement(new BinElement(new int[] { 70, 76, 78, 54, 68, 66 }, new int[] { 8, 2, 0, 24, 10, 12 }), 0, 0);
		model.setElement(new BinElement(new int[] { 39, 55, 66, 53, 62, 58 }, new int[] { 39, 23, 12, 25, 16, 20 }), 0, 1);
		BinProblem p = new BinProblem(model, 0, new int[][] { { 3, 4 } });

		BinCMRSolver solver = new BinCMRSolver();
		BinSolution soln = solver.solve(p);

		assertNotNull(soln.getXStar());

	}

	private BinModel makeBinModel() {
		BinModel model = new BinModel(18, 2);
		model.setElement(new BinElement(new int[] { 56, 54, 68, 35, 50, 56 }, new int[] { 22, 24, 10, 43, 28, 22 }), 0, 0);
		model.setElement(new BinElement(new int[] { 53, 52, 60, 48, 53, 59 }, new int[] { 25, 26, 18, 30, 25, 19 }), 0, 1);
		model.setElement(new BinElement(new int[] { 57, 47, 57, 46, 50, 51 }, new int[] { 21, 31, 21, 32, 28, 27 }), 1, 0);
		model.setElement(new BinElement(new int[] { 53, 51, 57, 47, 55, 56 }, new int[] { 25, 27, 21, 31, 23, 22 }), 1, 1);
		model.setElement(new BinElement(new int[] { 56, 54, 71, 47, 55, 52 }, new int[] { 22, 24, 7, 31, 23, 26 }), 2, 0);
		model.setElement(new BinElement(new int[] { 47, 52, 61, 45, 65, 62 }, new int[] { 31, 26, 17, 33, 13, 16 }), 2, 1);
		model.setElement(new BinElement(new int[] { 44, 56, 63, 44, 49, 48 }, new int[] { 34, 22, 15, 34, 29, 30 }), 3, 0);
		model.setElement(new BinElement(new int[] { 47, 54, 47, 53, 43, 55 }, new int[] { 31, 24, 31, 25, 35, 23 }), 3, 1);
		model.setElement(new BinElement(new int[] { 50, 50, 50, 45, 43, 52 }, new int[] { 28, 28, 28, 33, 35, 26 }), 4, 0);
		model.setElement(new BinElement(new int[] { 57, 54, 59, 43, 57, 64 }, new int[] { 21, 24, 19, 35, 21, 14 }), 4, 1);
		model.setElement(new BinElement(new int[] { 51, 61, 61, 49, 54, 50 }, new int[] { 27, 17, 17, 29, 24, 28 }), 5, 0);
		model.setElement(new BinElement(new int[] { 46, 47, 59, 43, 54, 57 }, new int[] { 32, 31, 19, 35, 24, 21 }), 5, 1);
		model.setElement(new BinElement(new int[] { 48, 53, 56, 41, 42, 48 }, new int[] { 30, 25, 22, 37, 36, 30 }), 6, 0);
		model.setElement(new BinElement(new int[] { 51, 56, 65, 53, 57, 60 }, new int[] { 27, 22, 13, 25, 21, 18 }), 6, 1);
		model.setElement(new BinElement(new int[] { 42, 42, 64, 52, 50, 54 }, new int[] { 36, 36, 14, 26, 28, 24 }), 7, 0);
		model.setElement(new BinElement(new int[] { 54, 53, 60, 48, 57, 66 }, new int[] { 24, 25, 18, 30, 21, 12 }), 7, 1);
		model.setElement(new BinElement(new int[] { 44, 51, 48, 49, 56, 52 }, new int[] { 34, 27, 30, 29, 22, 26 }), 8, 0);
		model.setElement(new BinElement(new int[] { 50, 51, 59, 43, 57, 67 }, new int[] { 28, 27, 19, 35, 21, 11 }), 8, 1);
		model.setElement(new BinElement(new int[] { 42, 53, 50, 50, 36, 48 }, new int[] { 36, 25, 28, 28, 42, 30 }), 9, 0);
		model.setElement(new BinElement(new int[] { 38, 47, 43, 43, 41, 51 }, new int[] { 40, 31, 35, 35, 37, 27 }), 9, 1);
		model.setElement(new BinElement(new int[] { 50, 57, 68, 33, 49, 58 }, new int[] { 28, 21, 10, 45, 29, 20 }), 10, 0);
		model.setElement(new BinElement(new int[] { 51, 56, 67, 43, 65, 63 }, new int[] { 27, 22, 11, 35, 13, 15 }), 10, 1);
		model.setElement(new BinElement(new int[] { 40, 64, 58, 40, 44, 52 }, new int[] { 38, 14, 20, 38, 34, 26 }), 11, 0);
		model.setElement(new BinElement(new int[] { 41, 50, 54, 57, 46, 62 }, new int[] { 37, 28, 24, 21, 32, 16 }), 11, 1);
		model.setElement(new BinElement(new int[] { 61, 61, 68, 41, 56, 54 }, new int[] { 17, 17, 10, 37, 22, 24 }), 12, 0);
		model.setElement(new BinElement(new int[] { 57, 60, 67, 58, 52, 66 }, new int[] { 21, 18, 11, 20, 26, 12 }), 12, 1);
		model.setElement(new BinElement(new int[] { 42, 56, 57, 49, 52, 44 }, new int[] { 36, 22, 21, 29, 26, 34 }), 13, 0);
		model.setElement(new BinElement(new int[] { 50, 47, 58, 39, 51, 50 }, new int[] { 28, 31, 20, 39, 27, 28 }), 13, 1);
		model.setElement(new BinElement(new int[] { 67, 74, 76, 53, 66, 66 }, new int[] { 11, 4, 2, 25, 12, 12 }), 14, 0);
		model.setElement(new BinElement(new int[] { 53, 59, 63, 51, 64, 62 }, new int[] { 25, 19, 15, 27, 14, 16 }), 14, 1);
		model.setElement(new BinElement(new int[] { 44, 53, 56, 51, 51, 53 }, new int[] { 34, 25, 22, 27, 27, 25 }), 15, 0);
		model.setElement(new BinElement(new int[] { 58, 58, 63, 54, 58, 71 }, new int[] { 20, 20, 15, 24, 20, 7 }), 15, 1);
		model.setElement(new BinElement(new int[] { 50, 47, 54, 45, 49, 57 }, new int[] { 28, 31, 24, 33, 29, 21 }), 16, 0);
		model.setElement(new BinElement(new int[] { 48, 46, 57, 44, 59, 53 }, new int[] { 30, 32, 21, 34, 19, 25 }), 16, 1);
		model.setElement(new BinElement(new int[] { 50, 56, 61, 39, 48, 47 }, new int[] { 28, 22, 17, 39, 30, 31 }), 17, 0);
		model.setElement(new BinElement(new int[] { 47, 52, 55, 54, 51, 60 }, new int[] { 31, 26, 23, 24, 27, 18 }), 17, 1);
		return model;
	}
}
