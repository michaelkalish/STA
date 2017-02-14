package au.edu.adelaide.fxmr.om;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class OMUtilTests {
	DoubleMatrix2D a = new DenseDoubleMatrix2D(new double[][] { { 1, 0, 1 }, { 0, 1, 1 } }).viewDice();
	DoubleMatrix2D b = new DenseDoubleMatrix2D(new double[][] { { 1, 0, -1 }, { 0, 1, -1 } }).viewDice();
	DoubleMatrix2D c = new DenseDoubleMatrix2D(new double[][] { { 1, 0, -1 }, { 0, 1, -1 }, { 1, 1, 1 } }).viewDice();
	DoubleMatrix2D d = new DenseDoubleMatrix2D(new double[][] { { 1, 0, -1 }, { 0, 1, -1 }, { -1, 0, 1 } }).viewDice();
	DoubleMatrix2D e = new DenseDoubleMatrix2D(
			new double[][] { { 1, 0, -1 }, { 0, 1, -1 }, { -1, 0, 1 }, { 2, 0, -2 } });

	@Test
	public void checkRankTest() {
		assertEquals(a, OMUtil.checkRank(a));
		assertEquals(b, OMUtil.checkRank(b));
		assertEquals(c, OMUtil.checkRank(c));

		DoubleMatrix2D d2 = OMUtil.checkRank(d);
		assertArrayEquals(new double[] { 1, 0, -1 }, d2.viewColumn(0).toArray(), 1e-15);
		assertArrayEquals(new double[] { 0, 1, -1 }, d2.viewColumn(1).toArray(), 1e-15);

		DoubleMatrix2D e2 = OMUtil.checkRank(e);
		assertArrayEquals(new double[] { 1, 0, -1, 2 }, e2.viewColumn(0).toArray(), 1e-15);
		assertArrayEquals(new double[] { 0, 1, 0, 0 }, e2.viewColumn(1).toArray(), 1e-15);
	}

	@Test
	public void minorsTest() {
		assertArrayEquals(new double[] { 1, 1, -1 }, OMUtil.minors(a).toArray(), 1e-15);
		assertArrayEquals(new double[] { 1, -1, 1 }, OMUtil.minors(b).toArray(), 1e-15);
		assertArrayEquals(new double[] { 3 }, OMUtil.minors(c).toArray(), 1e-15);
		assertArrayEquals(new double[] { 0, 0, 0, 0 }, OMUtil.minors(e).toArray(), 1e-15);
	}

	@Test
	public void orthocomTest() {
		assertArrayEquals(new double[] { -1, -1, 1 }, OMUtil.orthocom(a).viewRow(0).toArray(), 1e-15);
		assertArrayEquals(new double[] { 1, 1, 1 }, OMUtil.orthocom(b).viewRow(0).toArray(), 1e-15);
		assertNull(OMUtil.orthocom(c));
		assertArrayEquals(new double[] { -1, 0, 0, 0.5 }, OMUtil.orthocom(e).viewRow(0).toArray(), 1e-15);
		assertArrayEquals(new double[] { 0, 0, -1, -0.5 }, OMUtil.orthocom(e).viewRow(1).toArray(), 1e-15);
	}

	@Test
	public void circuitsTest() {
		assertArrayEquals(new int[] { 1, 1, -1 }, OMUtil.circuits(OMUtil.checkRank(a))[0]);
		assertArrayEquals(new int[] { 1, 0, 1, 0 }, OMUtil.circuits(OMUtil.checkRank(e))[0]);
		assertArrayEquals(new int[] { 1, 0, 0, -1 }, OMUtil.circuits(OMUtil.checkRank(e))[1]);
		assertArrayEquals(new int[] { 0, 0, 1, 1 }, OMUtil.circuits(OMUtil.checkRank(e))[2]);
	}

	@Test
	public void zEncodeTest() {
		int[][] d = { { 0, 0, 1 }, { 0, 1, 0 }, { 0, 0, 0 }, { 0, -1, -1 }, { 1, 0, -1 }, { 1, -1, 0 }, { 0, 0, -1 },
				{ 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, -1, 1 },
				{ 1, 0, 1 }, { 1, -1, 0 }, { 0, 0, 1 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 },
				{ 1, 0, 0 }, { 0, 0, 0 }, { 0, 1, 1 }, { 1, 0, 1 }, { 1, 1, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, { 1, 0, 0 },
				{ 0, 0, 0 }, { 0, 1, 1 }, { 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, 1 }, { 1, 1, 1 } };

		assertArrayEquals(new int[] { 1, 3, 0, -4, 8, 6, -1, -3, 9, 0, -3, 9, 0, -2, 10, 6, 1, -3, 9, 0, 1, 9, 0, 4, 10,
				12, 1, 3, 9, 0, 4, 5, 6, 7, 10, 13 }, OMUtil.zEncode(d));
	}

	@Test
	public void zDecodeTest() {
		int[][] zs = OMUtil.zDecode(new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13 }, 3);

		assertArrayEquals(new int[] { 0, 0, 1 }, zs[0]);
		assertArrayEquals(new int[] { 0, 1, -1 }, zs[1]);
		assertArrayEquals(new int[] { 0, 1, 0 }, zs[2]);
		assertArrayEquals(new int[] { 0, 1, 1 }, zs[3]);
		assertArrayEquals(new int[] { 1, -1, -1 }, zs[4]);
		assertArrayEquals(new int[] { 1, -1, 0 }, zs[5]);
		assertArrayEquals(new int[] { 1, -1, 1 }, zs[6]);
		assertArrayEquals(new int[] { 1, 0, -1 }, zs[7]);
		assertArrayEquals(new int[] { 1, 0, 0 }, zs[8]);
		assertArrayEquals(new int[] { 1, 0, 1 }, zs[9]);
		assertArrayEquals(new int[] { 1, 1, 0 }, zs[10]);
		assertArrayEquals(new int[] { 1, 1, 1 }, zs[11]);
	}

	@Test
	public void allSignVectorsTest() {
		int[][] exp = { { 0, 0, 1 }, { 0, 1, -1 }, { 0, 1, 0 }, { 0, 1, 1 }, { 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 },
				{ 1, 0, -1 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, -1 }, { 1, 1, 0 }, { 1, 1, 1 } };
		int[][] s3 = OMUtil.allSignVectors(3);
		for (int i = 0; i < exp.length; i++) {
			assertArrayEquals(exp[i], s3[i]);
		}
	}

	@Test
	public void signconTest() {
		boolean[][] exp = { { false, false, true }, { true, true, false }, { true, false, true }, { true, false, true },
				{ false, true, false }, { false, false, false }, { false, false, false }, { false, true, false },
				{ false, false, false }, { false, false, true }, { true, true, false }, { true, false, true },
				{ true, false, true }, { false, true, false }, { false, true, false }, { false, false, false },
				{ true, true, false }, { true, true, false }, { true, false, true }, { true, true, false },
				{ true, true, true }, { true, false, true }, { false, true, false }, { false, true, false },
				{ false, false, false }, { true, true, false }, { true, true, false }, { true, false, true },
				{ true, true, false }, { true, true, true }, { true, false, true }, { false, true, false },
				{ false, true, false }, { false, false, false }, { true, true, false }, { true, true, false },
				{ true, false, true }, { true, true, false }, { true, true, true }, { true, false, true } };
		int[][] z = { { 1, 0, 1, 0 }, { 1, 0, 0, -1 }, { 0, 0, 1, 1 } };
		int[][] x = OMUtil.allSignVectors(4);

		boolean[][] actual = OMUtil.signCon(x, z);
		for (int i = 0; i < exp.length; i++) {
			assertArrayEquals(exp[i], actual[i]);
		}
	}

	@Test
	public void covectorsTest() {
		int[][] exp = { { 0, 1, 1 }, { 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, 1 }, { 1, 1, 1 } };
		int[][] actual = OMUtil.covectors(a);
		for (int i = 0; i < exp.length; i++)
			assertArrayEquals(exp[i], actual[i]);

		int[][] exp2 = { { 0, 1, 0, 0 }, { 1, -1, -1, 1 }, { 1, 0, -1, 1 }, { 1, 1, -1, 1 } };
		int[][] actual2 = OMUtil.covectors(OMUtil.checkRank(e));
		for (int i = 0; i < exp2.length; i++)
			assertArrayEquals(exp2[i], actual2[i]);
	}

	@Test
	public void signClosureTest() {
		int[][] exp = { { 0, 0, 1 }, { 0, 1, -1 }, { 0, 1, 0 }, { 0, 1, 1 }, { 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 },
				{ 1, 0, -1 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, 0 }, { 1, 1, 1 } };
		int[][] actual = OMUtil.signClosure(OMUtil.covectors(a));
		for (int i = 0; i < exp.length; i++)
			assertArrayEquals(exp[i], actual[i]);

		int[][] exp2 = { { 0, 0, 0, 1 }, { 0, 0, 1, -1 }, { 0, 0, 1, 0 }, { 0, 1, -1, 0 }, { 0, 1, -1, 1 },
				{ 0, 1, 0, -1 }, { 0, 1, 0, 0 }, { 0, 1, 0, 1 }, { 0, 1, 1, -1 }, { 0, 1, 1, 0 }, { 1, -1, -1, 0 },
				{ 1, -1, -1, 1 }, { 1, -1, 0, 0 }, { 1, -1, 0, 1 }, { 1, 0, -1, 0 }, { 1, 0, -1, 1 }, { 1, 0, 0, 0 },
				{ 1, 0, 0, 1 }, { 1, 1, -1, 0 }, { 1, 1, -1, 1 }, { 1, 1, 0, 0 }, { 1, 1, 0, 1 } };
		int[][] actual2 = OMUtil.signClosure(OMUtil.covectors(OMUtil.checkRank(e)));
		for (int i = 0; i < exp2.length; i++)
			assertArrayEquals(exp2[i], actual2[i]);
	}

	@Test
	public void vectorsTest() {
		int[] exp = { 1, 1, -1 };
		int[][] actual = OMUtil.vectors(a);
		assertArrayEquals(exp, actual[0]);

		int[][] actual2 = OMUtil.vectors(c);
		assertNull(actual2);
	}

	@Test
	public void signAugmentTest() {
		int[][] exp = OMUtil.vectors(a);
		int[][] actual = OMUtil.signAugment(OMUtil.vectors(a));
		for (int i = 0; i < exp.length; i++)
			assertArrayEquals(exp[i], actual[i]);

		int[][] exp2 = { { 0, 1, 1 }, { 1, -1, -1 }, { 1, -1, 1 }, { 1, 0, 1 }, { 1, 1, -1 }, { 1, 1, 1 } };
		int[][] actual2 = OMUtil.signAugment(new int[][] { { 1, 0, 1 }, { 0, 1, 1 }, { -1, -1, 1 } });
		for (int i = 0; i < exp2.length; i++)
			assertArrayEquals(exp2[i], actual2[i]);

		int[][] exp3 = { { 1, 1, -1 }, { 1, 1, 0 }, { 1, 1, 1 } };
		int[][] actual3 = OMUtil.signAugment(new int[][] { { 1, 1, 0 } });
		for (int i = 0; i < exp3.length; i++)
			assertArrayEquals(exp3[i], actual3[i]);
	}

	@Test
	public void dimTest() {
		int[] exp = { 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 };
		int[] actual = OMUtil.dim(OMUtil.allSignVectors(3));
		assertArrayEquals(exp, actual);
	}
}
