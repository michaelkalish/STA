package au.edu.adelaide.fxmr.eqoptimiser;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.adelaide.fxmr.joptimizer.functions.QuadraticMultivariateRealFunction;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class EqTests {
	@Test
	public void simpleEQTest() {
		SimpleLinearConstraint[] eqMat = new SimpleLinearConstraint[3];
		eqMat[0] = new SimpleLinearConstraint(3, 4);
		eqMat[1] = new SimpleLinearConstraint(2, 3);
		eqMat[2] = new SimpleLinearConstraint(0, 5);

		double[][] pMatrix = new double[6][6];
		for (int i = 0; i < 6; i++)
			pMatrix[i][i] = 78;
		QuadraticMultivariateRealFunction objectiveFunction = new QuadraticMultivariateRealFunction(pMatrix, new double[] { -57, -47, -57, -46, -50, -51 }, 0);
		EQSolver eqs = new EQSolver(objectiveFunction, eqMat);

		double[] expected = new double[] { 0.6923076923076923, 0.6025641025641025, 0.6538461538461539, 0.6538461538461539, 0.6538461538461539,
				0.6923076923076923 };
		double[] sol = eqs.solve();

		assertEquals(-101.56410256410258, objectiveFunction.value(new DenseDoubleMatrix1D(sol)), 1e-15);

		assertArrayEquals(expected, eqs.solve(), 1e-15);
	}
}
