package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.model.mr.MRProblem;
import cern.colt.Arrays;

public class FeasStartTests {

	@Test
	public void basicTest() {
		double[] y = { 1, 2, 1, 1, 1, 1 };

		MRProblem problem = new MRProblem(y, null);

		problem.addConstraint(1, 2, y);
		assertTrue(!problem.testFeas(y));

		double[] ss = problem.findSimpleStart(0.01, y);
		// System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));

		problem.addConstraint(2, 3, ss);
		ss = problem.findSimpleStart(0.01, ss);
		// System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));

		ss[3] -= 0.5;
		// System.out.println(Arrays.toString(ss));
		ss = problem.findSimpleStart(0.01, ss);
		// System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));

		// Now make a cycle
		problem.addConstraint(3, 1, y);
		ss = problem.findSimpleStart(0.01, ss);
		assertTrue(ss == null);
	}

	@Test
	public void realTest() {
		double[] start = { 65.6, 48.6, 71.9, 62.4, 88.9, 98.98808573618385, 78.0, 85.66824966221256, 80.4, 87.7, 98.98808573665954, 85.66824966286786, 95.7,
				96.4, 82.7, 76.6 };
		MRProblem problem = new MRProblem(start, null);
		problem.addConstraint(11, 7, start);
		problem.addConstraint(11, 5, start);
		problem.addConstraint(12, 11, start);

		double[] ss = problem.findSimpleStart(1, start);

		// System.out.println(Arrays.toString(ss));
		// System.out.println(ss[11] + "<" + ss[7]);
		// System.out.println(ss[11] + "<" + ss[5]);
		// System.out.println(ss[12] + "<" + ss[11]);

		assertTrue(problem.testUnstrictFeas(ss));
	}
}
