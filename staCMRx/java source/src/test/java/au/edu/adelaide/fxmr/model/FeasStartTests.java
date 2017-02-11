package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.model.mr.MRProblem;

public class FeasStartTests {

	@Test
	public void basicTest() {
		double[] y = { 1, 2, 1, 1, 1, 1 };

		MRProblem problem = new MRProblem(y, null);

		problem.addConstraint(1, 2, y);
		assertTrue(!problem.testFeas(y));

		double[] ss = problem.findSimpleStart(0.01, y);
		//System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));

		problem.addConstraint(2, 3, ss);
		ss = problem.findSimpleStart(0.01, ss);
		//System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));
		
		ss[3] -= 0.5;
		//System.out.println(Arrays.toString(ss));
		ss = problem.findSimpleStart(0.01, ss);
		//System.out.println(Arrays.toString(ss));
		assertTrue(problem.testFeas(ss));
		
		//Now make a cycle
		problem.addConstraint(3, 1, y);
		ss = problem.findSimpleStart(0.01, ss);
		assertNull(ss);
	}
}
