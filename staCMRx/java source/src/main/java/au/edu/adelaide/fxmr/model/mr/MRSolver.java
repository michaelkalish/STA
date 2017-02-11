package au.edu.adelaide.fxmr.model.mr;

import au.edu.adelaide.fxmr.model.CMRSolver;

public abstract class MRSolver {
	public static final double TOL1 = 1e-10;
	public static final double TOL2 = CMRSolver.ZERO_TOL;

	protected double tolInitFeas = TOL1;
	protected double tolInit = TOL1 * 10;
	protected double tol2Feas = TOL2;
	protected double tol2 = TOL2 * 10;
	protected int calls;

	public abstract MRSolution solve(MRProblem p);// throws
													// CatastrophicMRFailure;

	public int getCalls() {
		return calls;
	}

	/**
	 * Reset to the original tolerance
	 */
	public void resetTolerance() {
		tolInitFeas = TOL1;
		tolInit = TOL1 * 10;
		tol2Feas = TOL2;
		tol2 = TOL2 * 10;
	}

	/**
	 * Change the tolerance. Lower values mean more accuracy but take longer to
	 * solve.
	 * 
	 * @param tol1
	 *            The initial tolerance to use. The secondary tolerance is also
	 *            set based on this value
	 */
	public void setTolerance(double tol1) {
		setTolerance(tol1, tol1 * 10000);
	}

	/**
	 * Change the tolerance. If zero is passed in, reset the tolerance
	 * (convenience behaviour)
	 * 
	 * @param tol1
	 *            the initial tolerance to use
	 * @param tol2
	 *            when the first attempt fails, use this tolerance along with
	 *            other more conservative parameters
	 */
	public void setTolerance(double tol1, double tolSecond) {
		if (tol1 <= 0) {
			resetTolerance();
		} else {
			tolInitFeas = tol1;
			tolInit = tol1 * 10;
			tol2Feas = tolSecond;
			tol2 = tolSecond * 10;
		}
	}
}
