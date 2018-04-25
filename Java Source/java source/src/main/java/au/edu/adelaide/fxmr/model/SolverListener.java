package au.edu.adelaide.fxmr.model;

public interface SolverListener {

	public boolean updateStatus(String string);

	public void setFinished();

	/**
	 * Status update
	 * 
	 * @param fBar
	 * @param fBar 
	 * @param upperFloor 
	 * @param remaining
	 * @param nIterThread
	 * @param collisions
	 * @param fBarReductions 
	 * @param cyclicAvoided 
	 * 
	 * @return false if the user has manually cancelled, indicating to the solver
	 *         to stop ASAP
	 */
	public boolean updateStatus(double fFloor, double fBar, double upperFloor, int remaining, int[] nIterThread, int collisions, int fBarReductions, int cyclicAvoided);

}
