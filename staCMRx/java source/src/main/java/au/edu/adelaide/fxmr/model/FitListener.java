package au.edu.adelaide.fxmr.model;

public interface FitListener {

	public void setFinished();

	public boolean updateStatus(double[] fits, double dataFit);

}
