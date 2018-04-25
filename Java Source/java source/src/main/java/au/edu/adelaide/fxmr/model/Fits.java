package au.edu.adelaide.fxmr.model;

public interface Fits {

	double getP();

	double[] getFits();

	double getDataFit();

	Badness[] getBadnesses();
	
	double[] getTimes();
}
