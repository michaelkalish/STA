package au.edu.adelaide.fxmr.model;

/**
 * Simple class that is stored as a record of each CMR iteration
 * 
 * 
 */
public class CMRIter {

	private double floor;
	private double upper;
	private double upperFloor;
	private int remainig;

	public CMRIter(double floor, double upper, double upperFloor, int remaining) {
		this(floor, upper, remaining);
		this.upperFloor = upperFloor;
	}

	public CMRIter(double floor, double upper, int remaining) {
		this.floor = floor;
		this.upper = upper;
		this.remainig = remaining;
	}

	public double getFloor() {
		return floor;
	}

	public double getUpper() {
		return upper;
	}

	@Override
	public String toString() {
		return floor + "\t" + upper + (upperFloor == 0 ? "" : "(" + upperFloor + ")");
	}

	public int getRemaining() {
		return remainig;
	}

	public double getUpperFloor() {
		return upperFloor;
	}

}
