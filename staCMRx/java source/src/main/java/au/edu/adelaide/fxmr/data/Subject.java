package au.edu.adelaide.fxmr.data;

/**
 * Subject of an experiment for the general model
 */
public class Subject {
	private int subjectNumber;
	private int betweenSubjectsCondition;
	private int dependentVariable;
	private double[] cond;

	public Subject(int subjectNumber, int betweenSubjectsCondition, int dependentVariable, double[] cond) {
		this.subjectNumber = subjectNumber;
		this.betweenSubjectsCondition = betweenSubjectsCondition;
		this.dependentVariable = dependentVariable;
		this.cond = cond;
	}

	public Subject(double[] row) {
		subjectNumber = (int) row[0];
		betweenSubjectsCondition = (int) row[1];
		dependentVariable = (int) row[2];

		// Ignore NAs - this allows people to import CSV files with non-values to have varying length.
		int condL = 0;
		for (int i = 3; i < row.length; i++)
			if (!Double.isNaN(row[i]))
				condL++;
		cond = new double[condL];
		condL = 0;
		for (int i = 3; i < row.length; i++)
			if (!Double.isNaN(row[i]))
				cond[condL++] = row[i];

		System.arraycopy(row, 3, cond, 0, condL);
	}

	public int getSubjectNumber() {
		return subjectNumber;
	}

	public void setSubjectNumber(int subjectNumber) {
		this.subjectNumber = subjectNumber;
	}

	public int getBetweenSubjectsCondition() {
		return betweenSubjectsCondition;
	}

	public void setBetweenSubjectsCondition(int betweenSubjectsCondition) {
		this.betweenSubjectsCondition = betweenSubjectsCondition;
	}

	public int getDependentVariable() {
		return dependentVariable;
	}

	public void setDependentVariable(int dependentVariable) {
		this.dependentVariable = dependentVariable;
	}

	public double[] getCond() {
		return cond;
	}

	public void setCond(double[] cond) {
		this.cond = cond;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(subjectNumber);
		sb.append(',');
		sb.append(betweenSubjectsCondition);
		sb.append(',');
		sb.append(dependentVariable);
		for (double c : cond) {
			sb.append(',');
			sb.append(c);
		}
		return sb.toString();
	}
}
