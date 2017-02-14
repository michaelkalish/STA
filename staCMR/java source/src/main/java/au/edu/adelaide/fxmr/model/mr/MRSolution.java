package au.edu.adelaide.fxmr.model.mr;

import java.text.DecimalFormat;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.jet.math.Functions;

public class MRSolution {
	private static final DecimalFormat FMT = new DecimalFormat("0.000000");
	private double fVal;
	private double[] xVector;
	private double timeS;
	private MRProblem problem;

	public MRProblem getProblem() {
		return problem;
	}

	private int iterations;

	public int getIterations() {
		return iterations;
	}

	public MRSolution(MRProblem problem, double[] xVector, double timeS, int iterations) {
		super();
		this.problem = problem;
		this.xVector = xVector;
		this.timeS = timeS;
		this.iterations = iterations;

		// Calculate the f value for MR (different from the QP fVal)
		DoubleMatrix1D diff = new DenseDoubleMatrix1D(xVector);
		diff.assign(new DenseDoubleMatrix1D(problem.getY()), Functions.minus);

		DoubleMatrix1D lhs = problem.getWeights().zMult(diff, new DenseDoubleMatrix1D(problem.getN()));
		fVal = lhs.zDotProduct(diff);
	}

	public double getTimeS() {
		return timeS;
	}

	public double[] getxVector() {
		return xVector;
	}

	public double getfVal() {
		return fVal;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		appendString(sb);
		return sb.toString();
	}

	public void appendString(StringBuilder sb) {
		sb.append("f = " + fVal + "\t(solved in " + timeS + "s)\n");
		sb.append("n\ty\t\tx\n");

		double[] y = problem.getY();

		for (int i = 0; i < problem.getN(); i++) {
			sb.append(i + "\t" + FMT.format(y[i]) + "\t" + FMT.format(xVector[i]) + "\n");
		}
	}
}
