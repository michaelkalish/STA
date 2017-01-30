package au.edu.adelaide.fxmr.model.bin;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;

public class BinSolution {

	private double fStar;
	private double g2Star;
	private double[][] xStar;
	private int iter;
	private double seconds;
	private int mrCalls;
	private HashSet<SimpleLinearConstraint>[] adjStar;

	public double getSeconds() {
		return seconds;
	}

	public BinSolution(double fStar, double g2Star, double[][] xStar, int iter, double seconds,
			HashSet<SimpleLinearConstraint>[] adjStar, int mrCalls) {
		this.fStar = fStar;
		this.g2Star = g2Star;
		this.xStar = xStar;
		this.iter = iter;
		this.seconds = seconds;
		this.mrCalls = mrCalls;
		this.adjStar = adjStar;
	}

	public int getMRCalls() {
		return mrCalls;
	}

	public double[][] getXStar() {
		return xStar;
	}

	public double[] getXStar(int col) {
		return xStar[col];
	}

	public int getIter() {
		return iter;
	}

	public double getFStar() {
		return fStar;
	}

	public HashSet<SimpleLinearConstraint>[] getAdjs() {
		return adjStar;
	}

	public double[][] getAdjMatrix(int index) {
		if (adjStar == null)
			return null;
		HashSet<SimpleLinearConstraint> curSet = adjStar[index];
		if (curSet == null)
			return null;

		int nCond = xStar[0].length;
		double[][] ret = new double[nCond][nCond];

		for (SimpleLinearConstraint c : curSet)
			ret[c.getPosIndex()][c.getNegIndex()] = 1;

		return ret;
	}

	public double getG2Star() {
		return g2Star;
	}
}
