package au.edu.adelaide.fxmr.model;

import java.util.HashSet;
import java.util.List;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;

public class CMRSolution {

	private double fStar;
	private double[][] xStar;
	private List<CMRIter> iters;
	private double seconds;
	private int mrCalls;
	private HashSet<SimpleLinearConstraint>[] adjStar;
	private int fBarReductions;

	public double getSeconds() {
		return seconds;
	}

	public CMRSolution(double fStar, double[][] xStar, List<CMRIter> iters, double seconds,
			HashSet<SimpleLinearConstraint>[] adjStar, int mrCalls, int fBarReductions) {
		this.fStar = fStar;
		this.xStar = xStar;
		this.iters = iters;
		this.seconds = seconds;
		this.mrCalls = mrCalls;
		this.adjStar = adjStar;
		this.fBarReductions = fBarReductions;
	}

	public int getfBarReductions() {
		return fBarReductions;
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

	public List<CMRIter> getIters() {
		return iters;
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

	public double[] getFloors() {
		int n = iters.size();
		double[] ret = new double[n];
		for (int i = 0; i < n; i++)
			ret[i] = iters.get(i).getFloor();
		return ret;
	}

	public double[] getUppers() {
		int n = iters.size();
		double[] ret = new double[n];
		for (int i = 0; i < n; i++)
			ret[i] = iters.get(i).getUpper();
		return ret;
	}

	public double[] getFufs() {
		int n = iters.size();
		double[] ret = new double[n];
		for (int i = 0; i < n; i++)
			ret[i] = iters.get(i).getUpperFloor();
		return ret;
	}

	public int[] getRemainings() {
		int n = iters.size();
		int[] ret = new int[n];
		for (int i = 0; i < n; i++)
			ret[i] = iters.get(i).getRemaining();
		return ret;
	}
}
