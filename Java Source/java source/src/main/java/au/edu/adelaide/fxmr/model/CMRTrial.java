package au.edu.adelaide.fxmr.model;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import cern.colt.matrix.DoubleMatrix2D;

public class CMRTrial implements Comparable<CMRTrial> {
	private double[][] xPrime;
	private double f;
	private MRProblem[] problems;
	private MRSolution[] solutions;
	private int index;

	public CMRTrial(int[][] rangeSet, int nvar, DoubleMatrix2D[] weights, double[][] means) {
		xPrime = new double[nvar][];
		problems = new MRProblem[nvar];
		solutions = new MRSolution[nvar];
		f = Double.NEGATIVE_INFINITY;
		index = -1;
		for (int i = 0; i < nvar; i++)
			problems[i] = new MRProblem(means[i], weights[i], rangeSet);
	}

	private CMRTrial() {
	}

	@Override
	public int compareTo(CMRTrial o) {
		int compareF = Double.compare(f, o.f);
		if (compareF == 0)
			return Integer.compare(index, o.index);
		return compareF;
	}

	public double[][] solveMRs(MRSolver solver) {
		int n = problems.length;
		f = 0;
		for (int i = 0; i < n; i++) {
			if (solutions[i] == null)
				solutions[i] = solver.solve(problems[i]);
			if (solutions[i] == null)
				// Failed optimisation
				return null;
			xPrime[i] = solutions[i].getxVector();
			f += solutions[i].getfVal();
		}
		return xPrime;
	}

	public double getF() {
		return f;
	}

	public CMRTrial split(int negIndex, int posIndex, int index) {
		CMRTrial newTrial = new CMRTrial();
		int nvar = this.problems.length;

		newTrial.index = index;
		newTrial.xPrime = new double[nvar][];
		newTrial.problems = new MRProblem[nvar];
		newTrial.solutions = new MRSolution[nvar];
		newTrial.f = f;

		for (int i = 0; i < nvar; i++) {
			newTrial.problems[i] = (MRProblem) problems[i].clone();
			newTrial.problems[i].addConstraint(negIndex, posIndex, solutions[i] == null ? null : solutions[i].getxVector());
			// Check feasability of old solution
			double[] curXVector = solutions[i].getxVector();
			if (curXVector[negIndex] < curXVector[posIndex])
				newTrial.solutions[i] = solutions[i];
		}

		return newTrial;
	}

	public HashSet<SimpleLinearConstraint>[] getAdjs() {
		@SuppressWarnings("unchecked")
		HashSet<SimpleLinearConstraint>[] adjs = new HashSet[problems.length];
		for (int i = 0; i < problems.length; i++)
			adjs[i] = problems[i].getMatIneq();
		return adjs;
	}
}
