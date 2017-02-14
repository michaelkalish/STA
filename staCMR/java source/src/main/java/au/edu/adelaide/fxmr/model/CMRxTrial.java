package au.edu.adelaide.fxmr.model;

import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import cern.colt.matrix.DoubleMatrix2D;

public class CMRxTrial implements Comparable<CMRxTrial> {
	private double[][] xPrime;
	private double f;
	private MRProblem[] problems;
	private MRSolution[] solutions;
	private int index;
	private MRSolver solver;
	private int solvedNum;
	private HashableAdjSet hAdjSet;

	public CMRxTrial(MRSolver solver, HashSet<SimpleLinearConstraint>[] adj, int nvar, DoubleMatrix2D[] weights, double[][] means) {
		this.solver = solver;
		xPrime = new double[nvar][];
		problems = new MRProblem[nvar];
		solutions = new MRSolution[nvar];
		f = Double.NEGATIVE_INFINITY;
		index = -1;
		for (int i = 0; i < nvar; i++)
			problems[i] = new MRProblem(means[i], weights[i], adj[i]);
		hAdjSet = new HashableAdjSet(adj);
	}

	private CMRxTrial() {
	}

	@Override
	public int compareTo(CMRxTrial o) {
		int compareF = Double.compare(f, o.f);
		if (compareF == 0)
			compareF = o.solvedNum - solvedNum;
		if (compareF == 0) {
			return Integer.compare(index, o.index);
		}
		return compareF;
	}

	public void run() {
		if (xPrime == null)
			// We failed previously, give up!
			return;

		int n = problems.length;

		f = 0;
		for (int i = 0; i < n; i++) {
			if (solutions[i] == null) {
				solutions[i] = solver.solve(problems[i]);
				solvedNum++;
			}
			if (solutions[i] == null) {
				solutions[i] = solver.solve(problems[i]);
				
				// Failed optimisation
				f = Double.POSITIVE_INFINITY;
				xPrime = null;
				return;
			}
			f += solutions[i].getfVal();
			xPrime[i] = solutions[i].getxVector();
		}
	}

	public double[][] getxPrime() {
		return xPrime;
	}

	public double getF() {
		return f;
	}

	public CMRxTrial split(int index) {
		CMRxTrial newTrial = new CMRxTrial();
		int nvar = this.problems.length;

		newTrial.solver = solver;
		newTrial.index = index;
		newTrial.xPrime = new double[nvar][];
		newTrial.problems = new MRProblem[nvar];
		newTrial.solutions = new MRSolution[nvar];
		newTrial.solvedNum = solvedNum;
		newTrial.f = f;

		for (int i = 0; i < nvar; i++) {
			newTrial.problems[i] = problems[i];
			newTrial.solutions[i] = solutions[i];
		}

		return newTrial;
	}

	public boolean addConstraint(int index, int posIndex, int negIndex) {
		double[] oldSolution = solutions[index] == null ? null  : solutions[index].getxVector();
		solutions[index] = null;
		solvedNum--;
		problems[index] = (MRProblem) problems[index].clone();
		boolean ret = problems[index].addConstraint(posIndex, negIndex, oldSolution);
		hAdjSet = new HashableAdjSet(getAdjs());
		return ret;
	}

	public void setConstraintsFrom(CMRxTrial bestTrial) {
		int n = problems.length;
		for (int i = 0; i < n; i++) {
			solutions[i] = null;
			problems[i].setMatIneq(new HashSet<SimpleLinearConstraint>(bestTrial.problems[i].getMatIneq()));
		}
		solvedNum = 0;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(solvedNum);
		sb.append("-");
		sb.append(f);
		// sb.append('\n');
		// for (MRProblem p : problems) {
		// sb.append(p.getMatA().toString());
		// sb.append("\n");
		// }

		return sb.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		for (MRProblem p : problems)
			result = prime * result + p.getMatIneq().hashCode();
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		CMRxTrial other = (CMRxTrial) obj;

		int n = problems.length;
		if (n != other.problems.length)
			return false;
		for (int i = 0; i < n; i++)
			if (!problems[i].getMatIneq().equals(other.problems[i].getMatIneq()))
				return false;

		return true;
	}

	public HashSet<SimpleLinearConstraint>[] getAdjs() {
		@SuppressWarnings("unchecked")
		HashSet<SimpleLinearConstraint>[] adjs = new HashSet[problems.length];
		for (int i = 0; i < problems.length; i++)
			adjs[i] = problems[i].getMatIneq();
		return adjs;
	}

	public int getSolvedNum() {
		return solvedNum;
	}

	public HashableAdjSet gethAdjSet() {
		return hAdjSet;
	}
}
