package au.edu.adelaide.fxmr.model.bin;

import java.util.HashSet;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.HashableAdjSet;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;

public class BinTrial implements Comparable<BinTrial> {
	private double[][] xPrime;
	private double f;
	private double g2;
	private MRProblem[] problems;
	private MRSolution[] solutions;
	private double[] mleBNs;
	private int index;
	private BinElement[] binElements;
	private HashableAdjSet hAdjSet;

	public BinTrial(int[][] rangeSet, int nvar, BinElement[] subj) {
		xPrime = new double[nvar][];
		problems = new MRProblem[nvar];
		solutions = new MRSolution[nvar];
		mleBNs = new double[nvar];
		f = Double.NEGATIVE_INFINITY;
		g2 = Double.NEGATIVE_INFINITY;
		index = -1;
		binElements = subj;
		for (int i = 0; i < nvar; i++) {
			// n == weights
			problems[i] = new MRProblem(subj[i].getMeans(), subj[i].getN(), rangeSet);
		}
		hAdjSet = new HashableAdjSet(problems);
	}

	private BinTrial() {
	}

	@Override
	public int compareTo(BinTrial o) {
		int compareG = Double.compare(g2, o.g2);
		if (compareG == 0) {
			int compareIndex = Integer.compare(index, o.index);
			return compareIndex;
		}
		return compareG;
	}

	public double[][] solveMRs(MRSolver solver) {
		int n = problems.length;
		f = 0;
		g2 = 0;
		for (int i = 0; i < n; i++) {
			MRSolution solutioni = solutions[i];
			
			if (solutioni == null){
				solutions[i] = solutioni = solver.solve(problems[i]);
				if (solutioni == null)
					// Failed optimisation
					return null;
				clampZeroOne(solutioni.getxVector());
				mleBNs[i] =  mleBN(binElements[i], solutioni.getxVector());
			}
			
			xPrime[i] = solutioni.getxVector();
			g2 += mleBNs[i];
			f += solutioni.getfVal();
		}
		return xPrime;
	}

	public double getG2() {
		return g2;
	}

	public static double mleBN(BinElement binElement, double[] x) {
		int nn = binElement.getNSum();
		int n = x.length;
		double sumG2 = 0;

		for (int i = 0; i < n; i++) {
			int curN = binElement.getN()[i];
			double xv = x[i];
			if (xv <= 0)
				xv = 0.5 / nn;
			else if (xv >= 1)
				xv = (nn - 0.5) / nn;
			double aHit = xv * curN;
			double aMiss = (1.0 - xv) * curN;

			int h = binElement.getHits()[i];
			if (h != 0)
				sumG2 += 2.0 * h * Math.log(h / aHit);

			int m = binElement.getMisses()[i];
			if (m != 0)
				sumG2 += 2.0 * m * Math.log(m / aMiss);
		}

		return sumG2 < 0 ? 0 : sumG2;
	}

	public static void clampZeroOne(double[] xp) {
		int n = xp.length;
		for (int i = 0; i < n; i++)
			if (xp[i] < 0)
				xp[i] = 0;
			else if (xp[i] > 1)
				xp[i] = 1;
	}

	public double getF() {
		return f;
	}

	public String toString() {
		int n = problems.length;
		f = 0;
		StringBuilder sb = new StringBuilder();
		sb.append(f);
		sb.append(": [");
		for (int i = 0; i < n; i++) {
			sb.append(solutions[i] == null ? "N" : "S");
		}
		sb.append("]");
		return sb.toString();
	}

	public boolean addConstraint(int k, int posIndex, int negIndex) {
		// k variable
		double[] curXVector = solutions[k].getxVector();
		problems[k] = (MRProblem) problems[k].clone();
		if (!problems[k].addConstraint(posIndex, negIndex, curXVector))
			// The constraint already existed!
			return false;
		// Check feasability of current (old) solution
		if (curXVector[posIndex] >= curXVector[negIndex])
			solutions[k] = null;
		hAdjSet = new HashableAdjSet(problems);
		return true;
	}

	public BinTrial split(int k, int posIndex, int negIndex, int index) {
		BinTrial newTrial = split(index);

		if (!newTrial.addConstraint(k, posIndex, negIndex))
			return null;
		return newTrial;
	}

	public BinTrial split(int index) {
		BinTrial newTrial = new BinTrial();
		int nvar = this.problems.length;

		newTrial.index = index;
		newTrial.xPrime = new double[nvar][];
		newTrial.problems = new MRProblem[nvar];
		newTrial.solutions = new MRSolution[nvar];
		newTrial.mleBNs = new double[nvar];
		newTrial.binElements = binElements;
		newTrial.f = f;
		newTrial.g2 = g2;
		newTrial.hAdjSet = hAdjSet;

		for (int i = 0; i < nvar; i++) {
			newTrial.problems[i] = problems[i];
			newTrial.solutions[i] = solutions[i];
			newTrial.mleBNs[i] = mleBNs[i];
		}

		return newTrial;
	}

	public BinTrial split(int negIndex, int posIndex, int index) {
		BinTrial newTrial = new BinTrial();
		int nvar = this.problems.length;

		newTrial.index = index;
		newTrial.xPrime = new double[nvar][];
		newTrial.problems = new MRProblem[nvar];
		newTrial.solutions = new MRSolution[nvar];
		newTrial.mleBNs = new double[nvar];
		newTrial.binElements = binElements;
		newTrial.f = f;
		newTrial.g2 = g2;

		for (int i = 0; i < nvar; i++) {
			newTrial.mleBNs[i] = mleBNs[i];
			newTrial.problems[i] = (MRProblem) problems[i].clone();
			newTrial.problems[i].addConstraint(negIndex, posIndex, null);
			// Check feasability of old solution
			double[] curXVector = solutions[i].getxVector();
			if (curXVector[negIndex] < curXVector[posIndex])
				newTrial.solutions[i] = solutions[i];
		}
		hAdjSet = new HashableAdjSet(problems);

		return newTrial;
	}

	public HashSet<SimpleLinearConstraint>[] getAdjs() {
		@SuppressWarnings("unchecked")
		HashSet<SimpleLinearConstraint>[] adjs = new HashSet[problems.length];
		for (int i = 0; i < problems.length; i++)
			adjs[i] = problems[i].getMatIneq();
		return adjs;
	}

	public double[][] getxPrime() {
		return xPrime;
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
		BinTrial other = (BinTrial) obj;

		int n = problems.length;
		if (n != other.problems.length)
			return false;
		for (int i = 0; i < n; i++)
			if (!problems[i].getMatIneq().equals(other.problems[i].getMatIneq()))
				return false;

		return true;
	}

	public int countAdj() {
		int adjN = 0;
		for (MRProblem p : problems) {
			adjN += p.getMatIneq().size();
		}
		return adjN;
	}

	public HashableAdjSet gethAdjSet() {
		return hAdjSet;
	}

	public boolean check() {
		int n = problems.length;
		for (int i = 0; i < n; i++) {
			MRSolution s = solutions[i];
			MRProblem p = problems[i];

			double[] xv = s.getxVector();
			
			for (SimpleLinearConstraint ineq : p.getMatIneq()) {
				if (xv[ineq.getPosIndex()]  > xv[ineq.getNegIndex()])
					return false;
			}
		}

		return true;
	}
}
