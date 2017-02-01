package au.edu.adelaide.fxmr.model.mr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

/**
 * Defines the input for the MR problem
 */
public class MRProblem implements Cloneable {
	/**
	 * Input vector y
	 */
	private double[] y;
	/**
	 * number of elements in y
	 */
	private int n;

	/**
	 * Weights matrix
	 */
	private DoubleMatrix2D weights;

	/**
	 * Inequality Constraints
	 */
	private HashSet<SimpleLinearConstraint> matIneq = new HashSet<>();

	public void setMatIneq(HashSet<SimpleLinearConstraint> matIneq) {
		this.matIneq = matIneq;
	}

	public HashSet<SimpleLinearConstraint> getMatIneq() {
		return matIneq;
	}

	/**
	 * 
	 * @return
	 */

	public DoubleMatrix2D getWeights() {
		return weights;
	}

	/**
	 * Contructor
	 * 
	 * @param y
	 *            TODO
	 * @param weights
	 *            TODO
	 * @param rangeSet
	 *            TODO (make comment about zero indexing?)
	 */
	public MRProblem(double[] y, double[][] weights, int[][] rangeSet) {
		this.y = y;
		this.n = y.length;
		makeA(rangeSet);

		if (weights == null) {
			this.weights = DoubleFactory2D.sparse.identity(n);
		} else {
			this.weights = new SparseDoubleMatrix2D(weights);
		}
	}

	public MRProblem(double[] y, int[] weights, int[][] rangeSet) {
		this.y = y;
		this.n = y.length;
		makeA(rangeSet);

		if (weights == null) {
			this.weights = DoubleFactory2D.sparse.identity(n);
		} else {
			this.weights = DoubleFactory2D.sparse.identity(n);
			for (int i = 0; i < n; i++)
				this.weights.set(i, i, weights[i]);
		}
	}

	public MRProblem(double[] y, DoubleMatrix2D weights, int[][] rangeSet) {
		this.y = y;
		this.n = y.length;
		makeA(rangeSet);

		if (weights == null) {
			this.weights = DoubleFactory2D.sparse.identity(n);
		} else {
			this.weights = weights;
		}
	}

	public MRProblem(double[] y, int[][] rangeSet) {
		this.y = y;
		this.n = y.length;
		makeA(rangeSet);
		this.weights = DoubleFactory2D.sparse.identity(n);
	}

	public MRProblem(double[] y, DoubleMatrix2D weights, HashSet<SimpleLinearConstraint> matA) {
		this.y = y;
		this.n = y.length;
		this.matIneq = matA;

		if (weights == null) {
			this.weights = DoubleFactory2D.sparse.identity(n);
		} else {
			this.weights = weights;
		}
	}

	/**
	 * Convert a set of ranges into an adjacency matrix
	 * 
	 * @param rangeSet
	 * @return
	 */
	void makeA(int[][] rangeSet) {
		if (rangeSet == null)
			return;

		for (int[] range : rangeSet) {
			int rl = range.length;
			if (rl >= 2) {
				// New code, fewer constraints
				int rlm1 = rl - 1;
				for (int i = 0; i < rlm1; i++)
					matIneq.add(new SimpleLinearConstraint(n, range[i], range[i + 1]));
			}
		}
	}

	public double[] getY() {
		return y;
	}

	public int getN() {
		return n;
	}

	public Object clone() {
		try {
			MRProblem newProblem = (MRProblem) super.clone();
			newProblem.matIneq = new HashSet<>();
			newProblem.matIneq.addAll(matIneq);
			return newProblem;
		} catch (CloneNotSupportedException e) {
			return null;
		}
	}

	public boolean addConstraint(int posIndex, int negIndex) {
		return matIneq.add(new SimpleLinearConstraint(n, posIndex, negIndex));
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		// sb.append(weights.toString());
		// if (findFeasibleStart(1) == null)
		// sb.append("CYCLIC!");
		// sb.append(Arrays.toString(y));
		for (SimpleLinearConstraint c : matIneq) {
			sb.append("\n");
			sb.append(c);
		}

		return sb.toString();
	}

	public boolean testFeas(double[] initial) {
		for (SimpleLinearConstraint c : matIneq)
			if (c.value(initial) > 0)
				return false;

		return true;
	}

	/**
	 * Treat matIneq as a Directed Graph and create a feasible initial starting
	 * point.
	 * 
	 * IF there is a cycle in the graph, we need to add constraints to MatEq and
	 * remove them from matIneq (i.e., this method can change the underlying
	 * structure!)
	 * 
	 * @return
	 */
	public double[] findFeasibleStart(double distanceScale, HashSet<SimpleLinearConstraint> ineq,
			HashSet<SimpleLinearConstraint> eq) {
		double sum = 0;
		boolean[] used = new boolean[n];

		ineq.addAll(matIneq);
		for (SimpleLinearConstraint c : ineq)
			sum += Math.abs(c.value(y));

		double meand = sum / n;
		meand *= distanceScale;

		// Find all cycles
		boolean[] cycle = new boolean[n];
		ArrayList<ArrayList<Integer>> cycleSets = null;
		for (int i = 0; i < n; i++) {
			Arrays.fill(cycle, false);

			if (findCycle(i, i, cycle, used, ineq, false)) {
				if (cycleSets == null)
					// Cycles should be rare
					cycleSets = new ArrayList<>();
				// One big equality loop
				ArrayList<Integer> cycleIndices = new ArrayList<>();
				cycleSets.add(cycleIndices);
				for (int j = 0; j < n; j++)
					if (cycle[j]) {
						cycleIndices.add(j);
						used[j] = true;
					}

				// Create our own equality constraints (without a cycle)
				int nCycle = cycleIndices.size();
				for (int j = 1; j < nCycle; j++)
					eq.add(new SimpleLinearConstraint(n, cycleIndices.get(j), cycleIndices.get(j - 1)));

				// Set those constraints to equalities
				Iterator<SimpleLinearConstraint> iter = ineq.iterator();
				while (iter.hasNext()) {
					SimpleLinearConstraint c = iter.next();
					if (cycle[c.getPosIndex()] && cycle[c.getNegIndex()])
						iter.remove();
				}
			}
		}
		cycle = null;

		if (cycleSets == null) {
			// Simplest case
			return findAcyclicStart(meand, ineq, used, n);
		}

		// Now globalise cycle
		int[] cycleIndex = new int[n];
		Arrays.fill(cycleIndex, -1);
		int remainingStart = n;
		for (int si = 0; si < cycleSets.size(); si++) {
			ArrayList<Integer> ccs = cycleSets.get(si);
			used[ccs.get(0)] = false;
			remainingStart -= ccs.size() - 1;
			for (int i : ccs) {
				cycleIndex[i] = si;
			}
		}

		ArrayList<SimpleLinearConstraint> newIneqs = new ArrayList<>();

		Iterator<SimpleLinearConstraint> iter = ineq.iterator();
		while (iter.hasNext()) {
			SimpleLinearConstraint c = iter.next();

			int newPos = c.getPosIndex();
			if (cycleIndex[newPos] != -1)
				newPos = cycleSets.get(cycleIndex[newPos]).get(0);

			int newNeg = c.getNegIndex();
			if (cycleIndex[newNeg] != -1)
				newNeg = cycleSets.get(cycleIndex[newNeg]).get(0);

			if (newPos != c.getPosIndex() || newNeg != c.getNegIndex()) {
				newIneqs.add(new SimpleLinearConstraint(c.getDim(), newPos, newNeg));
				iter.remove();
			}
		}

		ineq.addAll(newIneqs);
		newIneqs = null;

		double[] start = findAcyclicStart(meand, ineq, used, remainingStart);

		for (int si = 0; si < cycleSets.size(); si++) {
			ArrayList<Integer> ccs = cycleSets.get(si);
			double base = start[ccs.get(0)];
			for (int i : ccs) {
				start[i] = base;
			}
		}

		return start;
	}

	private double[] findAcyclicStart(double meand, HashSet<SimpleLinearConstraint> ineq, boolean[] used, int remainingStart) {
		int remaining = remainingStart;
		double[] start = new double[n];

		while (remaining > 0) {
			// Look at each variable - at least one HAS to be a tail
			for (int i = n - 1; i >= 0; i--) {
				if (!used[i]) {
					boolean curHasDependent = false;
					double bestParentMin = Double.POSITIVE_INFINITY;

					for (SimpleLinearConstraint c : ineq) {
						if (c.getPosIndex() == i) {
							if (!used[c.getNegIndex()]) {
								curHasDependent = true;
								break;
							} else if (start[c.getNegIndex()] < bestParentMin) {
								bestParentMin = start[c.getNegIndex()];
							}
						}
					}

					if (!curHasDependent) {
						if (bestParentMin == Double.POSITIVE_INFINITY) {
							start[i] = y[i];
						} else {
							start[i] = bestParentMin - meand;
						}
						used[i] = true;
						remaining--;
					}
				}
			}
		}

		return start;
	}

	/**
	 * Find a cycle starting from start, recursively
	 * 
	 * TODO: There are better algorithms for this but it's so rarely called...
	 * 
	 * @param start
	 * @param cycle
	 * @param lynchpins
	 * @param eq
	 * @return
	 */
	private boolean findCycle(int initial, int current, boolean[] cycle, boolean[] used,
			HashSet<SimpleLinearConstraint> ineq, boolean subSearch) {
		if (cycle[initial] && initial == current)
			return true;

		boolean foundCycle = false;
		if (!cycle[current]) {
			cycle[current] = true;
			for (SimpleLinearConstraint c : ineq) {
				if (c.getPosIndex() == current) {
					int next = c.getNegIndex();
					if (!used[next])
						foundCycle |= findCycle(initial, next, cycle, used, ineq, subSearch);
				}
			}

			if (!foundCycle)
				cycle[current] = false;
		}

		if (!subSearch && foundCycle && initial == current) {
			boolean[] subCycle = new boolean[n];

			// Need to check for any other overlapping cycles
			for (int i = 0; i < n; i++) {
				if (!used[i] && !cycle[i]) {
					Arrays.fill(subCycle, false);
					if (findCycle(i, i, subCycle, used, ineq, true)) {
						boolean crossover = false;
						for (int j = 0; j < n; j++) {
							if (cycle[j] && subCycle[j]) {
								crossover = true;
								break;
							}
						}
						if (crossover) {
							for (int j = 0; j < n; j++)
								cycle[j] |= subCycle[j];
						}
					}
				}
			}
		}

		return foundCycle;
	}

	public String toMatlabString() {
		StringBuilder sb = new StringBuilder();
		sb.append("W = [");
		for (int row = 0; row < weights.rows(); row++) {
			sb.append(weights.get(row, 0));
			for (int col = 1; col < weights.columns(); col++) {
				sb.append(',');
				sb.append(weights.get(row, col));
			}
			if (row != weights.rows() - 1)
				sb.append(';');
		}
		sb.append("];\ny=");
		sb.append(Arrays.toString(y));
		sb.append(";\nE={");
		for (SimpleLinearConstraint c : matIneq) {
			sb.append("[");
			sb.append(c.getPosIndex() + 1);
			sb.append(",");
			sb.append(c.getNegIndex() + 1);
			sb.append("] ");
		}
		sb.append("};\n[x,fit]=MR(y,W,E)");
		return sb.toString();
	}

	public boolean constraintsAlreadySatisfied() {
		return testFeas(y);
	}

	public void forceSymetry() {
		DoubleMatrix2D check = weights.copy();
		check.assign(weights.viewDice(), Functions.minus);
		double normCheck = Algebra.ZERO.normInfinity(check);
		if (normCheck > 1e-14) {
			// Force symmetry
			check.assign(weights);
			check.assign(weights.viewDice(), Functions.plus);
			check.assign(Functions.div(2));
			weights.assign(check);
		}
	}
}
