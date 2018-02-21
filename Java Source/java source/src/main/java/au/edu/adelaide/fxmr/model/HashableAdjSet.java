package au.edu.adelaide.fxmr.model;

import java.util.Arrays;
import java.util.HashSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRProblem;

public class HashableAdjSet {
	private SimpleLinearConstraint[][] adjs;
	int hashCode;

	public HashableAdjSet(HashSet<SimpleLinearConstraint>[] adjs) {
		final int prime = 31;
		hashCode = 1;

		this.adjs = new SimpleLinearConstraint[adjs.length][];
		int i = 0;
		for (HashSet<SimpleLinearConstraint> a : adjs) {
			hashCode = prime * hashCode + a.hashCode();
			this.adjs[i] = a.toArray(new SimpleLinearConstraint[a.size()]);
			Arrays.sort(this.adjs[i++]);
		}
	}

	public HashableAdjSet(MRProblem[] problems) {
		final int prime = 31;
		hashCode = 1;

		adjs = new SimpleLinearConstraint[problems.length][];
		int i = 0;
		for (MRProblem p : problems) {
			HashSet<SimpleLinearConstraint> a = p.getMatIneq();
			hashCode = prime * hashCode + a.hashCode();
			adjs[i] = a.toArray(new SimpleLinearConstraint[a.size()]);
			Arrays.sort(adjs[i++]);
		}		
	}

	@Override
	public int hashCode() {
		return hashCode;
	}

	@Override
	public boolean equals(Object obj) {
		HashableAdjSet other = (HashableAdjSet) obj;

		int n = adjs.length;
		if (n != other.adjs.length)
			return false;
		
		for (int i = 0; i < n; i++) {
			int n2 = adjs[i].length;
			if (other.adjs[i].length != n2)
				return false;
			for (int j = 0; j < n2; j++)
				if (!adjs[i][j].equals(other.adjs[i][j]))
					return false;
		}
		return true;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (SimpleLinearConstraint[] a : adjs) {
			sb.append(Arrays.toString(a));
			sb.append(",");
		}

		return sb.toString();
	}
}