package au.edu.adelaide.fxmr.joptimizer.functions;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * This implements JOptimiser's TwiceDifferentiableMultivariateRealFunction to
 * be an extremely lean simple linear constraint (yp - yn <= 0).
 *
 */
public class SimpleLinearConstraint implements TwiceDifferentiableMultivariateRealFunction, Comparable<SimpleLinearConstraint> {
	/**
	 * Dimension of the function argument.
	 */
	protected int dim = -1;

	private int posIndex;

	private int negIndex;

	public SimpleLinearConstraint(int len, int posIndex, int negIndex) {
		this.posIndex = posIndex;
		this.negIndex = negIndex;

		this.dim = len;
	}

	@Override
	public double value(DoubleMatrix1D x) {
		return x.getQuick(posIndex) - x.getQuick(negIndex);
	}

	public double value(double[] x) {
		return x[posIndex] - x[negIndex];
	}

	@Override
	public DoubleMatrix1D gradient(DoubleMatrix1D x) {
		return null;
	}

	@Override
	public void gradient(DoubleMatrix1D x, DoubleMatrix1D result) {
		// result.assign(grad);
	}

	public int getDim() {
		return this.dim;
	}

	public int getPosIndex() {
		return posIndex;
	}

	public int getNegIndex() {
		return negIndex;
	}

	public String toString() {
		return "x" + posIndex + " - x" + negIndex;
	}

	@Override
	public DoubleMatrix2D hessian() {
		return null;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + dim;
		result = prime * result + negIndex;
		result = prime * result + posIndex;
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
		SimpleLinearConstraint other = (SimpleLinearConstraint) obj;
		if (dim != other.dim)
			return false;
		if (negIndex != other.negIndex)
			return false;
		if (posIndex != other.posIndex)
			return false;
		return true;
	}

	@Override
	public int compareTo(SimpleLinearConstraint o) {
		if (posIndex == o.posIndex)
			return Integer.compare(negIndex, o.negIndex);
		return Integer.compare(posIndex, o.posIndex);
	}
}
