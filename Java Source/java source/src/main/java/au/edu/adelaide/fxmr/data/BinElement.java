package au.edu.adelaide.fxmr.data;

import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;

public class BinElement {
	private int[] hits;
	private int[] misses;
	private int[] n;
	private double[] means;

	public double[] getMeans() {
		return means;
	}

	private BinElement() {
		// Do nothing - private
	}

	/**
	 * Initialise BinElement with zero valued arrays of nCond size.
	 * 
	 * @param nCond
	 */
	public BinElement(int nCond) {
		hits = new int[nCond];
		misses = new int[nCond];
		n = new int[nCond];
	}

	/**
	 * Initialise BinElement with given hits and misses data.
	 * 
	 * We assume lengths are the same (nCond = length of both arrays)
	 * 
	 * stats are calculated from input.
	 * 
	 * @param hits
	 * @param misses
	 */
	public BinElement(int[] hits, int[] misses) {
		this.hits = hits;
		this.misses = misses;

		calcStats();
	}

	/**
	 * Convenience method
	 * 
	 * @param values
	 */
	public BinElement(int[][] values) {
		if (values.length == 2) {
			hits = values[0];
			misses = values[1];
		} else if (values[0].length == 2) {
			hits = new int[values.length];
			misses = new int[values.length];
			for (int i = 0; i < values.length; i++) {
				hits[i] = values[i][0];
				misses[i] = values[i][1];
			}
		} else {
			// Nope
			throw new IllegalArgumentException("Matrix must be 2xN or Nx2");
		}

		calcStats();
	}

	/**
	 * Calculate stats on provided hit and miss data
	 */
	public void calcStats() {
		int nCond = hits.length;
		n = new int[nCond];
		means = new double[nCond];
		for (int i = 0; i < nCond; i++) {
			n[i] = hits[i] + misses[i];
			means[i] = (double) hits[i] / n[i];
		}
	}

	public int[] getHits() {
		return hits;
	}

	public int[] getMisses() {
		return misses;
	}

	public int[] getN() {
		return n;
	}

	public int getNSum() {
		int sum = 0;
		for (int nn : n)
			sum += nn;
		return sum;
	}

	public BinElement resample(RandomEngine binRand) {
		BinElement newElement = new BinElement();
		newElement.n = n;
		int nCond = n.length;

		newElement.hits = new int[nCond];
		newElement.misses = new int[nCond];
		newElement.means = new double[nCond];
		
		Binomial bn = new Binomial(1, 0.5, binRand);
		for (int i = 0; i < nCond; i++) {
			int ni = n[i];
			double p = means[i];
			int b;
			if (p == 1)
				b = ni;
			else if (p == 0)
				b = 0;
			else
				b = bn.nextInt(ni, p);
			newElement.hits[i] = b;
			newElement.misses[i] = ni - b;
			newElement.means[i] = (double) b / ni;
		}

		return newElement;
	}

	public BinElement resample(double[] altMeans,RandomEngine binRand) {
		BinElement newElement = new BinElement();
		newElement.n = n;
		int nCond = n.length;

		newElement.hits = new int[nCond];
		newElement.misses = new int[nCond];
		newElement.means = new double[nCond];

		Binomial bn = new Binomial(1, 0.5, binRand);
		
		for (int i = 0; i < nCond; i++) {
			int ni = n[i];
			double p = altMeans[i];
			int b;
			if (p == 1)
				b = ni;
			else if (p == 0)
				b = 0;
			else
				b = bn.nextInt(ni, p);
			newElement.hits[i] = b;
			newElement.misses[i] = ni - b;
			newElement.means[i] = (double) b / ni;
		}

		return newElement;
	}
}
