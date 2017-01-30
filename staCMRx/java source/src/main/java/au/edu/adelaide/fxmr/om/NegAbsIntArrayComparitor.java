package au.edu.adelaide.fxmr.om;

import java.util.Comparator;

/**
 * Allows comparison of magnitude of int arrays, largest to smallest
 * 
 * Assumes arrays passed in are of the same size!
 */
public class NegAbsIntArrayComparitor implements Comparator<int[]> {

	@Override
	public int compare(int[] o1, int[] o2) {
		int n = o1.length;
		for (int i = 0; i < n; i++) {
			int a1 = Math.abs(o1[i]);
			int a2 = Math.abs(o2[i]);
			if (a1 < a2)
				return 1;
			else if (a1 > a2)
				return -1;
		}
		return 0;
	}
}
