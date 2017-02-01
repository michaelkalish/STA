package au.edu.adelaide.fxmr.data;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import au.edu.adelaide.fxmr.model.Badness;
import au.edu.adelaide.fxmr.model.CombinedStatsSTA;
import au.edu.adelaide.fxmr.model.StatsSTA;
import au.edu.adelaide.fxmr.model.mr.MRProblem;
import au.edu.adelaide.fxmr.model.mr.MRSolution;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import gnu.trove.set.hash.TIntHashSet;

/**
 * Class that stores a complete picture of a general model.
 */
public class GeneralModel {
	private ArrayList<Subject> data = new ArrayList<>();
	private int[] depNVar;
	private int[] nBSCond;
	private Badness badness = Badness.NONE;

	public MRSolution[] calcStaMR(CombinedStatsSTA[] stats, int[][] rangeSet, double shrinkage) {
		return this.calcStaMR(stats, rangeSet, shrinkage, false, false);
	}

	public MRSolution[] calcStaMR(CombinedStatsSTA[] stats, int[][] rangeSet) {
		return this.calcStaMR(stats, rangeSet, -1, true, false);
	}

	private MRSolution[] calcStaMR(CombinedStatsSTA[] stats, int[][] rangeSet, double shrinkage, boolean autoShrink,
			boolean flagForceShrink) {
		if (stats == null)
			stats = calcStats(shrinkage, autoShrink);

		if (rangeSet == null) {
			int n = stats[0].getMeans().size();
			int[] rangeSet0 = new int[n];
			for (int i = 0; i < n; i++)
				rangeSet0[i] = i;
			rangeSet = new int[][] { rangeSet0 };
		}

		// Note default MR tolerance used.
		MRSolver solver = new MRSolverAJOptimiser();
		int n = stats.length;
		MRSolution[] solutions = new MRSolution[n];
		for (int i = 0; i < n; i++) {
			CombinedStatsSTA curStats = stats[i];
			MRProblem problem = new MRProblem(curStats.getMeans().toArray(), curStats.getWeights(), rangeSet);
			problem.forceSymetry();
			solutions[i] = solver.solve(problem);
		}

		return solutions;
	}

	public CombinedStatsSTA[] calcStats() {
		return calcStats(-1, true);
	}

	public CombinedStatsSTA[] calcStats(double shrinkage) {
		return calcStats(shrinkage, false);
	}

	private CombinedStatsSTA[] calcStats(double shrinkage, boolean autoShrink) {
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;
		CombinedStatsSTA[] ret = new CombinedStatsSTA[nVar];

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];
			StatsSTA[] ssta = new StatsSTA[nCond];

			for (int j = 0; j < nCond; j++) {
				int n = 0;
				int curCond = cond[j];

				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond)
						n++;

				double[][] within = new double[n][];
				int w = 0;
				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond)
						within[w++] = s.getCond();

				DoubleMatrix2D dataWithin = new DenseDoubleMatrix2D(within);

				if (autoShrink)
					ssta[j] = new StatsSTA(dataWithin);
				else
					ssta[j] = new StatsSTA(dataWithin, shrinkage);

				if (ssta[j].getBad().ordinal() > badness.ordinal())
					badness = ssta[j].getBad();
			}

			ret[i] = new CombinedStatsSTA(ssta, var[i]);

			double wCond = new SingularValueDecomposition(ret[i].getWeights()).cond();
			if (wCond > 1e16)
				// Very bad scaling!
				throw new IllegalArgumentException("Bad condition number for weights matrix " + wCond);
		}

		return ret;
	}

	/**
	 * Get the number of between subjects conditions.
	 * 
	 * @return
	 */
	public int[] getNBSConds() {
		if (nBSCond == null) {
			TIntHashSet set = new TIntHashSet();
			for (Subject s : data)
				set.add(s.getBetweenSubjectsCondition());
			nBSCond = set.toArray();
			Arrays.sort(nBSCond);
		}
		return nBSCond;
	}

	/**
	 * Get the number of dependent variables
	 * 
	 * @return
	 */
	public int[] getNDepVar() {
		if (depNVar == null) {
			TIntHashSet set = new TIntHashSet();
			for (Subject s : data)
				set.add(s.getDependentVariable());
			depNVar = set.toArray();
			Arrays.sort(depNVar);
		}
		return depNVar;
	}

	public Subject addData(int subjectNumber, int betweenSubjects, int dependentVariable, double[] cond) {
		if (cond == null || cond.length == 0)
			throw new InvalidParameterException("Condition must be defined");

		Subject subject = new Subject(subjectNumber, betweenSubjects, dependentVariable, cond);
		data.add(subject);
		depNVar = nBSCond = null;
		return subject;
	}

	public Subject addData(Subject subject) {
		data.add(subject);
		depNVar = nBSCond = null;
		return subject;
	}

	public ArrayList<Subject> getData() {
		return data;
	}

	public double[][] getBlock(int var, int cond) {
		int count = 0;
		for (Subject s : data)
			if (s.getDependentVariable() == var && s.getBetweenSubjectsCondition() == cond)
				count++;

		double[][] ret = new double[count][];
		count = 0;
		for (Subject s : data)
			if (s.getDependentVariable() == var && s.getBetweenSubjectsCondition() == cond)
				ret[count++] = s.getCond();

		return ret;
	}

	/**
	 * Make a new general model, bootstrapped from this model
	 * 
	 * @return
	 */
	public GeneralModel bootStrap(Random rand) {
		ArrayList<Subject> newData = new ArrayList<>(data.size());
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;

		ArrayList<Subject> tmpData = new ArrayList<>();

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];

			for (int j = 0; j < nCond; j++) {
				int n = 0;
				int curCond = cond[j];

				tmpData.clear();

				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond) {
						n++;
						tmpData.add(s);
					}

				for (int k = 0; k < n; k++) {
					newData.add(tmpData.get(rand.nextInt(n)));
				}
			}
		}

		GeneralModel ret = new GeneralModel();
		ret.data = newData;
		return ret;
	}

	/**
	 * "movedata" = original data - means + tmpSolutionX
	 * 
	 * @param means
	 *            pull means from cached calcStats called on this
	 * @param xStar
	 * @return
	 */
	public GeneralModel moveFrom(double[][] means, double[][] xStar) {
		ArrayList<Subject> newData = new ArrayList<>(data.size());
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];
			double[] curMeans = means[i];
			double[] curXStar = xStar[i];

			int adder = 0;
			for (int j = 0; j < nCond; j++) {
				int curCond = cond[j];

				int curNCond = -1;
				for (Subject s : data) {
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond) {
						curNCond = s.getCond().length;
						double[] newCond = s.getCond().clone();
						for (int c = 0; c < curNCond; c++)
							newCond[c] += -curMeans[adder + c] + curXStar[adder + c];

						Subject newSubject = new Subject(s.getSubjectNumber(), curCond, curVar, newCond);
						newData.add(newSubject);
					}
				}
				adder += curNCond;
			}
		}

		GeneralModel ret = new GeneralModel();
		ret.data = newData;
		return ret;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Subject s : data) {
			sb.append(s.toString());
			sb.append('\n');
		}

		return sb.toString();
	}

	public String toMATLABString() {
		StringBuilder sb = new StringBuilder();
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;

		sb.append("data = cell(" + nCond + "," + nVar + ");\n");

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];
			for (int j = 0; j < nCond; j++) {
				sb.append("data{" + (j + 1) + "," + (i + 1) + "} = [");
				int curCond = cond[j];
				for (Subject s : data) {
					if (s.getBetweenSubjectsCondition() == curCond && s.getDependentVariable() == curVar) {
						double[] c = s.getCond();
						sb.append(c[0]);
						for (int k = 1; k < c.length; k++) {
							sb.append(",");
							sb.append(c[k]);
						}
						sb.append('\n');
					}
				}
				sb.append("];\n");
			}
		}

		return sb.toString();
	}

	public int getMaxCond() {
		int max = 0;
		for (Subject s : data)
			if (s.getCond().length > max)
				max = s.getCond().length;
		return max;
	}

	public int getNBetween(int var) {
		int max = 0;
		for (Subject s : data)
			if (s.getDependentVariable() == var && s.getCond().length > max)
				max = s.getCond().length;
		return max;
	}

	/**
	 * Check that this general model has:
	 * 
	 * Some data in all within conditions
	 * 
	 * No columns contain all NaNs
	 * 
	 * All data in the within condition is the same length
	 * 
	 * @return
	 */
	public boolean sanityCheck() {
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];

			for (int j = 0; j < nCond; j++) {
				int n = 0;
				int curCond = cond[j];

				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond)
						n++;

				if (n == 0)
					// Should be at least one sample
					return false;

				double[][] within = new double[n][];
				int w = 0;
				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond)
						within[w++] = s.getCond();

				int expectedLength = within[0].length;
				boolean[] set = new boolean[expectedLength];
				for (double[] curW : within) {
					if (curW.length != expectedLength)
						// Should all be same length within same within
						return false;
					for (w = 0; w < expectedLength; w++)
						if (!set[w] && !Double.isNaN(curW[w]))
							set[w] = true;
				}

				for (w = 0; w < expectedLength; w++)
					if (!set[w])
						// A column of all NaNs
						return false;
			}

		}

		return true;
	}

	public int calcGlobalNCond() {
		int[] cond = getNBSConds();
		int[] var = getNDepVar();
		int nVar = var.length;
		int nCond = cond.length;
		int maxCondSum = 0;

		for (int i = 0; i < nVar; i++) {
			int curVar = var[i];
			int curCondSum = 0;
			for (int j = 0; j < nCond; j++) {
				int curCond = cond[j];
				int curCondLength = 0;

				for (Subject s : data)
					if (s.getDependentVariable() == curVar && s.getBetweenSubjectsCondition() == curCond) {
						if (s.getCond().length > curCondLength)
							curCondLength = s.getCond().length;
					}
				curCondSum += curCondLength;
			}

			if (curCondSum > maxCondSum)
				maxCondSum = curCondSum;
		}

		return maxCondSum;
	}
}
