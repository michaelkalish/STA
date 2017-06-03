package au.edu.adelaide.fxmr.model.bin;

import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicBoolean;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.CMRxSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.om.OMUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class BinCMRxSolver {
	protected boolean onlyFeas = false;

	public BinCMRxSolver() {
	}

	public BinSolution[] solve(BinBaseProblem problem, AtomicBoolean running) {
		int n = problem.getModel().getnSubj();
		BinSolution[] solutions = new BinSolution[n];
		for (int i = 0; i < n; i++)
			solutions[i] = solve(new BinProblem(problem.getModel(), i, problem.getRangeSet()), running);

		return solutions;
	}

	public BinSolution[] solve(BinBaseProblem problem) {
		int n = problem.getModel().getnSubj();
		BinSolution[] solutions = new BinSolution[n];
		for (int i = 0; i < n; i++)
			solutions[i] = solve(new BinProblem(problem.getModel(), i, problem.getRangeSet()));

		return solutions;
	}

	public BinSolution[] solve(BinProblem[] problems) {
		int n = problems.length;
		BinSolution[] solutions = new BinSolution[n];
		for (int i = 0; i < n; i++)
			solutions[i] = solve(problems[i]);

		return solutions;
	}

	public BinSolution solve(BinProblem problem) {
		return solve(problem, null);
	}

	public BinSolution solve(BinProblem problem, AtomicBoolean running) {
		long start = System.nanoTime();

		BinModel model = problem.getModel();
		int nvar = model.getnVar();
		double gFloor = Double.NEGATIVE_INFINITY;
		double gBar = Double.POSITIVE_INFINITY;
		double fBar = Double.POSITIVE_INFINITY;
		double[][] xBar = null;

		int[][] cv = problem.getCv();
		TIntHashSet infeasZones = problem.getInfeasZones();

		BinElement[] data = model.getData()[problem.getSubjectIndex()];
		int ncond = data[0].getMeans().length;
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();
		int[] tmpCVSet = new int[nvar];
		int newTrialID = 0;

		HashSet<SimpleLinearConstraint>[] adjBar = null;
		MRSolver mrSolver = new MRSolverAJOptimiser();
		int iter = 0;
		// Get a reasonable upper bound
		BinTrial bestGreedy = getFeasible6BN(problem, mrSolver, zDecodeCache, tmpVolumes, tmpZoneNumbers);

		if (bestGreedy != null) {
			xBar = bestGreedy.getxPrime();
			fBar = bestGreedy.getF();
			gBar = bestGreedy.getG2();
			adjBar = bestGreedy.getAdjs();
			// For GC
			bestGreedy = null;

			if (onlyFeas) {
				// Only want best guess approximation
				return new BinSolution(fBar, gBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
						mrSolver.getCalls());
			}
		}

		BinVisitedSet visited = new BinVisitedSet();
		TreeSet<BinTrial> remaining = new TreeSet<>();
		// int saved = 0;

		// Add first CMRTrial
		remaining.add(new BinTrial(problem.getRangeSet(), nvar, data));

		while (!remaining.isEmpty() && gFloor < gBar && (running == null || running.get())) {
			BinTrial current = remaining.pollFirst();
			gFloor = current.getG2();
			iter++;

			// if (iter % 10000 == 0) {
			// int maxNV = 0;
			// int minRV = Integer.MAX_VALUE;
			// for (BinTrial r : remaining) {
			// int rv = r.countAdj();
			// if (rv < minRV)
			// minRV = rv;
			// }
			//
			// int nCullable = 0;
			// Iterator<BinTrial> iterRemove = visited.iterator();
			// while (iterRemove.hasNext()) {
			// int rv = iterRemove.next().countAdj();
			// if (rv > maxNV)
			// maxNV = rv;
			// if (rv < minRV) {
			// iterRemove.remove();
			// nCullable++;
			// }
			// }
			//
			// System.out.println(remaining.size() + ", " + visited.size() + ","
			// + saved + ", " + minRV + ", " + maxNV + ", " + nCullable);
			// }

			if (gFloor < gBar) {
				double[][] xPrimes = current.solveMRs(mrSolver);
				if (xPrimes == null) {
					// failed optimisation
					return null;
				}

				double gFit = current.getG2();

				int[] infeas = CMRxSolver.isFeasible3n(xPrimes, infeasZones, tmpVolumes, tmpZoneNumbers);

				if (gFit < gBar) {
					// Solution is better than current best
					if (infeas == null) {
						// Solution is feasible
						gBar = gFit;
						fBar = current.getF();
						xBar = xPrimes;
						adjBar = current.getAdjs();
						Iterator<BinTrial> iterR = remaining.descendingIterator();
						while (iterR.hasNext()) {
							if (iterR.next().getG2() > gFit)
								iterR.remove();
							else
								break;
						}
					} else {
						// Solution not feasible - make it so!
						int zone = infeas[2];
						int[] signVector = zDecodeCache.get(zone);
						if (signVector == null) {
							signVector = OMUtil.zDecode(zone, nvar);
							zDecodeCache.put(zone, signVector);
						}

						int negIndex = infeas[0];
						int posIndex = infeas[1];

						for (int[] covector : cv) {
							int nS = 0;
							for (int i = 0; i < nvar; i++)
								if (covector[i] != signVector[i] && signVector[i] != 0)
									tmpCVSet[nS++] = i;

							if (nS > 0) {
								BinTrial newTrial = null;
								for (int i = 0; i < nS; i++) {
									int k = tmpCVSet[i];
									if (covector[k] > 0)
										newTrial = current.split(k, posIndex, negIndex, newTrialID++);
									if (covector[k] < 0)
										newTrial = current.split(k, negIndex, posIndex, newTrialID++);
								}

								if (newTrial != null)
									if (!visited.contains(newTrial)) {
										remaining.add(newTrial);
										visited.add(newTrial);
									}
								// else {
								// saved++;
								// }
							} else {
							}
						}
					}
				}
			}
		}

		return new BinSolution(fBar, gBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls());
	}

	/**
	 * Do a depth first search of a greedy feasible solution in order to obtain
	 * a reasonable estimate of the upper bound.
	 * 
	 * @param zDecodeCache
	 * @param tmpZoneNumbers
	 * @param tmpVolumes2
	 * 
	 * @return
	 */
	public static BinTrial getFeasible6BN(BinProblem problem, MRSolver solver, TIntObjectHashMap<int[]> zDecodeCache, DoubleMatrix2D tmpVolumes,
			int[][] tmpZoneNumbers) {
		BinModel model = problem.getModel();
		int nvar = model.getnVar();
		TIntHashSet infeasZones = problem.getInfeasZones();
		int[][] cv = problem.getCv();
		BinElement[] data = model.getData()[problem.getSubjectIndex()];
		int[] tmpCVSet = new int[nvar];
		int trialIndex = 0;

		BinTrial curBest = new BinTrial(problem.getRangeSet(), nvar, data);
		curBest.solveMRs(solver);
		double[][] xPrime = curBest.getxPrime();
		int[] infeas = CMRxSolver.isFeasible3n(xPrime, infeasZones, tmpVolumes, tmpZoneNumbers);
		while (infeas != null) {
			int negIndex = infeas[0];
			int posIndex = infeas[1];

			int zone = infeas[2];
			int[] signVector = zDecodeCache.get(zone);
			if (signVector == null) {
				signVector = OMUtil.zDecode(zone, nvar);
				zDecodeCache.put(zone, signVector);
			}

			BinTrial newBest = null;
			for (int[] covector : cv) {
				int nS = 0;
				for (int i = 0; i < nvar; i++)
					if (covector[i] != signVector[i] && signVector[i] != 0)
						tmpCVSet[nS++] = i;

				if (nS > 0) {
					BinTrial newTrial = curBest.split(trialIndex++);

					for (int i = 0; i < nS; i++) {
						int k = tmpCVSet[i];
						if (covector[k] > 0)
							newTrial.addConstraint(k, posIndex, negIndex);
						else if (covector[k] < 0)
							newTrial.addConstraint(k, negIndex, posIndex);
					}

					newTrial.solveMRs(solver);
					if (newTrial.getxPrime() != null && (newBest == null || newTrial.getF() < newBest.getF()))
						newBest = newTrial;
				}
			}

			if (newBest == null) {
				// This would mean we failed very early on...
				return null;
			}

			curBest = newBest;
			infeas = CMRxSolver.isFeasible3n(curBest.getxPrime(), infeasZones, tmpVolumes, tmpZoneNumbers);
		}

		return curBest;
	}

	public boolean isOnlyFeas() {
		return onlyFeas;
	}

	public void setOnlyFeas(boolean onlyFeas) {
		this.onlyFeas = onlyFeas;
	}
}
