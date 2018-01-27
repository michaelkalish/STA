package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRSolver;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.ui.StatusFrame;
import au.edu.adelaide.fxmr.om.OMUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class CMRxSolver {
	/**
	 * If a call to MR fails, abort!
	 */
	private boolean easyFail;
	private boolean allowCyclic = true;
	private double tolerance = 0;
	protected double mrTolerance1 = 0;
	protected double mrTolerance2 = 0;
	private boolean onlyFeas = false;

	public CMRxSolver() {
	}

	public CMRSolution solve(CMRxProblem problem) {
		return this.solve(problem, false);
	}

	public CMRSolution solve(CMRxProblem problem, boolean showWindow) {
		SolverListener sl = null;
		if (showWindow)
			sl = new StatusFrame();
		return solve(problem, sl);
	}

	public CMRSolution solve(CMRxProblem problem, SolverListener sl) {
		return solve(problem, sl, false, -1);
	}

	public CMRSolution solve(CMRxProblem problem, SolverListener sl, boolean targetSet, double target) {
		long start = System.nanoTime();
		double finalTarget = targetSet ? target - problem.getfAddWeightedMeans() : 0;
		ArrayList<CMRIter> iter = new ArrayList<>();
		int nvar = problem.getNVar();
		double[][] means = problem.getMeans();
		DoubleMatrix2D[] weights = problem.getWeights();
		double fFloor = Double.NEGATIVE_INFINITY;
		double fBar = Double.POSITIVE_INFINITY;
		double[][] xBar = null;
		int[][] cv = problem.getCv();
		int[] tmpCVSet = new int[nvar];
		int ncond = problem.getNCond();
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		int newTrialID = 0;
		TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();
		int fBarReductions = 0;
		VisitedSet visited = new VisitedSet();

		if (sl != null)
			sl.updateStatus("CMRx checking feasability");

		TIntHashSet infeasZones = problem.getInfeasZones();
		int[] infeas = isFeasible3n(problem.getMeans(), infeasZones, tmpVolumes, tmpZoneNumbers);
		HashSet<SimpleLinearConstraint>[] adjBar = null;
		int cyclicAvoided = 0;

		if (infeas == null) {
			if (sl != null)
				sl.setFinished();
			return problem.createSolution(0, problem.getMeans(), iter,
					(double) (System.nanoTime() - start) / 1_000_000_000, null, 0, 0);
		}

		MRSolverAJOptimiser mrSolver = new MRSolverAJOptimiser();
		mrSolver.setTolerance(mrTolerance1, mrTolerance2);

		if (sl != null)
			sl.updateStatus("CMRx getting upper bound");

		// Get a reasonable upper bound
		CMRxTrial bestGreedy = getFeasibleX(problem, mrSolver, zDecodeCache);
		// CMRxTrial bestGreedy8 = getFeasible8(problem, mrSolver,
		// zDecodeCache);
		// if (Math.abs(bestGreedy8.getF() - bestGreedy.getF()) > 1e-6) {
		// System.out.println(bestGreedy8.getF() + ", " + bestGreedy.getF());
		// }

		if (bestGreedy != null) {
			xBar = bestGreedy.getxPrime();
			fBar = bestGreedy.getF();
			adjBar = bestGreedy.getAdjs();
			// For GC
			bestGreedy = null;
			fBarReductions++;

			if (onlyFeas || targetSet && fBar < finalTarget) {
				// upper bound is less than target!
				// System.out.println(fFloor + ", " + fBar + " <= " + target);
				if (sl != null)
					sl.setFinished();
				return problem.createSolution(fBar, xBar, null, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
						mrSolver.getCalls(), fBarReductions);
			}
		}

		if (sl != null)
			sl.updateStatus("Running CMRx");

		TreeSet<CMRxTrial> remaining = new TreeSet<>();

		// mrSolver.setFailQuietly(easyFail);
		mrSolver.setAllowCyclicProblems(allowCyclic);

		// Add first CMRxTrial
		remaining.add(new CMRxTrial(mrSolver, problem.getAdj(), nvar, weights, means));
		double tolerancem1 = 1.0 - tolerance;
		while (!remaining.isEmpty() && fFloor < fBar * tolerancem1) {
			CMRxTrial current = remaining.pollFirst();
			fFloor = current.getF();

			double upperFloor = remaining.isEmpty() ? fFloor : remaining.last().getF();
			if (sl != null && iter.size() % 100 == 0) {
				if (!sl.updateStatus(fFloor, fBar, upperFloor, remaining.size(), new int[] { iter.size() }, 0, fBarReductions,
						cyclicAvoided))
					break;
			}

			if (targetSet && (fBar < finalTarget || (fFloor >= finalTarget && fBar >= finalTarget))) {
				// upper bound < target or lower bound > target
				if (sl != null)
					sl.setFinished();
				return problem.createSolution(fBar, xBar, null, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
						mrSolver.getCalls(), fBarReductions);
			}

			iter.add(new CMRIter(fFloor, fBar, upperFloor, remaining.size()));

			if (fFloor < fBar * tolerancem1) {
				current.run();
				// current.run();
				double[][] xPrimes = current.getxPrime();
				if (xPrimes == null && easyFail) {
					if (sl != null)
						sl.setFinished();
					return null;
				} else if (xPrimes == null) {
					if (allowCyclic) {
						// This was caused by a failed MR call (the warning
						// would have been printed elsewhere)
						// System.err.println("Warning: ..." );
					} else {
						cyclicAvoided++;
					}
				} else {
					double fit = current.getF();
					infeas = isFeasible3n(xPrimes, infeasZones, tmpVolumes, tmpZoneNumbers);

					if (infeas == null) {
						// Solution is feasible
						if (fit < fBar) {
							// Solution is better than current best
							fBar = fit;
							fBarReductions++;
							xBar = xPrimes;
							adjBar = current.getAdjs();
							Iterator<CMRxTrial> iterR = remaining.descendingIterator();
							while (iterR.hasNext()) {
								if (iterR.next().getF() > fit)
									iterR.remove();
								else
									break;
							}
						}
					} else if (fit < fBar) {
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
								boolean useful = true;
								CMRxTrial newTrial = current.split(newTrialID++);

								for (int i = 0; i < nS; i++) {
									int k = tmpCVSet[i];
									if (covector[k] > 0)
										useful |= newTrial.addConstraint(k, posIndex, negIndex);
									if (covector[k] < 0)
										useful |= newTrial.addConstraint(k, negIndex, posIndex);
								}

								if (useful && !visited.contains(newTrial)) {
									remaining.add(newTrial);
									visited.add(newTrial);
								}
							} else {
							}
						}
					}
				}
			}
		}

		if (sl != null)
			sl.setFinished();

		// One last iter for plotting
		iter.add(new CMRIter(fBar, fBar, fBar, remaining.size()));

		return problem.createSolution(fBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls(), fBarReductions);
	}

	/**
	 * Choose the best getFeas to use
	 * 
	 * @param problem
	 * @param mrSolver
	 * @param zDecodeCache
	 * @return
	 */
	public static CMRxTrial getFeasibleX(CMRxProblem problem, MRSolverAJOptimiser solver, TIntObjectHashMap<int[]> zDecodeCache) {
		if (problem.getCv().length > 2)
			return getFeasible7(problem, solver, zDecodeCache);
		return getFeasible6(problem, solver, zDecodeCache);
	}

	/**
	 * modified version of isFeasible3n - given arbitrary ordered y vectors,
	 * will determine the largest non coupled monotonic pair.
	 * 
	 * @param y
	 *            2d array of means from StatsSTA
	 * @param infeasZones
	 *            zone numbers corresponding to infeasible points. For example,
	 *            [+1 -1] corresponds to zone number 2
	 * 
	 * @return null if feasible or the largest inversion if not (ZERO INDEXED!).
	 *         Zone number is also tacked on the end (Java implementation issue)
	 */
	public static int[] isFeasible3n(double[][] y, TIntHashSet infeasZones, DoubleMatrix2D volumes, int[][] zoneNumbers) {
		int nvar = y.length;
		int ny = y[0].length;
		for (int i = 0; i < ny; i++)
			Arrays.fill(zoneNumbers[i], 0);

		volumes.assign(1);

		for (int i = 0; i < nvar; i++) {
			double[] curY = y[i];
			long curPow = CMRUtil.POW_3[nvar - i - 1];

			for (int row = ny; --row >= 0;)
				for (int column = ny; --column >= row;) {
					double diff = curY[row] - curY[column];
					if (Math.abs(diff) > CMRSolver.ZERO_TOL) {
						volumes.setQuick(row, column, volumes.getQuick(row, column) * diff);
						zoneNumbers[row][column] += (int) (Math.signum(diff)) * curPow;
					}
				}
		}

		double maxVol = 0;
		int maxRow = -1;
		int maxColumn = -1;

		for (int row = ny; --row >= 0;) {
			for (int column = ny; --column >= row;) {
				int curZone = Math.abs(zoneNumbers[row][column]);
				if (infeasZones.contains(curZone)) {
					double curVol = Math.abs(volumes.getQuick(row, column));
					if (curVol > maxVol) {
						maxVol = curVol;
						maxRow = row;
						maxColumn = column;
					}
				}
			}
		}

		if (maxRow != -1)
			return new int[] { maxRow, maxColumn, zoneNumbers[maxRow][maxColumn] };
		return null;
	}

	/**
	 * Do a depth first search of a greedy feasible solution in order to obtain
	 * a reasonable estimate of the upper bound.
	 * 
	 * 
	 * @return
	 */
	// public static CMRxTrial getFeasible5(CMRxProblem problem, MRSolver
	// solver) {
	// int nvar = problem.getNVar();
	// int nvec = problem.getCv().length;
	// TIntHashSet infeasZones = problem.getInfeasZones();
	// int[][] cv = problem.getCv();
	// int ncond = problem.getNCond();
	// DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
	// int[][] tmpZoneNumbers = new int[ncond][ncond];
	//
	// double[][] xPrime = problem.getMeans();
	//
	// CMRxTrial[] curBests = new CMRxTrial[nvec];
	// for (int i = 0; i < nvec; i++)
	// curBests[i] = new CMRxTrial(solver, problem.getAdj(), nvar,
	// problem.getWeights(), problem.getMeans());
	//
	// int[] infeas = isFeasible3n(xPrime, infeasZones, tmpVolumes,
	// tmpZoneNumbers);
	// while (infeas != null) {
	// int negIndex = infeas[0];
	// int posIndex = infeas[1];
	//
	// CMRxTrial bestTrial = null;
	//
	// for (int ivec = 0; ivec < nvec; ivec++) {
	// for (int ivar = 0; ivar < nvar; ivar++) {
	// if (cv[ivec][ivar] >= 0) {
	// curBests[ivec].addConstraint(ivar, posIndex, negIndex);
	// } else if (cv[ivec][ivar] <= 0) {
	// curBests[ivec].addConstraint(ivar, negIndex, posIndex);
	// }
	// }
	// curBests[ivec].run();
	// if (bestTrial == null || curBests[ivec].getF() < bestTrial.getF())
	// bestTrial = curBests[ivec];
	// }
	//
	// if (bestTrial.getxPrime() == null)
	// // Hit circular constraints from every covector
	// return null;
	//
	// infeas = isFeasible3n(bestTrial.getxPrime(), infeasZones, tmpVolumes,
	// tmpZoneNumbers);
	// if (infeas == null) {
	// return bestTrial;
	// } else {
	// // Copy adjacency matries to all curBests
	// for (int i = 0; i < nvec; i++)
	// if (curBests[i] != bestTrial)
	// curBests[i].setConstraintsFrom(bestTrial);
	// }
	// }
	// // Something went wrong?!
	// return null;
	// }

	/**
	 * Do a depth first search of a greedy feasible solution in order to obtain
	 * a reasonable estimate of the upper bound.
	 * 
	 * @param zDecodeCache
	 * 
	 * @return
	 */
	public static CMRxTrial getFeasible6(CMRxProblem problem, MRSolver solver, TIntObjectHashMap<int[]> zDecodeCache) {
		int nvar = problem.getNVar();
		TIntHashSet infeasZones = problem.getInfeasZones();
		int[][] cv = problem.getCv();
		int ncond = problem.getNCond();
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		int[] tmpCVSet = new int[nvar];
		int trialIndex = 0;

		double[][] xPrime = problem.getMeans();

		CMRxTrial curBest = new CMRxTrial(solver, problem.getAdj(), nvar, problem.getWeights(), problem.getMeans());

		int[] infeas = isFeasible3n(xPrime, infeasZones, tmpVolumes, tmpZoneNumbers);
		while (infeas != null) {
			int negIndex = infeas[0];
			int posIndex = infeas[1];

			int zone = infeas[2];
			int[] signVector = zDecodeCache.get(zone);
			if (signVector == null) {
				signVector = OMUtil.zDecode(zone, nvar);
				zDecodeCache.put(zone, signVector);
			}

			CMRxTrial newBest = null;
			for (int[] covector : cv) {
				int nS = 0;
				for (int i = 0; i < nvar; i++)
					if (covector[i] != signVector[i] && signVector[i] != 0)
						tmpCVSet[nS++] = i;

				if (nS > 0) {
					CMRxTrial newTrial = curBest.split(trialIndex++);

					boolean isNew = false;

					for (int i = 0; i < nS; i++) {
						int k = tmpCVSet[i];
						if (covector[k] > 0)
							isNew |= newTrial.addConstraint(k, posIndex, negIndex);
						else if (covector[k] < 0)
							isNew |= newTrial.addConstraint(k, negIndex, posIndex);
					}

					if (isNew) {
						newTrial.run();
						if (newTrial.getxPrime() != null && (newBest == null || newTrial.getF() < newBest.getF()))
							newBest = newTrial;
					}
				}
			}

			if (newBest == null) {
				// This would mean we failed very early on...
				return null;
			}

			curBest = newBest;
			
			infeas = isFeasible3n(curBest.getxPrime(), infeasZones, tmpVolumes, tmpZoneNumbers);
		}

		return curBest;
	}

	/**
	 * Do a depth first search of a greedy feasible solution in order to obtain
	 * a reasonable estimate of the upper bound.
	 * 
	 * Difference with getFeas6 is that this will cache some MR calls.
	 * 
	 * @param zDecodeCache
	 * 
	 * @return
	 */
	public static CMRxTrial getFeasible7(CMRxProblem problem, MRSolver solver, TIntObjectHashMap<int[]> zDecodeCache) {
		int nvar = problem.getNVar();
		TIntHashSet infeasZones = problem.getInfeasZones();
		int[][] cv = problem.getCv();
		int ncond = problem.getNCond();
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		int[] tmpCVSet = new int[nvar];
		int trialIndex = 0;

		double[][] xPrime = problem.getMeans();

		CMRxTrial curBest = new CMRxTrial(solver, problem.getAdj(), nvar, problem.getWeights(), problem.getMeans());

		TrialCache trialCache = new TrialCache(nvar);

		int[] infeas = isFeasible3n(xPrime, infeasZones, tmpVolumes, tmpZoneNumbers);
		while (infeas != null) {
			int negIndex = infeas[0];
			int posIndex = infeas[1];
			trialCache.clear();

			int zone = infeas[2];
			int[] signVector = zDecodeCache.get(zone);
			if (signVector == null) {
				signVector = OMUtil.zDecode(zone, nvar);
				zDecodeCache.put(zone, signVector);
			}

			CMRxTrial newBest = null;
			for (int[] covector : cv) {
				int nS = 0;
				for (int i = 0; i < nvar; i++)
					if (covector[i] != signVector[i] && signVector[i] != 0)
						tmpCVSet[nS++] = i;

				if (nS > 0) {
					CMRxTrial newTrial = curBest.split(trialIndex++);

					boolean isNew = false;

					for (int i = 0; i < nS; i++) {
						int k = tmpCVSet[i];
						if (covector[k] > 0) {
							isNew |= newTrial.addConstraint(k, posIndex, negIndex);
						} else if (covector[k] < 0) {
							isNew |= newTrial.addConstraint(k, negIndex, posIndex);
						}
					}

					if (isNew) {
						newTrial.run(trialCache);
						if (newTrial.getxPrime() != null && (newBest == null || newTrial.getF() < newBest.getF()))
							newBest = newTrial;
					}
				}
			}

			if (newBest == null) {
				// This would mean we failed very early on...
				return null;
			}

			curBest = newBest;
			infeas = isFeasible3n(curBest.getxPrime(), infeasZones, tmpVolumes, tmpZoneNumbers);
		}

		return curBest;
	}

	/**
	 * generates a feasible solution by resolving each infeasible pair, one at a
	 * time based on Oleg's new algorithm from 12 Dec 2016 this simply resolves
	 * each infeasible pair and then does a final MR.
	 * 
	 * @param zDecodeCache
	 * 
	 * @return
	 */
	// public static CMRxTrial getFeasible8(CMRxProblem problem, MRSolver
	// solver, TIntObjectHashMap<int[]> zDecodeCache) {
	// int nvar = problem.getNVar();
	// int ncond = problem.getNCond();
	//
	// TIntHashSet infeasZones = problem.getInfeasZones();
	// int[][] cv = problem.getCv();
	// DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
	// int[][] tmpZoneNumbers = new int[ncond][ncond];
	// double[] diffs = new double[nvar];
	// double[] f = new double[nvar];
	// int nVec = cv.length;
	// double[] fit = new double[nVec];
	// double[][] xPrime = problem.getMeans();
	//
	// @SuppressWarnings("unchecked")
	// HashSet<SimpleLinearConstraint>[] curAdj = new HashSet[nvar];
	// for (int i = 0; i < nvar; i++)
	// curAdj[i] = new HashSet<>(problem.getAdj()[i]);
	//
	// isFeasible3n(xPrime, infeasZones, tmpVolumes, tmpZoneNumbers);
	//
	// DoubleMatrix2D[] w = problem.getWeights();
	//
	// TDoubleArrayList fits = new TDoubleArrayList();
	//
	// int[] infeasArray = infeasZones.toArray();
	// for (int inf : infeasArray)
	// for (int i = 0; i < ncond; i++)
	// for (int j = 0; j < ncond; j++)
	// if (Math.abs(tmpZoneNumbers[i][j]) == inf) {
	// for (int iVar = 0; iVar < nvar; iVar++) {
	// double xi = xPrime[iVar][i];
	// double xj = xPrime[iVar][j];
	//
	// diffs[iVar] = xi - xj;
	// double wi = w[iVar].get(i, i);
	// double wj = w[iVar].get(j, j);
	// double sumW = wi + wj;
	// double m = (wi * xi + wj * xj) / sumW;
	//
	// xi -= m;
	// xj -= m;
	// // cost of changing ivar
	// f[iVar] = wi * xi * xi + wj * xj * xj;
	// }
	// int mincvi = -1;
	// double minFit = Double.MAX_VALUE;
	// for (int cvi = 0; cvi < nVec; cvi++) {
	// double sum = 0;
	// int[] cvI = cv[cvi];
	// for (int iVar = 0; iVar < nvar; iVar++) {
	// if (signum(diffs[iVar]) != cvI[iVar]) {
	// sum += f[iVar];
	// }
	// }
	//
	// fit[cvi] = sum < 0 ? 0 : sum;
	// if (fit[cvi] < minFit) {
	// mincvi = cvi;
	// minFit = fit[cvi];
	// }
	// }
	//
	// fits.add(minFit);
	// for (int iVar = 0; iVar < nvar; iVar++) {
	// if (cv[mincvi][iVar] <= 0) {
	// curAdj[iVar].add(new SimpleLinearConstraint(i, j));
	// } else if (cv[mincvi][iVar] >= 0) {
	// curAdj[iVar].add(new SimpleLinearConstraint(j, i));
	// }
	// }
	// }
	//
	// CMRxTrial ret = new CMRxTrial(solver, curAdj, nvar, problem.getWeights(),
	// problem.getMeans());
	// ret.run();
	// return ret;
	// }

	// private static int signum(double d) {
	// if (d == 0.0)
	// return 0;
	// else if (d < 0)
	// return -1;
	// return 1;
	// }

	public void setEasyFail(boolean easyFail) {
		this.easyFail = easyFail;
	}

	public double getTolerance() {
		return tolerance;
	}

	public void setTolerance(double tolerance) {
		this.tolerance = tolerance;
	}

	public boolean isAllowCyclic() {
		return allowCyclic;
	}

	public void setAllowCyclic(boolean allowCyclic) {
		this.allowCyclic = allowCyclic;
	}

	public double getMrTolerance1() {
		return mrTolerance1;
	}

	public void setMrTolerance1(double mrTolerance1) {
		this.mrTolerance1 = mrTolerance1;
	}

	public double getMrTolerance2() {
		return mrTolerance2;
	}

	public void setMrTolerance2(double mrTolerance2) {
		this.mrTolerance2 = mrTolerance2;
	}

	public boolean isOnlyFeas() {
		return onlyFeas;
	}

	public void setOnlyFeas(boolean onlyFeas) {
		this.onlyFeas = onlyFeas;
	}
}
