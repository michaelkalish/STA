package au.edu.adelaide.fxmr.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.ui.StatusFrame;
import au.edu.adelaide.fxmr.om.OMUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

/**
 * Parallel implementation of CMRx
 * 
 * 
 */
public class ParCMRxSolver {
	private double fBar;
	private NavigableSet<CMRxTrial> remaining;
	private VisitedSet visited;
	private ArrayList<CMRIter> iter;
	private AtomicInteger newTrialID;
	private AtomicInteger waiting;
	private int nThread;
	private TIntHashSet infeasZones;
	private int nvar;
	private int ncond;
	private double[][] xBar;
	private int[][] cv;
	private boolean running;
	private int[] nIterThread;
	private int collisions;
	private AtomicLong nextUpdate;
	private HashSet<SimpleLinearConstraint>[] adjBar;
	private long updateFreq = 250;
	private double tolerance;
	private int fBarReductions;
	private boolean allowCyclic = true;
	private int grabEnd = -1;
	private int cyclesAvoided;
	private boolean catastrophicFailure;
	protected double mrTolerance1 = 0;
	protected double mrTolerance2 = 0;

	public ParCMRxSolver() {
	}

	public CMRSolution solve(CMRxProblem problem) {
		return this.solve(problem, false, Runtime.getRuntime().availableProcessors());
	}

	public CMRSolution solve(CMRxProblem problem, boolean showWindow) {
		return this.solve(problem, showWindow, Runtime.getRuntime().availableProcessors());
	}

	public CMRSolution solve(CMRxProblem problem, boolean showWindow, int nThread) {
		return solve(problem, showWindow ? new StatusFrame() : null, nThread);
	}

	public CMRSolution solve(CMRxProblem problem, SolverListener sl, int nThread) {
		if (nThread == -1)
			nThread = Runtime.getRuntime().availableProcessors();
		this.nThread = nThread;
		nIterThread = new int[nThread];
		long start = System.nanoTime();
		nvar = problem.getNVar();
		double[][] means = problem.getMeans();
		DoubleMatrix2D[] weights = problem.getWeights();
		fBar = Double.POSITIVE_INFINITY;
		cv = problem.getCv();
		iter = new ArrayList<>();
		ncond = problem.getNCond();
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		collisions = 0;

		newTrialID = new AtomicInteger();
		waiting = new AtomicInteger();
		nextUpdate = new AtomicLong();

		infeasZones = problem.getInfeasZones();
		int[] infeas = CMRxSolver.isFeasible3n(problem.getMeans(), infeasZones, tmpVolumes, tmpZoneNumbers);

		if (infeas == null)
			return problem.createSolution(0, problem.getMeans(), iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar, 0, 0);

		remaining = new TreeSet<CMRxTrial>();
		visited = new VisitedSet();

		MRSolverAJOptimiser mrSolver = new MRSolverAJOptimiser();
		mrSolver.setTolerance(mrTolerance1, mrTolerance2);
		mrSolver.setAllowCyclicProblems(allowCyclic);

		// Add first CMRxTrial
		CMRxTrial first = new CMRxTrial(mrSolver, problem.getAdj(), nvar, weights, means);
		remaining.add(first);
		visited.add(first);
		first = null;

		if (sl != null)
			sl.updateStatus("Running Parallel CMRx, nThread = " + nThread);

		// Get a reasonable upper bound
		TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();
		CMRxTrial bestGreedy = CMRxSolver.getFeasibleX(problem, mrSolver, zDecodeCache);

		if (bestGreedy != null) {
			xBar = bestGreedy.getxPrime();
			fBar = bestGreedy.getF();
			adjBar = bestGreedy.getAdjs();
			fBarReductions++;
			// For GC
			bestGreedy = null;
		}

		running = true;
		ExecutorService executor = Executors.newFixedThreadPool(nThread < 1 ? Runtime.getRuntime().availableProcessors() : nThread);
		ArrayList<CMRxTrialSolver> pool = new ArrayList<>();
		for (int i = 0; i < nThread; i++) {
			pool.add(new CMRxTrialSolver(i, sl));
		}
		try {
			executor.invokeAll(pool);
			executor.shutdown();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		if (sl != null)
			sl.setFinished();

		if (catastrophicFailure)
			// Return null since we can't trust the output
			return null;

		// Add one last item to iter to make it pretty
		iter.add(new CMRIter(fBar, fBar, fBar, remaining.size()));

		return problem.createSolution(fBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar, mrSolver.getCalls(),
				fBarReductions);
	}

	public double getTolerance() {
		return tolerance;
	}

	public void setTolerance(double tolerance) {
		this.tolerance = tolerance;
	}

	/**
	 * Parallel threads, grabbing items from the queue, solving and adding more
	 * 
	 * 
	 *
	 */
	private class CMRxTrialSolver implements Callable<Boolean> {
		private int threadNum;
		private SolverListener sl;

		public CMRxTrialSolver(int threadNum, SolverListener sl) {
			this.threadNum = threadNum;
			this.sl = sl;
		}

		@Override
		public Boolean call() {
			try {
				double tolerancem1 = 1.0 - tolerance;
				int[] tmpCVSet = new int[nvar];
				double fFloor = Double.NEGATIVE_INFINITY;
				double upperFloor = Double.NEGATIVE_INFINITY;
				DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
				int[][] tmpZoneNumbers = new int[ncond][ncond];
				int[] infeas;
				int[][] cv = ParCMRxSolver.this.cv;
				TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();

				ArrayList<CMRxTrial> tmpTrials = new ArrayList<>(cv.length);
				while (running) {
					CMRxTrial current = null;
					current = null;

					waiting.incrementAndGet();
					while (current == null) {
						synchronized (remaining) {
							current = grabEnd != -1 && remaining.size() > grabEnd ? remaining.pollLast() : remaining.pollFirst();
							if (current != null) {
								waiting.decrementAndGet();
								break;
							} else if (waiting.intValue() == nThread || !running) {
								running = false;
								return true;
							}
						}
						try {
							Thread.sleep(1);
						} catch (InterruptedException e) {
							// meh
						}
					}

					fFloor = current.getF();
					nIterThread[threadNum]++;
					synchronized (iter) {
						synchronized (remaining) {
							if (!remaining.isEmpty())
								upperFloor = remaining.last().getF();
							else
								upperFloor = fFloor;// This is the one we just
							// grabbed!
						}
						iter.add(new CMRIter(fFloor, fBar, upperFloor, remaining.size()));
					}

					if (fFloor < fBar * tolerancem1) {
						current.run();
						double[][] xPrimes = current.getxPrime();
						if (xPrimes == null) {
							// failed MR optimisation
							if (allowCyclic) {
								// Bad - there will be a warning message already
								catastrophicFailure = true;
								running = false;
							} else {
								synchronized (ParCMRxSolver.this) {
									cyclesAvoided++;
								}
							}
						} else {
							double fit = current.getF();
							infeas = CMRxSolver.isFeasible3n(xPrimes, infeasZones, tmpVolumes, tmpZoneNumbers);

							if (infeas == null) {
								// Solution is feasible
								synchronized (ParCMRxSolver.this) {
									if (fit < fBar) {
										// Solution is better than current best
										fBar = fit;
										xBar = xPrimes;
										adjBar = current.getAdjs();
										fBarReductions++;
										synchronized (remaining) {
											Iterator<CMRxTrial> iter = remaining.descendingIterator();
											while (iter.hasNext()) {
												if (iter.next().getF() > fit)
													iter.remove();
												else
													break;
											}
										}
									}
								}
							} else if (fit < fBar * tolerancem1) {
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
										CMRxTrial newTrial = current.split(newTrialID.getAndIncrement());

										boolean useful = false;
										for (int i = 0; i < nS; i++) {
											int k = tmpCVSet[i];
											if (covector[k] > 0)
												useful |= newTrial.addConstraint(k, posIndex, negIndex);
											if (covector[k] < 0)
												useful |= newTrial.addConstraint(k, negIndex, posIndex);
										}

										if (useful) {
											synchronized (visited) {
												if (visited.contains(newTrial))
													collisions++;
												else
													tmpTrials.add(newTrial);
											}
										}
									}
								}
								synchronized (visited) {
									visited.addAll(tmpTrials);
								}
								synchronized (remaining) {
									remaining.addAll(tmpTrials);
								}
								tmpTrials.clear();
							}
						}
					}

					if (sl != null && System.currentTimeMillis() > nextUpdate.get()) {
						nextUpdate.set(System.currentTimeMillis() + updateFreq);

						running = sl.updateStatus(fFloor, fBar, upperFloor, remaining.size(), nIterThread, collisions, fBarReductions,
								cyclesAvoided);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
			return false;
		}
	}

	public void setAllowCyclic(boolean allowCyclic) {
		this.allowCyclic = allowCyclic;
	}

	public int getGrabEnd() {
		return grabEnd;
	}

	public void setGrabEnd(int grabEnd) {
		this.grabEnd = grabEnd;
	}

	public long getUpdateFreq() {
		return updateFreq;
	}

	public void setUpdateFreq(long updateFreq) {
		this.updateFreq = updateFreq;
	}

	public boolean isCatastrophicFailure() {
		return catastrophicFailure;
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
}
