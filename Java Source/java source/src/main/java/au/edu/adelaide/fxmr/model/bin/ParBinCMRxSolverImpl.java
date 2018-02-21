package au.edu.adelaide.fxmr.model.bin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import au.edu.adelaide.fxmr.data.BinElement;
import au.edu.adelaide.fxmr.data.BinModel;
import au.edu.adelaide.fxmr.joptimizer.functions.SimpleLinearConstraint;
import au.edu.adelaide.fxmr.model.CMRxSolver;
import au.edu.adelaide.fxmr.model.CMRxTrial;
import au.edu.adelaide.fxmr.model.VisitedSet;
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
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
public class ParBinCMRxSolverImpl {
	private double fBar;
	private NavigableSet<BinTrial> remaining;
	private BinVisitedSet visited;
	private AtomicInteger newTrialID;
	private AtomicInteger waiting;
	private int nThread;
	private TIntHashSet infeasZones;
	private int nvar;
	private int ncond;
	private double[][] xBar;
	private int[][] cv;
	private boolean running;
	private int collisions;
	private HashSet<SimpleLinearConstraint>[] adjBar;
	private int grabEnd = -1;
	private boolean catastrophicFailure;
	private double gBar;
	private int[] nIterThread;
	private MRSolverAJOptimiser mrSolver;

	public ParBinCMRxSolverImpl() {
	}

	public BinSolution solve(BinProblem problem, int nThread) {
		if (nThread < 1)
			nThread = Runtime.getRuntime().availableProcessors();
		this.nThread = nThread;
		nIterThread = new int[nThread];
		long start = System.nanoTime();

		BinModel model = problem.getModel();
		nvar = model.getnVar();
		gBar = Double.POSITIVE_INFINITY;
		fBar = Double.POSITIVE_INFINITY;
		xBar = null;
		newTrialID = new AtomicInteger();
		waiting = new AtomicInteger();

		BinElement[] data = model.getData()[problem.getSubjectIndex()];
		ncond = data[0].getMeans().length;
		int[][] tmpZoneNumbers = new int[ncond][ncond];
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
		TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();

		cv = problem.getCv();
		infeasZones = problem.getInfeasZones();

		mrSolver = new MRSolverAJOptimiser();

		// Get a reasonable upper bound
		BinTrial bestGreedy = BinCMRxSolver.getFeasible6BN(problem, mrSolver, zDecodeCache, tmpVolumes, tmpZoneNumbers);

		if (bestGreedy != null) {
			xBar = bestGreedy.getxPrime();
			fBar = bestGreedy.getF();
			gBar = bestGreedy.getG2();
			adjBar = bestGreedy.getAdjs();
			// For GC
			bestGreedy = null;
		}

		remaining = new TreeSet<>();
		visited = new BinVisitedSet();
		
		// Add first CMRTrial
		remaining.add(new BinTrial(problem.getRangeSet(), nvar, data));

		running = true;
		ExecutorService executor = Executors.newFixedThreadPool(nThread);
		ArrayList<BinTrialSolver> pool = new ArrayList<>();
		for (int i = 0; i < nThread; i++) {
			pool.add(new BinTrialSolver(i));
		}
		try {
			executor.invokeAll(pool);
			executor.shutdown();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		if (catastrophicFailure)
			// Return null since we can't trust the output
			return null;

		int iter = 0;
		for (int i = 0; i < nThread; i++)
			iter += nIterThread[i];

		return new BinSolution(fBar, gBar, xBar, iter, (double) (System.nanoTime() - start) / 1_000_000_000, adjBar,
				mrSolver.getCalls());
	}

	/**
	 * Parallel threads, grabbing items from the queue, solving and adding more
	 * 
	 * 
	 *
	 */
	private class BinTrialSolver implements Callable<Boolean> {
		private int threadNum;

		public BinTrialSolver(int threadNum) {
			this.threadNum = threadNum;
		}

		@Override
		public Boolean call() {
			try {
				int[] tmpCVSet = new int[nvar];
				double gFloor = Double.NEGATIVE_INFINITY;
				DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(ncond, ncond);
				int[][] tmpZoneNumbers = new int[ncond][ncond];
				int[] infeas;
				int[][] cv = ParBinCMRxSolverImpl.this.cv;
				TIntObjectHashMap<int[]> zDecodeCache = new TIntObjectHashMap<>();

				ArrayList<BinTrial> tmpTrials = new ArrayList<>(cv.length);
				while (running) {
					BinTrial current = null;

					waiting.incrementAndGet();
					while (current == null) {
						synchronized (remaining) {
							current = grabEnd != -1 && remaining.size() > grabEnd ? remaining.pollLast() : remaining.pollFirst();
						}
						if (current != null) {
							waiting.decrementAndGet();
							break;
						} else if (waiting.intValue() == nThread || !running) {
							running = false;
							return true;
						}
						Thread.yield();
					}

					gFloor = current.getG2();
					nIterThread[threadNum]++;

					if (gFloor < gBar) {
						double[][] xPrimes = current.solveMRs(mrSolver);

						if (xPrimes == null) {
							// failed MR optimisation
							catastrophicFailure = true;
							running = false;
						} else {
							double gFit = current.getG2();
							infeas = CMRxSolver.isFeasible3n(xPrimes, infeasZones, tmpVolumes, tmpZoneNumbers);

							if (infeas == null) {
								// Solution is feasible
								synchronized (ParBinCMRxSolverImpl.this) {
									if (gFit < gBar) {
										// Solution is better than current best
										gBar = gFit;
										fBar = current.getF();
										xBar = xPrimes;
										adjBar = current.getAdjs();
										// fBarReductions++;
										synchronized (remaining) {
											Iterator<BinTrial> iter = remaining.descendingIterator();
											while (iter.hasNext()) {
												if (iter.next().getG2() > gFit)
													iter.remove();
												else
													break;
											}
										}
									}
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
												newTrial = current.split(k, posIndex, negIndex, newTrialID.getAndIncrement());
											if (covector[k] < 0)
												newTrial = current.split(k, negIndex, posIndex, newTrialID.getAndIncrement());
										}

										if (newTrial != null) {
											tmpTrials.add(newTrial);
										}
									}
								}
								synchronized (visited) {
									Iterator<BinTrial> iter = tmpTrials.iterator();
									while (iter.hasNext()) {
										BinTrial t = iter.next();
										if (visited.contains(t)) {
											collisions++;
											iter.remove();
										} else {
											visited.add(t);
										}
									}
								}
								synchronized (remaining) {
									remaining.addAll(tmpTrials);
								}
								tmpTrials.clear();
							}
						}
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
			return false;
		}
	}

	public int getGrabEnd() {
		return grabEnd;
	}

	public void setGrabEnd(int grabEnd) {
		this.grabEnd = grabEnd;
	}

	public boolean isCatastrophicFailure() {
		return catastrophicFailure;
	}

	public int getCollisions() {
		return collisions;
	}
}
