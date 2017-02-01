package au.edu.adelaide.fxmr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import au.edu.adelaide.fxmr.model.CMRSolution;
import au.edu.adelaide.fxmr.model.CMRxProblem;
import au.edu.adelaide.fxmr.model.CMRxSolver;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

/**
 * Simple test of asymptotic performance
 *
 */
public class DopeyTest {
	private double corCov;
	private DoubleMatrix2D[] weights;
	private CMRxSolver solver;

	public DopeyTest(int startIter, int endIter, int gapIter, double corCov, int samples, File outFile, double tolerance) {
		this.corCov = corCov;
		this.solver = new CMRxSolver();
		solver.setTolerance(tolerance);

		try (PrintWriter out = new PrintWriter(new FileWriter(outFile))) {
			out.println("Run, sample, k, r, f*,seconds,iter, tolerance");
			int run = 1;

			for (int i = startIter; i < endIter; i += gapIter) {
				weights = new DoubleMatrix2D[2];
				weights[0] = DoubleFactory2D.dense.identity(i).assign(Functions.mult(100));
				weights[1] = DoubleFactory2D.dense.identity(i).assign(Functions.mult(100));

				for (int s = 0; s < samples; s++) {
					CMRSolution r = runSample(i);
					out.println(run++ + "," + (s + 1) + "," + i + "," + corCov + "," + r.getFStar() + "," + r.getSeconds() + ","
							+ r.getIters().size() + "," + tolerance);

					out.flush();
					System.out.print(".");
				}
				System.out.println();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println("Complete");
	}

	private CMRSolution runSample(int k) {
		double[][] means = new double[2][k];

		for (int ki = 0; ki < k; ki++) {
			means[0][ki] = Math.random();
			means[1][ki] = corCov * means[0][ki] + Math.sqrt(1 - corCov * corCov) * Math.random();
		}

		CMRxProblem p = new CMRxProblem(means, weights, null, null);

		return solver.solve(p);
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			System.out
					.println("Usage: java DopeyTest <start-iter> <end-iter> <gap-iter> <cor-coef> <n-samples> <output csv> <tolerance=0>");
			return;
		}

		int startIter = Integer.parseInt(args[0]);
		int endIter = Integer.parseInt(args[1]);
		int gapIter = Integer.parseInt(args[2]);
		double corCov = Double.parseDouble(args[3]);
		int samples = Integer.parseInt(args[4]);
		File outFile = new File(args[5]);
		double tolerance = args.length > 6 ? Double.parseDouble(args[6]) : 0;

		new DopeyTest(startIter, endIter, gapIter, corCov, samples, outFile, tolerance);
	}

}
