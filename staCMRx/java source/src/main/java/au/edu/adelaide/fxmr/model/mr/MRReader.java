package au.edu.adelaide.fxmr.model.mr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class MRReader {
	public static MRProblem read(File file) {
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			double[] y = readY(br);
			int n = y.length;
			double[][] weights = readWeights(br, n);
			int[][] rangeSet = readRangeSet(br, n);
			return new MRProblem(y, weights, rangeSet);
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	private static int[][] readRangeSet(BufferedReader br, int n) throws NumberFormatException, IOException {
		ArrayList<Object> ranges = new ArrayList<>();

		String line;
		while ((line = br.readLine()) != null) {
			String[] split = line.split("\\s*,\\s*");
			int[] curRange = new int[split.length];
			for (int c = 0; c < split.length; c++)
				//-1 Is to convert into zero indexing
				curRange[c] = Integer.parseInt(split[c]) - 1;
			ranges.add(curRange);
		}

		int[][] rangeSet;
		if (ranges.size() == 0) {
			rangeSet = new int[1][n];
			for (int i = 0; i < n; i++)
				rangeSet[0][i] = i;
		} else {
			rangeSet = new int[ranges.size()][n];
			for (int i = 0; i < ranges.size(); i++)
				rangeSet[i] = (int[]) ranges.get(i);
		}

		return rangeSet;
	}

	private static double[][] readWeights(BufferedReader br, int n) throws IOException {
		double[][] w = null;

		for (int r = 0; r < n; r++) {
			String line = br.readLine();
			if (line == null)
				return null;
			String[] split = line.split("\\s*,\\s*");
			if (w == null)
				w = new double[n][n];
			if (split.length == n) {
				// Full matrix defined
				for (int c = 0; c < n; c++)
					w[r][c] = Double.parseDouble(split[c]);
			} else {
				// Diagonal defined
				w[r][r] = Double.parseDouble(split[0]);
			}
		}

		return w;
	}

	private static double[] readY(BufferedReader br) throws IOException {
		String[] ytxt = br.readLine().split("\\s*,\\s*");
		int n = ytxt.length;
		double[] y = new double[n];
		int[][] rangeSet = new int[1][n];
		for (int i = 0; i < n; i++) {
			y[i] = Double.parseDouble(ytxt[i]);
			rangeSet[0][i] = i;
		}

		return y;
	}
}
