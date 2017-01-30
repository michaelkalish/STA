package au.edu.adelaide.fxmr.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class GeneralModelReader {

	public static GeneralModel read(File file) {
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			GeneralModel gm = new GeneralModel();

			while ((line = br.readLine()) != null) {
				String[] split = line.split("\\s*,\\s*");
				double[] cond = new double[split.length - 3];
				int nCond = cond.length;
				for (int i = 0; i < nCond; i++) {
					cond[i] = Double.parseDouble(split[i + 3]);
				}
				gm.addData(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Integer.parseInt(split[2]), cond);
			}
			return gm;
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

}
