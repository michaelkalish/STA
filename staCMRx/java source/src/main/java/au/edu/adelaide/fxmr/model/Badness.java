package au.edu.adelaide.fxmr.model;

public enum Badness {
	NONE, DIAGONALISED, DIAGONALS_HAD_ZEROS;

	public static Badness worst(Badness a, Badness b) {
		return values()[Math.max(a.ordinal(), b.ordinal())];
	}
}
