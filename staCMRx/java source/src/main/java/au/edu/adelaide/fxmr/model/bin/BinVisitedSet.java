package au.edu.adelaide.fxmr.model.bin;

import java.util.ArrayList;
import java.util.HashSet;
import au.edu.adelaide.fxmr.model.HashableAdjSet;

public class BinVisitedSet {
	private HashSet<HashableAdjSet> store = new HashSet<>();

	public boolean add(BinTrial current) {
		return store.add(current.gethAdjSet());
	}

	public boolean contains(BinTrial newTrial) {
		return store.contains(newTrial.gethAdjSet());
	}

	public int size() {
		return store.size();
	}

	public void addAll(ArrayList<BinTrial> tmpTrials) {
		for (BinTrial t : tmpTrials)
			add(t);
	}
}
