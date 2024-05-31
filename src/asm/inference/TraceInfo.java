package asm.inference;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.evolution.tree.Tree;

/** collection of data from pairs of trace logs and pairs of tree logs **/
public class TraceInfo {
	
	public static DecimalFormat f = new DecimalFormat("#.###");
	public static DecimalFormat f1 = new DecimalFormat("#.#");

    /** tables of logs, one for each thread + one for the total**/
	List<Double>[][] logLines;
	
	/** column labels of logLines **/
	String [] labels;
	
    /** tables of trees, one for each thread + one for the total **/
	List<Tree> [] trees;

	DistanceMatrixCache distances;
	
	public List<Tree>[] getTrees() {
		return trees;
	}
	
	public TraceInfo(int chainCount) {
		f.setMinimumFractionDigits(3);
		f1.setMinimumFractionDigits(1);

		trees = new List[chainCount];
		for (int i = 0; i < chainCount; i++) {
			trees[i] = new ArrayList<>();
		}

		logLines = new List[chainCount][];
	}

	
	// maps "posterior", "likelihood" and "prior" columns to logLine columns
	private int [] map;
	private static int [] DEFAULT_MAP = {1,2,3};

	public int[] getMap() {
		if (map != null) {
			return map;
		}
		if (labels != null) {
			setUpMap(tracesString);
//			map = new int[3];
//			map[0] = indexOf(labels, "posterior", 1);
//			map[1] = indexOf(labels, "likelihood", 2);
//			map[2] = indexOf(labels, "prior", 3);
			return map;
		}
		return DEFAULT_MAP;
	}

	private int indexOf(String[] labels, String string, int defaultIndex) {
		for (int i = 0; i < labels.length; i++) {
			if (labels[i].equals(string)) {
				return i;
			}
		}
		return defaultIndex;
	}

	public int chainCount() {
		return trees.length;
	}

	public void setTrees(List<List<Tree>> treesList) {
		for (int i = 0; i < treesList.size(); i++) {
			trees[i] = treesList.get(i);
		}
	}

	
	String tracesString = null;
	
	public void setUpMap(String tracesString) {
		if (labels == null) {
			this.tracesString = tracesString;
			return;
		}
		String [] traces = tracesString.split(",");

		map = new int[traces.length];
		int k = 0;
		for (String trace : traces) {
			map[k] = indexOf(labels, trace.trim(), -1);
			if (map[k] == -1) {
				throw new IllegalArgumentException("Could not find label " + trace + " in trace log. "
						+ "Use one of " + Arrays.toString(labels));
			}
			k++;
		}
	}

}
