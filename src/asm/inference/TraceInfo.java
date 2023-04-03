package asm.inference;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Tree;

public class TraceInfo {
    /** tables of logs, one for each thread + one for the total**/
	List<Double>[][] logLines;
	
    /** tables of trees, one for each thread + one for the total **/
	List<Tree> [] trees;

	DistanceMatrixCache distances;
	
	public TraceInfo(int chainCount) {
		trees = new List[chainCount];
		for (int i = 0; i < chainCount; i++) {
			trees[i] = new ArrayList<>();
		}

		logLines = new List[chainCount][];
	}
	
}
