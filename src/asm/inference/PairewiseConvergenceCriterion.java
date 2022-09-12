package asm.inference;


import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;

@Description("Determine convergence based on two running chains")
public interface PairewiseConvergenceCriterion {

	/**
	 * @return true if the pair of chains converged
	 */
	boolean converged();

	/**
	 * initialise memory based on number of chains
	 * @param nChains number of chains
	 * @param logLines TODO
	 * @param trees TODO
	 */
	public void setup(int nChains, List<Double>[][] logLines, List<Node>[] trees);
	
	
	// public void process(int chainNr, Double [] logLine, Node root);
}
