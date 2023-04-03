package asm.inference;


import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.tree.Tree;

@Description("Determine convergence based on two running chains")
public interface PairewiseConvergenceCriterion {

	/**
	 * @return true if the pair of chains converged
	 */
	boolean converged();

	/**
	 * initialise memory based on number of chains
	 * @param nChains number of chains
	 */
	public void setup(int nChains, TraceInfo traceInfo);
	
	
	// public void process(int chainNr, Double [] logLine, Node root);
}
