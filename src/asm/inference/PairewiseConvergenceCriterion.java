package asm.inference;


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
	 */
	public void setup(int nChains);
	
	
	public void process(int chainNr, Double [] logLine, Node root);
}
