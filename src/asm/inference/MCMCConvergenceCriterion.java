package asm.inference;



import beast.base.core.Description;

@Description("Determine convergence based on one or more running chains")
public interface MCMCConvergenceCriterion {

	/**
	 * @param burnin TODO
	 * @return true if the pair of chains converged
	 */
	boolean converged(int[] burnin, int end);

	/**
	 * initialise memory based on number of chains
	 * @param nChains number of chains
	 */
	public void setup(int nChains, TraceInfo traceInfo);
	
	
	/**
	 * @return column header name of the convergence criterion for screen logging
	 */
	default public String header() {
		return getClass().getSimpleName();
	}

	
	// public void process(int chainNr, Double [] logLine, Node root);
}
