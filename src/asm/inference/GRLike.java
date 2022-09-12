package asm.inference;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;

@Description("Gelman-Rubin like criterion for convergence based on trees alone")
public class GRLike extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean converged() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void setup(int nChains) {
		// TODO Auto-generated method stub
		
	}
	

	@Override
	public void process(int chainNr, Double[] logLine, Node root) {
		// TODO Auto-generated method stub
		
	}

}
