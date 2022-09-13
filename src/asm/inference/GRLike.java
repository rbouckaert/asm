package asm.inference;


import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;

@Description("Gelman-Rubin like criterion for convergence based on trees alone")
public class GRLike extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);

	List<Node>[] trees;
	int nChains;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public boolean converged() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void setup(int nChains, List<Double>[][] logLines, List<Node>[] trees) {
		this.trees = trees;
		this.nChains = nChains;
		if (nChains != 2) {
			throw new IllegalArgumentException("Only 2 chains can be handled by " + this.getClass().getName() + ", not " + nChains);
		}
	}

}
