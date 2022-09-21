package asm.inference;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.util.ESS;

@Description("Gelman-Rubin like criterion for convergence based on trees alone")
public class GRLike extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
	public Input<Double> smoothingInput = new Input<>("smooting", "smoothing factor, which determines how many trees to disregard", 0.5);
	public Input<Double> bInput = new Input<>("b", "threshold determining acceptance bounds for PSFR like statistic", 0.05);

	List<Node>[] trees;
	int nChains;

	// number of consecutive trees added where the 
	// PSRF criterion passes, but the pseudo ESS criterion fails
	int consecutive = 0;
	int targetESS;
	double smoothing, upper, lower;
	
	@Override
	public void initAndValidate() {
		smoothing = smoothingInput.get();
		targetESS = targetESSInput.get();
		double b = bInput.get();
		upper = 1.0 + b;
		lower = 1.0 - b;
	}

	@Override
	public boolean converged() {
		int end = Integer.MAX_VALUE;
		for (List<Node> t : trees) {
			end = Math.min(end, t.size());
		}
		int start = (int)(end * smoothing);
		if (end - start < targetESS) {
			return false;
		}
		
		List<Double> psrf1 = new ArrayList<>();
		List<Double> psrf2 = new ArrayList<>();
		for (int x = start; x < end; x++) {
			psrf1.add(calcPSRF(trees[0], trees[1], x, start, end));
			psrf2.add(calcPSRF(trees[1], trees[0], x, start, end));
		}
		double psrf1median = median(psrf1);
		double psrf2median = median(psrf1);
		
		if (lower < psrf1median && psrf1median < upper && lower < psrf2median && psrf2median < upper) {
			consecutive++;
			if (consecutive >= targetESS) {
                int cutStart = end - consecutive + 1;
                int cutEnd = end;
                if (pseudoEss(trees[0], cutStart, cutEnd) >= targetESS && 
                	pseudoEss(trees[1], cutStart, cutEnd) >= targetESS) {
                	return true;
                }
			}
		} else {
			consecutive = 0;
		}
		
		
		return false;
	}

	private double pseudoEss(List<Node> trees, int cutStart, int cutEnd) {
		List<Double> trace = new ArrayList<>(cutEnd - cutStart);
		for (int i = cutStart; i < cutEnd; i++) {
			trace.add(distance(trees.get(i), trees.get(cutEnd)));
		}
		double ess = ESS.calcESS(trace);
		return ess;
	}

	private double median(List<Double> psrf) {
		Collections.sort(psrf);
		return psrf.get(psrf.size()/2);
	}

	private Double calcPSRF(List<Node> trees1, List<Node> trees2, int k, int start, int end) {
		double varIn = 0;
		for (int i = start; i < end; i++) {
			double d = distance(trees1.get(k), trees1.get(i));
			varIn += d * d;
		}
		double varBetween = 0;
		for (int i = start; i < end; i++) {
			double d = distance(trees1.get(k), trees2.get(i));
			varBetween += d * d;
		}
		double psrf = Math.sqrt(varIn/varBetween);
		return psrf;
	}

	private double distance(Node tree1, Node tree2) {
		// TODO implement
		return 0;
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
