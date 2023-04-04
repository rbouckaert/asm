package asm.inference;


import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.util.ESS;
import beastlabs.evolution.tree.RNNIMetric;

@Description("Tree ESS criterion for convergence based on trees alone")
public class TreeESS extends BEASTObject implements MCMCConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
	public Input<Double> smoothingInput = new Input<>("smoothing", "smoothing factor, which determines how proportion of trees to disregard: "
			+ "larger smoothing means more trees included in test", 0.9);

	public Input<Integer> cacheLimitInput = new Input<>("cacheLimit", 
			"Maximum size of the tree distance cache (default 1024). "
			+ "When limit is reached, half of the cache is purged", 1024);
	public Input<Integer> ESSSampleSizeInput = new Input<>("sampleSize",
			"number of trees used to calculated psuedo ESS (default 10)", 10);

	protected TraceInfo traceInfo;
	protected List<Tree>[] trees;
	protected int nChains;

	protected int targetESS;
	protected double smoothing;

	protected int cacheLimit;
	// delta = gap between sampled trees due to cache pruning
	protected int delta = 1;
	protected int N;
	
	
	@Override
	public void initAndValidate() {
		smoothing = smoothingInput.get();
		if (smoothing < 0 || smoothing >= 1) {
			throw new IllegalArgumentException("smoothing should be between 0 and 1, not " + smoothing);
		}

		targetESS = targetESSInput.get();
		if (targetESS <= 0 ) {
			throw new IllegalArgumentException("targetESS should be positive, not " + targetESS);
		}
		
		cacheLimit = cacheLimitInput.get();
		if (cacheLimit <= 0 ) {
			throw new IllegalArgumentException("cacheLimit should be positive, not " + cacheLimit);
		}
		if (cacheLimit <= targetESS ) {
			throw new IllegalArgumentException("cacheLimit (" + cacheLimit+ ") should be at least as big as targetESS, "
					+ " (" + targetESS +"), preferrably at least twice as large");
		}
		if (cacheLimit < 2*targetESS ) {
			Log.warning("cacheLimit (" + cacheLimit+ ") should preferrably be twice as large as targetESS, "
					+ " (" + targetESS +") to prevent late stopping");
		}
		
		N = ESSSampleSizeInput.get();
		if (N <= 0 ) {
			throw new IllegalArgumentException("sampleSize should be positive, not " + N);
		}
		if (N >= targetESS ) {
			throw new IllegalArgumentException("sampleSize should be less than targetESS ("+targetESS+"), not " + N);
		}
	}

	@Override
	public boolean converged(int[] burnin, int end) {
//		int end = Integer.MAX_VALUE;
//		for (List<Tree> t : trees) {
//			end = Math.min(end, t.size()-1);
//		}
		end = end - end % delta;
		
//		if (end % delta != 0) {
//			// only check when end us divisible by delta
//			Log.info("TreeESS direct fail");
//			return false;
//		}
		
		if (end/delta > cacheLimit) {
			delta *= 2;
//			if (indices != null) {
//				for (int i = 0; i < N; i++) {
//					indices[i] = indices[i] - indices[i] % delta;
//				}
//			}
			Log.warning("TreeESS Delta=" + delta);	
			return converged(burnin, end);
		}

		int start = (int)(end * (1.0-smoothing));
		start = start - start % delta;

		for (int i = 0; i < nChains; i++) {
			if (pseudoESS(i, start, end) < targetESS) {
				return false;
			}
		}
		return true;
	}

	
	int [] indices;
	// collects ESSs for N trees and returns the mean (should be median?)
	protected double pseudoESS(int treeSet, int cutStart, int cutEnd) {
		
		// subsample N trees in range [cutStart, cutEnd]
		if (indices == null) {
			indices = new int[N];
			for (int i = 0; i < N; i++) {
				int offset = (cutEnd-cutStart) * i / N;
				offset = offset - offset % delta;
				indices[i] = cutStart + offset;
			}
		} else {
			while (indices[0] < cutStart) {
				System.arraycopy(indices, 1, indices, 0, N-1);
				indices[N-1] = cutEnd - cutEnd % delta;
			}
		}

		
		System.out.println(Arrays.toString(indices));
		// calc sum of distances to the trees with index from `indices`
		Double [][] trace = new Double[N][(cutEnd-cutStart)/delta];
		int [] k = new int[N];
		for (int i = cutStart; i < cutEnd; i += delta) {
			for (int j = 0; j < N; j++) {
				if (i != indices[j]) {
					// only include distance to other trees
					trace[j][k[j]++] = (double) distancePlusOne(treeSet, i, treeSet, indices[j]);
				} else {
					trace[j][k[j]++] = 0.0;
				}
			}
		}

		// calc ESS for each trace
		double [] ess = new double[N];
		// TODO: use incremental ESS calculation
		// instead of calculating it from scratch every time
		for (int j = 0; j < N; j += 1) {
			ess[j] = ESS.calcESS(trace[j], 1);
		}
		// calculate mean ESS
		double sum = 0;
		for (double d : ess) {
			sum +=d;
		}
		double meanESS = sum / N;
		
		Log.info.print("pseudoESS = " + traceInfo.f1.format(meanESS) + " ");
		return meanESS;
	}

	protected float distancePlusOne(int treeSet1, int index1, int treeSet2, int index2) {
		if (treeSet1 == treeSet2 && index1 == index2) {
			return 0+1; // +1 so that we can use 0 to detect whether the distance is in the cache
		}
		float d = traceInfo.distances.getDistance(treeSet1, index1, treeSet2, index2);
		if (d > 0) {
			return d;
		}
		
		TreeInterface tree1 = trees[treeSet1].get(index1);
		TreeInterface tree2 = trees[treeSet2].get(index2);
		RNNIMetric m = new RNNIMetric();
		d = (float) m.distance(tree1, tree2) + 1; // +1 so that we can use 0 to detect whether the distance is in the cache
		traceInfo.distances.setDistance(treeSet1, index1, treeSet2, index2, d);
		
		// System.err.print(treeSet1 + "x" + treeSet2 + "[" + index1 + "," + index2+"] ");
		// System.out.print(".");
		return d;
	}

	@Override
	public void setup(int nChains, TraceInfo traceInfo) {
		this.traceInfo = traceInfo;
		if (traceInfo.distances == null) {
			traceInfo.distances = new DistanceMatrixCache(1024);
		}
		this.trees = traceInfo.trees;
		this.nChains = nChains;
		if (nChains != 2) {
			throw new IllegalArgumentException("Only 2 chains can be handled by " + this.getClass().getName() + ", not " + nChains);
		}
	}




}
