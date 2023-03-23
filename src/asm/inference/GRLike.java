package asm.inference;


import java.util.ArrayList;
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

@Description("Gelman-Rubin like criterion for convergence based on trees alone")
public class GRLike extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
	public Input<Double> smoothingInput = new Input<>("smooting", "smoothing factor, which determines how many trees to disregard", 0.5);
	public Input<Double> bInput = new Input<>("b", "threshold determining acceptance bounds for PSFR like statistic", 0.05);

	public Input<Integer> cacheLimitInput = new Input<>("cacheLimit", 
			"Maximum size of the tree distance cache (default 200). "
			+ "When limit is reached, half of the cache is purged", 200);
	public Input<Integer> ESSSampleSizeInput = new Input<>("sampleSize",
			"number of trees used to calculated psuedo ESS (default 10)", 10);

	List<Tree>[] trees;
	int nChains;

	// number of consecutive trees added where the 
	// PSRF criterion passes, but the pseudo ESS criterion fails
	private int consecutive = 0;
	private int targetESS;
	private double smoothing, upper, lower;

	private int cacheLimit;
	// delta = gap between sampled trees due to cache pruning
	private int delta = 1;
	private int N;

	
	@Override
	public void initAndValidate() {
		smoothing = smoothingInput.get();
		targetESS = targetESSInput.get();
		cacheLimit = cacheLimitInput.get();
		N = ESSSampleSizeInput.get();
		double b = bInput.get();
		upper = 1.0 + b;
		lower = 1.0 - b;
	}

	@Override
	public boolean converged() {
		int end = Integer.MAX_VALUE;
		for (List<Tree> t : trees) {
			end = Math.min(end, t.size());
		}
		if (end % delta != 0) {
			// only check when end us divisible by delta
			return false;
		}
		int start = (int)(end * smoothing);
		start = start - start % delta;
		if ((end - start) * delta < targetESS) {
			return false;
		}
		
		if (end/delta > cacheLimit) {
			delta *= 2;
			return converged();
		}
		
		List<Double> psrf1 = new ArrayList<>();
		List<Double> psrf2 = new ArrayList<>();
		for (int x = start; x < end; x += delta) {
			psrf1.add(calcPSRF(0, 1, x, start, end));
			psrf2.add(calcPSRF(1, 0, x, start, end));
		}
		double psrf1mean = mean(psrf1);
		double psrf2mean = mean(psrf1);
		
		if (lower < psrf1mean && psrf1mean < upper && lower < psrf2mean && psrf2mean < upper) {
			consecutive += delta;
			if (consecutive >= targetESS) {
                int cutStart = end - consecutive + 1;
                int cutEnd = end;
                if (pseudoESS(0, cutStart, cutEnd) >= targetESS && 
                	pseudoESS(1, cutStart, cutEnd) >= targetESS) {
                	return true;
                }
			}
		} else {
			consecutive = 0;
		}
		Log.info.print(consecutive + " ");
		return false;
	}

	

	// collects ESSs for N trees and returns the mean (should be median?)
	private double pseudoESS(int treeSet, int cutStart, int cutEnd) {
		
		// subsample N trees in range [cutStart, cutEnd]
		int [] indices = new int[N];
		for (int i = 0; i < N; i++) {
			int offset = (cutEnd-cutStart) * i / N;
			offset = offset - offset % delta;
			indices[i] = cutStart + offset;
		}

		// calc sum of distances to the trees with index from `indices`
		Double [][] trace = new Double[N][(cutEnd-cutStart)/delta];
		for (int i = cutStart; i < cutEnd; i += delta) {
			double d = 0;
			for (int j = 0; j < N; j++) {
				if (i != indices[j]) {
					// only include average distance to other trees
					d += distancePlusOne(treeSet, i, treeSet, indices[j]);
				}
				trace[j][i/delta] = d;
			}
		}

		// calc ESS for each trace
		double [] ess = new double[N];
		for (int j = 0; j < N; j += 1) {
			ess[j] = ESS.calcESS(trace[j], 1);
		}
		// calculate mean ESS
		double sum = 0;
		for (double d : ess) {
			sum +=d;
		}
		double meanESS = sum / N;
		return meanESS * delta;
	}

	private double mean(List<Double> psrf) {
		double sum = 0;
		for (double d : psrf) {
			sum += d;
		}
		double mean = sum / psrf.size();
		return mean;
	}

	private Double calcPSRF(int treeSet1, int treeSet2, int k, int start, int end) {
		double varIn = 0;
		for (int i = start; i < end; i += delta) {
			double d = distancePlusOne(treeSet1, k, treeSet1, i);
			varIn += d * d;
		}
		double varBetween = 0;
		for (int i = start; i < end; i += delta) {
			double d = distancePlusOne(treeSet1, k, treeSet2, i);
			varBetween += d * d;
		}
		double psrf = Math.sqrt(varBetween/varIn);
		return psrf;
	}

	
	class DistanceMatrixCache {
		// symmetric 2d distance matrix for trees 1
		float [] cache11;
		// symmetric 2d distance matrix for trees 2
		float [] cache22;
		// asymmetric 2d distance matrix between trees 1 and trees 2
		// represented by lower triangle and upper triangle half-matrices
		float [] cache12, cache21, diagonal;
		// matrix size
		int size;
		
		DistanceMatrixCache(int n) {
			cache11 = new float[n*(n-1)/2];
			cache22 = new float[n*(n-1)/2];
			cache12 = new float[n*(n-1)/2];
			cache21 = new float[n*(n-1)/2];
			diagonal = new float[n];
		}
		
		float getDistance(int treeSet1, int index1, int treeSet2, int index2) {
			if (index1 >= size || index2 >= size) {
				// resize
				size += 1024;
				cache11 = Arrays.copyOf(cache11, size*(size-1)/2);
				cache22 = Arrays.copyOf(cache22, size*(size-1)/2);
				cache12 = Arrays.copyOf(cache12, size*(size-1)/2);
				cache21 = Arrays.copyOf(cache21, size*(size-1)/2);
				diagonal = Arrays.copyOf(diagonal, size);
			}
			
			if (treeSet1 == 0) {
				if (treeSet2 == 0) {
					int i = index1 > index2 ?  
						index1 * (index1 - 1)/2 + index2:
						index2 * (index2 - 1)/2 + index1;
					return cache11[i];
				} else {
					if (index1 > index2) {  
						int i = index1 * (index1 - 1)/2 + index2;
						return cache12[i];
					} else if (index1 < index2) {  
						int i = index2 * (index2 - 1)/2 + index1;
						return cache21[i];
					} else {
						return diagonal[index1];
					}
				}
			} else {
				if (treeSet2 == 0) {
					int i = index2 * size + index1;  
					return cache12[i];
				} else {
					if (index1 > index2) {  
						int i = index1 * (index1 - 1)/2 + index2;
						return cache21[i];
					} else if (index1 < index2) {
						int i = index2 * (index2 - 1)/2 + index1;
						return cache12[i];
					} else {
						return diagonal[index1];
					}
				}
			}			
		}

		void setDistance(int treeSet1, int index1, int treeSet2, int index2, float d) {
			if (treeSet1 == 0) {
				if (treeSet2 == 0) {
					int i = index1 > index2 ?  
						index1 * (index1 - 1)/2 + index2:
						index2 * (index2 - 1)/2 + index1;
					cache11[i] = d;
				} else {
					if (index1 > index2) {  
						int i = index1 * (index1 - 1)/2 + index2;
						cache12[i] = d;
					} else if (index1 < index2) {
						int i = index2 * (index2 - 1)/2 + index1;
						cache21[i] = d;
					} else {
						diagonal[index1] = d;
					}
				}
			} else {
				if (treeSet2 == 0) {
					if (index1 > index2) {  
						int i = index1 * (index1 - 1)/2 + index2;
						cache21[i] = d;
					} else if (index1 < index2) {
						int i = index2 * (index2 - 1)/2 + index1;
						cache12[i] = d;
					} else {
						diagonal[index1] = d;
					}
				} else {
					int i = index1 > index2 ?  
							index1 * (index1 - 1)/2 + index2:
							index2 * (index2 - 1)/2 + index1;
					cache22[i] = d;
				}
			}			
			
		}
	}
	
	DistanceMatrixCache cache = new DistanceMatrixCache(1024);
	
	private float distancePlusOne(int treeSet1, int index1, int treeSet2, int index2) {
		if (treeSet1 == treeSet2 && index1 == index2) {
			return 0+1; // +1 so that we can use 0 to detect whether the distance is in the cache
		}
		float d = cache.getDistance(treeSet1, index1, treeSet2, index2);
		if (d > 0) {
			return d;
		}
		
		TreeInterface tree1 = trees[treeSet1].get(index1);
		TreeInterface tree2 = trees[treeSet2].get(index2);
		RNNIMetric m = new RNNIMetric();
		d = (float) m.distance(tree1, tree2) + 1; // +1 so that we can use 0 to detect whether the distance is in the cache
		cache.setDistance(treeSet1, index1, treeSet2, index2, d);
		return d;
	}

	@Override
	public void setup(int nChains, List<Double>[][] logLines, List<Tree>[] trees) {
		this.trees = trees;
		this.nChains = nChains;
		if (nChains != 2) {
			throw new IllegalArgumentException("Only 2 chains can be handled by " + this.getClass().getName() + ", not " + nChains);
		}
	}




}
