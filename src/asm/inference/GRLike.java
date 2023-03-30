package asm.inference;


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
	public Input<Double> smoothingInput = new Input<>("smooting", "smoothing factor, which determines how proportion of trees to disregard: "
			+ "larger smoothing means more trees included in test", 0.6);
	public Input<Double> bInput = new Input<>("b", "threshold determining acceptance bounds for PSFR like statistic", 0.05);

	public Input<Integer> cacheLimitInput = new Input<>("cacheLimit", 
			"Maximum size of the tree distance cache (default 256). "
			+ "When limit is reached, half of the cache is purged", 256);
	public Input<Integer> ESSSampleSizeInput = new Input<>("sampleSize",
			"number of trees used to calculated psuedo ESS (default 10)", 10);

	public Input<Boolean> twoSidedInput = new Input<>("twoSided", "if true (default) psrfs and pseudo ESSs for both chains are used,"
			+ "otherwise only one psrfs and pseudo ESS is calculated, which takes less computation but can be less robust", true);

	
	List<Tree>[] trees;
	int nChains;

	// number of consecutive trees added where the 
	// PSRF criterion passes, but the pseudo ESS criterion fails
	private int start0 = -1;
	private int targetESS;
	private double smoothing, upper, lower;

	private int cacheLimit;
	// delta = gap between sampled trees due to cache pruning
	private int delta = 1;
	private int N;
	
	private boolean twoSided;

	private DistanceMatrixCache cache = new DistanceMatrixCache(1024);

	
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
		
		twoSided = twoSidedInput.get();
		
		double b = bInput.get();
		upper = 1.0 + b;
		lower = 1.0 - b;
		
	}

	@Override
	public boolean converged() {
		if (twoSided) {
			return converged2sided();
			// return converged1sided(0) && converged1sided(1); 
		} else {
			return converged1sided(0);
		}
	}

	public boolean converged2sided() {
		try { 
			int end = Integer.MAX_VALUE;
			for (List<Tree> t : trees) {
				end = Math.min(end, t.size());
			}
			if (end % delta != 0) {
				// only check when end us divisible by delta
				return false;
			}
			int start = (int)(end * (1.0-smoothing));
			start = start - start % delta;
	
			if (end/delta > cacheLimit) {
				delta *= 2;
				Log.warning("Delta=" + delta);	
				return converged();
			}

			// This following commented out code follows the original pseudo code
			// However, only calculating psrf2mean if psrf1mean passes the test
			// makes things a bit faster.
			
//			List<Double> psrf1 = new ArrayList<>();
//			List<Double> psrf2 = new ArrayList<>();
//			for (int x = start; x < end; x += delta) {
//				psrf1.add(calcPSRF(0, 1, x, start, end));
//				psrf2.add(calcPSRF(1, 0, x, start, end));
//			}
//			double psrf1mean = mean(psrf1);
//			double psrf2mean = mean(psrf1);
//			
//			Log.info("psrf1mean = " + psrf1mean + " psrf2mean = " + psrf2mean);
//			if (lower < psrf1mean && psrf1mean < upper && lower < psrf2mean && psrf2mean < upper) {
			double [] psrf = new double[(end-start)/delta]; 
			for (int x = start; x < end; x += delta) {
				psrf[(x-start)/delta] = calcPSRF(0, 1, x, start, end);
			}
			double psrf1mean = mean(psrf);
			double psrf2mean = 0;
			
			Log.info("\npsrf1mean = " + psrf1mean);
			if (lower < psrf1mean && psrf1mean < upper) {
				for (int x = start; x < end; x += delta) {
					psrf[(x-start)/delta] = calcPSRF(1, 0, x, start, end);
				}
				psrf2mean = mean(psrf);
				Log.info("psrf2mean = " + psrf2mean);
				
				if (lower < psrf2mean && psrf2mean < upper) {
					if (start0 < 0) {
						start0 = start;
					}
	    			// Next condition has "<=" instead of "<", since the pseudoESS later on is based on
	    			// (end-start)/delta - 1 trees (the tree to compare with is removed from the sequence), 
	    			// so we need at least one more tree to make the pseudoESS >= targetESS
	    			if ((end - start0) <= targetESS) {
	    				return false;
	    			}
	                int cutStart = start0; 
	                int cutEnd = end;
	                cutStart = cutStart - cutStart % delta;
	                
	                if (pseudoESS(0, cutStart, cutEnd) >= targetESS && 
	                	pseudoESS(1, cutStart, cutEnd) >= targetESS) {
	                	return true;
	                }
				}
			}
			if (!(lower < psrf1mean && psrf1mean < upper && lower < psrf2mean && psrf2mean < upper)) {
				start0 = -1;
			}
		} catch (Throwable e) {
			// converged();
			Log.warning("Programmer error: Something is WRONG and needs fixing:");
			e.printStackTrace();
			start0 = -1;
			return false;
		}
		Log.info.print(start0 + " ");
		return false;
	}

	
	public boolean converged1sided(int side) {
		try { 
			int end = Integer.MAX_VALUE;
			for (List<Tree> t : trees) {
				end = Math.min(end, t.size());
			}
			if (end % delta != 0) {
				// only check when end us divisible by delta
				return false;
			}
			int start = (int)(end * (1.0-smoothing));
			start = start - start % delta;
			
			if (end/delta > cacheLimit) {
				delta *= 2;
				Log.warning("Delta=" + delta);	
				return converged();
			}

			double [] psrf = new double[(end-start)/delta]; 
			for (int x = start; x < end; x += delta) {
				psrf[(x-start)/delta] = calcPSRF(side-0, 1-side, x, start, end);
			}
			double psrf1mean = mean(psrf);
			
			Log.info("\npsrf1mean = " + psrf1mean);
			if (lower < psrf1mean && psrf1mean < upper) {
				
				
				if (start0 < 0) {
					start0 = start;
				}
				// Next condition has "<=" instead of "<", since the pseudoESS later on is based on
				// (end-start)/delta - 1 trees (the tree to compare with is removed from the sequence), 
				// so we need at least one more tree to make the pseudoESS >= targetESS
    			if ((end - start0) <= targetESS) {
    				return false;
    			}
                int cutStart = start0; 
                cutStart = cutStart - cutStart % delta;
                int cutEnd = end;
                if (pseudoESS(side-0, cutStart, cutEnd) >= targetESS) {
                	return true;
                }
			}
			if (!(lower < psrf1mean && psrf1mean < upper)) {
				start0 = -1;
			}
		} catch (Throwable e) {
			// converged();
			Log.warning("Programmer error: Something is WRONG and needs fixing:");
			e.printStackTrace();
			start0 = -1;
			return false;
		}
		Log.info.print(start0 + " ");
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
		Double [][] trace = new Double[N][(cutEnd-cutStart)/delta - 1];
		int [] k = new int[N];
		for (int i = cutStart; i < cutEnd; i += delta) {
			for (int j = 0; j < N; j++) {
				if (i != indices[j]) {
					// only include distance to other trees
					trace[j][k[j]++] = (double) distancePlusOne(treeSet, i, treeSet, indices[j]);
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
		double meanESS = delta * sum / N;
		
		Log.info("pseudoESS = " + meanESS);
		return meanESS;
	}

	private double mean(double [] psrf) {
		double sum = 0;
		for (double d : psrf) {
			sum += d;
		}
		double mean = sum / psrf.length;
		return mean;
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
				Log.warning("Resizing cache from " + size + " to " + (size+1024));
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
					if (index1 > index2) {  
						int i = index1 * (index1 - 1)/2 + index2;
						return cache21[i];
					} else if (index1 < index2) {
						int i = index2 * (index2 - 1)/2 + index1;
						return cache12[i];
					} else {
						return diagonal[index1];
					}
				} else {
					int i = index1 > index2 ?  
							index1 * (index1 - 1)/2 + index2:
							index2 * (index2 - 1)/2 + index1;
					return cache22[i];
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
		
		// System.err.print(treeSet1 + "x" + treeSet2 + "[" + index1 + "," + index2+"] ");
		// System.out.print(".");
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
