package asm.inference;



import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

@Description("Gelman-Rubin like criterion for convergence based on trees alone")
public class TreePSRF extends TreeESS implements MCMCConvergenceCriterion {
	public Input<Double> bInput = new Input<>("b", "threshold determining acceptance tolerance for PSFR like statistic", 0.05);

	public Input<Boolean> twoSidedInput = new Input<>("twoSided", "if true (default) psrfs and pseudo ESSs for both chains are used,"
			+ "otherwise only one psrfs and pseudo ESS is calculated, which takes less computation but can be less robust", true);
	
	public Input<Boolean> checkESSInput = new Input<>("checkESS", "whether to check Tree ESS exceeds the targetESS", false);

	//	TraceInfo traceInfo;
//	List<Tree>[] trees;
//	int nChains;

	// number of consecutive trees added where the 
	// PSRF criterion passes, but the pseudo ESS criterion fails
	private int start0 = -1;
//	private int targetESS;
	private double smoothing, upper, lower;
	private boolean prev = false;

//	private int cacheLimit;
	// delta = gap between sampled trees due to cache pruning
//	private int delta = 1;
//	private int N;
	
	private boolean twoSided, checkESS;

//	private DecimalFormat f = new DecimalFormat("#.###");
//	private DecimalFormat f1 = new DecimalFormat("#.#");

	public double[] getGrValues() {
		return grValues;
	}

	public void setGrValues(double[] grValues) {
		this.grValues = grValues;
	}

	private double[] grValues;
	
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
		
		twoSided = twoSidedInput.get();
		checkESS = checkESSInput.get();
		
		double b = bInput.get();
		upper = 1.0 + b;
		lower = 1.0 - b;
		
	}

	@Override
	public boolean converged(int[] burnin, int end) {
		if (twoSided) {
			prev = converged2sided(burnin, end);
			return prev;
			// return converged1sided(0) && converged1sided(1); 
		} else {
			prev = converged1sided(0, burnin, end);
			return prev;
		}
	}

	public boolean converged2sided(int[] burnin, int end) {
		try { 
//			int end = Integer.MAX_VALUE;
//			for (List<Tree> t : trees) {
//				end = Math.min(end, t.size());
//			}
			if (end % delta != 0) {
				// only check when end is divisible by delta
				return prev;
			}
			int start = (int)(end * (1.0-smoothing));
			start = start - start % delta;
	
			if (end/delta > cacheLimit) {
				delta *= 2;
				Log.warning("Delta=" + delta);	
				return converged(burnin, end);
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
				psrf[(x-start)/delta] = calcPSRF(0, 1, x, burnin, end);
			}
			double psrf1mean = mean(psrf);
			this.grValues[0] = psrf1mean; // Setting the value for logging
			double psrf2mean = 0;
			
			Log.info.print("psrf1mean = " + traceInfo.f.format(psrf1mean) + " ");
			if (lower < psrf1mean && psrf1mean < upper) {
				if (start0 < 0) {
					start0 = start;
				}
				for (int x = start; x < end; x += delta) {
					psrf[(x-start)/delta] = calcPSRF(1, 0, x, burnin, end);
				}
				psrf2mean = mean(psrf);
				this.grValues[1] = psrf2mean; // Setting the value for logging
				Log.info.print("psrf2mean = " + traceInfo.f.format(psrf2mean)+ " ");
				
				if (lower < psrf2mean && psrf2mean < upper) {
					if (!checkESS) {
						return true;
					}
	    			// Next condition has "<=" instead of "<", since the pseudoESS later on is based on
	    			// (end-start)/delta - 1 trees (the tree to compare with is removed from the sequence), 
	    			// so we need at least one more tree to make the pseudoESS >= targetESS
	    			if ((end - start0) <= targetESS) {
	    				Log.info.print("not enough samples ");
	    				return false;
	    			}
	                int cutStart = start0; 
	                int cutEnd = end;
	                cutStart = cutStart - cutStart % delta;

	                if (checkESS || pseudoESS(0, cutStart, cutEnd) >= targetESS &&
	                	pseudoESS(1, cutStart, cutEnd) >= targetESS) {
	                	return true;
	                }
	        		Log.info.print(start0 + " ");
	                return false;
				}
			}
			if (lower <= psrf1mean || psrf1mean >= upper || lower <= psrf2mean || psrf2mean >= upper) {
				Log.info.print("reset start ");
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

	
	public boolean converged1sided(int side, int[] burnin, int end) {
		try { 
//			int end = Integer.MAX_VALUE;
//			for (List<Tree> t : trees) {
//				end = Math.min(end, t.size());
//			}
			if (end % delta != 0) {
				// only check when end us divisible by delta
				return prev;
			}
			int start = (int)(end * (1.0-smoothing));
			start = start - start % delta;
			
			if (end/delta > cacheLimit) {
				delta *= 2;
				Log.warning("Delta=" + delta);	
				return converged(burnin, end);
			}

			double [] psrf = new double[(end-start)/delta]; 
			for (int x = start; x < end; x += delta) {
				psrf[(x-start)/delta] = calcPSRF(side-0, 1-side, x, burnin, end);
			}
			double psrf1mean = mean(psrf);
			this.grValues[0] = psrf1mean; // Setting the value for logging
			
			Log.info("\npsrf1mean = " + psrf1mean);
			if (lower < psrf1mean && psrf1mean < upper) {
				if (!checkESS) {
					return true;
				}
				
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
                return false;
			}
			if (lower >= psrf1mean || psrf1mean >= upper) {
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


	private double mean(double [] psrf) {
		double sum = 0;
		for (double d : psrf) {
			sum += d;
		}
		double mean = sum / psrf.length;
		return mean;
	}
	

	private Double calcPSRF(int treeSet1, int treeSet2, int k, int [] burnin, int end) {
		double varIn = 0;
		int start = start(burnin[treeSet1]);
		for (int i = start; i < end; i += delta) {
			double d = distancePlusOne(treeSet1, k, treeSet1, i);
			varIn += d * d;
		}
		double varBetween = 0;
		start = start(burnin[treeSet2]);
		for (int i = start; i < end; i += delta) {
			double d = distancePlusOne(treeSet1, k, treeSet2, i);
			varBetween += d * d;
		}
		double psrf;
		if (varIn != 0) {
			psrf = Math.sqrt(varBetween / varIn);
		} else {
			psrf = Math.sqrt(varBetween);
		}
		return psrf;
	}

	private int start(int start) {
		if (start % delta != 0) {
			start = start + delta- start % delta;
		}
		return start;
	}

	@Override
	public void setup(int nChains, TraceInfo traceInfo) {
		this.traceInfo = traceInfo;
		if (traceInfo.distances == null) {
			traceInfo.distances = new DistanceMatrixCache(1024);
		}
		this.trees = traceInfo.trees;
		this.numChains = nChains;
		if (nChains != 2) {
			throw new IllegalArgumentException("Only 2 chains can be handled by " + this.getClass().getName() + ", not " + nChains);
		}
		// Initializing the grValues -- used for logging
		this.grValues = new double[nChains];
		// setting it to -2.0 indicating it has not been calculated for later postprocessing
		Arrays.fill(grValues, -2.0);
	}

	public Map getLog() {
		// Returns the log values that will be output to the ECCLogger file
		Map<String, Double> logValues = new HashMap<>();

		// Todo add the tree ESS values here? or somewhere else?
		// Todo check whether the grValues each have been updated, if not no need to relog it?
		//  this is because currently it only checks the first, if not within boundary it just goes on to the next
		//  more efficient but harder to log sensibly
		for (int i = 0; i < this.numChains; i++) {
			logValues.put("GRT-" + i, this.grValues[i]);
		}

		return logValues;
	}



}
