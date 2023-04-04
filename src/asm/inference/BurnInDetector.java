package asm.inference;

import java.util.List;

import beast.base.core.Description;

@Description("Detects burnin as maximum of burnin of posterior, likelihood and prior based on "
		+ "Gelman-Rubin statistic of first and second half of trace after burnin")
public class BurnInDetector {
	TraceInfo traceInfo;
	
	double b = 0.05;
	
	public BurnInDetector(TraceInfo traceInfo) {
		this.traceInfo = traceInfo;
	}
	
	public int [] burnIn(int end) {
		int chainCount = traceInfo.chainCount();
		int [] burnin = new int[chainCount];
		for (int i = 0; i < chainCount; i++) {
			burnin[i] = burnIn(i, end);
		}
		return burnin;
	}
	
	private int burnIn(int chainNr, int end) {
		int maxBurnin = 0;
		int [] map = traceInfo.getMap();
		for (int j = 0; j < map.length; j++) {
			int burnin = burnIn(traceInfo.logLines[chainNr][map[j]], end);
			maxBurnin = Math.max(burnin,  maxBurnin);
		}
		return maxBurnin;
	}
	
	private int burnIn(List<Double> trace, int end) {
		int iStart = -1;
		for (int i = 0; i < end - 10; i+= 2) {
			int sampleCount = (end-i)/2;
//			double gr = calcGRStat(i, sampleCount, i + sampleCount, trace);
//			if (1.0-b < gr && gr <1.0+b) {
//				return i;
//			}
//			double p = calcKDEStat(i, sampleCount, i+sampleCount, trace);
//			if (p < 0.2) {
//				return i;
//			}
			
			double overlap = rangeOverlap(i, sampleCount, i+sampleCount, trace);
			if (overlap >= 0.75) {
				iStart = i;
				break;
			}
		}
		
		if (iStart > 0) {
			// try to see if higher ESS can be achieved by adding more burnin 
			int maxI = iStart;
			double maxESS = 0;
			int delta = Math.max(1,(end/2-iStart)/3); 
			for (int i = iStart; i < end/2; i+= delta) {
				double ess = TraceESS.calcESS(trace, i, end);
				if (ess > maxESS) {
					maxESS = ess;
					maxI = iStart;
				}
			}		
			return maxI;
		}
		return trace.size()-1;
	}
	
	private double rangeOverlap(int from1, int sampleCount, int from2, List<Double> trace) {
		
		double min1 = trace.get(from1), min2 = trace.get(from2), max1 = min1, max2 = min2;
		for (int i = from1; i < from1 + sampleCount; i++) {
			double d = trace.get(i);
			min1 = Math.min(min1, d);
			max1 = Math.max(max1, d);
		}
		for (int i = from2; i < from2 + sampleCount; i++) {
			double d = trace.get(i);
			min2 = Math.min(min2, d);
			max2 = Math.max(max2, d);
		}
		
		double overlap12, overlap21;
		
		if (min1 < min2) {
			if (max1 < max2) {
				overlap12 = (max1-min2) / (max1-min1);
				overlap21 = (max1-min2) / (max2-min2);
			} else {
				overlap12 = (max2-min2) / (max1-min1);
				overlap21 = (max2-min2) / (max2-min2);
			}
		} else {
			if (max1 < max2) {
				overlap12 = (max1-min1) / (max1-min1);
				overlap21 = (max1-min1) / (max2-min2);				
			} else {
				overlap12 = (max2-min1) / (max1-min1);
				overlap21 = (max2-min1) / (max2-min2);
			}
		}
		
		return Math.min(overlap12, overlap21);
	}

	/** original Gelman Rubin statistic for 2 chains **/	
	private double calcGRStat(final int from1, final int sampleCount, final int from2, List<Double> trace) {
		// calc means and squared means
		double mean1 = 0, mean2 = 0, sumsq1 = 0, sumsq2 = 0;
		for (int i = from1; i < from1 + sampleCount; i++) {
			Double d = trace.get(i);
			mean1 += d;
			sumsq1 += d * d;
		}
		mean1 /= sampleCount;
		for (int i = from2; i < from2 + sampleCount; i++) {
			Double d = trace.get(i);
			mean2 += d;
			sumsq2 += d * d;
		}
		mean2 /= sampleCount;

		// calculate variances for both chains
		double var1 = (sumsq1 - mean1 * mean1 * sampleCount)/(sampleCount - 1);
		double var2 = (sumsq2 - mean2 * mean2 * sampleCount)/(sampleCount - 1);
		
		// average variance for this item
		double fW = (var1 + var2) / 2;
		if (fW == 0) {
			return 1;
		}

		// sum to get totals
		double totalMean = (mean1 + mean2) / 2;
		double totalSq = mean1*mean1 + mean2*mean2;
		
		// variance for joint
		double fB = (totalSq - totalMean * totalMean * 2);
		
		
		double varR = ((sampleCount - 1.0)/sampleCount) + (fB/fW)*(1.0/sampleCount);
		double R = Math.sqrt(varR);
		return R;
	}
	
	
	
	final static private int RANGE = 10;
	private double calcKDEStat(final int from1, final int sampleCount, final int from2, List<Double> trace) {
		double mean = 0, sumsq = 0, min = trace.get(from1), max = trace.get(from1);
		for (int i = from1; i < from1 + sampleCount; i++) {
			double d = trace.get(i);
			mean += d;
			sumsq += d*d;
			min = Math.min(min, d);
			max = Math.max(max, d);
		}
		for (int i = from2; i < from2 + sampleCount; i++) {
			double d = trace.get(i);
			mean += d;
			sumsq += d*d;
			min = Math.min(min, d);
			max = Math.max(max, d);
		}
		if (max == min) {
			return 0;
		}
		
		int n = sampleCount * 2;
		
		mean /= n;
		double stdev0 = Math.sqrt((sumsq * sumsq - n * mean * mean) / (n-1));
		double stdev = 0.9 * stdev0 / Math.pow(n, 0.2);
		
		// pre-calculate contribution to plots
		double delta = (max - min) / RANGE;
		int tail = (int)(2*stdev / delta + 0.5);
		tail = 20;
		stdev = Math.sqrt(10.0*delta);
		
		double [] kernel = new double[tail*2 + 1];
		kernel[tail] = N(0.0, stdev);
		for (int i = 0; i < tail; i++) {
			kernel[i] = N((tail-i)*delta, stdev);
			kernel[tail*2-i] = kernel[i];
		}
		
		
		// create plots
		double [] plot1 = new double[RANGE];
		double [] plot2 = new double[RANGE];
		createPlot(trace, from1, sampleCount, plot1, min, kernel, delta, tail);
		createPlot(trace, from2, sampleCount, plot2, min, kernel, delta, tail);
		
		// calculate difference in plots
		double diff12  = 0, diff21 = 0;
		for (int i = 0; i < RANGE; i++) {
			double diff = plot1[i] - plot2[i];
			if (diff > 0) {
				diff12 += diff;
			} else {
				diff21 -= diff;
			}
		}
		// diff12 should be equal to diff21
		
		return diff12 + diff21;
	}
	
	private void createPlot(List<Double> trace, int from, int sampleCount, double[] plot1, double min, double[] kernel, double delta, int tail) {
		for (int j = from; j < from + sampleCount; j++) {
			double d = trace.get(j);
			int centre = (int)((d-min) / delta + 0.5);
			int lower = Math.max(0,  centre-tail);
			int upper = Math.min(RANGE-1, centre + tail);
			for (int i = lower; i < upper; i++) {
				plot1[i] += kernel[i-centre+tail];
			}
		}
		
		// normalise
		double sum = 0;
		for (double d : plot1) {
			sum += d;
		}
		for (int i = 0; i < plot1.length; i++) {
			plot1[i] /= sum;
		}
	}

	private double N(double x, double stdev) {
		double f = 1.0/Math.sqrt(2.0*Math.PI*stdev*stdev) * Math.exp(-x*x/(2.0*stdev*stdev));
		return f;
	}
	}
