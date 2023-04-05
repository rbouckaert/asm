package asm.inference;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Log;

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
		Log.info.print(" ");
		return maxBurnin;
	}
	
	private int burnIn(List<Double> trace, int end) {
		int iStart = -1;
		for (int i = 1; i < end - 10; i+= 2) {
			int sampleCount = (end-i)/2;
			
			double overlap = rangeOverlap(i, sampleCount, i+sampleCount, trace);
			if (overlap >= 0.75) {
				iStart = i;
				break;
			}
		}
		
		if (iStart > 0) {
			// try to see if higher ESS can be achieved by adding more burnin
			if (iStart + (end/2-iStart)/10 < end) {
				double maxESS = TraceESS.calcESS(trace, iStart, end);
				double ess = TraceESS.calcESS(trace, iStart + (end/2-iStart)/10, end);
				if (ess > maxESS) {
					Log.info.print("*");
					return iStart + (end/2-iStart)/10;
				}
			}
			return iStart;
//			int maxI = iStart;
//			int delta = Math.max(1,(end/2-iStart)/3); 
//			for (int i = iStart; i < end/2; i+= delta) {
//				double ess = TraceESS.calcESS(trace, i, end);
//				if (ess > maxESS) {
//					maxESS = ess;
//					maxI = i;
//					if (maxI != iStart) {
//						Log.info.print("*");
//					}
//				}
//			}	
//			if (maxI != iStart) {
//				Log.info.print(" ");
//			}
//			return maxI;
		}
		return trace.size()-1;
	}
	
	private double rangeOverlap(int from1, int sampleCount, int from2, List<Double> trace) {
		
		double min1 = trace.get(from1), min2 = trace.get(from2), max1 = min1, max2 = min2;
		for (int i = from1; i < from1 + sampleCount/2; i++) {
			double d = trace.get(i);
			min1 = Math.min(min1, d);
			max1 = Math.max(max1, d);
		}
		for (int i = from2; i < from2 + sampleCount/2; i++) {
			double d = trace.get(i);
			min2 = Math.min(min2, d);
			max2 = Math.max(max2, d);
		}
		double min3 = trace.get(from2 + sampleCount/2), max3 = min3;
		for (int i = from2 + sampleCount/2; i < from2 + sampleCount; i++) {
			double d = trace.get(i);
			min3 = Math.min(min3, d);
			max3 = Math.max(max3, d);
		}
		if (sampleCount <= 2) {
			return 0;
		}
//		double min1 = trace.get(from1)+trace.get(from1+1)+trace.get(from1+2), 
//			   min2 = trace.get(from2)+trace.get(from2+1)+trace.get(from2+2), 
//			   max1 = min1, max2 = min2;
//		for (int i = from1; i < from1 + sampleCount-2; i++) {
//			double d = trace.get(i)+trace.get(i+1)+trace.get(i+2);
//			min1 = Math.min(min1, d);
//			max1 = Math.max(max1, d);
//		}
//		for (int i = from2; i < from2 + sampleCount-2; i++) {
//			double d = trace.get(i)+trace.get(i+1)+trace.get(i+2);
//			min2 = Math.min(min2, d);
//			max2 = Math.max(max2, d);
//		}
		
		
//		double m1 = 0, m2 = 0, sq1 = 0, sq2 = 0;
//		for (int i = from1; i < from1 + sampleCount; i++) {
//			double d = trace.get(i);
//			m1 += d;
//			sq1 += d*d;
//		}
//		double stdev1 = Math.sqrt((m1 * m1 - sq1) / (sampleCount-1));
//		for (int i = from2; i < from2 + sampleCount; i++) {
//			double d = trace.get(i);
//			m2 += d;
//			sq2 += d*d;
//		}
//		double stdev2 = Math.sqrt((m2 * m2 - sq2) / (sampleCount-1));
//		
//		double min1 = m1 - stdev1;
//		double min2 = m2 - stdev2;
//		double max1 = m1 + stdev1;
//		double max2 = m2 + stdev2;
		
		double overlap12, overlap21, overlap13, overlap31;
		
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

		if (min1 < min2) {
			if (max1 < max2) {
				overlap13 = (max1-min2) / (max1-min1);
				overlap31 = (max1-min2) / (max2-min2);
			} else {
				overlap13 = (max2-min2) / (max1-min1);
				overlap31 = (max2-min2) / (max2-min2);
			}
		} else {
			if (max1 < max2) {
				overlap13 = (max1-min1) / (max1-min1);
				overlap31 = (max1-min1) / (max2-min2);				
			} else {
				overlap13 = (max2-min1) / (max1-min1);
				overlap31 = (max2-min1) / (max2-min2);
			}
		}

		return Math.min(Math.max(overlap12, overlap13), Math.max(overlap21, overlap31));
	}

}
