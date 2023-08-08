package asm.inference;


import java.util.List;

import beast.base.core.Description;
import beast.base.core.Log;

@Description("Detects burnin as maximum of burnin of posterior, likelihood and prior based on "
		+ "when running average starts to be inside mean +/- one stdev estimated on last quarter of trace")
public class BurnInDetector {
	TraceInfo traceInfo;
		
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
	
	final static int WINDOW_SIZE = 10;
	
	private int burnIn(List<Double> trace, int end) {
		// calc mean and stdev of last 25%
		double m2 = 0, sq2 = 0;
		int lb = 3*end/4;
		for (int i = lb; i < end; i++) {
			double d = trace.get(i);
			m2 += d;
		}
		m2 = m2 / (end - lb);
		for (int i = lb; i < end; i++) {
			double d = trace.get(i);
			sq2 += (d - m2) * (d - m2);
		}
		double stdev2 = Math.sqrt(sq2/(end-lb - 1));
	
		// pick first moving average that fits inside the m2 +/- stdev2 range
		for (int i = 0; i + WINDOW_SIZE < lb; i++) {
			double sum = 0;
			for (int j = 0; j < WINDOW_SIZE; j++) {
				sum += trace.get(i+j);
			}
			sum /= WINDOW_SIZE;
			if (m2 - stdev2 < sum && sum < m2 + stdev2) {
				return i + WINDOW_SIZE;
			}
		}
		return lb;
	}
	

}
