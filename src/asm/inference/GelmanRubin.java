package asm.inference;


import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;

public class GelmanRubin extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Double> acceptedThresholdInput = new Input<>("threshold", "level at which the biggest GR value is still acceptable", 1.05);

    /** tables of logs, one for each thread + one for the total**/
	List<Double>[][] m_logTables;
	
	/** pre-calculated sum of itmes and sum of itmes squared for all threads and all items */
	double [][] m_fSums;
	double [][] m_fSquaredSums;

	int nChains;
	double acceptedThreshold;

	@Override
	public void initAndValidate() {
		acceptedThreshold = acceptedThresholdInput.get();
	}

	@Override
	public boolean converged() {
		int available = m_logTables[0].length;
		for (List<Double>[] d : m_logTables) {
			available = Math.min(available, d.length);
		}
		
		// check all items for all pairs
		int nItems = m_logTables[0].length;
		for (int i = 0; i < nChains; i++) {
			for (int j = i+1; j < nChains; j++) {
				for (int k = 0; k < nItems; k++) {
					double GRstat = calcGRStat(available, m_logTables[i][k], m_logTables[j][k]);
					if (GRstat > acceptedThreshold) {
						return false;
					}
				}
			}
		}
		return true;
	}

	
	/** original Gelman Rubin statistic for 2 chains **/	
	private double calcGRStat(int sampleCount, List<Double> trace1, List<Double> trace2) {
		if (sampleCount >= trace1.size() || sampleCount > trace2.size()) {
			throw new IllegalArgumentException("Expected traces of sufficient length");
		}
		
		// calc means and squared means
		double mean1 = 0, mean2 = 0, sumsq1 = 0, sumsq2 = 0;
		for (int i = 0; i < sampleCount; i++) {
			Double d = trace1.get(i);
			mean1 += d;
			sumsq1 += d * d;
		}
		mean1 /= sampleCount;
		for (int i = 0; i < sampleCount; i++) {
			Double d = trace2.get(i);
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


	
	@Override
	public void setup(int nChains, List<Double>[][] logLines, List<Node>[] trees) {
		this.nChains = nChains;
		m_logTables = logLines;
	}
	
}
