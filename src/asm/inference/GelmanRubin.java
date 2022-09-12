package asm.inference;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;

public class GelmanRubin extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Double> acceptedThresholdInput = new Input<>("threshold", "level at which the biggest GR value is still acceptable", 1.05);

    /** tables of logs, one for each thread + one for the total**/
	List<Double[]>[] m_logTables;
	
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
		int available = m_logTables[0].size();
		for (List<Double[]> d : m_logTables) {
			available = Math.min(available, d.size());
		}
		return calcGRStats(available);
	}


//	http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/coda/html/gelman.diag.html
//	Brooks, SP. and Gelman, A. (1997) 
//	General methods for monitoring convergence of iterative simulations. 
//	Journal of Computational and Graphical Statistics, 
//	7, 
//	434-455. 
//
//  m = # threads
//	n = # samples
//	B = variance within chain
//	W = variance among chains
//	R=(m+1/m)(W(n-1)/n + B/n + B/(mn))/W - (n-1)/nm 
//	=>
//	R=(m+1/m)((n-1)/n + B/Wn + B/(Wmn)) - (n-1)/nm
//	=>
//	R=(m+1/m)((n-1)/n + B/W(1/n + 1/mn) - (n-1)/nm
//	=>
//	R=(m+1/m)((n-1)/n + B/W((m+1)/nm)) - (n-1)/nm
	/** This calculates the Gelman Rubin statistic from scratch (using 10% burn in)
	 * and reports the log of the first chain, annotated with the R statistic.
	 * This number approaches 1 on convergence, so during the run of the chain
	 * you can check how well the chain converges.
	 *
	 * Exploit potential for efficiency by storing means and squared means
	 * NB: when the start of the chain changes, this needs to be taken in account.
	 */
	boolean calcGRStats(int nCurrentSample) {
		int nLogItems = m_logTables[0].get(0).length;
		int nThreads = nChains;
		
		// calculate means and variance, use 10% burn in
		int nSamples = nCurrentSample - nCurrentSample/10;
		
		// the Gelman Rubin statistic for each log item 
		double [] fR = new double [nLogItems];
		if (nSamples > 5) {
			if (m_fSums == null) {
				m_fSums = new double[(nThreads+1)][nLogItems];
				m_fSquaredSums = new double[(nThreads+1)][nLogItems];
			}

			int nStartSample = nCurrentSample/10;
			int nOldStartSample = (nCurrentSample-1)/10;
			if (nStartSample != nOldStartSample) {
				// we need to remove log line from means
				// calc means and squared means
				int iSample = nOldStartSample;
				for (int iThread2 = 0; iThread2 < nThreads; iThread2++) {
					Double[] fLine = m_logTables[iThread2].get(iSample);
					for (int iItem = 1; iItem < nLogItems; iItem++) {
						m_fSums[iThread2][iItem] -= fLine[iItem];
						m_fSquaredSums[iThread2][iItem] -= fLine[iItem] * fLine[iItem];
					}
				}
				
				// sum to get totals
				for (int iItem = 1; iItem < nLogItems; iItem++) {
					double fMean = 0;
					for (int iThread2 = 0; iThread2 < nThreads; iThread2++) {
						fMean += m_logTables[iThread2].get(iSample)[iItem];
					}
					fMean /= nThreads;
					m_fSums[nThreads][iItem] -= fMean;
					m_fSquaredSums[nThreads][iItem] -= fMean * fMean;
				}
			}

			// calc means and squared means
			int iSample = nCurrentSample;
			for (int iThread2 = 0; iThread2 < nThreads; iThread2++) {
				Double[] fLine = m_logTables[iThread2].get(iSample);
				for (int iItem = 1; iItem < nLogItems; iItem++) {
					m_fSums[iThread2][iItem] += fLine[iItem];
					m_fSquaredSums[iThread2][iItem] += fLine[iItem] * fLine[iItem];
				}
			}
			
			// sum to get totals
			for (int iItem = 1; iItem < nLogItems; iItem++) {
				double fMean = 0;
				for (int iThread2 = 0; iThread2 < nThreads; iThread2++) {
					fMean += m_logTables[iThread2].get(iSample)[iItem];
				}
				fMean /= nThreads;
				m_fSums[nThreads][iItem] += fMean;
				m_fSquaredSums[nThreads][iItem] += fMean * fMean;
			}

			// calculate variances for all (including total counts)
			double [][] fVars = new double[(nThreads+1)][nLogItems];
			for (int iThread2 = 0; iThread2 < nThreads + 1; iThread2++) {
				for (int iItem = 1; iItem < nLogItems; iItem++) {
					double fMean = m_fSums[iThread2][iItem];
					double fMean2 = m_fSquaredSums[iThread2][iItem];
					fVars[iThread2][iItem] = (fMean2 - fMean * fMean);
				}
			}
			
			for (int iItem = 1; iItem < nLogItems; iItem++) {
				// average variance for this item
				double fW = 0;
				for (int i = 0 ; i < nThreads; i++ ){
					fW += fVars[i][iItem];
				}
				fW /= (nThreads*(nSamples -1));
				// variance for joint
				double fB = fVars[nThreads][iItem]/((nThreads-1) * nSamples);
				fR[iItem] = ((nThreads + 1.0)/nThreads) * ((nSamples-1.0) / nSamples + fB/fW * (nThreads+1)/(nSamples * nThreads)) - (nSamples-1.0)/(nSamples * nThreads); 
			}
		}

		for (double f : fR) {
			if (f > acceptedThreshold) {
				return false;
			}
		}
		return true;
	} // calcGRStats

	
	@Override
	public void setup(int nChains) {
		this.nChains = nChains;

		// start a thread to tail all log files
		m_logTables = new List[nChains + 1];
		for (int i = 0; i < nChains+1; i++) {
			m_logTables[i] = new ArrayList<Double[]>();
		}
	}
	
	@Override
	public void process(int chainNr, Double [] logLine, Node root) {
		processLogLine(chainNr, logLine);
	}
	
	
	/** Add log line to log table, and check whether other threads are up to date
	 * enough to report on a sample. 
	 */
	synchronized void processLogLine(int iThread, Double[] logLine) {
		m_logTables[iThread].add(logLine);
	}
	
	
	
	
}
