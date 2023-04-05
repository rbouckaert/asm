package asm.inference;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;

public class TraceESS extends BEASTObject implements MCMCConvergenceCriterion {
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
	public Input<Double> smoothingInput = new Input<>("smoothing", "smoothing factor, which determines how proportion of trees to disregard: "
			+ "larger smoothing means more trees included in test", 0.9);

	protected TraceInfo traceInfo;
	protected List<Double>[][] logLines;
	protected int nChains;
    final static int MAX_LAG = 2000;


	protected int targetESS;
	protected double smoothing;

//	private IncrementalESS[][] esss;
//	private double[][] ess;
	
	@Override
	public void initAndValidate() {
		smoothing = smoothingInput.get();
		targetESS = targetESSInput.get();
	}

	
	double [] trace = new double[1024];
	
	@Override
	public boolean converged(int[] burnin, int end) {
		// calculate space requirement:
		int total = 0;
		for (int i = 0; i < nChains; i++) {
			total += end - burnin[i];
		}
		while (trace.length < total) {
			trace = new double[trace.length + 1024];
		}
		trace = new double[total];
		
		// calculate minimal ESSs of combined logs
		double minESS = Double.POSITIVE_INFINITY;
		int [] map = traceInfo.getMap();
		for (int j = 0; j < map.length; j++) {
			int k = 0;
			for (int i = 0; i < nChains; i++) {
				for (int x = burnin[i]; x < end; x++) {
					trace[k++] = logLines[i][map[j]].get(x);
				}
			}
			double ess = calcESS(trace, 0, total);
			minESS = Math.min(ess,  minESS);
			Log.info.print(TraceInfo.f1.format(ess) + " ");
		}
		Log.info.print(":" + TraceInfo.f1.format(minESS) + "\t");
		return minESS >= targetESS * nChains;
		
		
		
//		// get ESSs
//		for (int i = 0; i < nChains; i++) {
//			for (int j = 0; j < esss[0].length; j++) {
//				// double value = logLines[i][map[j]].get(end);
//				// double ess = esss[i][j].log(value);
//				ess[i][j] = calcESS(logLines[i][map[j]], burnin[i], end);
//				Log.info.print(TraceInfo.f1.format(ess[i][j]) + " ");
//			}
//		}
//		
//		// calculate minimal sum of ESSs
//		double minESS = Double.POSITIVE_INFINITY;
//		for (int j = 0; j < esss[0].length; j++) {
//			double sum = 0;
//			for (int i = 0; i < nChains; i++) {
//				sum += ess[i][j];
//			}
//			minESS = Math.min(sum,  minESS);
//		}
//		Log.info.print(":" + TraceInfo.f1.format(minESS) + "\t");
//		return minESS >= targetESS * nChains;
	}
	
	
	
//	class IncrementalESS {
//	    /**
//	     * values from which the ESS is calculated *
//	     */
//	    protected List<Double> trace;
//	    /**
//	     * sum of trace, excluding burn-in  *
//	     */
//	    protected double sum = 0;
//	    /**
//	     * keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  *
//	     */
//	    protected List<Double> squareLaggedSums;
//
//	    public IncrementalESS() {
//	        trace = new ArrayList<>();
//	        squareLaggedSums = new ArrayList<>();
//	    }
//
//		
//	    public double log(double newValue) {
//	        sum += newValue;
//	        trace.add(newValue);	
//	        final int totalSamples = trace.size();
//	
//	        // take 10% burn in
//	        final int start = totalSamples / 10;
//	        if (start != ((totalSamples - 1) / 10)) {
//	            // compensate for 10% burnin
//	            sum -= trace.get((totalSamples - 1) / 10);
//	        }
//	        final int sampleCount = totalSamples - start;
//	        final int maxLag = Math.min(sampleCount, MAX_LAG);
//	
//	        // calculate mean
//	        final double mean = sum / sampleCount;
//	
//	        if (start != ((totalSamples - 1) / 10)) {
//	            // compensate for 10% burnin
//	            int traceIndex = ((totalSamples - 1) / 10);
//	            for (int lagIndex = 0; lagIndex < squareLaggedSums.size(); lagIndex++) {
//	                squareLaggedSums.set(lagIndex, squareLaggedSums.get(lagIndex) - trace.get(traceIndex) * trace.get(traceIndex + lagIndex));
//	            }
//	        }
//	
//	        while (squareLaggedSums.size() < maxLag) {
//	            squareLaggedSums.add(0.0);
//	        }
//	
//	        // calculate auto correlation for selected lag times
//	        double[] autoCorrelation = new double[maxLag];
//	        // sum1 = \sum_{start ... totalSamples-lagIndex-1} trace
//	        double sum1 = sum;
//	        // sum2 = \sum_{start+lagIndex ... totalSamples-1} trace
//	        double sum2 = sum;
//	        for (int lag = 0; lag < maxLag; lag++) {
//	            squareLaggedSums.set(lag, squareLaggedSums.get(lag) + trace.get(totalSamples - lag - 1) * trace.get(totalSamples - 1));
//	            // The following line is the same approximation as in Tracer 
//	            // (valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
//	            // though a more accurate estimate would be
//	            // autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2
//	            autoCorrelation[lag] = squareLaggedSums.get(lag) - (sum1 + sum2) * mean + mean * mean * (sampleCount - lag);
//	            autoCorrelation[lag] /= (sampleCount - lag);
//	            sum1 -= trace.get(totalSamples - 1 - lag);
//	            sum2 -= trace.get(start + lag);
//	        }
//	
//	        double integralOfACFunctionTimes2 = 0.0;
//	        for (int lagIndex = 0; lagIndex < maxLag; lagIndex++) {
//	            if (lagIndex == 0) {
//	                integralOfACFunctionTimes2 = autoCorrelation[0];
//	            } else if (lagIndex % 2 == 0) {
//	                // fancy stopping criterion - see main comment
//	                if (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex] > 0) {
//	                    integralOfACFunctionTimes2 += 2.0 * (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex]);
//	                } else {
//	                    // stop
//	                    break;
//	                }
//	            }
//	        }
//	
//	        // auto correlation time
//	        final double act = integralOfACFunctionTimes2 / autoCorrelation[0];
//	
//	        // effective sample size
//	        final double ess = sampleCount / act;
//	        return ess;
//	    } // log
//	}
//
	
	

    public static double calcESS(List<Double> trace, int start, int end) {
        return (end - start) / (ACT(trace, start, end));
    }

    public static double ACT(List<Double> trace, int start, int end) {
        /** sum of trace, excluding burn-in **/
        double sum = 0.0;
        /** keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  **/
        double[] squareLaggedSums = new double[MAX_LAG];
        double[] autoCorrelation = new double[MAX_LAG];
        for (int i = 0; i < end-start; i++) {
        	double traceI = trace.get(i+start);
            sum += traceI;
            // calculate mean
            final double mean = sum / (i + 1);

            // calculate auto correlation for selected lag times
            // sum1 = \sum_{start ... totalSamples-lag-1} trace
            double sum1 = sum;
            // sum2 = \sum_{start+lag ... totalSamples-1} trace
            double sum2 = sum;
            for (int lagIndex = 0; lagIndex < Math.min(i + 1, MAX_LAG); lagIndex++) {
            	double traceIMinusLagIndex = trace.get(i+start - lagIndex);
                squareLaggedSums[lagIndex] = squareLaggedSums[lagIndex] + traceIMinusLagIndex * traceI;
                // The following line is the same approximation as in Tracer
                // (valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
                // though a more accurate estimate would be
                // autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2
                autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex);
                autoCorrelation[lagIndex] /= (i+start + 1 - lagIndex);
                sum1 -= traceIMinusLagIndex;
                sum2 -= trace.get(start+lagIndex);
            }
        }

        final int maxLag = Math.min(end - start, MAX_LAG);
        double integralOfACFunctionTimes2 = 0.0;
        for (int lagIndex = 0; lagIndex < maxLag; lagIndex++) //{
            if (lagIndex == 0) //{
                integralOfACFunctionTimes2 = autoCorrelation[0];
            else if (lagIndex % 2 == 0)
                // fancy stopping criterion - see main comment in Tracer code of BEAST 1
                if (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex] > 0) //{
                    integralOfACFunctionTimes2 += 2.0 * (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex]);
                else
                    // stop
                    break;
        //}
        //}
        //}

        // auto correlation time
        return integralOfACFunctionTimes2 / autoCorrelation[0];
    }
	
    private double calcESS(double [] trace, int start, int end) {
        return (end - start) / (ACT(trace, start, end));
    }

    private double ACT(double [] trace, int start, int end) {
        /** sum of trace, excluding burn-in **/
        double sum = 0.0;
        /** keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  **/
        double[] squareLaggedSums = new double[MAX_LAG];
        double[] autoCorrelation = new double[MAX_LAG];
        for (int i = 0; i < end-start; i++) {
        	double traceI = trace[i+start];
            sum += traceI;
            // calculate mean
            final double mean = sum / (i + 1);

            // calculate auto correlation for selected lag times
            // sum1 = \sum_{start ... totalSamples-lag-1} trace
            double sum1 = sum;
            // sum2 = \sum_{start+lag ... totalSamples-1} trace
            double sum2 = sum;
            for (int lagIndex = 0; lagIndex < Math.min(i + 1, MAX_LAG); lagIndex++) {
            	double traceIMinusLagIndex = trace[i+start - lagIndex];
                squareLaggedSums[lagIndex] = squareLaggedSums[lagIndex] + traceIMinusLagIndex * traceI;
                // The following line is the same approximation as in Tracer
                // (valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
                // though a more accurate estimate would be
                // autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2
                autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex);
                autoCorrelation[lagIndex] /= (i+start + 1 - lagIndex);
                sum1 -= traceIMinusLagIndex;
                sum2 -= trace[start+lagIndex];
            }
        }

        final int maxLag = Math.min(end - start, MAX_LAG);
        double integralOfACFunctionTimes2 = 0.0;
        for (int lagIndex = 0; lagIndex < maxLag; lagIndex++) //{
            if (lagIndex == 0) //{
                integralOfACFunctionTimes2 = autoCorrelation[0];
            else if (lagIndex % 2 == 0)
                // fancy stopping criterion - see main comment in Tracer code of BEAST 1
                if (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex] > 0) //{
                    integralOfACFunctionTimes2 += 2.0 * (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex]);
                else
                    // stop
                    break;
        //}
        //}
        //}

        // auto correlation time
        return integralOfACFunctionTimes2 / autoCorrelation[0];
    }

    @Override
	public void setup(int nChains, TraceInfo traceInfo) {
		this.traceInfo = traceInfo;
		this.logLines = traceInfo.logLines;
		this.nChains = nChains;
		
//		esss = new IncrementalESS[nChains][3];
//		for (int i = 0; i < nChains; i++) {
//			for (int j = 0; j < esss[0].length; j++) {
//				esss[i][j] = new IncrementalESS();
//			}
//		}
		// ess = new double[nChains][3];
	}


}
