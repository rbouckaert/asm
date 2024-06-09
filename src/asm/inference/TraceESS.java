package asm.inference;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

@Description("Stopping criterion based on ESS of selected items from trace")
public class TraceESS extends BEASTObject implements MCMCConvergenceCriterion {
    public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
    public Input<String> tracesInput = new Input<>("traces", "comma separated string of trace entries to be tracked", "posterior,prior,likelihood");

    protected TraceInfo traceInfo;
    /** Tables of logs, one for each thread; [chainIndex][attribute]{values} */
    protected List<Double>[][] logLines;
    protected int nChains;
    final static int MAX_LAG = 2000;

    protected int targetESS;
    private double[] currentESSs;

    @Override
    public void initAndValidate() {
        targetESS = targetESSInput.get();
    }

    double[] trace = new double[1024];

    @Override
    public boolean converged(int[] burnin, int end) {
        return converged(burnin, end, true);
    }

    public boolean converged(int[] burnin, int end, boolean useMapping) {
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
        int[] map = traceInfo.getMap();
        this.currentESSs = new double[map.length];
        for (int j = 0; j < map.length; j++) {
            int traceIndex = useMapping ? map[j] : j;
            int k = 0;
            for (int i = 0; i < nChains; i++) {
                for (int x = burnin[i]; x < end; x++) {
                    trace[k++] = logLines[i][traceIndex].get(x);
                }
            }
            double ess = calcESS(trace, 0, total);
            currentESSs[j] = ess;
            minESS = Math.min(ess, minESS);
            Log.info.print(TraceInfo.f1.format(ess) + " ");
        }
        Log.info.print(":" + TraceInfo.f1.format(minESS) + "\t");


        return minESS >= targetESS * nChains;
    }

    public static double calcESS(List<Double> trace, int start, int end) {
        return (end - start) / (ACT(trace, start, end));
    }

    public static double ACT(List<Double> trace, int start, int end) {
        /** sum of trace, excluding burn-in **/
        double sum = 0.0;
        /** keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  **/
        double[] squareLaggedSums = new double[MAX_LAG];
        double[] autoCorrelation = new double[MAX_LAG];
        for (int i = 0; i < end - start; i++) {
            double traceI = trace.get(i + start);
            sum += traceI;
            // calculate mean
            final double mean = sum / (i + 1);

            // calculate auto correlation for selected lag times
            // sum1 = \sum_{start ... totalSamples-lag-1} trace
            double sum1 = sum;
            // sum2 = \sum_{start+lag ... totalSamples-1} trace
            double sum2 = sum;
            for (int lagIndex = 0; lagIndex < Math.min(i + 1, MAX_LAG); lagIndex++) {
                double traceIMinusLagIndex = trace.get(i + start - lagIndex);
                squareLaggedSums[lagIndex] = squareLaggedSums[lagIndex] + traceIMinusLagIndex * traceI;
                // The following line is the same approximation as in Tracer
                // (valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
                // though a more accurate estimate would be
                // autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2
                autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex);
                autoCorrelation[lagIndex] /= (i + start + 1 - lagIndex);
                sum1 -= traceIMinusLagIndex;
                sum2 -= trace.get(start + lagIndex);
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

    private double calcESS(double[] trace, int start, int end) {
        return (end - start) / (ACT(trace, start, end));
    }

    private double ACT(double[] trace, int start, int end) {
        /** sum of trace, excluding burn-in **/
        double sum = 0.0;
        /** keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  **/
        double[] squareLaggedSums = new double[MAX_LAG];
        double[] autoCorrelation = new double[MAX_LAG];
        for (int i = 0; i < (end - start); i++) {
            double traceI = trace[i + start];
            sum += traceI;
            // calculate mean
            final double mean = sum / (i + 1);

            // calculate auto correlation for selected lag times
            // sum1 = \sum_{start ... totalSamples-lag-1} trace
            double sum1 = sum;
            // sum2 = \sum_{start+lag ... totalSamples-1} trace
            double sum2 = sum;
            for (int lagIndex = 0; lagIndex < Math.min(i + 1, MAX_LAG); lagIndex++) {
                double traceIMinusLagIndex = trace[i + start - lagIndex];
                squareLaggedSums[lagIndex] = squareLaggedSums[lagIndex] + traceIMinusLagIndex * traceI;
                // The following line is the same approximation as in Tracer
                // (valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
                // though a more accurate estimate would be
                // autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2
                autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex);
                autoCorrelation[lagIndex] /= (i + start + 1 - lagIndex);
                sum1 -= traceIMinusLagIndex;
                sum2 -= trace[start + lagIndex];
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

        traceInfo.setUpMap(tracesInput.get());
    }

    /** @return a map of values to be logged */
    public Map<String, Double> getLogMap() {
        String[] traces = traceInfo.getTraceLabels();
        Map<String, Double> logValues = new HashMap<>(traces.length);
        for (int i = 0; i < traces.length; i++) {
            logValues.put(traces[i], currentESSs[i]);
        }
        return logValues;
    }

    public String[] getTraceLabels() {
        return traceInfo.getTraceLabels();
    }

}
