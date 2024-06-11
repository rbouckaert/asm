package asm.inference;

import beast.base.evolution.tree.Tree;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The collection of data from a set of trace logs and set of tree logs.
 *
 * @author Remco Bouckaert
 */
public class TraceInfo {

    public static DecimalFormat f = new DecimalFormat("#.###");
    public static DecimalFormat f1 = new DecimalFormat("#.#");

    /** tables of logs, one for each thread + one for the total; [chainIndex][attribute]{values} */
    List<Double>[][] logLines;

    /** column labels of logLines */
    String[] labels;

    /** tables of trees, one for each thread + one for the total */
    List<Tree>[] trees;

    DistanceMatrixCache distances;

    public TraceInfo(int chainCount) {
        f.setMinimumFractionDigits(3);
        f1.setMinimumFractionDigits(1);

        trees = new List[chainCount];
        for (int i = 0; i < chainCount; i++) {
            trees[i] = new ArrayList<>();
        }

        logLines = new List[chainCount][];
    }


    // maps "posterior", "likelihood" and "prior" columns to logLine columns
    private int[] map;
    private static int[] DEFAULT_MAP = {1, 2, 3};

    private String[] traces;

    String tracesString = null;

    public int[] getMap() {
        if (map != null) {
            return map;
        }
        if (labels != null) {
            setUpMap(tracesString);
            // map = new int[3];
            // map[0] = indexOf(labels, "posterior", 1);
            // map[1] = indexOf(labels, "likelihood", 2);
            // map[2] = indexOf(labels, "prior", 3);
            return map;
        }

        // default case
        this.traces = new String[]{"posterior", "likelihood", "prior"};
        return DEFAULT_MAP;
    }

    public void setUpMap(String tracesString) {
        if (labels == null) {
            this.tracesString = tracesString;
            return;
        }

        this.traces = tracesString.split(",");

        map = new int[traces.length];
        int k = 0;
        for (String trace : traces) {
            map[k] = indexOf(labels, trace.trim(), -1);
            if (map[k] == -1) {
                throw new IllegalArgumentException("Could not find label " + trace + " in trace log. "
                        + "Use one of " + Arrays.toString(labels));
            }
            k++;
        }
    }

    public String[] getTraceLabels() {
        if (this.traces == null) {
            getMap();
        }

        return this.traces;
    }

    private int indexOf(String[] labels, String string, int defaultIndex) {
        for (int i = 0; i < labels.length; i++) {
            if (labels[i].equals(string)) {
                return i;
            }
        }
        return defaultIndex;
    }

    public int chainCount() {
        return trees.length;
    }

    public List<Tree>[] getTrees() {
        return trees;
    }

    public void setTrees(List<List<Tree>> treesList) {
        for (int i = 0; i < treesList.size(); i++) {
            trees[i] = treesList.get(i);
        }
    }

    public void setLogs(List<Double>[][] logLines) {
        this.logLines = logLines;
    }

    public void setLabels(String[] traceLabels) {
        this.labels = traceLabels;
    }
}
