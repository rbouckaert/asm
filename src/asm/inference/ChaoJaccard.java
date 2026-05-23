package asm.inference;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.*;

@Description("Convergence criterion based on incidence-based Chao-Jaccard similarity " +
        "of features between chains (Chao et al. 2005, Ecology Letters 8:148-159). " +
        "Treats each tree as a sampling unit and each clade (or split) as a species. " +
        "Estimates true Jaccard similarity accounting for unseen shared features.")
public class ChaoJaccard extends BEASTObject implements MCMCConvergenceCriterion {

    // For logging and accessing logs
    public static final class Keys {

        private Keys() {
        }

        public static final String CHAIN_MIN = "Chao-Chain-min";
        public static final String PAIR = "Chao-pair";
        public static final String GLOBAL_MIN = "Chao-global-min";
    }


    public Input<Double> thresholdInput = new Input<>("threshold",
            "Chao-Jaccard similarity above which chains are considered converged", 0.99);

    public Input<String> featureLevelInput = new Input<>("featureLevel",
            "Feature level: 'clade' (CCD0), 'split' (CCD1)", "clade");

    private double threshold;
    private int nChains;
    private List<Tree>[] trees;

    /**
     * feature incidence maps: feature string -> number of trees containing it
     */
    private Map<String, Integer>[] featureMaps;

    /**
     * number of trees processed per chain (T_A, T_B)
     */
    private int[] numTrees;

    /**
     * total feature incidences per chain (n_+, m_+)
     */
    private int[] totalIncidence;

    private boolean useSplits;
    private int current = 0;

    // Logging tools
    private double[][] pairwiseSimilarity;
    private double[] chainMinSimilarity;
    private double globalMinSimilarity;
    public FeatureLevel featureLevel;

    public enum FeatureLevel {
        CLADE,
        SPLIT
    }

    @Override
    public void initAndValidate() {
        threshold = thresholdInput.get();
        String levelInput = featureLevelInput.get().toUpperCase();

        try {
            this.featureLevel = FeatureLevel.valueOf(levelInput);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(
                    "Unkown featureLevel: " + featureLevelInput.get() +
                            ". Expected clade or split"
            );
        }

        useSplits = (featureLevel == FeatureLevel.SPLIT);
    }

    @Override
    public void setup(int nChains, TraceInfo traceInfo) {
        this.nChains = nChains;

        pairwiseSimilarity = new double[nChains][nChains];
        chainMinSimilarity = new double[nChains];
        globalMinSimilarity = Double.POSITIVE_INFINITY;

        this.trees = traceInfo.trees;
        this.featureMaps = new Map[nChains];
        this.numTrees = new int[nChains];
        this.totalIncidence = new int[nChains];
        for (int i = 0; i < nChains; i++) {
            featureMaps[i] = new HashMap<>();
        }
    }

    @Override
    public boolean converged(int[] burnin, int available) {
        // process new trees incrementally
        for (int i = current; i < available; i++) {
            for (int j = 0; j < nChains; j++) {
                List<String> features = new ArrayList<>();
                traverse(trees[j].get(i).getRoot(), features);
                Map<String, Integer> map = featureMaps[j];
                for (String feature : features) {
                    map.merge(feature, 1, Integer::sum);
                }
                numTrees[j]++;
                totalIncidence[j] += features.size();
            }
        }
        current = available;

        // compute minimum pairwise Chao-Jaccard similarity
        // 1. compute pairwise
        for (int i = 0; i < nChains; i++) {
            for (int j = i + 1; j < nChains; j++) {

                double cj = chaoJaccardIncidence(
                        featureMaps[i], totalIncidence[i], numTrees[i],
                        featureMaps[j], totalIncidence[j], numTrees[j]);

                pairwiseSimilarity[i][j] = cj;
                pairwiseSimilarity[j][i] = cj;

                globalMinSimilarity = Math.min(globalMinSimilarity, cj);
            }
        }

        // 2. compute per-chain minima
        for (int i = 0; i < nChains; i++) {
            chainMinSimilarity[i] = Double.POSITIVE_INFINITY;

            for (int j = 0; j < nChains; j++) {
                if (i != j) {
                    chainMinSimilarity[i] =
                            Math.min(chainMinSimilarity[i], pairwiseSimilarity[i][j]);
                }
            }
        }
        return globalMinSimilarity >= threshold;
    }

    // TODO implement logger for when running xmls, WIP, see dissonance for example
//    public Map<String, Double> getLogMap() {
//        Map<String, Double> log = new HashMap<>();
//
//        // per-chain minima (like entropy per chain)
//        for (int i = 0; i < nChains; i++) {
//            log.put(Keys.CHAIN_MIN + i, chainMinSimilarity[i]);
//        }
//
//        // pairwise values (optional but powerful)
//        for (int i = 0; i < nChains; i++) {
//            for (int j = i + 1; j < nChains; j++) {
//                log.put(Keys.PAIR + i + "-" + j, pairwiseSimilarity[i][j]);
//            }
//        }
//
//        // global summary
//        log.put(Keys.GLOBAL_MIN, globalMinSimilarity);
//
//        return log;
//    }

    public double getGlobalMinSimilarity() {
        return globalMinSimilarity;
    }

    public double getChainMinSimilarity(int i) {
        return chainMinSimilarity[i];
    }

    /**
     * Compute the incidence-based Chao-Jaccard similarity estimator
     * (Chao et al. 2005, Ecology Letters 8:148-159, equations 11-13).
     * <p>
     * Each tree is a sampling unit, each feature (clade or split) is a
     * "species". X_i is the number of trees in run A containing feature i.
     *
     * @param mapA  feature incidence counts for run A
     * @param nPlus total feature incidences in A (= T_A * features per tree)
     * @param tA    number of trees in run A
     * @param mapB  feature incidence counts for run B
     * @param mPlus total feature incidences in B (= T_B * features per tree)
     * @param tB    number of trees in run B
     * @return estimated Jaccard similarity
     */
    static double chaoJaccardIncidence(Map<String, Integer> mapA, int nPlus, int tA,
                                       Map<String, Integer> mapB, int mPlus, int tB) {
        if (nPlus == 0 || mPlus == 0) {
            return 0.0;
        }

        // find shared features (D_12)
        Set<String> shared = new HashSet<>(mapA.keySet());
        shared.retainAll(mapB.keySet());

        if (shared.isEmpty()) {
            return 0.0;
        }

        // observed relative incidence sums for shared features
        double uObs = 0.0;
        double vObs = 0.0;

        // singleton/doubleton counts among shared features
        int f1plus = 0; // shared features with X_i = 1
        int f2plus = 0; // shared features with X_i = 2
        int fPlus1 = 0; // shared features with Y_i = 1
        int fPlus2 = 0; // shared features with Y_i = 2

        // conditional incidence sums
        double sumXgivenY1 = 0.0; // sum of X_i/n_+ for shared features where Y_i = 1
        double sumYgivenX1 = 0.0; // sum of Y_i/m_+ for shared features where X_i = 1

        for (String feature : shared) {
            int x = mapA.get(feature);
            int y = mapB.get(feature);

            uObs += (double) x / nPlus;
            vObs += (double) y / mPlus;

            if (x == 1) f1plus++;
            if (x == 2) f2plus++;
            if (y == 1) fPlus1++;
            if (y == 2) fPlus2++;

            if (y == 1) sumXgivenY1 += (double) x / nPlus;
            if (x == 1) sumYgivenX1 += (double) y / mPlus;
        }

        // bias correction: if no doubletons, use 1 to avoid division by zero
        if (f2plus == 0) f2plus = 1;
        if (fPlus2 == 0) fPlus2 = 1;

        // incidence-based estimators (equations 11-12), capped at 1
        double uHat = Math.min(1.0,
                uObs + ((double) (tB - 1) / tB) * ((double) fPlus1 / (2 * fPlus2)) * sumXgivenY1);
        double vHat = Math.min(1.0,
                vObs + ((double) (tA - 1) / tA) * ((double) f1plus / (2 * f2plus)) * sumYgivenX1);

        // Chao-Jaccard similarity (equation 13)
        double denom = uHat + vHat - uHat * vHat;
        if (denom == 0.0) {
            return 0.0;
        }
        return (uHat * vHat) / denom;
    }

    /**
     * Extract features from a tree. At the clade level, each internal node
     * defines a clade (sorted set of descendant taxa). At the split level,
     * each internal node defines a parent-clade -> child-clade split.
     */
    private int[] traverse(Node node, List<String> features) {
        if (node.isLeaf()) {
            return new int[]{node.getNr()};
        }

        int[] left = traverse(node.getLeft(), features);
        int[] right = traverse(node.getRight(), features);

        // merge-sort the two child arrays
        int[] merged = new int[left.length + right.length];
        int i = 0, l = 0, r = 0;
        while (l < left.length && r < right.length) {
            if (left[l] < right[r]) {
                merged[i++] = left[l++];
            } else {
                merged[i++] = right[r++];
            }
        }
        while (l < left.length) merged[i++] = left[l++];
        while (r < right.length) merged[i++] = right[r++];

        if (useSplits) {
            // split feature: "parent|child" where child is the smaller partition
            String parent = arrayToString(merged);
            String leftStr = arrayToString(left);
            String rightStr = arrayToString(right);
            // use the smaller child as the split identifier (canonical form)
            String child = (left.length <= right.length) ? leftStr : rightStr;
            features.add(parent + "|" + child);
        } else {
            // clade feature
            features.add(arrayToString(merged));
        }

        return merged;
    }

    private static String arrayToString(int[] arr) {
        StringBuilder sb = new StringBuilder();
        for (int k = 0; k < arr.length; k++) {
            if (k > 0) sb.append(',');
            sb.append(arr[k]);
        }
        return sb.toString();
    }
}
