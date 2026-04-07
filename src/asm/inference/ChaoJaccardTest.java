package asm.inference;

import java.util.HashMap;
import java.util.Map;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class ChaoJaccardTest {

    private static final double EPS = 1e-10;

    @Test
    void identicalChains() {
        // two chains with identical feature counts -> J = 1.0
        Map<String, Integer> map = Map.of("0,1", 50, "0,1,2", 50);
        // each tree has 2 features, 50 trees -> nPlus = 100
        double j = ChaoJaccard.chaoJaccardIncidence(map, 100, 50, map, 100, 50);
        assertEquals(1.0, j, EPS);
    }

    @Test
    void disjointChains() {
        // no shared features -> J = 0.0
        Map<String, Integer> mapA = Map.of("0,1", 10);
        Map<String, Integer> mapB = Map.of("2,3", 10);
        double j = ChaoJaccard.chaoJaccardIncidence(mapA, 10, 10, mapB, 10, 10);
        assertEquals(0.0, j, EPS);
    }

    @Test
    void emptyChain() {
        Map<String, Integer> mapA = Map.of("0,1", 10);
        Map<String, Integer> mapB = new HashMap<>();
        assertEquals(0.0, ChaoJaccard.chaoJaccardIncidence(mapA, 10, 10, mapB, 0, 0), EPS);
        assertEquals(0.0, ChaoJaccard.chaoJaccardIncidence(mapB, 0, 0, mapA, 10, 10), EPS);
    }

    @Test
    void symmetry() {
        // J(A,B) == J(B,A)
        Map<String, Integer> mapA = Map.of("0,1", 10, "0,2", 5, "1,2", 1);
        Map<String, Integer> mapB = Map.of("0,1", 8, "0,2", 1, "1,3", 3);
        double jAB = ChaoJaccard.chaoJaccardIncidence(mapA, 16, 10, mapB, 12, 8);
        double jBA = ChaoJaccard.chaoJaccardIncidence(mapB, 12, 8, mapA, 16, 10);
        assertEquals(jAB, jBA, EPS);
    }

    @Test
    void partialOverlap() {
        // hand-computed example
        // mapA: "0,1"->10, "0,2"->5, "1,2"->1; tA=10, nPlus=16
        // mapB: "0,1"->8,  "0,2"->1, "1,3"->3;  tB=8,  mPlus=12
        //
        // shared = {"0,1", "0,2"}
        // uObs = 10/16 + 5/16 = 15/16
        // vObs = 8/12 + 1/12 = 9/12 = 3/4
        //
        // f1plus=0, f2plus=0->1, fPlus1=1 ("0,2" has Y=1), fPlus2=0->1
        // sumXgivenY1 = 5/16 (only "0,2" has Y=1)
        // sumYgivenX1 = 0 (no shared feature has X=1)
        //
        // uHat = min(1, 15/16 + (7/8)*(1/2)*(5/16)) = min(1, 1.074..) = 1.0
        // vHat = 3/4 + 0 = 0.75
        //
        // J = 1.0 * 0.75 / (1.0 + 0.75 - 0.75) = 0.75
        Map<String, Integer> mapA = Map.of("0,1", 10, "0,2", 5, "1,2", 1);
        Map<String, Integer> mapB = Map.of("0,1", 8, "0,2", 1, "1,3", 3);
        double j = ChaoJaccard.chaoJaccardIncidence(mapA, 16, 10, mapB, 12, 8);
        assertEquals(0.75, j, EPS);
    }

    @Test
    void biasCorrection() {
        // singletons in shared set should trigger upward correction
        // Chain A: 100 trees, features "a"->100, "b"->100, "c"->1
        // Chain B: 100 trees, features "a"->100, "b"->1, "c"->100
        //
        // tA=tB=100, nPlus=mPlus=201
        // shared = {"a", "b", "c"}
        //
        // uObs = (100+100+1)/201 = 201/201 = 1.0
        // vObs = (100+1+100)/201 = 201/201 = 1.0
        //
        // No correction needed since uObs = vObs = 1.0 already (all features shared)
        // J = 1.0
        Map<String, Integer> mapA = Map.of("a", 100, "b", 100, "c", 1);
        Map<String, Integer> mapB = Map.of("a", 100, "b", 1, "c", 100);
        double j = ChaoJaccard.chaoJaccardIncidence(mapA, 201, 100, mapB, 201, 100);
        assertEquals(1.0, j, EPS);
    }

    @Test
    void biasCorrectionWithUnsharedFeatures() {
        // shared features have singletons; unshared features exist
        // Chain A: "a"->50, "b"->1, "x"->10; tA=50, nPlus=61
        // Chain B: "a"->50, "b"->1, "y"->10; tB=50, mPlus=61
        //
        // shared = {"a", "b"}
        // uObs = (50+1)/61
        // vObs = (50+1)/61
        //
        // f1plus=1 ("b"), f2plus=0->1
        // fPlus1=1 ("b"), fPlus2=0->1
        // sumXgivenY1 = 1/61 ("b" has Y=1, X=1)
        // sumYgivenX1 = 1/61 ("b" has X=1, Y=1)
        //
        // uHat = min(1, 51/61 + (49/50)*(1/2)*(1/61))
        //      = min(1, 0.83607 + 0.00803) = 0.84410
        // vHat = same by symmetry = 0.84410
        //
        // J = 0.84410^2 / (2*0.84410 - 0.84410^2)
        Map<String, Integer> mapA = Map.of("a", 50, "b", 1, "x", 10);
        Map<String, Integer> mapB = Map.of("a", 50, "b", 1, "y", 10);
        double j = ChaoJaccard.chaoJaccardIncidence(mapA, 61, 50, mapB, 61, 50);

        // verify symmetry holds for this case
        double uObs = 51.0 / 61;
        double correction = (49.0 / 50) * (1.0 / 2) * (1.0 / 61);
        double uHat = uObs + correction;
        double expected = (uHat * uHat) / (2 * uHat - uHat * uHat);
        assertEquals(expected, j, EPS);

        // J should be less than 1 because of unshared features
        assertTrue(j < 1.0);
        // but greater than naive Jaccard (|shared|/|union| = 2/4 = 0.5)
        assertTrue(j > 0.5);
    }

    @Test
    void convergesWithMoreData() {
        // as chains accumulate more samples of the same features, J should increase
        // simulate two chains with mostly shared features and one unique each
        double jSmall = computeWithSampleSize(10);
        double jMedium = computeWithSampleSize(100);
        double jLarge = computeWithSampleSize(1000);

        assertTrue(jSmall <= jMedium, "J should not decrease with more data");
        assertTrue(jMedium <= jLarge, "J should not decrease with more data");
    }

    private double computeWithSampleSize(int n) {
        // shared features scale with n, unique features stay at 1
        Map<String, Integer> mapA = new HashMap<>();
        mapA.put("0,1", n);
        mapA.put("0,2", n);
        mapA.put("1,2", n);
        mapA.put("only_a", 1);

        Map<String, Integer> mapB = new HashMap<>();
        mapB.put("0,1", n);
        mapB.put("0,2", n);
        mapB.put("1,2", n);
        mapB.put("only_b", 1);

        int nPlus = 3 * n + 1;
        return ChaoJaccard.chaoJaccardIncidence(mapA, nPlus, n, mapB, nPlus, n);
    }

    @Test
    void perfectOverlapDifferentCounts() {
        // same feature set, different counts -> J = 1.0
        // (all features shared, so J should be 1 regardless of count differences)
        Map<String, Integer> mapA = Map.of("0,1", 90, "0,2", 10);
        Map<String, Integer> mapB = Map.of("0,1", 10, "0,2", 90);
        double j = ChaoJaccard.chaoJaccardIncidence(mapA, 100, 50, mapB, 100, 50);
        assertEquals(1.0, j, EPS);
    }
}
