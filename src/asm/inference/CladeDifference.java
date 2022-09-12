package asm.inference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Pairwise convergence criterion based on the difference in clade support between the chains")
public class CladeDifference extends BEASTObject implements PairewiseConvergenceCriterion {
	public Input<Double> acceptedThresholdInput = new Input<>("threshold", "level at which the biggest clade support difference is still acceptable", 0.25);

	/** maximum difference of clade probabilities for chain 1 & 2 **/
	List<Double> m_fMaxCladeProbDiffs;
	
	/** for each thread, keeps track of the frequency of clades **/
	Map<String, Integer> [] m_cladeMaps;
	
	/** total nr of clades in a tree **/
	int m_nClades = 1;
	
	int nChains;
	double acceptedThreshold;
	
	@Override
	public void initAndValidate() {
		acceptedThreshold = acceptedThresholdInput.get();

	}

	@Override
	public boolean converged() {
		double fMaxCladeProbDiff = 0;
		for (int i = 0; i < nChains; i++) {
			for (int k = 0; k < nChains; k++) {
				fMaxCladeProbDiff = Math.max(fMaxCladeProbDiff, calcMaxCladeDifference(i, k));
			}
		}
		m_fMaxCladeProbDiffs.add(fMaxCladeProbDiff);
		return fMaxCladeProbDiff < acceptedThreshold;
	}

	@Override
	public void setup(int nChains) {
		this.nChains = nChains;
		m_fMaxCladeProbDiffs = new ArrayList<>();
		m_cladeMaps = new Map[nChains];
		for (int i = 0; i < nChains; i++) {
			m_cladeMaps[i] = new HashMap<>();
		}

	}
	
	/** calculate maximum difference of clade probabilities 
	 * of 2 threads. It uses only the clades in thread1, so
	 * to it requires two calls to make sure no clades in thread2
	 * are missed, i.e. use max(calcMaxCladeDifference(iThread1, iThread2), calcMaxCladeDifference(iThread2, iThread2)) **/
	/** TODO: can be done incrementally?!? **/
	double calcMaxCladeDifference(int iThread1, int iThread2) {
		if (iThread1 == iThread2) {
			return 0;
		}
		int nTotal = 0;
		int nMax = 0;
		Map<String, Integer> map1 = m_cladeMaps[iThread1];
		Map<String, Integer> map2 = m_cladeMaps[iThread2];
		for (String sClade : map1.keySet()) {
			int i1 = map1.get(sClade);
			int i2 = 0;
			if (map2.containsKey(sClade)) {
				i2 = map2.get(sClade);
			}
			nTotal += i1;
			nMax = Math.max(nMax, Math.abs(i1 - i2));
		}
		return nMax / (double) (nTotal/m_nClades);
	} // calcMaxCladeDifference


	/** get clades from tree and store them in a list in String format **/
	int [] traverse(Node node, List<String> sClades) {
		int [] clade = null;
		if (node.isLeaf()) {
			clade = new int[1];
			clade[0] = node.getNr();
		} else {
			int [] leftClade = traverse(node.getLeft(), sClades);
			int [] rightClade = traverse(node.getRight(), sClades);
			
			// merge clade with rightClade
			clade = new int[leftClade.length + rightClade.length];
			int i = 0, iLeft = 0, iRight = 0;
			while (i < clade.length) {
				if (leftClade[iLeft] < rightClade[iRight]) {
					clade[i] = leftClade[iLeft++];
					if (iLeft == leftClade.length) {
						while (iRight < rightClade.length) {
							clade[++i] = rightClade[iRight++];
						}
					}
				} else {
					clade[i] = rightClade[iRight++]; 
					if (iRight == rightClade.length) {
						while (iLeft < leftClade.length) {
							clade[++i] = leftClade[iLeft++];
						}
					}
				}
				i++;
			}
			String sClade = "";
			for (i = 0; i < clade.length; i++) {
				sClade += clade[i] + ",";
			}
			sClades.add(sClade);
		}
		return clade;
	}
	

	@Override
	public void process(int chainNr, Double [] logLine, Node root) {
		List<String> sClades = new ArrayList<String>();
		m_nClades = traverse(root, sClades).length;
		Map<String, Integer> cladeMap = m_cladeMaps[chainNr];
		for (String sClade : sClades) {
			if (cladeMap.containsKey(sClade)) {
				cladeMap.put(sClade, cladeMap.get(sClade) + 1);
			} else {
				cladeMap.put(sClade, 1);
			}
		}

	}
}
