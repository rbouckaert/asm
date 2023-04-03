package asm.inference;

import java.util.Arrays;

import beast.base.core.Log;

public class DistanceMatrixCache {
	// symmetric 2d distance matrix for trees 1
	private float [] cache11;
	// symmetric 2d distance matrix for trees 2
	private float [] cache22;
	// asymmetric 2d distance matrix between trees 1 and trees 2
	// represented by lower triangle and upper triangle half-matrices
	private float [] cache12, cache21, diagonal;
	// matrix size
	private int size;
	
	public DistanceMatrixCache(int n) {
		cache11 = new float[n*(n-1)/2];
		cache22 = new float[n*(n-1)/2];
		cache12 = new float[n*(n-1)/2];
		cache21 = new float[n*(n-1)/2];
		diagonal = new float[n];
	}
	
	float getDistance(int treeSet1, int index1, int treeSet2, int index2) {
		if (index1 >= size || index2 >= size) {
			// resize
			Log.warning("Resizing cache from " + size + " to " + (size+1024));
			size += 1024;
			cache11 = Arrays.copyOf(cache11, size*(size-1)/2);
			cache22 = Arrays.copyOf(cache22, size*(size-1)/2);
			cache12 = Arrays.copyOf(cache12, size*(size-1)/2);
			cache21 = Arrays.copyOf(cache21, size*(size-1)/2);
			diagonal = Arrays.copyOf(diagonal, size);
		}
		
		if (treeSet1 == 0) {
			if (treeSet2 == 0) {
				int i = index1 > index2 ?  
					index1 * (index1 - 1)/2 + index2:
					index2 * (index2 - 1)/2 + index1;
				return cache11[i];
			} else {
				if (index1 > index2) {  
					int i = index1 * (index1 - 1)/2 + index2;
					return cache12[i];
				} else if (index1 < index2) {  
					int i = index2 * (index2 - 1)/2 + index1;
					return cache21[i];
				} else {
					return diagonal[index1];
				}
			}
		} else {
			if (treeSet2 == 0) {
				if (index1 > index2) {  
					int i = index1 * (index1 - 1)/2 + index2;
					return cache21[i];
				} else if (index1 < index2) {
					int i = index2 * (index2 - 1)/2 + index1;
					return cache12[i];
				} else {
					return diagonal[index1];
				}
			} else {
				int i = index1 > index2 ?  
						index1 * (index1 - 1)/2 + index2:
						index2 * (index2 - 1)/2 + index1;
				return cache22[i];
			}
		}			
	}

	void setDistance(int treeSet1, int index1, int treeSet2, int index2, float d) {
		if (treeSet1 == 0) {
			if (treeSet2 == 0) {
				int i = index1 > index2 ?  
					index1 * (index1 - 1)/2 + index2:
					index2 * (index2 - 1)/2 + index1;
				cache11[i] = d;
			} else {
				if (index1 > index2) {  
					int i = index1 * (index1 - 1)/2 + index2;
					cache12[i] = d;
				} else if (index1 < index2) {
					int i = index2 * (index2 - 1)/2 + index1;
					cache21[i] = d;
				} else {
					diagonal[index1] = d;
				}
			}
		} else {
			if (treeSet2 == 0) {
				if (index1 > index2) {  
					int i = index1 * (index1 - 1)/2 + index2;
					cache21[i] = d;
				} else if (index1 < index2) {
					int i = index2 * (index2 - 1)/2 + index1;
					cache12[i] = d;
				} else {
					diagonal[index1] = d;
				}
			} else {
				int i = index1 > index2 ?  
						index1 * (index1 - 1)/2 + index2:
						index2 * (index2 - 1)/2 + index1;
				cache22[i] = d;
			}
		}			
		
	}
}