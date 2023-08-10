package asm.tools;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import asm.inference.TreePSRF;
import asm.inference.TraceInfo;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.FastTreeSet;
import beastfx.app.util.TreeFile;

@Description("Calculate pairwise Tree PSRF statistics for tree sets")
public class TreeConvergenceCheck extends Runnable {

	final public Input<List<TreeFile>> treesInput = new Input<>("tree", "tree files to be compared. Need at least two.", new ArrayList<>());
	public Input<Double> smoothingInput = new Input<>("smoothing", "smoothing factor, which determines how proportion of trees to disregard: "
			+ "larger smoothing means more trees included in test", 0.6);
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	
	
	
	public Input<Double> bInput = new Input<>("b", "threshold determining acceptance bounds for PSFR like statistic", 0.05);
	public Input<Boolean> twoSidedInput = new Input<>("twoSided", "if true (default) psrfs and pseudo ESSs for both chains are used,"
			+ "otherwise only one psrfs and pseudo ESS is calculated, which takes less computation but can be less robust", true);
	public Input<Boolean> checkESSInput = new Input<>("checkESS", "whether to check Tree ESS exceeds the targetESS", false);
	public Input<Integer> targetESSInput = new Input<>("targetESS", "target effective sample size per chain (default 100)", 100);
	
	
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		List<TreeFile> treeFiles = treesInput.get();
		int n = treeFiles.size();
		if (n < 2) {
			Log.warning("At least two trees should be specified");
			return;
		}

		
		for (int i = 0; i < n-1; i++) {
			List<Tree> trees1 = getTrees(i);
			for (int j = i + 1; j < n; j++) {
				List<Tree> trees2 = getTrees(j);
				TraceInfo traceInfo = new TraceInfo(2);
				traceInfo.setTrees(trees1, trees2);
				TreePSRF treePSRF = new TreePSRF();
				treePSRF.initByName(
						"smoothing", smoothingInput.get(), 
						"b", bInput.get(), 
						"twoSided", twoSidedInput.get(), 
						"checkESS", checkESSInput.get(),
						"targetESS", targetESSInput.get());
				treePSRF.setup(2, traceInfo);
				int end = Math.min(trees1.size(), trees2.size());
				int start = burnInPercentageInput.get() * end / 100; 
				Log.info("");
				Log.info("Comparing " + treeFiles.get(i).getName() + " and " + treeFiles.get(j).getName());
				boolean converged = treePSRF.converged(new int[] {start, start}, end);
				Log.info("");
				Log.info((converged ? "Converged" : "Not converged" ) + " according to the GRLike criterion");
			}
		}
		
		Log.warning("Done");
	}

	private PrintStream nullstream = new PrintStream(new OutputStream() {
		@Override
		public void write(int b) throws IOException {
		}
		@Override
		public void write(byte[] b) throws IOException {
		}
		@Override
		public void write(byte[] b, int off, int len) throws IOException {
		}
    });

    private List<Tree> getTrees(int i) throws IOException {
		PrintStream originalErr = System.err;
	    System.setErr(nullstream);
	    
	    try {
			String treeFile1 = treesInput.get().get(i).getPath();
			FastTreeSet treeset1 = new TreeAnnotator().new FastTreeSet(treeFile1, 0);
			List<Tree>  trees1 = new ArrayList<>();
			while (treeset1.hasNext()) {
				trees1.add(treeset1.next());
			}
		    System.setErr(originalErr);
			return trees1;
		} catch (IOException e) {
		    System.setErr(originalErr);
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}
	}

	public static void main(String[] args) throws Exception {
		new Application(new TreeConvergenceCheck(), "Tree Convergence Check", 512, 700, args);
	}

}
