package asm.inference;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;



@Description("Runs multiple MCMC chains in parallel and reports "
		+ "1. Rubin-Gelman statistic while running the chain for each item in the first log file, "
		+ "2. The maximum difference in clade probability for every pair of chains, "
		+ "3. Tree distance stopping criterion. " +
		"Note that the log and tree log should have the same sample frequency.")
public class AutoStopMCMC extends MCMC {
	public Input<Integer> nrOfChainsInput = new Input<>("chains", "number of chains to run in parallel (default 2)", 2);
	public Input<List<MCMCConvergenceCriterion>> stoppingCriterionInput = new Input<>("stoppingCriterion", "one or more stopping criterion for tracking progress of the chains", new ArrayList<>());
	
	/** plugins representing MCMC with model, loggers, etc **/
	MCMCChain [] m_chains;
	
	/** threads for running MCMC chains **/
	Thread [] m_threads;
	
	/** keep track of time taken between logs to estimate speed **/
    long m_nStartLogTime;
	

    TraceInfo traceInfo;
    BurnInDetector burnInDetector;

	
	
	/** last line for which log is reported for all chains */
	int m_nLastReported = 0;
	
	
	
	List<MCMCConvergenceCriterion> stoppingCriteria;
	

	/** index of log and tree log among the MCMC loggers**/
	int m_iTreeLog = 0;
	int m_iLog = 0;

	
	private boolean slept = false;

	@Override
	public void initAndValidate() {
		m_chains = new MCMCChain[nrOfChainsInput.get()];
		stoppingCriteria = new ArrayList<>();
		stoppingCriteria.addAll(stoppingCriterionInput.get());
		if (stoppingCriteria.size() == 0) {
			Log.warning("=========================");
			Log.warning("=========================");
			Log.warning("No stoppingCriteria specified -- will run chain for chainLength steps without checking for convergence");
			Log.warning("=========================");
			Log.warning("=========================");
		}
		stoppingCriterionInput.get().clear();
		
		// the difference between the various chains is
		// 1. it runs an MCMC, not an AutoStopMCMC
		// 2. remove chains attribute
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("chains=[^ /]*", "");
		sXML = sXML.replaceAll("targetESS=[^ /]*", "");
		
		String sMultiMCMC = this.getClass().getName();
		while (sMultiMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b"+sMultiMCMC+"\\b", MCMCChain.class.getName());
			if (sMultiMCMC.indexOf('.') >= 0) {
				sMultiMCMC = sMultiMCMC.substring(sMultiMCMC.indexOf('.')+1);
			} else {
				sMultiMCMC = "";
			}
		}
		
		// create new chains
		XMLParser parser = new XMLParser();
		for (int i = 0; i < m_chains.length; i++) {
			String sXML2 = sXML;
			sXML2 = sXML2.replaceAll("fileName=\"", "fileName=\"chain" + i+ "-");
			if (sXML2.equals(sXML)) {
				// Uh oh, no loggers to file
				throw new IllegalArgumentException("Use file loggers, otherwise there are no trace a tree logs to track");
			}
			try {
				m_chains[i] = (MCMCChain) parser.parseFragment(sXML2, true);
			} catch (XMLParserException e) {
				throw new IllegalArgumentException(e);
			}
			// remove log to stdout, if any
			for (int iLogger = m_chains[i].loggersInput.get().size()-1; iLogger >= 0; iLogger--) {
				if (m_chains[i].loggersInput.get().get(iLogger).fileNameInput.get() == null) {
					m_chains[i].loggersInput.get().remove(iLogger);
				}
			}
		}
	
		// collect indices for tree log file names
		while (m_chains[0].loggersInput.get().get(m_iTreeLog).mode != Logger.LOGMODE.tree) {
			m_iTreeLog++;
		}
		while (m_chains[0].loggersInput.get().get(m_iLog).mode != Logger.LOGMODE.compound) {
			m_iLog++;
		}
		int nEveryLog = m_chains[0].loggersInput.get().get(m_iLog).everyInput.get();
		int nEveryTree = m_chains[0].loggersInput.get().get(m_iTreeLog).everyInput.get();
		if (nEveryLog != nEveryTree) {
			throw new IllegalArgumentException("log frequencey and tree log frequencey should be the same.");
		}
	} // initAndValidate
	
	@SuppressWarnings("unchecked")
	@Override 
	public void run() throws IOException {
		long start = System.currentTimeMillis();
		// memory for trees & tree distances
		traceInfo = new TraceInfo(m_chains.length);
		burnInDetector = new BurnInDetector(traceInfo);
		
		for (MCMCConvergenceCriterion stoppingCriterium : stoppingCriteria) {
			stoppingCriterium.setup(m_chains.length, traceInfo);
		}

		// start threads with individual chains here.
		m_threads = new Thread[m_chains.length];
		int k = 0;
		for (final MCMC mcmc : m_chains) {
			mcmc.setStateFile(stateFileName + "." +k, restoreFromFile);
			// need this to keep regression testing time reasonable
			mcmc.chainLengthInput.setValue(chainLengthInput.get(), this);
			m_threads[k] = new Thread() {
				public void run() {
					try {
						mcmc.run();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			};
			m_threads[k].start();
			k++;
		}

		new LogWatcherThread().start();
		
		// wait for the chains to finish
        m_nStartLogTime = System.currentTimeMillis();
		for (Thread thread : m_threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		// wait 2 seconds for the log to complete and close before exiting the process
		Log.warning("Wait for the log to complete");
		try {
			Thread.sleep(2000);
		} catch (InterruptedException e) {
			// ingore
		}
		
		combinelogs(burnInDetector.burnIn(m_nLastReported));
		
		
		long end = System.currentTimeMillis();
		Log.warning("All done in " + (end-start)/1000.0 + " seconds");
	} // run
	
	private void combinelogs(int[] burnIn) throws IOException {
		for (Logger logger: loggersInput.get()) {
			if (logger.fileNameInput.get() != null && logger.fileNameInput.get().length() > 0) {

				long sampleNr = 0;
				long sampleDelta = (long)(traceInfo.logLines[0][0].get(1) - traceInfo.logLines[0][0].get(0));
				logger.init();
				PrintStream out = logger.getM_out();
				if (logger.mode == Logger.LOGMODE.tree) {
					
					for (int i = 0; i < nrOfChainsInput.get(); i++) {
						// tree log
				        BufferedReader fin = new BufferedReader(new FileReader("chain" + i + "-" + logger.fileNameInput.get()));
				        String str = null;
				        
				        // skip burnin
				        int skipped = 0;
				        while (fin.ready() && skipped < burnIn[i]) {
				            str = fin.readLine();
				            if (str.matches("^tree STATE.*")) {
				            	skipped++;
				            }
				        }

				        // append trees to log
				        while (fin.ready()) {
				            str = fin.readLine();
				            if (str.matches("^tree STATE.*")) {
			                    str = str.replaceAll("^tree STATE_[^\\s]*", "");
			                	out.print("tree STATE_" + sampleNr + str);
			                	out.println();
								sampleNr += sampleDelta;
				            }
				        }
				        
				        fin.close();
					}					
				} else {
					
					// trace log
					for (int i = 0; i < nrOfChainsInput.get(); i++) {
						List<Double> [] logLines = traceInfo.logLines[i];
						int j = burnIn[i];
						while (j < logLines[0].size()) {
							out.print(sampleNr);
							out.print('\t');
							for (int k = 1; k < logLines.length; k++) {
								out.print(logLines[k].get(j));
								out.print('\t');
							}
							out.print('\n');
							sampleNr += sampleDelta;
							j++;
						}
					}
				}
				logger.close();
			}
		}
	}

	/** Represents class that tails all logs files and tree log files.
	 * When a new line is added, this is processed */ 
	class LogWatcherThread extends Thread {
		@Override
		public void run() {
			try {
				int nThreads = m_chains.length;
				/* file handle pairs; two for each thread, 
				 * even numbered files are for the log file, odd numbers for the tree file */
				BufferedReader [] fin = new BufferedReader[nThreads*2];
				int nFilesOpened = 0;

				// wait a seconds, the log file should be available
				sleep(1000);
				// open files
				while (nFilesOpened < nThreads*2) {
					for (int i = 0; i < nThreads*2; i++) {
						if (fin[i] == null) {
							String sFileName = m_chains[i/2].loggersInput.get().get(i%2==0 ? m_iLog : m_iTreeLog).fileNameInput.get(); 
							fin[i] = new BufferedReader(new FileReader(sFileName));
							if (fin[i] != null) {
								nFilesOpened++;
							}
						}
					}
					if (nFilesOpened < nThreads*2) {
						sleep(1000);
					}
				}

				// keep polling the tree logs file every second
				while (true) {
					int nLinesRead = 0;
					// grab a tree from every thread
					slept = false;
					while (nLinesRead < nThreads*2) {
						boolean [] bDone = new boolean[nThreads*2];
						for (int i = 0; i < nThreads*2; i++) {
							if (!bDone[i]) {
								boolean bLineRead = (i%2==0?readLogLines(i/2, fin[i]) : readTreeLogLine(i/2, fin[i]));
								if (bLineRead) {
									nLinesRead++;
									bDone[i] = true;
								}
							}
						}
						if (nLinesRead< nThreads*2) {
							// wait a second before seeing if there is more
							sleep(1000);
							slept = true;
						}
					}
					
					if (!slept && m_nLastReported > 1) {
						// Log.warning("Checks cannot keep up with log files: expect ESSs higher than targetESS");
					}
					if (check()) {
						return;
					}
					m_nLastReported++;	
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		private boolean check() {
			Log.info.print("Check " + m_nLastReported + " ");
			
			boolean converged = true;
			long start = System.currentTimeMillis();
			int [] burnin = burnInDetector.burnIn(m_nLastReported);
			Log.info.print("burnin:" + Arrays.toString(burnin));
			for (MCMCConvergenceCriterion crit : stoppingCriteria) {
				if (!crit.converged(burnin, m_nLastReported)) {
					converged = false;
					// break;
				}
			}
			long end = System.currentTimeMillis();					
			//Log.info.println((converged?"Succss!":"failed") + " in " + (end-start) + " mseconds");
			Log.info.println(" in " + (end-start) + " mseconds");
			if (converged) {
				// stop all threads
				for (MCMCChain t : m_chains) {
					t.terminate();
				}
				return true;
			}
			return false;
		}	
	} // class LogWatcherThread
	
	/** read a line from the log, return true if successfull */
	boolean readLogLines(int iThread, BufferedReader fin) {
		try {
			String sStr = null;
			String [] sStrs = null;;
			do {
				sStr = fin.readLine();
				if (sStr == null) {
					return false;
				}
				sStrs = sStr.split("\\s+");
			} while (sStr.startsWith(("#")) || sStrs.length == 1); // ignore comment lines
			int nItems = sStrs.length;
			
			if (sStr.startsWith("Sample")) {
				traceInfo.labels = sStr.split("\\s+");
				return false;
			}

			if (traceInfo.logLines[0] == null) {
				for (int i = 0; i < m_chains.length; i++) {
					traceInfo.logLines[i] = new List[nItems];
				}
				for (int i = 0; i < m_chains.length; i++) {
					for (int j = 0; j < nItems; j++) {
						traceInfo.logLines[i][j] = new ArrayList<>();
					}
				}
			}
			
			try {
				for (int i = 0; i < nItems; i++) {
					try {
						traceInfo.logLines[iThread][i].add(Double.parseDouble(sStrs[i]));
					} catch (NumberFormatException e) {
						traceInfo.logLines[iThread][i].add(0.0);
					}
				}
			} catch (Exception e) {
				//ignore, probably a parse errors
				if (iThread == 0) {
					System.out.println(sStr);
				}
				return false;
			}
			return true;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	} // readLogLine



	/** read a single tree from the tree log file, return true if successful **/
	boolean readTreeLogLine(int iThread, BufferedReader fin) {
		String sStr = null;
		do {
			try {
				sStr = fin.readLine();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			if (sStr == null) {
				try {
					// wait 2.5 seconds till tree becomes available
					Thread.sleep(1000);
					slept = true;
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		} while (sStr == null || !sStr.matches("tree STATE.*")); // ignore non-tree lines

		sStr = sStr.substring(sStr.indexOf("("));
		TreeParser parser = new TreeParser();
		parser.offsetInput.setValue(1, parser);
		Node root = parser.parseNewick(sStr);
		parser.setRoot(root);
		Tree tree = parser.copy();
		traceInfo.trees[iThread].add(tree);
		return true;
	} // readTreeLogLine


	

	
} // class AutoStopMCMC



