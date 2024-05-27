package asm.tools;

import asm.inference.AutoStopMCMC;
import asm.inference.MCMCChain;
import beast.base.core.Description;
import beast.base.inference.Logger;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;

@Description("Works as AutoStopMCMC but with emulated MCMC chains")
public class EmulatedAutoStopMCMC extends AutoStopMCMC {

	private String inputLogDir, outputLogDir;

	@Override
	protected void initChains() {
		// the difference between the various chains is
		// 1. it runs an MCMC, not an AutoStopMCMC
		// 2. remove chains attribute
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("chains=[^ /]*", "");
		sXML = sXML.replaceAll("targetESS=[^ /]*", "");
		sXML = sXML.replaceAll("burnInStrat=[^ /]*", "");
		sXML = sXML.replaceAll("burnInPercent=[^ /]*", "");

		String sMultiMCMC = this.getClass().getName();
		while (sMultiMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b"+sMultiMCMC+"\\b", EmulatedMCMCChain.class.getName());
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
			((EmulatedMCMCChain)m_chains[i]).setLogDirs(inputLogDir, outputLogDir);
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
			throw new IllegalArgumentException("log frequency and tree log frequency should be the same.");
		}
	}

	protected void setLogDirs(String in, String out) {
		inputLogDir = in;
		outputLogDir = out;
	}
}
