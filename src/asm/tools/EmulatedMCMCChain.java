package asm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import asm.inference.MCMCChain;
import beast.base.core.Description;
import beast.base.inference.Logger;

@Description("Emulated chain that does not actually run MCMC, but gets copies from the log files")
public class EmulatedMCMCChain extends MCMCChain {

	private List<BufferedReader> inputLoggers;
	private List<PrintStream> outputLoggers;
	private String inputLogDir, outputLogDir;
	private int every;
	
	protected void setLogDirs(String in, String out) {
		inputLogDir = in;
		outputLogDir = out;
	}
	
	
	@Override
	public void run() throws IOException, SAXException, ParserConfigurationException {
        if (restoreFromFile) {
        	throw new RuntimeException("Cannot use 'resume' with emulator");
        }

        // initialise loggers
        loggers = loggersInput.get();
        inputLoggers = new ArrayList<>();
        outputLoggers = new ArrayList<>();
        every = -1;
        for (final Logger log : loggers) {
        	if (every < 0) {
        		every = log.everyInput.get();
        	} else {
        		if (every != log.everyInput.get()) {
        			throw new IllegalArgumentException("loggers must have the same log frequency");
        		}
        	}
        	
        	String fn = log.fileNameInput.get();
        	PrintStream out = new PrintStream(outputLogDir + "/" + fn);
        	outputLoggers.add(out);
            
            // get header info from previous run
            BufferedReader fin = new BufferedReader(new FileReader(inputLogDir + "/" + fn));
            String str = fin.readLine();
            out.println(str);
            if (str.startsWith("#NEXUS")) {
            	// tree file
                while (fin.ready() && (str == null || !str.startsWith("tree"))) {
                    str = fin.readLine();
                    out.println(str);
                }
            } else {
                while (fin.ready() && (str == null || str.startsWith("#"))) {
                    str = fin.readLine();
                    out.println(str);
                }
            }
            inputLoggers.add(fin);
        }

        doLoop();
	}
	
	@Override
	protected void doLoop() throws IOException {
        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
        	log(sampleNr);
        	try {
        		Thread.sleep(500);
        	} catch (Throwable e) {
        		// ignore
        	}
        }
	}
	
	@Override
	public void log(long sampleNr) {
        if ((sampleNr < 0) || (sampleNr % every > 0)) {
            return;
        }
        for (int i = 0; i < inputLoggers.size(); i++) {
        	try {
        		String str = inputLoggers.get(i).readLine();
        		outputLoggers.get(i).println(str);
        	} catch (IOException e) {
        		// should not get here
        		e.printStackTrace();
        	}
        }
	}
	
}
