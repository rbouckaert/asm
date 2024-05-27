package asm.tools;

import java.io.File;

import asm.inference.AutoStopMCMC;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Runnable;
import beast.base.parser.XMLParser;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.tools.Application;
import beastfx.app.util.XMLFile;

@Description("Emulate ASM analysis -- useful for comparing different stopping criteria on pre-run analyses")
public class ASMEmulator extends Runnable {
	final public Input<XMLFile> xmlInput = new Input<>("xml", "BEAST XML file containing an ASM analysis", Validate.REQUIRED);
	final public Input<File> inputLogDirInput = new Input<>("inputLogDir", "directory containing ASM trace and tree logs already run", Validate.REQUIRED);
	final public Input<File> outputLogDirInput = new Input<>("outputLogDir", "directory where emulator creates trace and tree logs", Validate.REQUIRED);

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		String sXML = BeautiDoc.load(xmlInput.get());
		sXML = sXML.replaceAll("AutoStopMCMC", "EmulatedAutoStopMCMC");
		
		XMLParser parser = new XMLParser();
		Object o = parser.parseFragment(sXML, true);
		if (!(o instanceof AutoStopMCMC)) {
			throw new IllegalArgumentException("XML file should contain an AutoStopMCMC analysis");
		}
		EmulatedAutoStopMCMC asm = (EmulatedAutoStopMCMC) o;
		asm.setLogDirs(inputLogDirInput.get().getPath(), outputLogDirInput.get().getPath());
		asm.run();
		
	}

	public static void main(String[] args) throws Exception {
		new Application(new ASMEmulator(), "ASMEmulator", args);
	}
}
