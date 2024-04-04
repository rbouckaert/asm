package asm.inference;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;

import java.io.PrintStream;
import java.util.Map;

public class GRTLogger extends BEASTObject implements Loggable {

    final public Input<TreePSRF> psrfInput = new Input<>("grt", "Tree GR to be logged", Input.Validate.REQUIRED);

    TreePSRF grt;  // todo rename the grt to something appropriate

    @Override
    public void initAndValidate() {grt = psrfInput.get();}

    @Override
    public void init(PrintStream out) {
        // Todo currently only works for two chains -- like the rest of the criterion atm
        out.print("chain0-GRT" + "\t");
        out.print("chain1-GRT" + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        // Todo add the pseudo ESS to this log too
        Map logValues = grt.getLog();

        int numChains = grt.getNumChains();

        for (int i = 0; i < numChains; i++) {
            if (Double.isNaN((Double) logValues.get("GRT-" + i))) {
                // todo this only happens in the check for sample 0 -- not sure why there is NaN
                out.print("-2.0\t");
            } else {
                out.print(logValues.get("GRT-" + i) + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
