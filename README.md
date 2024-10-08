# ASM: Auto Stopping MCMC/MC3

ASM package for [BEAST 2](http://beast2.org)

ASM is a package that allows MCMC to stop automatically instead of having run till an arbitrary number a samples have been obtained.
This makes it more convenient to set up and run an analysis.

# Installing the ASM package

To install ASM, start BEAUti
* choose menu `File => Manage packages`
* Select ASM from the list of packages by clicking on it in the list. 
	If ASM is not in the list, add the package-extra-2.7 repository,
	which you can do in BEAUti through button `Package Repositories` 
	then click `Add URL` 
	then add 
	`https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml`
	 in the entry and click `OK`. 
	 Then click `Close` and the package should appear in the list (together with a few other experimental packages).
* click the `Install` button.

Restart BEAUti before using the package.

# Using the ASM package

Set up the analysis as per usual for an MCMC analysis, then go the the MCMC tab.
From the drop-down box select `Automatic Stopping MCMC`, and a number MCMC options disappear and are replaced by ASM options.

<img width="916" alt="asm" src="https://raw.githubusercontent.com/rbouckaert/asm/main/doc/asm-select.png">

You can set parameters by clicking the TreePSRF button for the Gelman-Rubin-like criterion or TraceESS button for the trace-ESS criterion.

<img width="916" alt="asm" src="https://raw.githubusercontent.com/rbouckaert/asm/main/doc/asm-options.png">


# Interpreting the output

Screen output may look something like this:

```
Check 19   burnin:[14, 14]psrf1mean = 1.472 reset start -1 10.0 5.1 10.0 :5.1  in 5 mseconds
Check 20   burnin:[15, 15]psrf1mean = 1.399 reset start -1 7.2 5.3 10.0 :5.3   in 5 mseconds
Check 21   burnin:[15, 15]psrf1mean = 1.426 reset start -1 10.9 7.2 12.0 :7.2  in 5 mseconds
Check 22   burnin:[16, 16]psrf1mean = 1.466 reset start -1 12.0 12.0 12.0 :12.0      in 3 mseconds
Check 23   burnin:[17, 17]psrf1mean = 1.525 reset start -1 8.8 12.0 8.9 :8.8   in 7 mseconds
Check 24   burnin:[18, 18]psrf1mean = 1.504 reset start -1 7.6 12.0 6.9 :6.9   in 18 mseconds
```

Let’s have a look at the first entry:

```
Check 19   burnin:[14, 14]psrf1mean = 1.472 reset start -1 10.0 5.1 10.0 :5.1  in 5 mseconds
```

**Check 19**   check number

**burnin:[14, 14]** log items removed as burnin, first number for chain1, second for chain2

**psrf1mean = 1.472**  potential scale reduction factor — this statistic represents a Gelman-Rubin like statistic for the trees. It should go towards 1 at convergence.

**reset start -1** whether the start of the psrf is recalculated from the start

**10.0 5.1 10.0** ESSs of the statistics specified in the `asm.inference.TraceESS` statistic

**:5.1** minimum of the ESSs of the `asm.inference.TraceESS statistic`. This should go to the threshold specified in the `targetESS` attribute.

**in 5 mseconds** time used to calculate the statistics




# Paper 

If you publish a paper using this package, please cite<br>
Lars Berling, Remco Bouckaert, and Alex Gavryushkin<br>
**Automated convergence diagnostic for phylogenetic MCMC analyses**<br>
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2024
[doi: 10.1109/TCBB.2024.3457875]([https://doi.org/10.1109/TCBB.2024.3457875)<br>
Preprint on *BioRxiv 2023*<br>
[https://doi.org/10.1101/2023.08.10.552869](https://doi.org/10.1101/2023.08.10.552869)<br>
