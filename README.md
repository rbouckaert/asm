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
You can set parameters by clicking the GRLike button for the GR-like criterion or TraceESS button for the trace-ESS criterion.


<img width="916" alt="asm" src="https://github.com/rbouckaert/asm/assets/3612959/cfe9e9f4-0fb2-4f13-9d68-781dcde8b3ac">



# Paper 

Lars Berling, Remco Bouckaert, Alex Gavryushkin.
Automated convergence diagnostic for phylogenetic MCMC analyses
BioRxiv 2023
