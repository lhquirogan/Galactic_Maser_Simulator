# GaMa: Galactic Methanol Maser distribution.
## This simulation generates samples of 6.7 GHz methanol masers in the entire Galaxy and then fit those maser to a Galactic model in order to obtain the Fundamental Galactic Parameters. A detailed description of this model can be found in Quiroga-Nuñez et al. 2017.
## This codes were develped by: Mark Reid (CfA), van Langevelde (JIVE/Leiden) and Quiroga-Nuñez (Leiden/JIVE)

## Introduction: Basic Structure.
The code is splitted in three folders:
	
	1. galaxy_generator: run_several_gala_v22.py

	This code generates the distribution of 6.7GHz galactic methanol masers
	based on the parameter file called: para_v22_...txt. For each maser,
	position, velocities and luminosity in different frames (equatorial,
	galactic and LSR) are simulated.
	
	2. galaxy_fitting: fit_gala_parallel_v22.py
	
	This codes take a whole galactic maser simulation or a subsample created
	by galaxy_generator part and fit it to find the best galactic parameters.

	3. galaxy_digest: digest_Reid_v22.py

	This codes digest the galactic parameters produced by the galaxy_fitting
	part. It produces the pearson coefficents and readable tables to compare
	with initial parameters intrduced in galaxy_generator part, especifically
	in the the parameter file called: para_v22_...txt.

# Setting up:

	1. Compile the fortran code in /galaxy_generator/BeSSeL_codes called
	fit_galaxy_bayesian_15.f as an executable called "a.out".
	Note: if there was a previous a.out please erase it before compiling.

	2. Set the paramerter values for your simulation in the file called para_vXX.txt in /galaxy_generator/

	3. Rename para_vxx.txt as para_vXX_XXX.txt, i.e., XXX sould be an 
	alfanumeric combination that help the user to remember the aim of that
	specific simulation. 
	Note: It is important to add some XXX, because the outputs will appear with this
	name

	4. Set the file control_file_bayesian in /galaxy_generator/BeSSeL_codes/ with
	the proper values to fit the simulations after the generation. Use the
	para_v22_XXX.txt settled before to "guide" the fitting.

	5. Check the file called control_file

# Running
	1. Creating Galaxies

		Run the run_several_gala_v22.py. This file takes GaMe_LHQN_v22.py and 
		run it several times, one for each galaxy desired. Therefore, the code will
		need the number of galaxies desired.
		Note:all of them will use the same para_v22_XXX.txt and control_file_bayesian.inp,
		so each galaxy contains the same distribution but different indivudal physical 
		paramaters per masers due to come from several distributions (gaussians,
		power laws and Monte Carlo)
		
		============================================================================
		Ouput: 
			A folder called output will be created if it the first time running this
			code. Inside that folder a new folder is created called si22_XXX which contains
			all the galaxies organized in folders. The names of the galaxies are selected
			using the date and time at the moment of the simulation.

			Inside each galaxy-folder several subfolders are created (depending of the samples
			desired in para_v22_XXX.txt). An extra folder of plots are also created (if plots
			were selected in para_v22_XXX.txt). Inside the subfolders, the fortran codes are
			added and each maser is saved as an independent file following this structure:
	2. Fitting Galaxies

		Run fit_gala_parallel_v22.py inside "/galaxy_fitting/". It will ask you
		for the path where the galaxies that you want to fit are. Typically, they
		are saved in: "/galaxy_generator/output/si22_XXX/". This codes is 
		parallelized, so it will ask you for the number of cores that you want to
		use and it will suggest the max. number avaliable in the computer that you used.
		Note: The fitting of each galaxy (only one sample in) take around 2.5 hours.

		============================================================================
		Ouput:

			The output files (one per sample) has ".prt" extension and they saved in a folder
			called "prt_files_XXX". Each prt files gives a lot of information, well worth reading.
			The solved-for parameters have binned, marginalized, probability density functions
			listed, as well as peak and 68% and 95% confidence ranges.
			Detailed "data, model, residuals" are listed for each source
			The parameter correlations are also listed.
			There are also 2 output files (inside each galaxy folder but not in prt_files_XXX):
			-fort.99   Lists sources dropped by distance limits and near/far issues
			-fort.7    Gives the "thinned" McMC trial values.  These can be used to re-calculate
			probability distributions as well as to check for convergence
			(by plottting Value vs trial #)
	# 3. Digesting Results

		Run digest_reid_v22.py inside "/galaxy_digest/". It will ask you
		for the path where the prt files are: "/galaxy_fitting/prt_files/prt_files_XXX"
		as well as if you want to make plots of the results distribution and pearson
		coefficent values. The information of the prt files will be classify and it makes
		some statistics of the results for each galaxy and each sample.
	
			============================================================================
			Ouput:

			It will produce a folder for each simulation as output_2_XXX. There a html file
			which contains the initial and final results is going to be saved. Besides, plots
			will be also saved here 
