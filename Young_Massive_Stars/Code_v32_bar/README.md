#Readme 
#Code Develped by:
#Mark Reid (CfA), van Langevelde (JIVE/Leiden) and Quiroga-Nunez (Leiden/JIVE)
#-----------------------------------------------------------------------------
#Basic Structure: The code is splitted in three folders:
	-galaxy_generator: run_several_gala_v22.py

	This code generates the distribution of 6.7GHz galactic methanol masers
	based on the parameter file called: para_v22_...txt. For each maser,
	position, velocities and luminosity in different frames (equatorial,
	galactic and LSR) are simulated.
	
	-galaxy_fitting: fit_gala_parallel_v22.py
	
	This codes take a whole galactic maser simulation or a subsample created
	by galaxy_generator part and fit it to find the best galactic parameters.

	-galaxy_digest: digest_Reid_v22.py

	This codes digest the galactic parameters produced by the galaxy_fitting
	part. It produces the pearson coefficents and readable tables to compare
	with initial parameters intrduced in galaxy_generator part, especifically
	in the the parameter file called: para_v22_...txt.
#-----------------------------------------------------------------------------
#Setting up:

	1. Compile the fortran code in /galaxy_generator/BeSSeL_codes called
	fit_galaxy_bayesian_15.f as an executable called "a.out".
	Note: if there was a previous a.out please erase it before compiling.

	2. Set the file called para_v22.txt in /galaxy_generator/: 

		#Parameters for Methanol Maser Simulation   
		#Developed by H. van Langevelde and L.H. Quiroga-Nunez 
		#Simulating BeSSeL data 
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Fundamental Galactic Parameters    
		#Distance from the galactic center to the sun (kpc)  
		r0=8.34
		#Circular rotational speed of the sun (km/s)
		v0t=240
		#Sun peculiar motion (km/s)
		usun=11.1
		vsun=14.6
		wsun=7.2
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Spatial Distribution
		#Standard deviation for the masers position wrt to closest arm (kpc)
		fat=0.35
		#Standard deviation of vertical position (kpc)
		hz=0.070  
		#Disk Scale Lenght
		halfmr=2.44
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Velocity Distribution
		#mean radial velocity (km/s)
		v0r=0 
		# Standard deviation for radial velocity (km/s)
		sigvr=5
		#Standard deviation for tangential velocity (km/s)
		sigva=15
		#Mean vertical velocity   
		v0z=0
		#Standard deviation for vertical velocity (km/s)     
		sigvz=5
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Mean source motion in U and V
		#-------------------------------------------------------------------------------
		#Mean source motion in U (km/s)
		v0r_us=9.0 
		# Standard deviation for mean source motion in U (km/s)
		sigvr_us=5.0
		#Counter Mean source motion in V (km/s) (warning:it contains a minus an inner minus)
		v0t_vs=-25.0
		# Standard deviation for mean source motion in V (km/s)
		sigva_vs=15.0
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Luminosity Distribution as a power law N(L)=k*L^index...min_l<L<max_l in (L_sun)
		index=-1.561
		min_l=0.00000001
		max_l=0.001
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Truncations
		# maximum distance or radius of the spiral arms final point (kpc)
		rtrunc=15 
		# Minimum distance or radius of the spiral arms final point (kpc)
		rmintrunc=3
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Fragmentation of spiral arms
		# number of semisegments that you to split the spiral arms in order to reproduce the segment
		nspirseg=40 
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Comparison with surveys
		# Sigma value of the MMB survey (1sigma=0.17 Jy) to compare results of the Luminosity Function
		MMB_sigma=4.0
		# Sigma value of the Arecibo survey (1sigma=0.27 Jy) to compare results of the Luminosity Function
		Arecibo_sigma=1.0 
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Extras
		#Number of maser/points per simulation (entire galaxy)
		nobj=1541
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Writing files for Reid's code
		#Type of sample selection
		#Sample selection (print out all sources=0) (1=yes, 0=no)
		sam_selec=1
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#if yes in sample selection:
		#Sample Selection limited by brightest sources and declination (1=yes, 0=no) (still will print the whole sample for comparison)
		sam_bri_dec=1
		#How many samples?
		many=2
		#Number of brightest sources, if they are more than one sample then write and array
		N_sour=100,200
		#declination limit down (up>delta>down), if they are more than one sample then write and array
		down_limit=-20,-20
		#declination limit up (up>delta>down), if they are more than one sample then write and array
		up_limit=91,91
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#if yes in sample selection:
		#Random Sample Selection limited by a defined number of sources (1=yes, 0=no)
		samp_rand=0
		#How many samples?
		many_rand=3
		#If yes, then Number of sources
		N_sour_rand=100,200,400
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Generate plots of the model (1=yes, 0=no)
		plots=1
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------
		#Assign measurement errors to the model (1=yes, 0=no----> all errors (parallax,v_lsr,mu_x and mu_y) are set to 0.0)
		Error_confirmation=0
		#-------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------

	3. Rename para_v22.txt as para_v22_XXX.txt, i.e., XXX sould be an 
	alfanumeric combination that help the user to remember the aim of that
	specific simulation. 
	Note: It is important to add XXX, because the outputs will appear with this
	name

	4. Set the file control_file_bayesian in /galaxy_generator/BeSSeL_codes/ with
	the proper values to fit the simulations after the generation. Use the
	para_v22_XXX.txt settled before to "guide" the fitting. Here is an explanation
	of the file:

		The program  fit_galaxy_bayesian_15.f  reads parallax and proper motion
		information for a number of sources and fits a model of the Galaxy,
		Solar Motion, and average source peculiar motions.   It uses a Bayesian
		Markov chain Monte Carlo technique to explore parameter probability
		distributions.

		-------------------------
		control_file_bayesian.inp

		! Control file for "fit_galaxy_bayesian" program
		! McMC control parameters first
		! Model parameters last
		! For model parameters, enter prior value, prior uncertainty (<0 for flat), est. posteriori uncertainty
		!
		      5                          ! number burn-in stages
		 100000                          ! number McMC iterations in a burn-in stage
		1000000      100                 ! number of final McMC iterations; store every nth for PDFs (max 10^6)
		      1    F T T T               ! error-tolerant fitting flag (0=least-squares, 1=error-tolerant); data types to use (par,mu_x,mu_y,Vlsr)
		     1.    -937895               ! initial step_fraction for McMC parameter adjustments; RAN1 seed (large negative integer)
		     5.0      0.0      4.0    0  ! "Virial" noise (km/s), dist_min (kpc), dist_GC_min (kpc), Rot.Curve
		  8.4        -1.       0.16      ! Ro (kpc)              
		240.0       -10.       6.0       ! To (km/s)             
		  0.0        -1.       0.5       ! dT/dr (km/s/kpc)      
		 11.1         2.1      1.0       ! Solar Motion Uo (km/s)   eg, SBD10 Uo = 11.1 +/- 2.0
		 12.2         2.       1.8       ! Solar Motion Vo (km/s)             Vo = 12.2 +/- 2.1
		  7.2         2.       0.5       ! Solar Motion Wo (km/s)             Wo =  7.2 +/- 2.0
		  5.        -10.       7.0       ! Ao (km/s) Avg pec motion counter to Gal Rot (km/s)
		  3.0       -10.       2.0       ! Bo   "     "   "    "    toward Gal Cen     (km/s)
		  0.          0.       0.        ! a1: Alternative rotation curve parameter 1
		  0.          0.       0.        ! a2: Alternative rotation curve parameter 2
		  0.          0.       0.        ! a3: Alternative rotation curve parameter 3

		Line 1: Number of "burn-in" stages (to go from initial parameter values to near optimal values and determine McMC step size)
		Line 2: Number of McMC iterations for each burn-in stage
		Line 3: Number of final McMC iterations;  Thinning parameter (store every N'th iteration)
		Line 4: Likelihood function code: 0 for Gaussian data errors (Least-squares), 1 for Sivia&Skilling "conservative formulation"
		             I use "1" to indicate outliers, 
		             then I spot those with stars "***" in the printout (eg, 1-star for 3-sigma deviation, 2-stars for 6-sig, etc) 
		             then remove the outliers from "source_files.inp", 
		             and refit with "0" to get an optimumL-S solution (which assumes data have Gaussian error distributions) 
		        Data types to fit, F (false) or T (true) for each data type (parallax, pm_x, pm_y, Vlsr)
		             I strongly recommend using "F T T T"
		Line 5: Expected 1-dimensional Virial "noise" in km/s units
		        Miminum distance from Sun to use source when fitting
		        Miminum distance from G.C. to use source when fitting
		        Rotation curve flag: see subroutine get_Tr for details:
		                             0 = linear (specified by Ro, To, dT/dR)
		                             1 = Clemens "10 kpc" model (rescaled to fitted Ro,To values)
		                             2 = Clemens "8.5 kpc" model (rescaled "    "     "       "
		                             3 = Brand-Blitz (specified by Ro, a1, a2, a3)
		                             4 = "Alternative"  "       "      "   "   "
		                             5 = "Universal"    "       "      "   "   "  (see Persic & Salucci 1997)
		Lines 6-16: Parameter lines, each with 3 values:
		        Initial value for parameter
		        Prior uncertainty (0=hold constant; negative=solve for with no prior info; positive=solve for with prior uncertainty constraint)
		        Posteriori uncertainty: estimate of final uncertainty (only used to control binning of parameter histograms)

		        For example, the Ro line:  8.4  -1. 0.16 
		                     means start at Ro=8.4kpc, solve for it with no constraints, and expect final uncertainty will be 0.16 kpc

				     the dT/dr line: 0.  0.  0.
		                     means hold dR/dR = 0 and don't solve for it

		                     the Vo line: 12.2   2.  1.8
		                     means start at Vo=12.2km/s, solve for it with a prior constriant of +/-2km/s, expecing final unc to be 1.8km/s

	5. Check the file called control_file
#-----------------------------------------------------------------------------
#Running
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
	3. Digesting Results

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