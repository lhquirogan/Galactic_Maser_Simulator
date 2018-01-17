      program fit_galaxy_bayesian

c     Bayesian fitting with "error-tolerant" data errors option

c     Fits a model of Galactic rotation and Solar Motion to 
c     parallax and proper motion to data

c     Started with fit_galaxy_bayesian_11.f ...

c     V11 Generalize to other rotation curves:
c         Make k_dist more robust to sources near center/anti-center
c         Add 3 rotation curve parameters: p(9), p(10), p(11) 
c         Add iRC flag: 
c           0 for Standard   Vr = Vo + (dV/dR)*rho    
c           1 Clemens-Ro=10  Vr empiracle formula; then scaled to Ro,To 
c           2 Clemens-Ro=8.5 Vr empiracle formula; then scaled to Ro,To
c           3 Brand-Blitz    Vr = a1(rho)^a2 + a3
c           4 Alternative    Vr = a1 + a2*(rho-1) + a3*(rho-1)**2
c           5 Universal      Vr = disk + halo model
c         Allow more iterations by thinning (new control parameter) 
c         stored values (keeping dim's as is)
c     V12 removed earlier change to re-calculate near/far code
c         in calc_model (based on trial Ro
c     V13 changed subroutine k_dist to k_dist_smooth, which removes
c         discontinuous D_K(Vlsr), but needs near/far as input
c     V14 allows any or all of the 4 data types (pi,mu_x,mu_y,Vlsr) to be fitted
c         Initialiazed random number seed explicitly
c     V15 added minimum percentage parallax uncertainty parameter
c     V16 de-bias distances that are based on parallaxes:
c         The mean distance from asymmetric distance PDF will be slightly
c         larger than it should.   Correct for this.

c----------------------------
c     Input files:
c        source_files.inp            ! contains names of input files; 
c                                    !   each file has data from one source
c        control_file_bayesian.inp   ! various controls and parameters
c----------------------------

      implicit real*8 (a-h,o-z)

      character*48  data_file, out_file, data_files(2000)
      character*12  src, source(2000)

      real*8        data(8000), model(8000), resids(8000), res_err(8000)

      real*8        params(11), prior_val(11), post_sigma(11) 
      real*8        prior_unc(11), post_unc(11)
      character*16  parnames(11)

c     Following for McMC trials; save memory with R*4's here
      real*4     params_stored(11,1000000)
      real*4     params_burnin(11,100000)

      real*8     params_best(11), delta_param(11)

      logical    err_tol, accepted, use_source, use_data_type(4)

      character*16 data_names(4)
      data       data_names/'Parallax (mas)  ','X-motion (mas/y)',
     +                      'Y-motion (mas/y)','V_Helio  (km/s) '/

      integer    data_type(8000), source_number(8000)
      real*8     ra(2000),dec(2000),parallaxes(2000),par_errors(2000),
     +           vhelios(2000),err_vhelios(2000),near_fars(2000)
      common /model_info/ source_number, data_type, ra, dec,
     +                    near_fars,parallaxes,par_errors,
     +                    vhelios,err_vhelios


c     Set some maximum control values for above array dimensions
      max_stored     = 1000000    ! maximum final McMC trials stored        
      max_burnin     =  100000    ! number trials for each burnin stage
      max_num_files  =    2000    ! one file per source used
      max_num_data   =    8000    ! allows for 4 data points per source 
      max_parameters =      11    ! maximum number parameters 

      lu_print     = 6
      lu_out       = 7
      lu_data      = 8
      lu_control   = 9

      write (lu_print,1000)
 1000 format(/' Program fit_galaxy_bayesian_16: 2016 Nov 02',//)

      pi = 4.d0*atan(1.d0)
      deg_to_rad = pi/180.d0

c     Random number seeds (idum now overriden with "call control_parameters")
      idum  = -725391
      xseed = 0.75381d0

c     =============================================================
c     Get input data file names...
      call get_source_file_names ( lu_control, lu_print, max_num_files,
     +                             data_files, n_files )

c     =============================================================
c     Get input control parameters...
      call control_parameters ( lu_control, lu_print,
     +      num_stages, iterations_per_stage, itermax, n_thin,
     +      err_tol, step_fraction, use_data_type, idum,
     +      num_params, params, prior_unc, post_unc, parnames,
     +      par_virial, dist_min, dist_GC_min, iRC, par_max_unc )

      if ( iterations_per_stage .gt. max_burnin ) then
         write (lu_print,1003) iterations_per_stage, max_burnin
 1003    format(' STOP: number iterations per stage',i9,' greater than',
     +          ' dimension limit',i9)
         STOP
      endif

      if ( num_params .gt. max_parameters ) then
         write (lu_print,1004) num_params, max_parameters
 1004    format(' STOP: number parameters',i3,' greater than',
     +          ' dimension limit',i3)
         STOP
      endif

      if ( itermax/n_thin .gt. max_stored ) then
         write (lu_print,1005) itermax,n_thin, max_stored
 1005    format(' STOP: requested',i9,' iterations divided by',
     +     ' "thinning" factor"',i4,' is greater',
     +     ' than dimension limit of 10^6 for stored trials:',i9)
         STOP
      endif

      if ( err_tol ) write (lu_print,1010)
 1010 format(/' Using "error-tolerant" data uncertainty formulation')

      write (lu_print,1020) step_fraction
 1020 format(/' McMC parameter steps: fraction of expected', 
     +        ' posteriori sigmas:',f9.6)

c     Store priors...
      do n_p = 1, num_params
         prior_val(n_p) = params(n_p)
      enddo

c     =============================================================
c     Read all data files...

      num_used    = 0
      i_d         = 0

      Ro = params(1) ! initial value used for checking  near/far code
      write (lu_print,1030) Ro
 1030 format(/' Assuming Ro =',f5.2,' kpc for checking near/far code')

      do i_s = 1, n_files

c        Input data for a source ...  
         data_file = data_files(i_s)
         call read_one_source ( lu_data, data_file, lu_print, 
     +        src,ra_hhmmss,dec_ddmmss, farnear, par, err_par, 
     +        pm_x, err_pm_x, pm_y, err_pm_y, vhelio, err_v )

c        Convert ra,dec to decimal; calculate Galactic coords
         call hmsrad ( ra_hhmmss,  ra_rad )
         call dmsrad ( dec_ddmmss,dec_rad )
         call radec_to_galactic ( ra_rad, dec_rad, 
     +                            ell_II, bee_II )

c        Decide if want to use this source...
         use_source = .true.             ! default is to use it

         dist_par   = 1.d0/par           ! kpc
         if ( dist_par .lt. dist_min ) then
            use_source = .false.
            write (lu_print,1050) src, dist_par, dist_min
c           write (99,1050)       src, dist_par, dist_min
 1050       format(' Discarding ',a12,': Solar D =',f5.1,
     +             ' < ',f5.1,' kpc')
         endif

c        Discard sources with unacceptably large fractional parallax uncertainty
         par_unc_percentage = 100.d0 * abs(err_par/par)    ! percentage
         if ( par_unc_percentage.gt.par_max_unc .or. 
     +        par_unc_percentage.gt.50.d0            ) then
            use_source = .false.
            write (lu_print,1052) src,par_unc_percentage,par_max_unc
c           write (99,1050)       src,par_unc_percentage,par_max_unc
 1052       format(' Discarding ',a12,': parallax unc (%)',f5.1,
     +             ' > ',f5.1,'%')
         endif

         cos_l   = cos( ell_II * deg_to_rad )
         r_sq = Ro**2 + dist_par**2 - 2.d0*Ro*dist_par*cos_l
         dist_GC = sqrt( r_sq )          ! kpc

         if ( dist_GC  .lt. dist_GC_min ) then
            use_source = .false.
            write (lu_print,1055) src, dist_GC, dist_GC_min
c          write (99,1055)       src, dist_GC, dist_GC_min
 1055       format(' Discarding ',a12,': Galactocentric D =',f5.1,
     +             ' < ',f5.1,' kpc')
         endif

c        Check, but don't change, near/far based on parallax distance and input Ro value
         dist_tp = Ro * cos_l            ! kpc
         farnear_calc  = 0.d0            ! start at "near"
         if ( dist_par.gt.dist_tp .and. 
     +        dist_tp.gt.0.d0           ) farnear_calc = 1.d0 
         if ( abs(farnear - farnear_calc) .gt. 0.1 ) then
c            write (99,9732) src, farnear, farnear_calc
c 9732       format(' NOTE: near/far code for ',a12,
c     +             '; input=',f4.1,'; calculated=',f4.1)
         endif

         if ( use_source ) then 

            num_used = num_used + 1

c           Store information on this source...
            source(num_used)    = src
            ra(num_used)        = ra_rad              ! radians
            dec(num_used)       = dec_rad             ! radians
            near_fars(num_used) = farnear
            parallaxes(num_used)= par                 ! mas
            par_errors(num_used)= err_par             ! mas

c           Calculate distance bias correction for asymmetric PDF
            call distance_bias ( par, err_par,
     +                           bias_factor )
            write (lu_print,1062) par, err_par, bias_factor
 1062       format(/' When converting parallax =',f7.3,' +/-',f6.3,
     +              ' to distance,'
     +             /4x,'will divide distance by bias',
     +              '-factor of',f7.4)

c           Reset "err_v" to include contribution from Virial noise
            err_v = sqrt( err_v**2 + par_virial**2 )  ! km/s
            vhelios(num_used)  = vhelio      ! store vhelio by source
            err_vhelios(num_used) = err_v    ! store error on vel by src

c           Transfer 4 types of data to 1-D arrays...
            i_d = i_d + 1
            data(i_d)          = par

c           Add uncertainty in modeling the parallax (from a kinematic
c           distance due to the error in measured velocity) to the  
c           measured parallax error.  
            call get_kdist ( params, iRC, ra_rad, dec_rad, farnear,
     +                       par, err_par, vhelio, err_v,
     +                       dist, unc_dist )
            err_kd_par_sq = ( unc_dist / dist**2 )**2
            res_err(i_d)  = sqrt( err_par**2 + err_kd_par_sq )
            write (lu_print,1060) err_par,res_err(i_d),dist,unc_dist
 1060       format(' Parallax unc: input',f6.3,'; inflated to',f6.3,
     +             ' mas,'
     +            /'    based on Dk =',f7.2,' +/-',f7.2,' kpc')

            data_type(i_d)     = 1           ! indicates Parallax datum
            source_number(i_d) = num_used    ! source number

c           Proper motions...
c           set error floors for proper motion data based on 
c           uncertainty in source radial velocity...
c           conversion from km/s/kpc to mas/yr...(0.211 = 1/4.74)
            convert   = 365.2422d0 * 86400.d0 / 1.49597871d8
            err_pm_floor = convert * ( err_v * par ) ! mas/yr
            write (lu_print,1065) err_pm_floor
 1065       format(' Using proper motion err-floor of',f7.3,' mas/yr')

            i_d = i_d + 1
            data(i_d)          = pm_x
            res_err(i_d)       = sqrt( err_pm_x**2 + err_pm_floor**2 )
            data_type(i_d)     = 2           ! indicates X-motion datum
            source_number(i_d) = num_used    ! source number

            i_d = i_d + 1
            data(i_d)          = pm_y
            res_err(i_d  )     = sqrt( err_pm_y**2 + err_pm_floor**2 )
            data_type(i_d)     = 3           ! indicates Y-motion datum
            source_number(i_d) = num_used    ! source number

            i_d = i_d + 1
            data(i_d)          = vhelio
            res_err(i_d)       = err_v       
            data_type(i_d)     = 4           ! indicates Vhelio datum
            source_number(i_d) = num_used    ! source number

         endif   ! edit sources by distances from Sun and GC
                                             
      enddo

      num_sources = num_used 
      num_data    = i_d

      num_removed = n_files - num_used
      write (lu_print,9127) num_used, num_removed
 9127 format(//'=============================================',
     +         '=============================================',
     +        /' Number of sources used=',i3,'; number removed=',i3)
c     write (99,9128) num_sources, num_data
 9128 format(' Total # sources used=',i3,'; # data points=',i3)

c     ============================================================
c     ============================================================
c     Markov chain Monte Carlo method of Bayesian fitting

c     Define characteristic steps in parameters space for McMC exploration
c     of the posteriori parameter probability density function; 
c     this is an "art" not a "science" and needs experimentation for 
c     each type of problem and number of parameters.  Aim for about 25%
c     acceptance rate for Metropolis-Hastings algorithm.

c     ============================================================
c     "Burn-in" trials...
      do n_p = 1, num_params
         post_sigma(n_p)  = post_unc(n_p)
         delta_param(n_p) = step_fraction * post_sigma(n_p)
      enddo

      write (lu_print,1099)
 1099 format(/' Starting burn-in trials...',/)

      do n_pass = 1, num_stages

         n_accepted = 0
         n_best     = 0
         sum_ln_P_old  = 0.d0
         sum_ln_P_best =-9.d99

         do iter = 1, iterations_per_stage

            call mcmc_exploration (num_params, params, iRC,
     +                 delta_param, prior_val, prior_unc, err_tol,
     +                 num_data, data, res_err, use_data_type,
     +                 idum,xseed, sum_ln_P_old, sum_ln_P_new, 
     +                 accepted)
 
c           Store all parameters from this trial...
            do n_p = 1, num_params
               params_burnin(n_p,iter) = params(n_p)
            enddo

            if ( accepted ) n_accepted = n_accepted + 1
         
c           Keep track of best trial parameters and iteration number...
            if ( sum_ln_P_new .gt. sum_ln_P_best ) then
                 sum_ln_P_best = sum_ln_P_new
                 do n_p = 1, num_params
                    params_best(n_p) = params(n_p)
                 enddo
                 n_best = iter
            endif

         enddo

c        -------------------------
c        Estimate parameter posterior 1-sigma values
c        and adjust parameter step-sizes accordingly...

c        Posteriori probability distribution for globals...
         call marginalize_burnin (iterations_per_stage, num_params, 
     +                   parnames, params_burnin, params_best,
     +                   prior_unc, post_sigma)

c        Reset parameter step-sizes to achieve 40% acceptance rate
c        for this problem which has an "intermediate" number of parameters
         acc_ratio = float(n_accepted)/float(iterations_per_stage)
         acceptance_percent = 100.d0*acc_ratio
         step_new = step_fraction*(acceptance_percent/40.d0)
c        But adjust slowly...
         step_fraction = ( step_new + step_fraction ) / 2.d0
         write (lu_print,1100) n_pass, num_stages,
     +                         acceptance_percent, step_fraction
 1100    format(' Burn-in stage',i3,' of',i3,': acceptance rate =',f6.1,
     +          '%; new step-fraction =',f8.5)
         do n_p = 1, num_params
            delta_param(n_p) = step_fraction * post_sigma(n_p)
         enddo

         write (lu_print,1115) (params_best(n_p), n_p=1,num_params)
 1115    format('         best param values:',11f10.3)
         write (lu_print,1120) (post_sigma(n_p), n_p=1,num_params)
 1120    format('         posteriori_sigmas:',11f10.3)

      enddo     ! burnin stages

c     ==============================================================
c     Now final McMC  trials...
      sum_ln_P_old  = 0.d0
      sum_ln_P_best =-9.d99
      n_best        = 0
      n_accepted    = 0
      n_stored      = 0

      do iter = 1, itermax 

         call mcmc_exploration (num_params, params, iRC,
     +                 delta_param, prior_val, prior_unc, err_tol,
     +                 num_data, data, res_err, use_data_type, 
     +                 idum,xseed, sum_ln_P_old, sum_ln_P_new, 
     +                 accepted)

         if ( mod(iter,n_thin) .eq. 0 ) then
c            Store the global parameters from this trial...
             n_stored = n_stored + 1
             do n_p = 1, num_params
                params_stored(n_p,n_stored) = params(n_p)
             enddo
         endif

         if ( accepted ) n_accepted = n_accepted + 1

c        Write to file for plotting
         if ( mod(iter,n_thin) .eq. 0 ) then
           write (lu_out,2222)iter,(params(n_p),n_p=1,num_params)
 2222      format(i8,2f10.4,f8.3,3f8.3,2f8.3,3f8.2)
         endif

c        Keep track of best trial parameters and iteration number...
         if ( sum_ln_P_new .gt. sum_ln_P_best ) then
              sum_ln_P_best = sum_ln_P_new
              do n_p = 1, num_params
                 params_best(n_p) = params(n_p)
              enddo
              n_best = iter
         endif

      enddo

      percent_accepted = 100.d0*n_accepted/float(itermax)
      write (lu_print,2230) percent_accepted 
 2230 format(//' McMC percent of trials accepted (aim for ~40%)',
     +         f6.1,'%')

c     ----------------------------------------------------------------
c     Best trial global parameter values...
      write (lu_print,2233) n_best
 2233 format(/' Best parameter values for trial',i9)
      call best_parameters ( params_best, iRC, prior_unc,
     +                  data, model, resids, res_err,
     +                  num_params, num_data, source,
     +                  use_data_type,
     +                  lu_print, lu_plot )

c     Marginalized posteriori probability distribution for parameters...
      call marginalize ( lu_print, n_stored, num_params, 
     +                   parnames, params_stored, params_best,
     +                   prior_unc, post_unc, post_sigma )

c     PDF for (To+Vsun)/Ro
      call est_angular_rotation(lu_print, n_stored, num_parameters, 
     +                   parnames, params_stored, params_best, 
     +                   prior_unc, post_unc, post_sigma_Om)
      write (lu_print,2300) post_sigma_Om
 2300 format(/18x,' Combined 1-sigma uncertainty for (To+Vsun)/Ro =',
     +        f7.2,' km/s/kpc')

c     ----------------------------------------------------------------
c     Residuals printout
      write (lu_print,3000) 
 3000 format (//17x,'         Data       Model    Residual      ',
     +               ' Error     Res/Err'/)
      n = 0
      do i_s = 1, num_sources
         write (lu_print,3050) source(i_s)
 3050    format(/' Data for ',a12)

         do nn = 1, 4                      ! 4 data points per source
            n = n + 1
            wtres = resids(n) / res_err(n)
            write (lu_print,3100) data_names(nn),
     +                            data(n), model(n), resids(n), 
     +                            res_err(n), wtres
 3100       format (1x,a16,1x,5f12.4)
         enddo

      enddo

c     ----------------------------------------------------------------
c     Correlation matrix...
      call correlation_matrix(params_stored,num_params,n_stored,
     +                        lu_print )


      stop
      end

c     ================================================================
c     ================================================================

c     ================================================================
      subroutine mcmc_exploration (num_params, params, iRC,
     +                 delta_param, prior_val, prior_unc, err_tol,
     +                 num_data, data, res_err, use_data_type, 
     +                 idum,xseed, sum_ln_P_old, sum_ln_P_new, 
     +                 accepted)

c     Explores parameter space via a Markov chain Monte Carlo method
c     Return new parameter values and sum_ln_P 

      implicit real*8 ( a-h,o-z )

      real*8  params(11), params_wiggled(11), delta_param(11)
      real*8  data(8000), model(8000), resids(8000), res_err(8000)

      real*8  prior_val(11), prior_unc(11)

      logical err_tol, accepted, use_data_type(4)

c     Set new trial parameters...
      call trial_param_values( idum,xseed, num_params, 
     +                         params, iRC, delta_param,
     +                         params_wiggled ) 
      
c     Compute Metropolis ratio:
c     Calculate data probability...
      call calc_residuals( params_wiggled, iRC, num_data,
     +                     data,
     +                     model, resids )

      call calc_ln_P_data (err_tol, num_data, 
     +                     resids, res_err, use_data_type,
     +                     sum_ln_P_data)

c     Add priors probability...
      call calc_ln_P_priors ( num_params, params_wiggled, 
     +                        prior_val, prior_unc,
     +                        sum_ln_P_priors )

      sum_ln_P_new = sum_ln_P_data + sum_ln_P_priors

      if ( sum_ln_P_old .ne. 0.d0 ) then
            ratio = exp( sum_ln_P_new - sum_ln_P_old )
         else
            ratio = 1.d0
      endif

c     Metropolis accept/decline decision...
c     call gauss_random (xseed, u1, u2, g1, g2 )
      call gauss_ran1   (idum,  u1, u2, g1, g2 )
      if ( u2 .le. ratio ) then
c           Accept trial
            accepted = .true.
            do n_p = 1, num_params
               params(n_p) = params_wiggled(n_p)
            enddo
c            print *,'DEBUG: accepted ln(Ps); ratio; U:',
c     +              sum_ln_P_old,sum_ln_P_new,ratio,u2
            sum_ln_P_old = sum_ln_P_new
         else
c           Decline trial
            accepted = .false.
c            print *,'DEBUG: declined ln(Ps)=',
c     +              sum_ln_P_old,sum_ln_P_new,ratio,u2
            sum_ln_P_new = sum_ln_P_old
      endif
        
      return
      end

c     =============================================================
      subroutine trial_param_values(idum, xseed, num_params,
     +                              params, iRC, delta_param,
     +                              params_wiggled ) 

      implicit real*8 (a-h,o-z)

      real*8  params(11), params_wiggled(11), delta_param(11)

c     Set a new set of trial parameters...
      do i_p = 1, num_params

         if ( delta_param(i_p) .gt. 0.d0 ) then

c              call gauss_random (xseed, u1, u2, g1, g2 )
               call gauss_ran1   (idum,  u1, u2, g1, g2 )
c              Gaussian distribution of parameter changes...
               delta = g1*delta_param(i_p)
               params_wiggled(i_p) = params(i_p) + delta

            else

c              Parameter not to be varied...
               params_wiggled(i_p) = params(i_p)

         endif 
  
      enddo

c     Check Rotation Curve models for which To can be
c     calculated from parameters Ro, a1, a2, and a3...
      if ( iRC .ge. 3 ) then
c        Setup to get To, so that we can get its PDF later
         Ro   = params_wiggled(1)
         To   = 0.d0                ! will be calcuated
         dTdr = 0.d0                ! not used
         a1   = params_wiggled(9)
         a2   = params_wiggled(10)
         a3   = params_wiggled(11)

         r    = Ro
         call get_Tr( Ro, To, dTdr, iRC, a1, a2, a3, r,
     +                Tr )

c        Insert calculated value of To in "params_wiggled" 
         params_wiggled(2) = Tr
      endif

c     Do 50% of trials anti-corellated for (Vsun,Ao) 
c     and directly-correlated for (Usun,Bo)
c     call gauss_random (xseed, u1, u2, g1, g2 )
      call gauss_ran1   (idum,  u1, u2, g1, g2 )
      if ( u2 .gt. 0.5d0 ) then

c        call gauss_random (xseed, u1, u2, g1, g2 )
         call gauss_ran1   (idum,  u1, u2, g1, g2 )

         if ( delta_param(5).gt.0.d0 .and.
     +        delta_param(7).gt.0.d0       ) then
c           Do (Vsun,-Ao) = params(5,-7); negatively correlated
            delta = g1 * delta_param(5)
            params_wiggled(5) = params(5) + delta
            delta =-g1 * delta_param(7)
            params_wiggled(7) = params(7) + delta
         endif

         if ( delta_param(4).gt.0.d0 .and.
     +        delta_param(8).gt.0.d0       ) then
c           Do (Usun,+Bo) = params(4,8); directly correlated
            delta = g2 * delta_param(4)
            params_wiggled(4) = params(4) + delta
            delta = g2 * delta_param(8)
            params_wiggled(8) = params(8) + delta
         endif

      endif

      return
      end

c     ============================================================
      subroutine calc_ln_P_data (err_tol, num_data, 
     +                           resids, res_err, use_data_type,
     +                           sum_ln_P)

      implicit real*8 (a-h,o-z)

      real*8     resids(8000), res_err(8000)

      logical    err_tol, use_data_type(4)

      sum_ln_P = 0.d0

      if ( err_tol ) then

c          Use error-tolerant data uncertainties
c          Calculate ln(probability) when pdf(sigma|sigma_o) = sigma_o/sigma^2

           do n_d = 1, num_data
              n_dt = mod((n_d-1),4) + 1
              if ( use_data_type(n_dt) ) then
                 R_2 = ( resids(n_d) / res_err(n_d) )**2
                 arg = -R_2/2.d0
                 sum_ln_P = sum_ln_P+log((1.d0 - exp(arg))/R_2) 
              endif
           enddo

        else

c          Normal Gaussian distribution of errors...
c          Calculate ln(probability) for Gaussian distribution of data 
           do n_d = 1, num_data
              n_dt = mod((n_d-1),4) + 1
              if ( use_data_type(n_dt) ) then
                 R_2 = ( resids(n_d) / res_err(n_d) )**2
                 sum_ln_P = sum_ln_P - R_2 / 2.d0
              endif
           enddo

      endif

      return
      end

c     ==============================================================
      subroutine calc_ln_P_priors ( num_params, params, 
     +                              prior_val, prior_unc,
     +                              sum_ln_P )

      implicit real*8 (a-h,o-z)

      real*8  params(11), prior_val(11), prior_unc(11)

c     Calculate ln(probability) for Gaussian distribution of 
c     parameters with priors

      sum_ln_P = 0.d0

      do n = 1, num_params

         if ( prior_unc(n) .gt. 0.d0 ) then

            resid = params(n) - prior_val(n) 
            R_2 = ( resid / prior_unc(n) )**2
            sum_ln_P = sum_ln_P - R_2 / 2.d0

         endif

      enddo

      return
      end

c     ==========================================================
      subroutine best_parameters ( params, iRC, prior_unc,
     +                  data, model, resids, res_err,
     +                  num_params, num_data, source,
     +                  use_data_type,
     +                  lu_print, lu_plot )

c     Best parameters' output...

      implicit real*8 (a-h,o-z)

      real*8     params(11), prior_unc(11)
      real*8     data(8000), model(8000), resids(8000), res_err(8000)
      real*8     wtres(4)

      logical    use_data_type(4)

      character*12 source(2000)
      character*5  stars(6)
      data  stars/'     ','*    ','**   ','***  ','**** ','*****'/


      write (lu_print,2244) (params(n_p), n_p=1,num_params) 
 2244 format(2x,11f10.3)

c     Calculate residuals for best parameters...
      call calc_residuals( params, iRC, num_data,
     +                     data,
     +                     model, resids )

c     Calculate chi-squared per degree of freedom for each data type...
      chisq_par = 0.d0
      chisq_mux = 0.d0
      chisq_muy = 0.d0
      chisq_vel = 0.d0
      n = 0
      num_masers = num_data/4
      do n_d = 1, num_masers
         n = n + 1
         chisq_par = chisq_par + (resids(n)/res_err(n))**2
         n = n + 1
         chisq_mux = chisq_mux + (resids(n)/res_err(n))**2
         n = n + 1
         chisq_muy = chisq_muy + (resids(n)/res_err(n))**2
         n = n + 1
         chisq_vel = chisq_vel + (resids(n)/res_err(n))**2
      enddo

c     Count number of parameters actually varied...
      n_varied_params = 0
      do n_p = 1, num_params
         if ( prior_unc(n_p) .ne. 0 ) then
            n_varied_params = n_varied_params + 1
         endif
      enddo

c     What to do for numbers of parameters in the d.o.f. calculation 
c     for one data type?  Will assume 1/4 of the total number of parameters 
c     for each of the 4 types of data
      deg_freedom_1D = num_masers - n_varied_params/4.d0
      if ( deg_freedom_1D .gt. 0 ) then
            reduced_chisq_par = chisq_par / deg_freedom_1D 
            reduced_chisq_mux = chisq_mux / deg_freedom_1D 
            reduced_chisq_muy = chisq_muy / deg_freedom_1D 
            reduced_chisq_vel = chisq_vel / deg_freedom_1D 
         else
            reduced_chisq_par = -1.d0
            reduced_chisq_mux = -1.d0
            reduced_chisq_muy = -1.d0
            reduced_chisq_vel = -1.d0
      endif

c     Now do true combined chi-squared for all types of data used...
      chisq  = 0.d0
      ntypes = 0
      if ( use_data_type(1) ) then 
         chisq = chisq + chisq_par
         ntypes = ntypes + 1
      endif
      if ( use_data_type(2) ) then 
         chisq = chisq + chisq_mux
         ntypes = ntypes + 1
      endif
      if ( use_data_type(3) ) then 
         chisq = chisq + chisq_muy
         ntypes = ntypes + 1
      endif
      if ( use_data_type(4) ) then 
         chisq = chisq + chisq_vel
         ntypes = ntypes + 1
      endif

      n_deg_freedom = ntypes*num_data/4.d0 - n_varied_params
      if ( n_deg_freedom .gt. 0 ) then
            reduced_chisq = chisq / float( n_deg_freedom )
         else
            reduced_chisq = -1.d0
      endif

      write (lu_print,3000) chisq, reduced_chisq, n_deg_freedom, 
     +                      reduced_chisq_par, reduced_chisq_mux,
     +                      reduced_chisq_muy, reduced_chisq_vel
 3000 format(/' Fit parameters & weighted residuals for best ',
     +   'trial: chisq =',f10.3,'; chisq_pdf =',f10.3,';',i5,' d.o.f.',
     +       /' Individual data-type "chisq_pdfs":',4f12.3,
     +      //' Source              ',
     +        ' chiPar    chiMuX    chiMuY    chiVel')

c     Print out weighted residuals
      i_d = 0
      do n_m = 1, num_masers

         max_sigma = 0
         do ii = 1, 4
            i_d = i_d + 1
            wtres(ii) = resids(i_d) / res_err(i_d)
            if (max_sigma.lt.abs(wtres(ii))) max_sigma=abs(wtres(ii))
         enddo

c        Flag with "*'s"; one * for multiple of 3-sigma
         n_stars = 1 + ( max_sigma/3.d0 )
         if ( n_stars .gt. 6 ) n_stars = 6

         write (lu_print,3010) source(n_m), wtres, stars(n_stars)
 3010    format(1x,a12,5x,4f10.3,5x,a5)

      enddo

      return
      end

c     =================================================================
      subroutine marginalize ( lu_print, itermax, num_parameters, 
     +           parnames, params_stored, params_best, 
     +           prior_unc, post_unc, post_sigma )

c     Determines posteriori probability distribution functions for parameters

      implicit real*8 (a-h,o-z)

      real*4        params_stored(11,1000000)
      real*8        prior_unc(11), post_unc(11), params_best(11), 
     +              post_sigma(11)
      character*16  parnames(11)

c     Bin the parallax values, which marginalizes (integrates) over the other params.  
      do n_p = 1, num_parameters

c        Check if "solved for" parameter: 
c        if prior uncertainty > 0 (Gaussian prior) or < 0 (for a flat prior)
c        But if prior uncertainty =0, held it constant
         post_sigma(n_p) = 0.d0
         if ( prior_unc(n_p) .ne. 0.d0 ) then

            n_bins    = 200
            par_best  = params_best(n_p)
            par_delta = 0.15d0 * post_unc(n_p)              ! 15% of posteriori uncertainty
            par_start = par_best - (n_bins/2)*par_delta - par_delta/2.d0

            if ( lu_print .gt. 0 ) then
               write (lu_print,1000) parnames(n_p) 
 1000          format(//' Binned probability for ',a16)
            endif

            p_max = 0
            n_max = 0

            integrate_up = 0
            integrate_dn = itermax + 1

c           "2-sigma" one-sided threshold (0.5*4.6%  = 2.3%) 
            n_2sig         = 0.023d0 * itermax + 0.5d0
            two_sigma_low  = 0
            two_sigma_high = 0

c           "1-sigma" one-sided threshold (0.5*31.8% = 15.9%) 
            n_1sig         = 0.159d0 * itermax + 0.5d0
            one_sigma_low  = 0
            one_sigma_high = 0

            
            do n_bin = 1, n_bins 

               p0 = par_start + par_delta*n_bin
               p1 = p0 + par_delta

               n_in_bin = 0
               do iter = 1, itermax

                  par_val = params_stored(n_p,iter)

                  if ( par_val.ge.p0 .and. par_val.lt.p1) then 
                     n_in_bin = n_in_bin + 1
                  endif

               enddo

               if ( n_in_bin .gt. 0 ) then 
                  p_mid = ( p0 + p1 ) / 2.d0
                  if ( lu_print .gt. 0 ) then
                     write (lu_print,2000) p_mid, n_in_bin
 2000                format(f15.4,i10)
                  endif
               endif
               
c              keep track of maximum in p.d.f.
               if ( n_in_bin .gt. n_max ) then
                  n_max = n_in_bin                         ! largest value of p.d.f.
                  p_max = p_mid                            ! parameter value of largest value
               endif

c              integrate p.d.f to get 68% and 95% confidence ranges
               integrate_up = integrate_up + n_in_bin      ! counts up from 0
               integrate_dn = integrate_dn - n_in_bin      ! counts down from itermax

c              check for crossing "2-sigma" thresholds...
               if (two_sigma_low.eq.0 .and.integrate_up.ge.n_2sig) then
                  two_sigma_low = p_mid
               endif
               if (two_sigma_high.eq.0.and.integrate_dn.le.n_2sig) then
                  two_sigma_high= p_mid
               endif

c              check for crossing "1-sigma" thresholds...
               if (one_sigma_low.eq.0 .and.integrate_up.ge.n_1sig) then
                  one_sigma_low = p_mid
               endif
               if (one_sigma_high.eq.0.and.integrate_dn.le.n_1sig) then
                  one_sigma_high= p_mid
               endif

            enddo

c           Calculate best estimate of posterior 1-sigma uncertainty
            one_sigma_est = abs( one_sigma_high - one_sigma_low ) / 2.d0
            two_sigma_est = abs( two_sigma_high - two_sigma_low ) / 4.d0
            post_sigma(n_p) = ( one_sigma_est + two_sigma_est ) / 2.d0

            if ( lu_print .gt. 0 ) then
c              Then printout pdf statistics...
               write (lu_print,3000) parnames(n_p), p_max, n_max
 3000          format(/1x,a16,': maximum of posteriori distrib ',
     +                'function:',f10.3,' (',i7,')')

               p_middle = 0.5d0 *( one_sigma_low + one_sigma_high )
               write (lu_print,4000) p_middle, one_sigma_low,
     +                               one_sigma_high,one_sigma_est
 4000          format(19x,'68% confidence range: center and edges:',
     +                3f10.3,' (1-sigma=',f10.3,')')

               p_middle = 0.5d0 *( two_sigma_low + two_sigma_high )
               write (lu_print,5000) p_middle, two_sigma_low,
     +                               two_sigma_high,two_sigma_est 
 5000          format(19x,'95% confidence range: center and edges:',
     +                3f10.3,' (1-sigma=',f10.3,')')
            endif

         endif

      enddo

      return
      end

c     =================================================================
      subroutine est_angular_rotation(lu_print, itermax, num_parameters, 
     +           parnames, params_stored, params_best, 
     +           prior_unc, post_unc, post_sigma_Om)

c     Determines posteriori probability distribution functions for 
c     the full angular rotation rate of the Sun: (To + Vsun)/Ro

      implicit real*8 (a-h,o-z)

      real*4        params_stored(11,1000000)
      real*8        prior_unc(11), post_unc(11), params_best(11)
      character*16  parnames(11)

c     Calculate and bin the (To + Vsun)/Ro values, marginalizing over all other params.  

      post_sigma(n_p) = 0.d0

      n_bins    = 200
      par_best  = (params_best(2) + params_best(5))/params_best(1)
      par_delta = 0.15d0 * 1.5d0*post_unc(2)/params_best(1) ! crude est of 15% of post unc.
      par_start = par_best - (n_bins/2)*par_delta - par_delta/2.d0

      if ( lu_print .gt. 0 ) then
         write (lu_print,1000) 
 1000    format(//' Binned probability for (To+Vsun)/Ro')
      endif

      p_max = 0
      n_max = 0

      integrate_up = 0
      integrate_dn = itermax + 1

c     "2-sigma" one-sided threshold (0.5*4.6%  = 2.3%) 
      n_2sig         = 0.023d0 * itermax + 0.5d0
      two_sigma_low  = 0
      two_sigma_high = 0

c     "1-sigma" one-sided threshold (0.5*31.8% = 15.9%) 
      n_1sig         = 0.159d0 * itermax + 0.5d0
      one_sigma_low  = 0
      one_sigma_high = 0

            
      do n_bin = 1, n_bins 

         p0 = par_start + par_delta*n_bin
         p1 = p0 + par_delta

         n_in_bin = 0
         do iter = 1, itermax
            To   = params_stored(2,iter)
            Vsun = params_stored(5,iter)
            Ro   = params_stored(1,iter)
            par_val = (To+Vsun)/Ro                   ! km/s/kpc
            if ( par_val.ge.p0 .and. par_val.lt.p1) then 
               n_in_bin = n_in_bin + 1
            endif
         enddo

         if ( n_in_bin .gt. 0 ) then 
            p_mid = ( p0 + p1 ) / 2.d0
            if ( lu_print .gt. 0 ) then
               write (lu_print,2000) p_mid, n_in_bin
 2000          format(f15.4,i10)
            endif
         endif
               
c        keep track of maximum in p.d.f.
         if ( n_in_bin .gt. n_max ) then
            n_max = n_in_bin                         ! largest value of p.d.f.
            p_max = p_mid                            ! parameter value of largest value
         endif

c        integrate p.d.f to get 68% and 95% confidence ranges
         integrate_up = integrate_up + n_in_bin      ! counts up from 0
         integrate_dn = integrate_dn - n_in_bin      ! counts down from itermax

c        check for crossing "2-sigma" thresholds...
         if (two_sigma_low.eq.0 .and.integrate_up.ge.n_2sig) then
             two_sigma_low = p_mid
         endif
         if (two_sigma_high.eq.0.and.integrate_dn.le.n_2sig) then
             two_sigma_high= p_mid
         endif

c        check for crossing "1-sigma" thresholds...
         if (one_sigma_low.eq.0 .and.integrate_up.ge.n_1sig) then
             one_sigma_low = p_mid
         endif
         if (one_sigma_high.eq.0.and.integrate_dn.le.n_1sig) then
             one_sigma_high= p_mid
         endif

      enddo

c     Calculate best estimate of posterior 1-sigma uncertainty
      one_sigma_est = abs( one_sigma_high - one_sigma_low ) / 2.d0
      two_sigma_est = abs( two_sigma_high - two_sigma_low ) / 4.d0
      post_sigma_Om = ( one_sigma_est + two_sigma_est ) / 2.d0

      if ( lu_print .gt. 0 ) then
c        Then printout pdf statistics...
         write (lu_print,3000) p_max, n_max
 3000    format(/1x,'Omega (km/s/kpc)',': maximum of posteriori ',
     +              'distrib function:',f10.3,' (',i7,')')

         p_middle = 0.5d0 *( one_sigma_low + one_sigma_high )
         write (lu_print,4000) p_middle, one_sigma_low,
     +                                   one_sigma_high
 4000    format(19x,'68% confidence range: center and edges:',
     +          3f10.3)

         p_middle = 0.5d0 *( two_sigma_low + two_sigma_high )
         write (lu_print,5000) p_middle, two_sigma_low,
     +                                   two_sigma_high 
 5000    format(19x,'95% confidence range: center and edges:',
     +                3f10.3)
      endif


      return
      end

c     =================================================================
      subroutine marginalize_burnin ( iterations, num_parameters, 
     +           parnames, params_stored, params_best, 
     +           prior_unc, post_sigma )

c     Determines posteriori probability distribution functions for parameters
c     Stripped down version of "marginalize" used only to find posterior
c     parameter sigmas.   Can't use full subroutine because of dimensioning
c     issues.

      implicit real*8 (a-h,o-z)

c     Note following dimension for "params_stored" for this routine!!
      real*4        params_stored(11,100000)
      real*8        prior_unc(11), params_best(11), post_sigma(11)
      character*16  parnames(11)

      n_bins    = 200
c     "1-sigma" one-sided threshold (0.5*31.8% = 15.9%) 
      n_1sig    = 0.159d0 * iterations + 0.5d0
c     "2-sigma" one-sided threshold (0.5*4.6%  = 2.3%) 
      n_2sig    = 0.023d0 * iterations + 0.5d0

c     Bin the parallax values, which marginalizes (integrates) over the other params.  
      do n_p = 1, num_parameters

c        Check if "solved for" parameter: 
c        if prior uncertainty > 0 (Gaussian prior) or < 0 (for a flat prior)
c        But if prior uncertainty =0, held it constant

         if ( prior_unc(n_p) .ne. 0.d0 ) then

            par_best  = params_best(n_p)
            par_delta = 0.15d0 * post_sigma(n_p)     ! 15% of current posteriori uncertainty

c           May need to change binning intervals to get enough bins with "hits"
c           to determine "sigmas"...
            n_bins_hit = 0
            n_cycles   = 0
            do while ( ( n_bins_hit.lt.10 .or. n_bins_hit.gt.100 ) .and. 
     +                   n_cycles.lt.5 )

               n_cycles   = n_cycles + 1

               p_max = 0
               n_max = 0
               integrate_up = 0
               integrate_dn = iterations + 1
               two_sigma_low  = 0
               two_sigma_high = 0
               one_sigma_low  = 0
               one_sigma_high = 0

c              Adjust binning intervals appropriately for 2nd thru 5th "cycles"...
               if ( n_cycles .gt. 1 ) then
                  if ( n_bins_hit .lt. 10 ) par_delta = par_delta/10.d0
                  if ( n_bins_hit .gt.100 ) par_delta = par_delta*5.d0
               endif
               par_start=par_best-(n_bins/2)*par_delta-par_delta/2.d0

c              Do binning...
               n_bins_hit = 0
               do n_bin = 1, n_bins 

                  p0 = par_start + par_delta*n_bin
                  p1 = p0 + par_delta

                  n_in_bin = 0
                  do iter = 1, iterations
                     par_val = params_stored(n_p,iter)
                     if ( par_val.ge.p0 .and. par_val.lt.p1 ) then 
                        n_in_bin = n_in_bin + 1
                     endif
                  enddo

                  if ( n_in_bin .gt. 0 ) then
                     n_bins_hit = n_bins_hit + 1
                     p_mid      = ( p0 + p1 ) / 2.d0
                  endif
               
c                 keep track of maximum in p.d.f.
                  if ( n_in_bin .gt. n_max ) then
                     n_max = n_in_bin ! largest value of p.d.f.
                     p_max = p_mid ! parameter value of largest value
                  endif

c                 integrate p.d.f to get 68% and 95% confidence ranges
                  integrate_up = integrate_up + n_in_bin ! counts up from 0
                  integrate_dn = integrate_dn - n_in_bin ! counts down from "iterations"

c                 check for crossing "2-sigma" thresholds...
                  if(two_sigma_low.eq.0.and.integrate_up.ge.n_2sig) then
                     two_sigma_low = p_mid
                  endif
                  if(two_sigma_high.eq.0.and.integrate_dn.le.n_2sig)then
                     two_sigma_high= p_mid
                  endif

c                 check for crossing "1-sigma" thresholds...
                  if(one_sigma_low.eq.0.and.integrate_up.ge.n_1sig) then
                     one_sigma_low = p_mid
                  endif
                  if(one_sigma_high.eq.0.and.integrate_dn.le.n_1sig)then
                     one_sigma_high= p_mid
                  endif

               enddo

               if ( n_cycles .gt. 1 ) then
                  write (6,2000) n_p, n_cycles, p_max, n_bins_hit
 2000             format(' Cycling: ',2i5,f15.4,i10)
               endif

            enddo       ! cycles

c           Check n_cycles didn't hit limit
            if ( n_cycles.eq.5 .and. n_bins_hit.lt.10 ) then
               print *,'STOP: Hit limit of 5 cycles adjusting',
     +                 ' binsize for parameter number ',n_p
               STOP
            endif

c           Calculate new estimate of posterior 1-sigma uncertainty
            one_sigma_est = abs( one_sigma_high - one_sigma_low ) / 2.d0
            two_sigma_est = abs( two_sigma_high - two_sigma_low ) / 4.d0
            post_sigma_est= ( one_sigma_est + two_sigma_est ) / 2.d0
c           Update post_sigma smoothly...
            post_sigma(n_p) = (post_sigma(n_p) + post_sigma_est)/2.d0

         endif

      enddo

      return
      end

c======================================================================
      subroutine open_ascii_file ( lu_in, ascii_file, lu_print )

c     Opens ascii file and reads and prints comments (if first character is "!") 
C     Leaves file open, ready to read data.

      character*48       ascii_file

      character*80       comment
      character*1        c_1
      equivalence       (c_1, comment)

      logical            first_comment

C     Open input ascii file...

      first_comment = .true.

      write (lu_print,1000) ascii_file
 1000 format(/' ====================================================',
     +       /' Opening input ascii file "',a32,'"')

      open (unit=lu_in, file=ascii_file, status='old')

C     Read all comment lines (but stop at data)
      write (6,1010) 
 1010 format(' Comment lines follow:')

      c_1 = '!'
      do while ( c_1 .eq. '!' )
	 read  (lu_in,1020) comment
 1020    format(a80)
	 if ( c_1 .eq. '!' ) then
               write (lu_print,1030) comment
 1030          format(1x,a80)
            else
               backspace (unit=lu_in)
	 endif
      enddo

      write (lu_print,1040)     ! print a blank line at end
 1040 format(' End of Comment lines',/1x)

      return
      end

c===================================================================
      subroutine read_one_source ( lu_data, data_file, lu_print, 
     +     src,ra_hhmmss,dec_ddmmss, farnear, par, err_par, 
     +     pm_x, err_pm_x, pm_y, err_pm_y, vhelio, err_v )

      implicit real*8 (a-h,o-z)

      character*48  data_file

      character*12   src, near_far

      call open_ascii_file ( lu_data, data_file, lu_print )

C     Now read in the source information and parallax and motion data...
      read(lu_data,*) src, farnear
      read(lu_data,*) ra_hhmmss
      read(lu_data,*) dec_ddmmss
      read(lu_data,*) par, err_par
      read(lu_data,*) pm_x, err_pm_x
      read(lu_data,*) pm_y, err_pm_y
      read(lu_data,*) vlsr, err_v

c     Archival printout...
      call hmsrad ( ra_hhmmss, ra_rad )
      call dmsrad (dec_ddmmss, dec_rad)
      call radec_to_galactic ( ra_rad, dec_rad,
     +                         ell, bee )

c     Convert VLSR to Heliocentric (in order to remove OLD values
c     of Solar Motion used to define VLSR)
      call standard_vlsr_to_helio ( ra_rad, dec_rad, V_proj )
      vhelio = vlsr + V_proj

      write (lu_print,1000)  src, ra_hhmmss,dec_ddmmss, ell,bee
 1000 format(2x,a12,2f10.1,2f10.3)
      write (lu_print,1010) par, err_par
 1010 format(2x,'Parallax',2f10.3)
      write (lu_print,1020) pm_x, err_pm_x 
 1020 format(2x,'X-motion',2f10.2)
      write (lu_print,1022) pm_y, err_pm_y 
 1022 format(2x,'Y-motion',2f10.2)
      write (lu_print,1030) vhelio, err_v, vlsr
 1030 format(2x,'V_Helio ',2f10.1,' (from V_LSR =',f7.1,')')

      near_far = 'near    '
      if ( farnear .ne. 0.d0 ) near_far = 'far     ' 
      write (lu_print,1040) near_far
 1040 format(2x,'Assumed at ',a5,' kinematic distance')

      close ( unit=lu_data )

      return
      end

c===================================================================
      subroutine control_parameters ( lu_control, lu_print,
     +      num_stages, iterations_per_stage, itermax, n_thin,
     +      err_tol, step_fraction, use_data_type, idum,
     +      num_params, params, prior_unc, post_unc, parnames,
     +      par_virial, dist_min, dist_GC_min, iRC, par_max_unc)

      implicit real*8 (a-h,o-z)

      real*8          params(11), prior_unc(11), post_unc(11)
      character*16    parnames(11), data_type(4)
      logical         err_tol, use_data_type(4), use_flag

      character*48    control_file
      control_file = 'control_file_bayesian.inp'   ! hardwired input control file name

      data_type(1) = 'Parallax data'
      data_type(2) = 'X-motion data'
      data_type(3) = 'Y-motion data'
      data_type(4) = 'V_lsr data'

      call open_ascii_file ( lu_control, control_file, lu_print )

c     Start with important program control parameters:
c     Number of burn-in stages
      read (lu_control,*)  num_stages
  
c     Number of McMC burnin-iterations per stage; beware of dimension limit
      read (lu_control,*)  iterations_per_stage

c     Number of McMC final iterations; beware of dimension limit
      read (lu_control,*)  itermax, n_thin

c     Use "conservative" error-tolerant data uncertainty formulation of 
c     Sivia in "A Bayesian Tutorial"?
      read (lu_control,*)  tolerate, use_data_type
      err_tol = .false.
      if ( tolerate .ne. 0.d0 ) err_tol = .true.

      use_flag = .false.
      do n_dt = 1, 4
         if ( use_data_type(n_dt) ) then
            use_flag = .true.
            write(lu_print,8100) data_type(n_dt)
 8100       format(' Fitting ',a16)   
         endif
      enddo
      if ( use_flag ) then
            write (lu_print,8101)
 8101       format(' ')
         else
            print *,'STOP: No data type specified to fit'
            stop
      endif

c     MCMC parameter step size specified as a fraction of the prior sigma
c     (see Gregory p.319 who cites Gilks, Richardson & Spiegelhalter 1996 
c     which suggests aiming for a ~25% parameter trial acceptance rate for 
c     "high-dimensional models" and ~50% for few model parameters).

      read (lu_control,*)  step_fraction, idum
      if ( idum .gt. 0 ) idum = -idum
      write (lu_print,8111) idum
 8111 format(' RAN1 seed (large negative integer) =',i12,/' ')

c     "Virial" parameter (to be added in quadrature with the Vhelio uncertainty)
c         used to better approximate true kinematic errors;
c     Minumum distance from Sun 
c     Minimum distance from G.G.
c     Rotation curve model flag
      read (lu_control,*)   par_virial, dist_min, dist_GC_min, iRC, 
     +                      par_max_unc
      write (lu_print,8120) par_virial, dist_min, dist_GC_min, iRC,
     +                      par_max_unc
 8120 format(' Virial velocity uncertainty (km/s):',f6.1,
     +      /' Minimum distance from Sun (kpc):   ',f6.1,
     +      /' Minimum distance from G.C. (kpc):  ',f6.1,
     +      /' Rotation curve model flag:         ',i4,
     +      /' Maximum parallax uncertainty (%):  ',f6.1 )
        
c     ============================================================
c     Physical parameters of Galaxy model...
c     Galactic structure parameters
      read (lu_control,*) params(1), prior_unc(1), post_unc(1)  ! Ro [kpc]
      read (lu_control,*) params(2), prior_unc(2), post_unc(2)  ! To [km/s]
      read (lu_control,*) params(3), prior_unc(3), post_unc(3)  ! dT/dr [km/s/kpc]
c     Solar Motion parameters
      read (lu_control,*) params(4), prior_unc(4), post_unc(4)  ! Uo [km/s]
      read (lu_control,*) params(5), prior_unc(5), post_unc(5)  ! Vo [km/s]
      read (lu_control,*) params(6), prior_unc(6), post_unc(6)  ! Wo [km/s]
c     Source peculiar motion in its local Galactocentric frame
      read (lu_control,*) params(7), prior_unc(7), post_unc(7)  ! Ao [km/s]
      read (lu_control,*) params(8), prior_unc(8), post_unc(8)  ! Bo [km/s]
c     Alternative Rotation Curve parameters
      read (lu_control,*) params(9), prior_unc(9), post_unc(9)  ! model dependent
      read (lu_control,*) params(10),prior_unc(10),post_unc(10) ! 
      read (lu_control,*) params(11),prior_unc(11),post_unc(11) ! 

      if ( iRC .ge. 1 ) then
c        Zero dT/dr parameters for all alternative RCs
         params(3)    =   0.d0
         prior_unc(3) =   0.d0
         post_unc(3)  =   0.d0

         if ( iRC .ge. 3 ) then
c           Set parameters for To so that it will have marginalized PDF evaluated
c           For these RC models, To will be updated automatically based on
c              Ro, a1, a2, and a3
c           So, will override user entered values just to be safe
            params(2)    =   1.d0               ! dummy value, not used
            prior_unc(2) = -10.d0               ! flat prior
            post_unc(2)  =  10.d0               ! expected unc for PDF binning
         endif
      endif

      if ( iRC .lt. 3 ) then
c        Zero a1, a2, a3 parameters for "standard RCs"
         do i_p = 9, 11
            params(i_p)    =  0.d0
            prior_unc(i_p) =  0.d0
            post_unc(i_p)  =  0.d0
         enddo
      endif

c     Set parameter names...
      parnames(1)    = 'Ro [kpc]        '
      parnames(2)    = 'To [km/s]       '
      parnames(3)    = 'dT/dr [km/s/kpc]'
      parnames(4)    = 'Uo [km/s]       '
      parnames(5)    = 'Vo [km/s]       '
      parnames(6)    = 'Wo [km/s]       '
      parnames(7)    = 'Ao [km/s]       '
      parnames(8)    = 'Bo [km/s]       '
      parnames(9)    = 'a1              '
      parnames(10)   = 'a2              '
      parnames(11)   = 'a3              '
        
c     Print out parameter values and aprioris...
      num_params = 11                      ! for the 11 parameters above

      write (lu_print,8005)
 8005 format(/' PRIORS...',/' Input parameters:')

      do n = 1, num_params
c        Don't allow negative posteriori uncertainties
         post_unc(n)  = abs( post_unc(n) )

c        Expected posteriori uncertainties must be .le. prior unc (if not flat priors)...
         if ( prior_unc(n) .ge. 0.d0 ) then
            post_unc(n) = min( prior_unc(n), post_unc(n) )
         endif

c        If posteriori unc = 0, set prior unc = 0 (since will not be varied)
         if ( post_unc(n) .eq. 0.d0 ) prior_unc(n) = 0.d0

         write (lu_print,8010) parnames(n), params(n), 
     +                         prior_unc(n), post_unc(n)
 8010    format(4x,a16,2x,'priors: ',f12.3,'  +/-',f9.3,
     +          '; expected posteriori uncertainty +/-',f9.3)
      enddo

      close(unit=lu_control)

      return
      end

c===================================================================
        subroutine get_source_file_names ( lu_control, lu_print, 
     +                  max_num_files,
     +                  data_files, n_files )


        implicit real*8 (a-h,o-z)

        character*48 control_file, file_name, blank48, data_files(2000)
        character*1  c1
        equivalence (file_name,c1)

        blank48 = '                                                '

        control_file = 'source_files.inp'    ! file containing list of source files

        call open_ascii_file ( lu_control, control_file, lu_print )

        n_files = 0
        ieof    = 0
        do while ( ieof.ge.0 .and. n_files.lt.max_num_files )

         read(lu_control,*,iostat=ieof) file_name

c        Check for end of file...
         if ( ieof.ge.0 .and. file_name.ne.blank48 .and.
     +                        c1       .ne. "!"  ) then

            n_files = n_files + 1
            data_files( n_files ) = file_name

         endif

        enddo
        
        print *,' Number of source file names =', n_files
        
        close(unit=lu_control)

        return
        end

c===================================================================
	subroutine calc_residuals ( params, iRC, num_data, 
     +                              data,
     +                              model, resids )
 
C	Calculates models and resduals for entire data set...

	implicit real*8 ( a-h,o-z )

	real*8     params(11)
	real*8	   data(8000), model(8000), resids(8000)
        integer    data_type(8000), source_number(8000)

        real*8  ra(2000),dec(2000),parallaxes(2000),par_errors(2000),
     +          vhelios(2000),err_vhelios(2000),near_fars(2000)
        common /model_info/ source_number, data_type, ra, dec,
     +                    near_fars,parallaxes,par_errors,
     +                    vhelios,err_vhelios
     

        do n = 1, num_data

           id_type      = data_type(n)
           i_src        = source_number(n)

           src_parallax = parallaxes(i_src)
           src_parallax_unc = par_errors(i_src)
           src_v_helio  = vhelios(i_src)
           err_v_helio  = err_vhelios(i_src)

           call calc_model ( params, iRC, id_type, i_src, 
     +          src_parallax, src_parallax_unc,
     +          src_v_helio, err_v_helio,
     +          value ) 

	   model(n)  = value
	   resids(n) = data(n) - model(n)

	enddo
	
	return
	end

c===================================================================
	subroutine calc_model ( params, iRC, id_type, i_src, 
     +             src_parallax, src_parallax_unc,
     +             src_v_helio, err_v_helio,
     +             model_value ) 

c       Calculates a model_value for a single datum:
c       if   id_type=1    parallax                  [mas]
c            id_type=2    Eastward proper motion    [mas/yr]
c            id_type=3    Northward proper motion   [mas/yr]
c            id_type=4    V_helio                   [km/s]

c       Parameter iRC specifies rotation curve model

	implicit real*8 ( a-h, o-z )

        real*8     model_value
        real*8     params(11)

        integer    data_type(8000), source_number(8000)

        real*8  ra(2000),dec(2000),parallaxes(2000),par_errors(2000),
     +          vhelios(2000),err_vhelios(2000),near_fars(2000)
        common /model_info/ source_number, data_type, ra, dec,
     +                    near_fars,parallaxes,par_errors,
     +                    vhelios,err_vhelios


c       Check that data type is 1 to 4
        if ( id_type.lt.1 .or. id_type.gt.4 ) then
           print *,' SUB CALC_MODEL: invalid data type ',id_type
           stop
        endif

c       Get source coordinates from common block
        ra_rad  = ra(  i_src )                     ! radians
        dec_rad = dec( i_src )                     ! radians
        farnear = near_fars( i_src )

        if ( id_type .eq. 1 ) then
c          Parallax datum... 
           call get_kdist ( params, iRC, ra_rad, dec_rad, farnear,
     +                      src_parallax, src_parallax_unc,
     +                      src_v_helio, err_v_helio,
     +                      dist, unc_dist )
           model_value = 1.d0 / dist             ! mas
        endif

        if ( id_type .eq. 2 ) then
c          Eastward proper motion datum... 
           call calc_xy_motions ( params, iRC, 
     +                            ra_rad, dec_rad, 
     +                            src_parallax, src_parallax_unc,
     +                            x_motion, y_motion )
           model_value = x_motion                  ! mas/yr
        endif

        if ( id_type .eq. 3 ) then
c          Northward proper motion datum... 
           call calc_xy_motions ( params, iRC,
     +                            ra_rad, dec_rad, 
     +                            src_parallax, src_parallax_unc,
     +                            x_motion, y_motion )
           model_value = y_motion                  ! mas/yr
        endif

        if ( id_type .eq. 4 ) then
c          V_heliocentric datum... 
           call calc_v_lsr ( params, iRC, 
     +                       ra_rad, dec_rad, 
     +                       src_parallax, src_parallax_unc,
     +                       v_lsr )
           Uo = params(4)      ! new solar motion toward G.C. [km/s]
           Vo = params(5)      ! new solar motion toward Gal. rotation
           Wo = params(6)      ! new solar motion toward NGP
           call new_lsr_to_helio ( ra_rad, dec_rad, Uo, Vo, Wo,
     +                             V_proj )

           model_value = v_lsr + V_proj             ! km/s
        endif

	return 
	end

c====================================================================
	subroutine get_kdist ( params, iRC, ra_rad, dec_rad, farnear,  
     +                         src_parallax, src_parallax_unc, 
     +                         v_helio, err_v,
     +                         dist, unc_dist )

c	Calculates source kinematic distance [kpc] from Galaxy model
c       parameters, source coordinates, Heliocentric velocity.

c       Also, it returns an uncertainty in the kinematic distance,
c       based on the uncertainty in the velocity.

        implicit real*8 (a-h, o-z)

        real*8     params(11)

        integer    data_type(8000), source_number(8000)

        real*8  ra(2000),dec(2000),parallaxes(2000),par_errors(2000),
     +          vhelios(2000),err_vhelios(2000),near_fars(2000)
        common /model_info/ source_number, data_type, ra, dec,
     +                    near_fars,parallaxes,par_errors,
     +                    vhelios,err_vhelios

        pi    = 4.d0*atan(1.d0)
        deg_to_rad = pi / 180.d0

        Ro   = params(1)                    ! kpc
        To   = params(2)                    ! km/s
        dTdr = params(3)                    ! km/s/kpc

        Uo = params(4)                      ! new solar motion toward G.C. [km/s]
        Vo = params(5)                      ! new solar motion toward Gal. rotation
        Wo = params(6)                      ! new solar motion toward NGP

        Ao   = params(7)                    ! km/s
        Bo   = params(8)                    ! km/s

        a1   = params(9)                    ! RC parameters; units
        a2   = params(10)                   ! depend on model
        a3   = params(11)

        unc_dist_min = 0.5d0                ! minimum D_k uncertainty [kpc]

        dist_max_in  = 20.d0                ! maximum D_k for outer Galaxy
        dist_max_out = 10.d0                ! maximum D_k for inner Galaxy

        unc_dist_min = 0.5d0                ! minimum D_k uncertainty [kpc]
        unc_dist_max = 10.d0                ! maximum D_k uncertainty

c       Calculate source's Calactic coordinates...
        call radec_to_galactic (ra_rad, dec_rad, 
     +                          gal_long, gal_lat)

        if ( gal_long.gt.359.9d0 .or. gal_long.lt.0.1d0 ) then 

c             Assume source is Sgr A* (special case)
              dist     = Ro             ! kpc
              unc_dist = unc_dist_min   ! kpc  Unsure what to put here

           else

c             Convert from heliocentric to new "lsr" using solved-for solar motion...
              call new_lsr_to_helio ( ra_rad, dec_rad, Uo, Vo, Wo,
     +                                V_proj )
              v_lsr = v_helio - V_proj

c             Correct for the "global" source peculiar motion parameters
              d    = 1.d0 / src_parallax                  ! kpc
c             Apply bias correction for asymmetric PDF
              call distance_bias ( src_parallax, src_parallax_unc,
     +                             bias_factor )
              d = d / bias_factor

              ell_rad = gal_long * deg_to_rad             ! radians
              sin_l   = sin( ell_rad )
              cos_l   = cos( ell_rad )
              r_sq = Ro**2 + d**2 - 2.d0*Ro*d*cos_l
              r    = sqrt( r_sq )                         ! kpc
              sin_beta =   d  * sin_l     / r
              cos_beta = ( Ro - d*cos_l ) / r
              beta     = atan2( sin_beta, cos_beta )      ! radians
              ellbeta = ell_rad + beta                    ! radians
              v_lsr = v_lsr - (-Ao*sin(ellbeta) + Bo*cos(ellbeta) )   ! km/s

c             Calculate D_k, given Galactic parameters, Vlsr, and far/near
c             code as input.
              call k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3,
     +                             v_lsr,  gal_long, farnear,
     +                             dist )

c             Now estimate kin. dist. uncertainty, from vel uncertainty...
              v_high = v_lsr + err_v
              call k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3,
     +                             v_high, gal_long, farnear, 
     +                             d_vhigh )

              v_low  = v_lsr - err_v
              call k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3,
     +                             v_low,  gal_long, farnear,
     +                             d_vlow )
              unc_dist1 = abs( d_vhigh - d_vlow ) / 2.d0
              unc_dist  = max( unc_dist1, unc_dist_min )

c             Deal with highly uncertain kinematic distances for 
c             sources within 15 degrees of center/anticenter.  
              if ( (gal_long.gt.345.d0 .or.  gal_long.lt.15.d0) .or.
     +             (gal_long.gt.165.d0 .and. gal_long.lt.195.d0) ) then
c                Get a "2nd-opinion" by going off 2-sigma in err_v
                 v_high2 = v_lsr + 2.d0*err_v
                 call k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3,
     +                                v_high2, gal_long, farnear, 
     +                                d_vhigh2 )

                 v_low2  = v_lsr - 2.d0*err_v
                 call k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3,
     +                                v_low2,  gal_long, farnear,
     +                                d_vlow2 )
                 unc_dist2 = abs( d_vhigh2 - d_vlow2 ) / 4.d0
c                Adopt maximum of 1-sigma and 2-sigma estimates
                 unc_dist  = max ( unc_dist1, unc_dist2 )
                 unc_dist  = max( unc_dist, unc_dist_min )
              endif

              if ( gal_long.gt.345.d0 .or.  gal_long.lt.15.d0 ) then
                 if ( dist .lt. dist_min    ) dist = dist_min      ! kpc
                 if ( dist .gt. dist_max_in ) dist = dist_max_in
                 if ( abs(sin_l) .gt. 1.d-06 ) then
cnotneeded          unc_dist = unc_dist / abs(sin_l)          ! inflate unc
                    unc_dist = min( unc_dist, unc_dist_max )  ! but cap unc
                 endif
              endif
              if ( gal_long.gt.165.d0 .and. gal_long.lt.195.d0) then
                 if ( dist .lt. dist_min    ) dist = dist_min      ! kpc
                 if ( dist .gt. dist_max_out) dist = dist_max_out
                 if ( abs(sin_l) .gt. 1.d-06 ) then
cnotneeded          unc_dist = unc_dist / abs(sin_l)          ! inflate unc
                    unc_dist = min( unc_dist, unc_dist_max )  ! but cap unc
                 endif
              endif

        endif
 
	return
	end

c====================================================================
        subroutine k_dist_smooth ( Ro, To, dTdr, iRC, a1, a2, a3, 
     +                             v_lsr, gal_long, farnear,
     +                             D_k )

c	Calculates source kinematic distance [kpc] from Galaxy model
c       parameters, source coordinates, and LSR velocity.

c       This version requires user to specify far/near distance
c       (farnear = 0 for near; =1 for far), in order to avoid
c       discontinuities in D_k(v_lsr)


        implicit real*8 (a-h, o-z)

c       Following appropriate for HMSFRs and current realistic measurements
        d_min = 0.5d0                    ! minimum allowed D_k [kpc]
        d_max = 20.d0                    ! maximum allowed D_k

        pi      = 4.d0*atan(1.d0)
        deg_rad = pi / 180.d0

        sinl   = sin( gal_long * deg_rad )
        cosl   = cos( gal_long * deg_rad )
        Rosinl = Ro*sinl                                   ! kpc
        Rocosl = Ro*cosl

        nf = 2                              ! far distance
        if ( farnear .eq. 0.d0 ) nf = 1     ! near distance

c       Need d to get r in order to get Tr
c       So, iterate until convergence...
        iter = 0
        r    = Ro               ! forces Tr=To to start off
        d_old = 1.d0
        d     = 9.d0
        max_iter = 100
        do while ( abs(d-d_old).gt.0.001 .and. iter.lt.max_iter )

           iter  = iter + 1
           d_old = d

           call get_Tr( Ro, To, dTdr, iRC, a1, a2, a3, r,
     +                  Tr )

           r = Rosinl * ( Tr / (v_lsr + To*sinl) )          ! kpc

c          Check for negative or wild large distances (usually for anti-center sources)
           if ( r .lt. 0.d0 )  r = 30.d0       ! needed to keep d vs v_lsr continuous 
           if ( r .gt. 30.d0 ) r = 30.d0

           rootterm = r**2 - Rosinl**2
           if ( rootterm .lt. 0.d0 ) rootterm = 0.d0

           if ( gal_long.ge.0.d0   .and. gal_long.lt.90.d0 ) then
              if ( nf .eq. 1 ) then
                 d = Rocosl - sqrt( rootterm )
                 if ( v_lsr .le. 0.d0 ) d = d_min  ! for Q1 & near dist, V <0 not allowed
              endif
              if ( nf .eq. 2 ) d = Rocosl + sqrt( rootterm )
           endif

           if ( gal_long.ge.90.d0  .and. gal_long.lt.180.d0 ) then
              d = Rocosl + sqrt( rootterm )
              if ( v_lsr .ge. 0.d0 ) d = d_min
           endif

           if ( gal_long.ge.180.d0 .and. gal_long.lt.270.d0 ) then
              d = Rocosl + sqrt( rootterm )
              if ( v_lsr .le. 0.d0 ) d = d_min
           endif

           if ( gal_long.ge.270.d0 .and. gal_long.le.360.d0 ) then
              if ( nf .eq. 1 ) then 
                 d = Rocosl - sqrt( rootterm )
                 if ( v_lsr .ge. 0.d0 ) d = d_min  ! for Q4 & near dist, V >0 not allowed
              endif
              if ( nf .eq. 2 ) d = Rocosl + sqrt( rootterm )
           endif
              
        enddo

        if ( iter .ge. max_iter ) then
c       Typical problem for ell~180 
c       print *,' K_DIST: max_iter for nf=',nf,': (ell=',gal_long,
c    +          '; D=',d
           d = 3.d0                        ! arbitrary realistic value
        endif

c       Deal with small, zero, or negative distance...
        if ( d .lt. d_min ) d = d_min

c       Check for wild large distances (usually for ell~180)
        if ( d .gt. d_max ) d = d_max

        D_k = d

        return
        end

c====================================================================
        subroutine calc_xy_motions ( params, iRC, ra_rad, dec_rad, 
     +                               src_parallax, src_parallax_unc,
     +                               x_motion, y_motion )

c       Calculates expected proper motion of a source from a model
c       of Galactic rotation

	implicit real*8 ( a-h,o-z )

	real*8     params(11)

        pi        = 4.0d0*atan(1.0d0)
        deg_rad   = pi/180.d0   

c       conversion from km/s/kpc to mas/yr...(0.211 = 1/4.74)
        convert   = 365.2422d0 * 86400.d0 / 1.496d8

c       Galactic rotation parameters
        Ro   = params(1)                 ! kpc
        To   = params(2)                 ! km/s
        dTdr = params(3)                 ! km/s/kpc

c       Solar Motion parameters
        Uo   = params(4)                 ! km/s
        Vo   = params(5)                 ! km/s
        Wo   = params(6)                 ! km/s

c       Source "global" peculiar motion parameters
        Ao   = params(7)                 ! km/s
        Bo   = params(8)                 ! km/s

c       Alternative rotation curve parameters (for iRC>=3)
        a1   = params(9)
        a2   = params(10)
        a3   = params(11)

c       Get Galactic coordintes of source
        call radec_to_galactic ( ra_rad, dec_rad,
     +                           ell, bee )

        ell_rad = ell * deg_rad
        sin_l   = sin( ell_rad )
        cos_l   = cos( ell_rad )

        bee_rad = bee * deg_rad
        sin_b   = sin( bee_rad )
        cos_b   = cos( bee_rad )

c       Calculate Theta(r) for source at its measured distance
        d      = 1.d0 / src_parallax          ! kpc
c       Apply bias correction for asymmetric PDF
        call distance_bias ( src_parallax, src_parallax_unc,
     +                       bias_factor )
        d = d / bias_factor 

c       Project into Galactic plane...
        d_proj = d * cos_b                    ! kpc

        r_sq   = Ro**2 + d_proj**2 - 2.d0*Ro*d_proj*cos_l
        r_proj = sqrt( r_sq )                 ! kpc

cold        Tr   = To + dTdr * ( r_proj - Ro )    ! km/s
        call get_Tr ( Ro, To, dTdr, iRC, a1, a2, a3, r_proj,
     +                Tr )

c       Special case of Sgr A*...can't use Tr calculation
        if ( (ell.gt.359.9d0 .and. ell.lt.0.1d0) .and.
     +        abs(bee).lt.0.1d0 )    then
c          Have Sgr A* case
           Tr = 0.d0
        endif

c       Calculate expected motion in Galactic coordinates
        sin_beta =   d_proj  * sin_l     / r_proj
        cos_beta = ( Ro - d_proj*cos_l ) / r_proj
        beta     = atan2( sin_beta, cos_beta )      ! radians
        ellbeta  = ell_rad + beta                   ! radians

        vl_source = Tr*cos_b * cos( ellbeta )       ! km/s
        vl_sun    = ( To + Vo )*cos_l - Uo*sin_l    ! km/s
        v_ell     = vl_source - vl_sun              ! km/s

c       Add expected "global" peculiar source motion
        v_ell_pec = -Ao*cos(ellbeta) - Bo*sin(ellbeta)  ! km/s
        v_ell_pec_proj = v_ell_pec * cos_b

        v_ell     = v_ell + v_ell_pec               ! km/s

c       Since we can't really predict pec vertical motion, assume zero
        v_bee     = -Wo                             ! km/s

        pm_ell_kmskpc = v_ell / d                   ! km/s/kpc
        pm_ell    = pm_ell_kmskpc * convert         ! mas/yr

        pm_bee_kmskpc = v_bee / d
        pm_bee    = pm_bee_kmskpc * convert         ! mas/yr

c       Add in motions (as offsets over 1000 years) 
c       and get new (ra,dec).
        ell_1 = ell + pm_ell / 3600.d0              ! deg
        bee_1 = bee + pm_bee / 3600.d0              ! deg

        call galactic_to_radec ( ell, bee,
     +                           ra_0, dec_0 )

        call galactic_to_radec ( ell_1, bee_1,
     +                           ra_1, dec_1 )

c       Difference (ra,dec)'s to get motion
        del_ra = ( ra_1 - ra_0 ) * (180.d0/pi)*3600.d0   ! arcsec/1000yr
        x_motion  = del_ra * cos(dec_rad)                ! mas/yr

        del_dec=( dec_1 - dec_0 )* (180.d0/pi)*3600.d0   ! arcsec/1000yr
        y_motion  = del_dec                              ! mas/yr

c        print *,' CALC_XY_MOT: (mu_x,mu_y) =', x_motion, y_motion

        return
        end

c====================================================================
        subroutine calc_v_lsr( params, iRC, 
     +                         ra_rad, dec_rad, 
     +                         src_parallax, src_parallax_unc,
     +                         v_lsr )

c       Calculates expected LSR velocity of a source from a model
c       of Galactic rotation.

	implicit real*8 ( a-h,o-z )

	real*8     params(11)

        pi        = 4.0d0*atan(1.0d0)
        deg_rad   = pi/180.d0   

c       Galactic rotation parameters
        Ro   = params(1)                 ! kpc
        To   = params(2)                 ! km/s
        dTdr = params(3)                 ! km/s/kpc

c       Solar Motion parameters
        Uo   = params(4)                 ! km/s
        Vo   = params(5)                 ! km/s
        Wo   = params(6)                 ! km/s

c       Global source peculiar motion parameters
        Ao   = params(7)                 ! km/s
        Bo   = params(8)                 ! km/s

c       Get rotation curve parameters (used if iRC>=3)
        a1   = params(9)
        a2   = params(10)
        a3   = params(11)

c       Get Galactic coordintes of source
        call radec_to_galactic ( ra_rad, dec_rad,
     +                           ell, bee )

        ell_rad = ell * deg_rad
        sin_l   = sin( ell_rad )
        cos_l   = cos( ell_rad )

        bee_rad = bee * deg_rad
        sin_b   = sin( bee_rad )
        cos_b   = cos( bee_rad )

c       Calculate Theta(r) for source at its measured distance
        d      = 1.d0 / src_parallax          ! kpc
c       Apply bias correction for asymmetric PDF
        call distance_bias ( src_parallax, src_parallax_unc,
     +                        bias_factor )
        d = d / bias_factor 

c       Project to Galactic Plane
        d_proj = d * cos_b

        r_sq   = Ro**2 + d_proj**2 - 2.d0*Ro*d_proj*cos(ell_rad)
        r_proj = sqrt( r_sq )                 ! kpc

        call get_Tr ( Ro, To, dTdr, iRC, a1, a2, a3, r_proj,
     +                Tr )

c       Special case of Sgr A*...can't use Tr calculation
        if ( (ell.gt.359.9d0 .and. ell.lt.0.1d0) .and.
     +        abs(bee).lt.0.1d0 )    then
c          Have Sgr A* case
           Tr = 0.d0
        endif

c       Calculate expected motion in Galactic coordinates
        sin_beta =   d_proj  * sin_l     / r_proj
        cos_beta = ( Ro - d_proj*cos_l ) / r_proj
        beta     = atan2( sin_beta, cos_beta )      ! radians
        ellbeta  = ell_rad + beta                   ! rad
        vr_source = Tr*cos_b * sin( ellbeta )       ! km/s

c       Add in expected "global" source peculiar motions
c       to model observed radial velocities
c       changed Bo sign here for version 6d
        vr_pec    = -Ao*sin(ellbeta) + Bo*cos(ellbeta) ! km/s
        vr_pec_proj = vr_pec * cos_b

        vr_source = vr_source + vr_pec_proj

c       Calculate circular motion of LSR, corrected for projection
c       of source to plane...
        vr_lsr    = To * cos_b * sin_l              ! km/s

        v_lsr     = vr_source - vr_lsr              ! km/s

        return
        end

c       ========================================================
        subroutine distance_bias(src_parallax,src_parallax_unc,
     +                           bias_factor )

c       Given a source's parallax and parallax uncertainty,
c       use a lookup table to calculate the bias in the 
c       mean distance relative to 1/parallax (the peak of the
c       distance PDF.   

c       Used program parallax_vs_distance_pdfs_3.f to generate
c       the table of biases.
       
	implicit real*8	(a-h,o-z)

        real*8 fract_unc(26), bias(26)

        n_table = 26

        fract_unc(1)  = 0.00d0
        bias(1)       = 1.0000d0
        fract_unc(2)  = 0.02d0
        bias(2)       = 1.0004d0
        fract_unc(3)  = 0.04d0
        bias(3)       = 1.0016d0
        fract_unc(4)  = 0.06d0
        bias(4)       = 1.0036d0
        fract_unc(5)  = 0.08d0
        bias(5)       = 1.0065d0
        fract_unc(6)  = 0.10d0
        bias(6)       = 1.0103d0
        fract_unc(7)  = 0.12d0
        bias(7)       = 1.0151d0
        fract_unc(8)  = 0.14d0
        bias(8)       = 1.0209d0
        fract_unc(9)  = 0.16d0
        bias(9)       = 1.0279d0
        fract_unc(10) = 0.18d0
        bias(10)      = 1.0362d0
        fract_unc(11) = 0.20d0
        bias(11)      = 1.0462d0
        fract_unc(12) = 0.22d0
        bias(12)      = 1.0582d0
        fract_unc(13) = 0.24d0
        bias(13)      = 1.0726d0
        fract_unc(14) = 0.26d0
        bias(14)      = 1.0898d0
        fract_unc(15) = 0.28d0
        bias(15)      = 1.1100d0
        fract_unc(16) = 0.30d0
        bias(16)      = 1.1328d0
        fract_unc(17) = 0.32d0
        bias(17)      = 1.1577d0
        fract_unc(18) = 0.34d0
        bias(18)      = 1.1840d0
        fract_unc(19) = 0.36d0
        bias(19)      = 1.2109d0
        fract_unc(20) = 0.38d0
        bias(20)      = 1.2374d0
        fract_unc(21) = 0.40d0
        bias(21)      = 1.2628d0
        fract_unc(22) = 0.42d0
        bias(22)      = 1.2866d0
        fract_unc(23) = 0.44d0
        bias(23)      = 1.3083d0
        fract_unc(24) = 0.46d0
        bias(24)      = 1.3276d0
        fract_unc(25) = 0.48d0
        bias(25)      = 1.3444d0
        fract_unc(26) = 0.50d0
        bias(26)      = 1.3592d0

c       Fractional uncertainty in parallax
        fraction = src_parallax_unc / src_parallax

c       Find table entries that straddle "fraction"
        n_lower = 0
        do n = 1, n_table
           if ( fraction .gt. fract_unc(n) ) n_lower = n
        enddo
        n_upper = n_lower + 1

c       Catch sources with fractional uncertainty > 0.5
        if ( n_upper .gt. n_table ) then
           print *,'STOP: cannot debias source with ',
     +             'fractional parallax uncertainty >0.5'
           STOP500
        endif

        if ( n_upper .le. 1 ) bias_factor = 1.d0

c       Interpolate for bias_factor
        if ( n_upper .ge. 2 ) then
           delta   = fraction - fract_unc(n_lower)
           f_diff  = fract_unc(n_upper) - fract_unc(n_lower)
           b_diff  = bias(n_upper) - bias(n_lower)
           bias_factor = bias(n_lower) + b_diff*(delta/f_diff)
        endif

c        write (99,1000) src_parallax,src_parallax_unc,bias_factor
c 1000   format(3f10.4)

        return
        end

c       ========================================================
        subroutine get_Tr( Ro, To, dTdr, iRC, a1, a2, a3, r,
     +                     Tr )

c       Given rotation curve info (various params &iRC) and  
c       projected Galactocentric radius (r), 
c       calculates circular rotation speed (Tr)

c       Models with iRC < 3, require Ro, To, dTdr;
c       those with  iRC >= 3 require Ro, a1, a2, a3

	implicit real*8	(a-h,o-z)

        Tr = 0.d0

c       -------------------------------------------------
c       Following models use Ro, To and dTdr to calculate Tr:

c       Standard parametrized model (iRC = 0)
        if ( iRC .eq. 0 ) then
           Tr   = To + dTdr * ( r - Ro )       ! km/s
        endif

c       Clemens-10kpc rotation curve (iRC = 1)
        if ( iRC .eq. 1 ) then
           r_scaled = r / (Ro/10.d0)
           call clemens_10_250 ( r_scaled, Tr )
c          Also, scale rotations speeds to parameter value of To
           Tr = Tr*(To/250.d0)
        endif

c       Clemens-8.5kpc rotation curve (iRC = 2)
c       Clemens Ro=8.5, To=220 curve, which is a bit flatter than Clemens-10 for R<Ro.  
        if ( iRC .eq. 2 ) then
           r_scaled = r / (Ro/8.5d0)
           call clemens_8p5_220 ( r_scaled, Tr )
c          Also, scale rotations speeds to parameter value of To
           Tr = Tr*(To/220.d0)
        endif

c       ------------------------------------------------
c       Following models use Ro, a1, a2, a2 to calculate Tr:
c       NB: calculate To in order to have McMC trials

c       Brand-Blitz formulation (iRC = 3 )
        if ( iRC .eq. 3 ) then
c          Brand-Blitz (1993, A&A, 275, 67) "like" rotation curve
c          parameterization.  Use km/s units for a1 and a3 
c          instead of dimensionless params multiplied by "To".
           rho = r/Ro
           To = a1 + a3                        ! km/s
           Tr = a1 * rho**a2 + a3              ! km/s
        endif

c       "Alternative" formulation (iRC = 4 )
        if ( iRC .eq. 4 ) then
           rho = (r/Ro) - 1.d0
           To  = a1                            ! km/s
           Tr  = a1 + a2*rho + a3*rho**2       ! km/s
        endif

c       "Universal" or disk-halo formulation (iRC = 5 )
        if ( iRC .eq. 5 ) then
c          Disk plus halo parameterization of rotation curve...
c          see http://ned.ipac.caltech.edu/level5/March01/Battaner/node7.html
c          referenced as following  Persic and Salucci (1997)
           beta   = 0.72d0      ! 0.72 + 0.44 ln(L/L*)
           V2_Ropt= a1**2       ! a1 = V(Ropt)  
           Ropt   = a2 * Ro     ! a2 = Ropt/Ro; (NB: Ropt=3.2Rd encloses 83% of light)
                                ! where Rd=scale length ~ 4.5 kpc (but often ~3 kpc used)
                                ! P.C. van der Kruit & K.C. Freeman ARAA, 2011, 49, 301
           a      = a3          ! "a"     = 1.5(L/L*)^(1/5)
c          So, for L = L*, beta=0.72; a1~To; a2~3.2*(3.0to4.5)/8.3 = 1.2to1.7; a3~1.5

c          Calculate Tr...
           rho    = r/Ropt
c          Exponential thin disk component...
           top = 1.97d0*rho**1.22
           bot = (rho**2 + 0.78d0**2)**1.43d0
           V2_disk = V2_Ropt * beta * top / bot
c          Halo component...
           top = rho**2
           bot = rho**2 + a**2
           V2_halo = V2_Ropt * (1.d0-beta) * (1.d0+a**2)*top/bot
c          Total...
           Tr = sqrt( V2_disk + V2_halo )      ! km/s 

c          Also, calculate To...
           rho    = Ro/Ropt
c          Exponential thin disk component...
           top = 1.97d0*rho**1.22
           bot = (rho**2 + 0.78d0**2)**1.43d0
           V2_disk = V2_Ropt * beta * top / bot
c          Halo component...
           top = rho**2
           bot = rho**2 + a**2
           V2_halo = V2_Ropt * (1.d0-beta) * (1.d0+a**2)*top/bot
c          Total...
           To = sqrt( V2_disk + V2_halo )      ! km/s 

        endif

c       Check that we have calculated a value for Tr
        if ( Tr .le. 0.d0 ) then

c          print *,'DEBUG: get_Tr: ',Ro,To,dTdr,r,Tr

c           print *,' Sub get_Tr: Tr not calculated; iRC=',iRC,': STOP'
c           stop
        endif

        return
        end

      subroutine clemens_10_250( R, Theta )

c     calculates Theta(R) from Clemens 1985 paper

      implicit real*8 (a-h,o-z)

      real*8  A_old(7), B_old(6), C_old(8)

c     Table 5 values for Ro=10, To=250  valid for  R/Ro < 0.09
      data    A_old /     0.0d0, 2618.86d0, -11486.7d0,
     +                27263.8d0,-36117.9d0,  24768.0d0,
     +                 -6819.d0 /

c     Table 5 values for Ro=10, To=250  valid for 0.09 < R/Ro < 0.45
      data    B_old / 319.8354d0, -194.3951d0, 155.30181d0,
     +               -63.00332d0, 12.142556d0, -0.869573d0 /

c     Table 5 values for Ro=10, To=250  valid for 0.45 < R/Ro < 1.6
      data    C_old / -1931.3363d0, 1768.47176d0, -606.461977d0,
     +                112.102138d0,-11.9698505d0,   0.7367828d0,
     +               -0.02423453d0, 0.00032952d0 /

c     Table 5 value  for Ro=10, To=250  valid for 1.6 < R/Ro 
      data    D_old / 264.76d0 /

      Ro = 10.0d0                             ! kpc

      theta = 0.d0
      ratio = R/Ro

      if ( ratio.lt.0.09d0 ) then
         do i = 1, 7
            theta = theta + A_old(i)*R**(i-1)
         enddo
      endif

      if ( ratio.ge.0.09d0 .and. ratio.lt.0.45d0 ) then
         do i = 1, 6
            theta = theta + B_old(i)*R**(i-1)
         enddo
      endif

      if ( ratio.ge.0.45d0 .and. ratio.le.1.6d0 ) then
         do i = 1, 8
            theta = theta + C_old(i)*R**(i-1)
         enddo
      endif

      if ( ratio.gt.1.6d0 ) theta = D_old

      return
      end

c     =================================================================
      subroutine clemens_8p5_220( R, Theta )

c     calculates Theta(R) from Clemens 1985 paper

      implicit real*8 (a-h,o-z)

      real*8  A_new(7), B_new(6), C_new(8)

c     Table 5 values for Ro=8.5, To=220  valid for  R/Ro < 0.09
      data    A_new /      0.d0, 3069.81d0, -15809.8d0,
     +                43980.1d0,-68287.3d0,  54904.0d0,
     +               -17731.0d0 /

c     Table 5 values for Ro=8.5, To=220  valid for 0.09 < R/Ro < 0.45
      data    B_new / 325.0912d0, -248.1467d0, 231.87099d0,
     +              -110.73531d0, 25.073006d0, -2.110625d0 /

c     Table 5 values for Ro=8.5, To=220  valid for 0.45 < R/Ro < 1.6
      data    C_new / -2342.6564d0, 2507.60391d0,-1024.068760d0,
     +                224.562732d0,-28.4080026d0,   2.0697271d0,
     +               -0.08050808d0, 0.00129348d0 /

      data    D_new / 234.88d0 /

      Ro =  8.5d0                             ! kpc

      theta = 0.d0
      ratio = R/Ro

      if ( ratio.lt.0.09d0 ) then
         do i = 1, 7
            theta = theta + A_new(i)*R**(i-1)
         enddo
      endif

      if ( ratio.ge.0.09d0 .and. ratio.lt.0.45d0 ) then
         do i = 1, 6
            theta = theta + B_new(i)*R**(i-1)
         enddo
      endif

      if ( ratio.ge.0.45d0 .and. ratio.le.1.6d0 ) then
         do i = 1, 8
            theta = theta + C_new(i)*R**(i-1)
         enddo
      endif

      if ( ratio.gt.1.6d0 ) theta = D_new

      return
      end


c       ======================================================
	subroutine calc_sigsq (num_data, num_solved_for,
     +			resids, res_err, 
     +			sigsq_old, sigsq, sigma_pdf)

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8		resids(8000), res_err(8000)

	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 1, num_data
c          Exclude parallax data...
           if ( mod(i_d,4) .ne. 1 ) then
              sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
           endif
	enddo
c       excluding parallax data
	num_deg_freedom = 3.d0*num_data/4.d0 - num_solved_for

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

	return
	end

c     ========================================================
      subroutine radec_to_galactic ( ra_rad, dec_rad, 
     +                               ell_II, bee_II )

c     inputs:    ra_rad     (radians)  in J2000  !!!
c                dec_rad    (radians)
c     outputs:   ell_II (galactic longitude in degrees)
c                bee_II (galactic latitude in degrees)

      implicit real*8 (a-h,o-z)

      pi        = 4.0d0*atan(1.0d0)
      deg_rad   = pi/180.d0     ! convert degrees to radians
      hr_deg    = 15.d0

c     Galactic plane-pole data for B1950...
c      g_alpha  = 282.25d0      ! RA(?) of Galactic pole (deg)
c      g_dec    = 62.6d0        ! Dec(?) of Galactic pole (deg)
c      g_33     = 33.0d0        ! l(II) - l(I) (?) (deg)
c     J2000 Galactic plane-pole data...
      g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0     ! hrs
      g_alpha   = (g_pole_ra + 6.d0) * hr_deg                 ! deg
                       ! = 282.859507
      g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0      ! deg
      g_dec     = 90.d0 - g_pole_dec                          ! deg
                       ! = 62.8716631
      g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                       ! = 32.932
                       ! Determined g_33 by requiring that the
                       ! B1950 defined G.C., when precessed to J2000,
                       ! still gives l=0.0 and b=0.0.  (works to 0.0001 deg)

c     Useful quantities...
      cos_a     = cos( ra_rad - g_alpha*deg_rad )
      sin_a     = sin( ra_rad - g_alpha*deg_rad )
      cos_d     = cos( dec_rad )
      sin_d     = sin( dec_rad )
      cos_gd    = cos( g_dec*deg_rad )
      sin_gd    = sin( g_dec*deg_rad )

c     Now calculate central galactic coordinates given ra,dec...
      sin_bII   = sin_d*cos_gd - cos_d*sin_a*sin_gd
      bee       = asin( sin_bII )
      cos_bII   = cos( bee )

      cos_ell   = cos_d*cos_a / cos_bII
      sin_ell   = ( cos_d*sin_a*cos_gd + sin_d*sin_gd ) / cos_bII
      ell       = atan2( sin_ell,cos_ell )

      ell_II    = ell/deg_rad + g_33
      if ( ell_II .lt.  0.d0 ) ell_II = ell_II + 360.d0
      if ( ell_II .gt.360.d0 ) ell_II = ell_II - 360.d0

      bee_II    = bee/deg_rad

      return
      end

      subroutine galactic_to_radec ( ell_deg, bee_deg,
     +                               ra_rad, dec_rad )
 
c       ========================================================
c       converts galactic to RA, Dec coordinates

c       enter galactic coords in degrees
c       outputs equatorial coordinates in radians

        implicit real*8 (a-h,o-z)

        pi      = 4.0*atan(1.0) 
        deg_rad = pi/180.        
        hr_deg  = 15.d0

c       J2000 Galactic plane-pole data...
        g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0 ! hrs
        g_alpha   = (g_pole_ra + 6.d0) * hr_deg ! deg
                      ! = 282.859507

        g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0 ! deg
        g_dec     = 90.d0 - g_pole_dec ! deg
                      ! = 62.8716631
        g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                      ! = 32.932
                      ! Determined g_33 by requiring that the
                      ! B1950 defined G.C., when precessed to J2000,
                      ! still gives l=0.0 and b=0.0.  (works to 0.0001 deg)

c        Now convert to RA, Dec...
         cos_b = cos( bee_deg * deg_rad )         
         sin_b = sin( bee_deg * deg_rad ) 

         cos_l = cos( ell_deg * deg_rad )
         sin_l = sin( ell_deg * deg_rad )

         sin_l33 = sin( (ell_deg - g_33) * deg_rad )
         cos_l33 = cos( (ell_deg - g_33) * deg_rad )

         sin_62  = sin( g_dec * deg_rad )
         cos_62  = cos( g_dec * deg_rad )

         sin_d = cos_b*sin_l33*sin_62 + sin_b*cos_62
         dec   = asin( sin_d )
         cos_d = cos( dec )   
         dec_deg = dec / deg_rad 

         temp = cos_b * sin_l33 - sin_d * sin_62
         sin_a282 = temp / ( cos_d * cos_62 )
         cos_a282 = cos_b * cos_l33 / cos_d

         a282_rad = atan2( sin_a282, cos_a282 )           ! radians
         a282_deg = a282_rad / deg_rad                    ! degrees

         ra_deg   = a282_deg + g_alpha                    ! degrees
         if ( ra_deg .lt. 0.d0   ) ra_deg = ra_deg + 360.d0
         if ( ra_deg .gt. 360.d0 ) ra_deg = ra_deg - 360.d0

         ra_rad = ra_deg * deg_rad

         dec_rad = dec_deg * deg_rad

         return
         end

c       ========================================================
      subroutine standard_vlsr_to_helio ( ra_rad, dec_rad, V_proj )

c     Enter (RA,Dec) = (ra_rad,dec_rad)  in radians

c     Gives V_proj   =  V(Heliocentric) - V(LSR)  in km/s
c         ie, add V_proj to V(LSR) to get V(Heliocentric)

c     Uses Standard Solar Motion (old values)

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )

c     LSR defined by removing peculiar Solar Motion of
c     20.0 km/s toward 18.0 hours, +30.0 degrees (RA,Dec)
c     (Actually defined in 1900.0 system, technically should
c     precess to J2000)

      Vsun     = 20.d0                  ! km/s

c     RA_Vsun  = 18.d0                  ! hours      1900.0
c     Dec_Vsun = 30.d0                  ! degrees
      RA_Vsun  = 18.0640d0              ! hours      J2000.0
      Dec_Vsun = 30.0024d0              ! degrees

      cos_RA = cos( RA_Vsun * pi / 12.d0 )
      sin_RA = sin( RA_Vsun * pi / 12.d0 )

      cos_Dec= cos( Dec_Vsun* pi /180.d0 )
      sin_Dec= sin( Dec_Vsun* pi /180.d0 )

c     In equatorial Cartesian frame the Solar Motion is

      Xo = Vsun * cos_RA * cos_Dec       ! toward RA=0h in equat plane
      Yo = Vsun * sin_RA * cos_Dec       ! toward RA=6h in equat plane
      Zo = Vsun * sin_Dec                ! toward Dec=90 (North)

c     Projection of Solar Motion toward a source at RA,Dec

      cos_alpha = cos( ra_rad )
      sin_alpha = sin( ra_rad )

      cos_delta = cos( dec_rad )
      sin_delta = sin( dec_rad )

      V_proj = -Xo*cos_alpha*cos_delta - Yo*sin_alpha*cos_delta -
     +          Zo*sin_delta

      return
      end

c       ========================================================
      subroutine new_lsr_to_helio ( ra_rad, dec_rad, Uo, Vo, Wo,
     +                              V_proj )

c     Enter ra_rad   R.A. [radians]
c           dec_rad  Dec. [radians]
c           Uo       solar motion inward toward G.C.       [km/s]
c           Vo       solar motion toward Galactic rotation [km/s]
c           Wo       solar motion toward N.G.P.            [km/s]

c     Gives V_proj   =  V(Heliocentric) - V(LSR)
c         ie, add V_proj to V(LSR) to get V(Heliocentric)

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )

      call radec_to_galactic ( ra_rad, dec_rad,
     +                         ell_deg, bee_deg )

      ell_rad = ell_deg * pi / 180.d0
      sin_l   = sin( ell_rad )
      cos_l   = cos( ell_rad )

      bee_rad = bee_deg * pi / 180.d0
      sin_b   = sin( bee_rad )
      cos_b   = cos( bee_rad )

      V_sun   =  (Vo*sin_l + Uo*cos_l)*cos_b + Wo*sin_b

      V_proj  = -V_sun

      return
      end


c       ========================================================
        subroutine hmsrad ( hms, hms_rad )

C       Converts hhmmss.sss format to radians

        implicit real*8 (a-h,o-z)

        pi = 4.d0 * atan(1.d0)

        xhms = abs(hms)

        ih = xhms/10000.d0 + 1.d-10
        im = dmod(xhms,10000.d0)/100.d0 + 1.d-10
        s  = dmod(xhms,100.d0)

        hms_rad = dfloat(ih) + dfloat(im)/60.d0 + s/3600.d0
        hms_rad = hms_rad * pi / 12.d0
        if ( hms .lt. 0.d0 )  hms_rad = -hms_rad

        return
        end
c       ========================================================
        subroutine dmsrad ( dms, dms_rad )

c       Converts hhmmss.sss format to radians

        implicit real*8 (a-h,o-z)

        pi = 4.d0 * atan(1.d0)
        deg_rad = pi / 180.d0

        xdms = abs(dms)

        id = xdms/10000.d0
        im = dmod(xdms,10000.d0)/100.d0
        s  = dmod(xdms,100.d0)

        dms_rad = dfloat(id) + dfloat(im)/60.d0 + s/3600.d0
        dms_rad = dms_rad * deg_rad
        if ( dms .lt. 0.d0 )  dms_rad = -dms_rad

        return
        end
c       =========================================================
        subroutine radians_to_hhmmss ( ra_rad, ra_hhmmss)

c       Input :  ra_rad       in radians
c       Output:  ra_hhmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

        pi = 4.d0 * atan(1.d0)
        rad_to_hr = 12.d0/ pi

        ra_hr = ra_rad * rad_to_hr

        ihr  = ra_hr
        imin = (ra_hr - ihr*1.d0) * 60.d0
        sec  = (ra_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0

        ra_hhmmss  = ihr*10000.d0 + imin*100.d0 + sec

        return
        end
c       =========================================================
        subroutine radians_to_ddmmss ( dec_rad, dec_ddmmss )

c       Input :  dec_rad       in radians
c       Output:  dec_ddmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

        pi = 4.d0 * atan(1.d0)
        rad_deg = 180 / pi 

        if ( dec_rad .ge. 0.d0 ) then
             dec_sign =  1.d0
          else
             dec_sign = -1.d0
        endif

        dec_deg = abs( dec_rad ) * rad_deg

        ideg = dec_deg
        imin = (dec_deg - ideg*1.d0) * 60.d0
        sec  = (dec_deg - ideg*1.0d0 - imin/60.0d0) * 3600.d0

        dec_ddmmss  = ideg*10000.d0 + imin*100.d0 + sec
        dec_ddmmss  = dec_ddmmss * dec_sign

        return
        end

c===============================================================
      subroutine gauss_ran1 ( idum, u1,u2, g1,g2 )

C     Routine to return two uniformly distributed random numbers (0-1)
C     called "u1" and "u2", and two Gaussian distributed random numbers
C     called "g1" and "g2" with zero mean and unity sigma.

c     Uses Numerical Recipies "RAN1" and the Box-Muller transformation

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )
      twopi = 2.d0 * pi

      u1 = ran1( idum )
      u2 = ran1( idum )

      amp = sqrt( -2.0d0*log(u1) )
      arg = twopi*u2
      g1 = amp*cos(arg)
      g2 = amp*sin(arg)
        
      return
      end

c     ===========================================================
      function ran1( idum )

c     RAN1 random number generator from Numerical Recipies

      implicit real*8 (a-h,o-z)

      parameter (ia=16807, im=2147483647, am=1.d0/im, iq=127773,
     +           ir=2836, ntab=32, ndiv=1+(im-1)/ntab, 
     +           eps=1.2d-07, rnmx=1.d0-eps )

      integer iv(ntab)

      save    iv, iy

      data    iv/ntab*0/, iy/0/

c     Initialize...
      if ( idum.le.0 .or. iy.eq.0 ) then
         idum = max(-idum,1)
         do j = ntab+8,1,-1
            k = idum / iq
            idum = ia*(idum - k*iq) - ir*k
            if ( idum .lt. 0    ) idum = idum + im
            if ( j    .le. ntab ) iv(j)= idum
         enddo
         iy = iv(1)
      endif

      k = idum / iq

      idum = ia*(idum - k*iq) - ir*k
      if ( idum .lt. 0 ) idum = idum + im

      j = 1 + iy/ndiv
      iy    = iv(j)
      iv(j) = idum

      ran1 = min( am*iy, rnmx )

      return
      end

c     ===========================================================
      subroutine correlation_matrix (p,num_global,itermax,lu_print)

      implicit real*8 (a-h,o-z)

      real*4     p(11,1000000)

      real*8     cor(11,11)

      do k = 1, num_global

        do j = k, num_global

          if ( k .ne. j ) then    

            n     = 0
            x_sum = 0.d0
            x_ssq = 0.d0
            y_sum = 0.d0
            y_ssq = 0.d0
            xy_sum= 0.d0

            do n = 1, itermax

               x = p(k,n)
               y = p(j,n)

               x_sum = x_sum + x
               x_ssq = x_ssq + x**2

               y_sum = y_sum + y
               y_ssq = y_ssq + y**2

               xy_sum= xy_sum + x*y
                     
            enddo

c           correlation
            top  =  n*xy_sum - x_sum*y_sum

            bot2 = (n*x_ssq - x_sum**2)*(n*y_ssq - y_sum**2)

            if ( bot2 .gt. 0.d0 ) then
               r_xy = top / sqrt(bot2)
                 else
               r_xy = 0.d0
            endif

            cor(k,j) = r_xy

           else

            cor(k,j) = 1.0d0

          endif

          cor(j,k) = cor(k,j)

        enddo

      enddo


      write (lu_print,1000) 
 1000 format(/'  Correlation Matrix:',/)

      do k = 1, num_global

            write (lu_print,2000) ( cor(k,j), j=1,num_global )
 2000       format(20f7.3)

      enddo

      return
      end
