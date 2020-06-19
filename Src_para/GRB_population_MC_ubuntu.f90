PROGRAM Core
  USE Physics
  USE Cosmology
  USE ParameterFile
  USE Stats
  USE MPI
  USE mtmod
  IMPLICIT None

  ! --------------------------------- !
  ! ----- Variables of the code ----- !
  ! --------------------------------- !
 
  LOGICAL                  :: Mode_test  = .FALSE.                            ! Waits for 'enter' to proceed to next GRB
  LOGICAL                  :: reprise = .FALSE.                               ! Reprise mode
  INTEGER                  :: run_mode = 0                                    ! Run mode : 0 (default) is one run mode, 1 is param search, etc...
  INTEGER, PARAMETER       :: one_run = 0                                     ! Value for one run mode
  INTEGER, PARAMETER       :: param_search = 1                                ! Value for param search mode
  INTEGER, PARAMETER       :: post_process = 2                                ! Value for post processing mode of previously-run good models
  INTEGER, PARAMETER       :: MCMC = 3                                        ! Value for MCMC mode
  INTEGER                  :: hist_flag = 0                                   ! Flag to tell MonteCarlo routine what to generate (0: just Chi2 hists, 1: also prop hists, 2: save everything)
  CHARACTER(*), PARAMETER  :: InitFile = '../Input_para/GRB_pop.init'         ! Name of input file
  CHARACTER(Len=255)       :: path                                            ! Path for output
  CHARACTER(Len=255)       :: param_search_path='../Input_para/'              ! Path for parameter search file
  INTEGER,      PARAMETER  :: verbose    = 0                                  ! Verbosity of the code (0, 1 or 2)
  LOGICAL,      PARAMETER  :: Save_all_GRB = .TRUE.                          ! Writes properties for every GRB in GRB_prop
  LOGICAL                  :: SHOW_READ  = .FALSE.                            ! Verbosity of the ReadInitFile subroutine
  CHARACTER(*), PARAMETER  :: Format_RIF = '(A,A13,A22,A14,ES12.5,A)'         ! Format for ReadInitFile verbosity
  LOGICAL                  :: NaNtest_prop = .FALSE.                          ! To test if any properties are NaN  !! only works during one_run mode
  LOGICAL                  :: NaNtest_hist = .FALSE.                          ! To test if any histograms are NaN  !! only works during one_run mode
  INTEGER                  :: Nb_lines = 0                                    ! Number of lines in reprise file
  INTEGER                  :: starting_i = 0                                  ! Indice at which to start the MCMC loop
  INTEGER                  :: Nb_GRB                                          ! Number of GRBs 
  INTEGER                  :: n_call = 0                                      ! If n_call == 0, will generate histograms, else will just reset them to 0
  INTEGER                  :: pp_flag = 0                                     ! Flag for post-processing mode to exit cleanly
 
  ! ---------- RNG variables ---------- !
  INTEGER, PARAMETER                     :: N_proc_max = 8                       ! Maximum number of processors (8 by default)
  INTEGER, DIMENSION(1:3,0:N_proc_max-1) :: TabInit_ijk_Kiss                     ! Table for initatlizing KISS generator
  INTEGER, DIMENSION(1:3,0:N_proc_max-1) :: TabSave_ijk_Kiss                     ! Table for saving KISS generator
  INTEGER, DIMENSION(0:N_proc_max-1)     :: TabInit_seed_MT                      ! Table for initializing MT19937 generator
  INTEGER                                :: iKiss, jKiss, kKiss                  ! v.a. uniforme : KISS generator
  INTEGER                                :: iKiss_GRB, jKiss_GRB, kKiss_GRB      ! v.a. uniforme : KISS generator
  INTEGER                                :: iKiss_save, jKiss_save, kKiss_save   ! to save v.a.
  INTEGER                                :: MT_seed                              ! MT19937 seed
  INTEGER                                :: Kiss_rng = 0, MT19937 = 1            ! Different cases for RNG
  INTEGER                                :: RNG = 999                            ! Random number generator (Kiss_rng or MT)

  ! ----------- MPI variables ---------- !
  INTEGER :: nb_procs, rank, code, name_length                                ! Variables used to store number of procs, rank of each proc etc..
  INTEGER, PARAMETER :: master_proc = 0                                       ! Rank of the master processor that does the Chi2 calculations
  CHARACTER(len=1)   :: str_rank                                              ! Characeter version of the rank (used to write the BAT6 redshift distribution)

  ! --------- GRB samples and constraints --------- !
  ! ----------------------------------------------- !

  ! --- Samples --- !
  INTEGER, PARAMETER                        :: N_Samples          = 10        ! Number of samples
  INTEGER                                   :: i_Sample                       ! indice for sample loops
  INTEGER, PARAMETER                        :: Sample_Intrinsic   = 0         ! indice for intrinsic sample    (=all GRBs)
  INTEGER, PARAMETER                        :: Sample_Kommers     = 1         ! indice for Kommers sample      (50-300 keV & trigger efficiency from Kommers et al. 2000)
  INTEGER, PARAMETER                        :: Sample_Preece      = 2         ! indice for Preece sample       (50-300 keV & Peak flux > 5     ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_Stern       = 3         ! indice for Stern sample        (50-300 keV)
  INTEGER, PARAMETER                        :: Sample_SWIFTweak   = 4         ! indice for SWIFT sample        (15-150 keV & Peak flux > 0.01  ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_SWIFT       = 5         ! indice for SWIFT sample        (15-150 keV & Peak flux > 0.2   ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_SWIFTbright = 6         ! indice for SWIFT sample        (15-150 keV & Peak flux > 1     ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_HETE2       = 7         ! indice for HETE2 sample        (2-10 keV & 30-400 keV & Peak flux > 1 ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_EpGBM       = 8         ! indice for GBM sample          (50-300 keV & Peak flux >= 0.9  ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_eBAT6       = 9         ! indice for eBAT6 sample        (15-150 keV & Peak flux >= 2.6  ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_SVOM        = 10        ! indice for SVOM sample         (4-150 keV & detailed threshold)

  ! INTEGER, PARAMETER                        :: Sample_HETE2FRE    = 6         ! indice for FREGATE sample      (30-400 keV & Peak flux > 1     ph/cm2/s)
  ! INTEGER, PARAMETER                        :: Sample_HETE2WXM    = 7         ! indice for WXM sample          (2-10   keV & Peak flux > 1     ph/cm2/s)
  ! INTEGER, PARAMETER                        :: Sample_SWIFTBand   = 8         ! formule de band a faire plus tard
  LOGICAL(8), DIMENSION(0:N_Samples)        :: Sample_Included = .FALSE.      ! Table to include sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: TabSample_name                 ! Table for output filename of each sample 
  REAL(8), DIMENSION(0:N_Samples)           :: Threshold                      ! Table for threshold for each sample
  REAL(8), DIMENSION(0:N_Samples)           :: Prob_det                       ! Table for detection probability of instrument for each sample
   
  ! --- Instruments --- !
  INTEGER, PARAMETER                            :: N_Instruments       = 5         ! Number of instruments
  INTEGER                                        :: i_Instrument                    ! indice for instrment loops
  CHARACTER(Len=150), DIMENSION(0:N_Instruments) :: TabInstrument_name              ! Table for output filename of each instrument 
  INTEGER, PARAMETER                            :: Instrument_BATSE    = 1         ! indice for BATSE   instrument (50-300 keV)
  INTEGER, PARAMETER                            :: Instrument_BAT      = 2         ! indice for SWIFT   instrument (15-150 keV)
  INTEGER, PARAMETER                            :: Instrument_FREGATE  = 3         ! indice for FREGATE instrument (30-400 keV)
  INTEGER, PARAMETER                            :: Instrument_WXM      = 4         ! indice for WXM     instrument ( 2-10  keV)
  INTEGER, PARAMETER                            :: Instrument_ECLAIRs  = 5         ! indice for ECLAIRs instrument ( 4-150 keV)
  LOGICAL(8), DIMENSION(0:N_Instruments)        :: Instrument_Included = .FALSE.   ! Table to include sample
  REAL(8), DIMENSION(1:N_Instruments)           :: TabEmin                         ! Table for lower energy limit of instrument 
  REAL(8), DIMENSION(1:N_Instruments)           :: TabEmax                         ! Table for higher energy limit of instrument 
  REAL(8), DIMENSION(1:N_Instruments)           :: Peakflux_Instrument             ! Table for Peak fluxes for each instrument 


  ! --- SVO:/ECLAIRs --- !
  CHARACTER(len=150)                   :: PathSVOM = './SVOM/'
  CHARACTER(len=15)                    :: ExtEclairs
  REAL(8), PARAMETER                   :: nsigmasECLAIRs = 6.5d0
  INTEGER                              :: NeffECLAIRS = 985
  INTEGER                              :: NbkgECLAIRS = 985
  REAL(8), DIMENSION(1:985)            :: TABeffECLAIRsE, TABeffECLAIRsA
  REAL(8), DIMENSION(1:985)            :: TABbkgECLAIRsE, TABbkgECLAIRsB
  INTEGER                              :: istartECLAIRs = 0, iendECLAIRs = 0
  REAL(8)                              :: bkgECLAIRsB1 = 0.d0
  INTEGER, PARAMETER                   :: NoffECLAIRs = 44
  REAL(8), DIMENSION(-NoffECLAIRs:NoffECLAIRs, -NoffECLAIRs:NoffECLAIRs)  ::TABoffECLAIRs,TABomegaECLAIRs
  REAL(8)                              :: omegaECLAIRs = 0.d0
  REAL(8)                              :: Delta_t_pflx = 1.d0
  
  ! --- Constraints --- !
  INTEGER, PARAMETER                     :: N_Constraints       = 6           ! Number of Constraints
  INTEGER                                :: i_Constraint                      ! indice for Constraint loops
  INTEGER, PARAMETER                     :: Constraint_Kommers  = 1           ! indice for Kommers Constraint (LogN LogP  Kommers et al. 2000, Table 2) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_Preece   = 2           ! indice for Preece  Constraint (LogN LogEp Preece  et al. ????, Table ?) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_Stern    = 3           ! indice for Stern   Constraint (LogN LogP  Stern   et al. 2001, Table ?) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_HETE2    = 4           ! indice for X-ray Flash fraction constraint (Frac_XRFHETE2, sigma_XRFHETE2) Requires HETE2FRE and HETE2WXM
  INTEGER, PARAMETER                     :: Constraint_EpGBM    = 5           ! indice for Ep constraint (from GBM catalog, Gruber at al. 2014) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_eBAT6    = 6           ! indice for extended BAT6 constraint (from eBAT6 catalog, Pescalli at al. 2016) Requires Swift
  LOGICAL(8), DIMENSION(0:N_Constraints) :: Constraint_Included = .FALSE.     ! Table to include constraint
  LOGICAL(8), DIMENSION(0:N_Constraints) :: Constraint_save     = .FALSE.     ! Table to save constraint
  CHARACTER(Len=150), DIMENSION(0:N_Constraints) :: TabConstraint_name         ! Table for output filename of each constraint
  
  ! --- Histogram limits --- !
  INTEGER, PARAMETER           :: N_prop         = 6                          ! Number of properties (used to generate the same number of histograms)
  REAL(8), DIMENSION(1:N_prop) :: TabHistlim_inf = 0.d0                       ! Table for inferior limits to the various histograms
  REAL(8), DIMENSION(1:N_prop) :: TabHistlim_sup = 0.d0                       ! Table for superior limits to the various histograms
  INTEGER, PARAMETER           :: Prop_LogL      = 1                          ! indice for luminosity property
  INTEGER, PARAMETER           :: Prop_L         = 1                          ! indice for luminosity property (same as LogL, used in different instances for my sanity of mind)
  INTEGER, PARAMETER           :: Prop_z         = 2                          ! indice for redshift property
  INTEGER, PARAMETER           :: Prop_Ep        = 3                          ! indice for Peak energy property
  INTEGER, PARAMETER           :: Prop_LogP      = 4                          ! indice for Peak flux property
  INTEGER, PARAMETER           :: Prop_alpha     = 5                          ! indice for alpha property 
  INTEGER, PARAMETER           :: Prop_beta      = 6                          ! indice for beta property  
  REAL(8), PARAMETER           :: z_maximum      = 20.d0                      ! maximum redshift of GRBs



  ! ---------------------------------------------------------- !
  ! --- Source properties / GRB properties in source frame --- !
  ! ---------------------------------------------------------- !
  
  ! --------- Luminosity [erg/s] --------- !
  ! -------------------------------------- !

  ! Assumptions

  INTEGER                   :: Model_Lum         = 0                          ! Luminosity Model indice
  INTEGER, PARAMETER        :: Model_LumFix      = 0                          ! Fixed luminosity.   Parameters : L0
  INTEGER, PARAMETER        :: Model_LumPL       = 1                          ! Power-law.          Parameters : Lmin, Lmax, slope
  INTEGER, PARAMETER        :: Model_LumBPL_evol = 2                          ! Broken power-law.   Parameters : Lmin, Lmax, Lbreak, slopeL, slopeH
  INTEGER, PARAMETER        :: Model_LumPL_evol  = 3                          ! Evolving power-law. Parameters : Lmin, Lmax, slope, k_evol
  INTEGER, PARAMETER        :: Model_LumSch      = 4                          ! Schechter function. Parameters : Lmin, Lstar, slope
  
  INTEGER, PARAMETER               :: NParam_Lum       = 10                   ! Maximum number of parameters for luminosity model
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum     = -999.d0              ! Table of Parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum_old = -999.d0              ! Table of new Parameters for the Luminosity Models (used in MCMC)
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum_min = -999.d0              ! Table of minima for parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum_max = -999.d0              ! Table of maxima for parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabBestParam_Lum = -999.d0              ! Table for the best parameters of the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: Step_Lum         = -999.d0              ! Table of step for the Luminosity Models (used in MCMC)
  INTEGER                          :: lum_explore = 0                         ! Flag for exploration of parameter space of Luminosity model
  
  
  INTEGER, PARAMETER        :: Param_Lum_L0      = 1                           ! index for L0
  INTEGER, PARAMETER        :: Param_Lum_Lmin    = 2                           ! index for Lmin
  INTEGER, PARAMETER        :: Param_Lum_Lmax    = 3                           ! index for Lmax
  INTEGER, PARAMETER        :: Param_Lum_Lbreak  = 4                           ! index for Lbreak
  INTEGER, PARAMETER        :: Param_Lum_slope   = 5                           ! index for slope
  INTEGER, PARAMETER        :: Param_Lum_slopeL  = 6                           ! index for slope(low lum)
  INTEGER, PARAMETER        :: Param_Lum_slopeH  = 7                           ! index for slope(high lum)
  INTEGER, PARAMETER        :: Param_Lum_k_evol = 8                           ! index for slope of Lum evolution : (1+z)**k_evol

  ! Histogram

  INTEGER, PARAMETER                    :: N_L = 50                           ! Number of bins in TabLogL
  REAL(8), DIMENSION(0:N_L)             :: TabLogL                            ! Table equally spaced in Log(L) 
  REAL(8), DIMENSION(0:N_Samples,1:N_L) :: TabHistLogL = 0.d0                 ! Table for Histogram of L
  REAL(8), DIMENSION(0:N_Samples,1:N_L) :: TabHistLogL_master = 0.d0          ! Table for master Histogram of L
  INTEGER                               :: iBinL                              ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Lsave = .FALSE.                ! Saves (TRUE) Luminosity histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: LFile                          ! Outputfile name for Luminosity for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: LErrorFile                     ! Outputfile name for Luminosity Error for each sample


  ! Random draw
  INTEGER, PARAMETER        :: M_L = 1000                                     ! Precision of Distribution Function integration (reference=1000)(only used if not analytical)
  REAL(8), DIMENSION(0:M_L) :: TabLogL_CDF                                    ! Table for Distribution Function of L
  REAL(8), DIMENSION(0:M_L) :: TabFctDistrLum                                 ! Table for Distribution Function of L

  ! ----------- Redshift ----------- !
  ! -------------------------------- !

  ! Assumptions

  INTEGER                   :: Model_z         = 0                            ! Redshift Model indice
  INTEGER, PARAMETER        :: Model_zFix      = 0                            ! Fixed redshift.       Parameters : z0
  INTEGER, PARAMETER        :: Model_zUniform  = 1                            ! Uniform rate.         Parameters : zmax 
  INTEGER, PARAMETER        :: Model_zSH       = 2                            ! Springel & Hernquist. Parameters : zmax, zm, a, b 
  INTEGER, PARAMETER        :: Model_zDaigne   = 3                            ! Daigne et al. 2006.   Parameters : a, b, c, d
  INTEGER, PARAMETER        :: Model_z_evol    = 4                            ! Redshift evolution    Parameters : zmax, zm, a, b zeta
  INTEGER, PARAMETER        :: Model_zLi       = 5                            ! Li 2008.              Parameters : a, b, c, d
  INTEGER, PARAMETER        :: Model_zPesc     = 6                            ! Pescalli et al. 2016  Parameters : None
  INTEGER, PARAMETER        :: Model_zBPL      = 7                            ! Broken Power Law      Parameters : zmax, zm, a, b
  INTEGER, PARAMETER        :: Model_zBExp     = 8                            ! Broken Exponential    Parameters : zmax, zm, a, b 

  INTEGER, PARAMETER             :: NParam_z       = 10                       ! Maximum number of parameters for redshift model  
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z     = -999.d0                  ! Table of Parameters for the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z_old = -999.d0                  ! Table of new Parameters for the Redshift Models (used in MCMC)
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z_min = -999.d0                  ! Table of minima for parameters for the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z_max = -999.d0                  ! Table of maxima for parameters for the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: TabBestParam_z = -999.d0                  ! Table for the best parameters of the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: Step_z         = -999.d0                  ! Table of step for Parameters for the Redshift Models (used in MCMC)
  INTEGER                        :: z_explore = 0                             ! Flag for exploration of parameter space of redshift
  
  INTEGER, PARAMETER        :: Param_z_z0      = 1                            ! index for z0
  INTEGER, PARAMETER        :: Param_z_zmax    = 2                            ! index for zmax
  INTEGER, PARAMETER        :: Param_z_zm      = 3                            ! index for zm   (S&H, _evol)
  INTEGER, PARAMETER        :: Param_z_a       = 4                            ! index for a    (S&H, Daigne, Li)
  INTEGER, PARAMETER        :: Param_z_b       = 5                            ! index for b    (S&H, Daigne, Li)
  INTEGER, PARAMETER        :: Param_z_c       = 6                            ! index for c    (Daigne, Li)
  INTEGER, PARAMETER        :: Param_z_d       = 7                            ! index for d    (Daigne, Li)
  INTEGER, PARAMETER        :: Param_z_zeta    = 8                            ! index for zeta (z_evol)

  ! Histogram 

  INTEGER, PARAMETER                    :: N_z = 40                           ! Number of bins in Tabz
  REAL(8), DIMENSION(0:N_z)             :: Tabz                               ! Table equally spaced in z
  REAL(8), DIMENSION(0:N_Samples,1:N_z) :: TabHistz = 0.d0                    ! Table for Histogram of z
  REAL(8), DIMENSION(0:N_Samples,1:N_z) :: TabHistz_master = 0.d0             ! Table for master Histogram of z
  INTEGER                               :: iBinz                              ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: zsave = .FALSE.                ! Saves (TRUE) redshift histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: zFile                          ! Outputfile name for Redshift for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: zFile_cumul                    ! Outputfile name for cumulative Redshift distribution for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: zErrorFile                     ! Outputfile name for Redshift Error for each sample


  ! Random draw

  REAL(8), DIMENSION(0:INT(SIZE(Tabprecisez)-1)) :: TabFctDistrz              ! Table for Distribution Function of z (Tabprecisez comes from cosmology.f90)
  INTEGER                                        :: izmax                     ! Index for zmax

  ! ----------- Spectral shape ----------- !
  ! -------------------------------------- !

  ! Assumptions

  INTEGER                   :: Model_Spec        = 0                          ! Spectral Model indice
  INTEGER, PARAMETER        :: Model_SpecBPLFix  = 0                          ! Fixed BPL.  Parameters : alpha, beta
  INTEGER, PARAMETER        :: Model_SpecBPLK    = 1                          ! BPL with alpha,beta = Kaneko distribution (involves random draw)
  INTEGER, PARAMETER        :: Model_SpecBandFix = 2                          ! Fixed Band. Parameters : alpha, beta
  INTEGER, PARAMETER        :: Model_SpecBandK   = 3                          ! Band with alpha ,beta = Kaneko distribution (involves random draw)
  INTEGER, PARAMETER        :: Model_SpecBandD   = 4                          ! Band with alpha, beta = Daigne distribution (involves random draw)
  INTEGER, PARAMETER        :: Model_SpecBandGBM = 5                          ! Band with alpha, beta from GBM catalog histograms (involves random draw)

  INTEGER, PARAMETER                :: NParam_spec   = 6                      ! Maximum number of parameters for spectral shape model  
  REAL(8), DIMENSION(1:NParam_spec) :: TabParam_Spec = -999.d0                ! Table of parameters for the Spectral shape model
  REAL(8), DIMENSION(1:NParam_spec) :: TabParam_Spec_min = -999.d0            ! Table of minima for parameters for the Spectral shape model
  REAL(8), DIMENSION(1:NParam_spec) :: TabParam_Spec_max = -999.d0            ! Table of maxima for parameters for the Spectral shape model
  INTEGER                           :: spec_explore = 0                       ! Flag for exploration of parameter space of spectral model
  
  INTEGER, PARAMETER        :: Param_spec_alpha = 1                           ! Index for alpha
  INTEGER, PARAMETER        :: Param_spec_beta  = 2                           ! Index for beta

  
  ! Histogram
  
  INTEGER, PARAMETER                :: N_GBM_alpha = 116                         ! Number of bins in GBM alpha histogram
  INTEGER, PARAMETER                :: N_GBM_beta  = 200                         ! Number of bins in GBM beta histogram
  REAL(8), DIMENSION(0:N_GBM_alpha) :: TabFctDistrGBM_alpha = 0.d0               ! Table for Distribution Function of alpha, used for random draws.
  REAL(8), DIMENSION(0:N_GBM_beta)  :: TabFctDistrGBM_beta = 0.d0                ! Table for Distribution Function of beta, used for random draws.
  REAL(8), DIMENSION(0:N_GBM_alpha) :: TabGBM_alpha = 0.d0                       ! Table for values of alpha, used for random draws.
  REAL(8), DIMENSION(0:N_GBM_beta)  :: TabGBM_beta = 0.d0                        ! Table for values of beta, used for random draws.
  INTEGER, PARAMETER                :: N_spec_a = 30                             ! Number of bins in TabSpec_a (for alpha)
  INTEGER, PARAMETER                :: N_spec_b = 30                             ! Number of bins in TabSpec_b (for  beta)
  REAL(8), DIMENSION(0:N_spec_a)    :: TabSpec_a                                 ! Table equally spaced in alpha
  REAL(8), DIMENSION(0:N_spec_b)    :: TabSpec_b                                 ! Table equally spaced in beta
  REAL(8), DIMENSION(0:N_Samples,1:N_spec_a) :: TabHistalpha = 0.d0              ! Table for Histogram of alpha
  REAL(8), DIMENSION(0:N_Samples,1:N_spec_a) :: TabHistalpha_master = 0.d0       ! Table for master Histogram of alpha
  REAL(8), DIMENSION(0:N_Samples,1:N_spec_b) :: TabHistbeta = 0.d0               ! Table for Histogram of beta
  REAL(8), DIMENSION(0:N_Samples,1:N_spec_b) :: TabHistbeta_master = 0.d0        ! Table for master Histogram of beta
  INTEGER                           :: iBina                                     ! Indice for the value drawn
  INTEGER                           :: iBinb                                     ! Indice for the value drawn

  
  LOGICAL,           DIMENSION(0:N_Samples) :: Specsave = .FALSE.             ! Saves (TRUE) Spec histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: SpecFile_a                     ! Outputfile name for alpha for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: SpecFile_b                     ! Outputfile name for beta for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: SpecErrorFile_a                ! Outputfile name for error on alpha for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: SpecErrorFile_b                ! Outputfile name for error on beta for each sample
  



  ! --------- Peak Energy [keV] --------- !
  ! ------------------------------------- !

  ! Assumptions

  INTEGER                    :: Model_Ep          = 0                         ! Peak Energy Model indice
  INTEGER, PARAMETER         :: Model_EpFix       = 0                         ! Fixed Peak Energy (source frame).       Parameters : Ep0
  INTEGER, PARAMETER         :: Model_EpLogNormal = 1                         ! Log-normal distribution (source frame). Parameters : Ep0, sigmaLog
  INTEGER, PARAMETER         :: Model_EpY         = 2                         ! Yonetoku relation TO BE IMPLEMENTED
  INTEGER, PARAMETER         :: Model_EpAmati     = 3                         ! Amati-like relation (source frame)      Parameters : Ep0, sigmaLog, alpha_amati, (Ep = Ep0 * (L/L0)**alpha_amati

  INTEGER, PARAMETER              :: NParam_Ep        = 6                     ! Maximum number of parameters for model of Ep
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep      = -999.d0               ! Table of Parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep_old  = -999.d0               ! Table of new Parameters for the Peak Energy Models (used in MCMC)
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep_min  = -999.d0               ! Table of minima for parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep_max  = -999.d0               ! Table of maxima for parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabBestParam_Ep  = -999.d0               ! Table for the best parameters of the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Ep) :: Step_Ep          = -999.d0               ! Table of step for Parameters for the Peak Energy Models (used in MCMC)
  INTEGER                         :: Ep_explore = 0                           ! Flag for exploration of parameter space of Peak energy
  
  INTEGER, PARAMETER         :: Param_Ep_Ep0         = 1                      ! index for Ep0
  INTEGER, PARAMETER         :: Param_Ep_sigmaLog    = 2                      ! index for sigmaLog
  INTEGER, PARAMETER         :: Param_Ep_L0          = 3                      ! index for L0 (for Amati)
  INTEGER, PARAMETER         :: Param_Ep_alpha_amati = 4                      ! index for alpha_amati 


  ! Histogram 

  INTEGER, PARAMETER                     :: N_Ep = 20                         ! Number of bins in TabLogEp
  REAL(8), DIMENSION(0:N_Ep)             :: TabLogEp                          ! Table equally spaced in Log(Ep)
  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEp = 0.d0               ! Table for Histogram of Log(Ep)
  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEp_master = 0.d0        ! Table for master Histogram of Log(Ep)
  INTEGER                                :: iBinEp                            ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Epsave = .FALSE.               ! Saves (TRUE) Peak Energy histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: EpFile                         ! Outputfile name for Peak Energy for each sample
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: EpErrorFile                    ! Outputfile name for Peak Energy Error for each sample


  ! ----------------------------------------- !
  ! ---------- Observed properties ---------- !
  ! ----------------------------------------- !



  ! ------ Observed peak energy [keV] ------ !
  ! ---------------------------------------- !

  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEpobs = 0.d0             ! Table for Histogram of Log(Epobs)
  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEpobs_master = 0.d0      ! Table for master Histogram of Log(Epobs)
  INTEGER                                :: iBinEpobs                          ! Indice for the value drawn

  ! ---- Observed peak flux [ph/cm^2/s] ---- !
  ! ---------------------------------------- !

  REAL(8), DIMENSION(1:N_Samples)       :: Peakflux                           ! Table for Peak fluxes for each sample                 

  ! Histogram

  INTEGER, PARAMETER                    :: N_P = 40                           ! Number of bins in TabP
  REAL(8), DIMENSION(0:N_P)             :: TabLogP                            ! Table equally spaced in Log(P) 
  REAL(8), DIMENSION(1:N_Samples,1:N_P) :: TabHistLogP = 0.d0                 ! Table for Histogram of peak flux P for each sample
  REAL(8), DIMENSION(1:N_Samples,1:N_P) :: TabHistLogP_master = 0.d0          ! Table for Histogram of peak flux P for each sample
  INTEGER, DIMENSION(1:N_Samples)       :: iBinPeakflux                       ! Indices for values drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Psave = .FALSE.                ! Saves (TRUE) Peak Flux histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: PFile                          ! Outputfile name for Peak Flux for each sample  
  CHARACTER(Len=150), DIMENSION(0:N_Samples) :: PErrorFile                     ! Outputfile name for Peak Flux Error for each sample


  ! --------------------------------------- !
  ! ------ Observational constraints ------ !
  ! --------------------------------------- !

  ! --- BATSE23 : peak flux [ph/cm2/s] distribution Kommers et al. 2000 --- !
  ! ----------------------------------------------------------------------- !

  ! Histogram

  INTEGER, PARAMETER           :: N_Komm                 = 25                     ! Number of bins in TabP23
  REAL(8), DIMENSION(1:N_Komm) :: TabHistKomm_P23        = 0.d0                   ! Table for Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Komm) :: TabHistKomm_P23_master = 0.d0                   ! Table for Master Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Komm) :: TabHistKomm_DRDP       = 0.d0                   ! Table for Histogram of Rate ( = TabHistP23(i)/(Nb_GRB * (P23(i)-P23(i-1))) )
  REAL(8), DIMENSION(1:N_Komm) :: TabHistKomm_DRDPerr    = 0.d0                   ! Table for error of Histogram of Rate ( = TabHistP23(i) / (Nb_GRB * (P23(i)-P23(i-1))) )
  REAL(8)                      :: k_Kommers              = 0.d0                   ! Normalization coefficient 
  INTEGER                      :: iBinKomm                                        ! Index for histogram

  CHARACTER(Len=150)            :: KommFile       = "Kommers_constraint"           ! Outputfile name for Kommers
  CHARACTER(Len=150)            :: KommErrorFile  = "Kommers_constrainterr"        ! Outputfile name for Kommers Error

  ! Observational data

  REAL(8), DIMENSION(0:N_Komm) :: TabKomm_P23, TabKomm_LogP23              ! Table 2 of Kommers et al. 2000 
  REAL(8), DIMENSION(1:N_Komm) :: TabHistKomm_P23obs = 0.d0                ! Table for observed BATSE23 Histogram of P23 (data from Kommers et al 2000)
  REAL(8), DIMENSION(1:N_Komm) :: TabKomm_DRobs                            ! Table for observed BATSE23 Rate [GRB/yr/sr] (Kommers et al 2000)
  REAL(8), DIMENSION(1:N_Komm) :: TabKomm_DRDPobs                          ! Table for observed BATSE23 Rate per Peak flux [(GRB/yr/sr) / (ph/sm2/s)] (Kommers et al 2000)
  REAL(8), DIMENSION(1:N_Komm) :: TabKomm_DRobserr, TabKomm_DRDPobserr     ! Table for error on observed BATSE23 rate (Kommers et al 2000)


  ! --- BATSE23 : peak energy [keV] distribution Preece --- !
  ! ------------------------------------------------------- !
  
  ! Histogram

  INTEGER, PARAMETER             :: N_Preece                = 10                             ! Number of bins in TabPreece_Ep
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Ep        = 0.d0                           ! Table for Histogram of Peak energy
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Ep_master = 0.d0                           ! Table for Histogram of Peak energy of master proc
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Eperr     = 0.d0                           ! Table for Histogram of Peak energy error
  REAL(8)                        :: k_Preece                = 0.d0                           ! Normalization coefficient
  INTEGER                        :: iBinPreece                                               ! Index for histogram

  CHARACTER(Len=150)              :: PreeceFile       = "Preece_constraint"                   ! Outputfile name for generated Ep histogram
  CHARACTER(Len=150)              :: PreeceErrorFile  = "Preece_constrainterr"                ! Outputfile name for Error on generated Ep histogram 

  ! Observational data

  REAL(8), DIMENSION(0:N_Preece) :: TabPreece_Ep, TabPreece_LogEp                            ! Table from Preece
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Epobs    = 0.d0                            ! Table for observed BATSE23 Histogram of Ep (data from Preece)
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Epobserr = 0.d0                            ! Table for observed BATSE23 Histogram error of Ep (data from Preece)

  ! --- BATSE23 : LogN LogP distribution Stern --- !
  ! ---------------------------------------------- !
  
  ! Histogram

  INTEGER, PARAMETER            :: N_Stern                 = 27                              ! Number of bins in TabStern_P23
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23        = 0.d0                            ! Table for Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23_master = 0.d0                            ! Table for Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23err     = 0.d0                            ! Table for the normalized residuals (Model-Obs)/Error
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23_forlnL = 0.d0                            ! Table for Histogram of Peak flux used in lnL (not corrected for useful time and dlogP)
  REAL(8)                       :: k_Stern = 0.d0                                            ! Normalization coefficient
  REAL(8)                       :: k_Stern_forlnL = 0.d0                                     ! Normalization coefficient for lnL calculation
  INTEGER                       :: iBinStern                                                 ! Index for histogram

  CHARACTER(Len=150)              :: SternFile       = "Stern_constraint"                     ! Outputfile name for generated P23 
  CHARACTER(Len=150)              :: SternErrorFile  = "Stern_constrainterr"                  ! Outputfile name for Error on generated P23 

 ! Observational data

  REAL(8), DIMENSION(0:N_Stern) :: TabStern_P23, TabStern_LogP23                                 ! Table from Stern
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obs        = 0.d0                             ! Table for LogN LogP from BATSE23 (data from Stern)
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obs_forlnL = 0.d0                             ! Table for LogN LogP from BATSE23 (data from Stern)
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obserr     = 0.d0                             ! Table for LogN LogP error from BATSE23 (data from Stern)
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obserr_forlnL = 0.d0                          ! Table for LogN LogP error from BATSE23 (data from Stern)
  REAL(8)                       :: Delta_t_Stern              = 9.1d0                            ! Duration of the observation of the Stern sample (in years)
  REAL(8)                       :: Omega_div_4Pi              = 0.67d0                           ! Fraction of the solid angle of sky observed by BATSE (for the Stern sample)
  
  ! --- HETE2 : X-ray Flash fraction --- !
  ! ------------------------------------ !
  
  ! Fraction
  
  REAL(8) :: Frac_XRFHETE2
  REAL(8) :: sigma_XRFHETE2
  REAL(8) :: NGRB_XRFHETE2=0.d0, NGRB_HETE2=0.d0
  CHARACTER(Len=150) :: XRFHETE2File      = "XRFHETE2_constraint"                             ! Outputfile name for generated XRF fraction 
  CHARACTER(Len=150) :: XRFHETE2ErrorFile = "XRFHETE2_constrainterr"                          ! Outputfile name for error on generated XRF fraction
  
 ! Observational data

  REAL(8), PARAMETER :: Frac_XRFHETE2obs  = 0.35d0                                           ! X-ray flash fraction detected by HETE2
  REAL(8), PARAMETER :: sigma_XRFHETE2obs = 0.15d0                                           ! Error on X-ray flash fraction detected by HETE2


  ! --- GBM : Ep distribution --- !
  ! ----------------------------- !
  
  INTEGER, PARAMETER            :: N_EpGBM                   = 9                      ! Number of bins in TabEpGBM
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobs        = 0.d0                   ! Table for Histogram of Peak energy
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobs_master = 0.d0                   ! Table for Master Histogram of Peak energy
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobserr     = 0.d0                   ! Table for Histogram of error on peak energy
  REAL(8)                       :: k_EpGBM                   = 0.d0                   ! Normalization coefficient
  INTEGER                       :: iBinEpGBM                                          ! Index for histogram
  REAL(8)                       :: Delta_t_EpGBM             = 8.28d0                 ! Duration of the observation of the EpGBM sample (in years)

  CHARACTER(Len=150)             :: EpGBMFile       = "EpGBM_constraint"               ! Outputfile name for EpGBM
  CHARACTER(Len=150)             :: EpGBMErrorFile  = "EpGBM_constrainterr"            ! Outputfile name for EpGBM Error

  ! Observational data

  REAL(8), DIMENSION(0:N_EpGBM) :: TabEpGBM_Epobs, TabEpGBM_LogEpobs                  ! Extracted histogram from GBM catalog
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobsobs    = 0.d0                    ! Table for observed GBM Histogram of Epobs (data from Gruber et al 2014 catalog)
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobsobserr = 0.d0                    ! Table for observed error of GBM Histogram of Epobs (data from Gruber et al 2014 catalog)

  
  ! --- eBAT6 : Ep-L plane --- !
  ! -------------------------- !
  
  INTEGER, PARAMETER                              :: N_eBAT6_EpL             = 97                         ! Number of bins in TabeBAT6_EpL
  REAL(8), DIMENSION(0:1,0:N_eBAT6_EpL)           :: TabeBAT6_EpL            = 0.d0                       ! Table for bins of EpL distribution
  REAL(8), DIMENSION(1:N_eBAT6_EpL,1:N_eBAT6_EpL) :: TabHisteBAT6_EpL        = 0.d0                       ! Table for Histogram of EpL distribution
  REAL(8), DIMENSION(1:N_eBAT6_EpL,1:N_eBAT6_EpL) :: TabHisteBAT6_EpL_master = 0.d0                       ! Table for Master Histogram of EpL distribution
  INTEGER                                         :: iBineBAT6_Ep, iBineBAT6_L                            ! Indexes for histogram

  INTEGER, PARAMETER                    :: Indice_Ep = 0
  INTEGER, PARAMETER                    :: Indice_L  = 1

  CHARACTER(Len=150)             :: eBAT6_EpLFile       = "eBAT6_EpL"               ! Outputfile name for eBAT6 Ep L plane
  CHARACTER(Len=150)             :: eBAT6_EpLErrorFile  = "eBAT6_EpLerr"            ! Outputfile name for eBAT6 Error on Ep L plane


  ! --- eBAT6 : redshift distribution --- !
  ! ------------------------------------- !
  
  INTEGER, PARAMETER            :: N_eBAT6               = 15                         ! Number of bins in TabeBAT6
  REAL(8), DIMENSION(1:N_eBAT6) :: TabHisteBAT6_z        = 0.d0                       ! Table for Histogram of redshift distribution
  REAL(8), DIMENSION(1:N_eBAT6) :: TabHisteBAT6_z_master = 0.d0                       ! Table for Master Histogram of redshift distribution
  REAL(8), DIMENSION(1:N_eBAT6) :: TabHisteBAT6_zerr     = 0.d0                       ! Table for Histogram of error on redshift distribution
  REAL(8)                       :: norm_eBAT6            = 0.d0                       ! Normalization coefficient
  INTEGER                       :: iBineBAT6                                          ! Index for histogram
  REAL(8)                       :: k_eBAT6               = 0.d0                       ! Normalization coefficient
  REAL(8)                       :: Delta_t_eBAT6         = 9.3d0                      ! Duration of the observation of the eBAT6 sample (in years, date first burst - date last burst)
  
  CHARACTER(Len=150)             :: eBAT6File       = "eBAT6_constraint"               ! Outputfile name for eBAT6
  CHARACTER(Len=150)             :: eBAT6ErrorFile  = "eBAT6_constrainterr"            ! Outputfile name for eBAT6 Error

  ! Observational data

  REAL(8), DIMENSION(0:N_eBAT6) :: TabeBAT6_z                                         ! Extracted histogram bins from extended BAT6 catalog
  REAL(8), DIMENSION(1:N_eBAT6) :: TabHisteBAT6_zobs    = 0.d0                        ! Table for observed eBAT6 Histogram of redshift (data from Pescalli et al. 2016)
  REAL(8), DIMENSION(1:N_eBAT6) :: TabHisteBAT6_zobserr = 0.d0                        ! Table for observed error of eBAT6 Histogram of redshift (data from Pescalli et al. 2016)



  ! -------- Chi2 -------- !
  ! ---------------------- !

  REAL(8), DIMENSION(0:N_Constraints) :: Chi2                                                ! Table for Chi2 adjusment of the model for each constraint
  REAL(8), DIMENSION(0:N_Constraints) :: Chi2_old                                            ! Table for Chi2 adjusment of the model for each constraint (used in MCMC)
  INTEGER :: dof                                                                             ! Numbers of degrees of freedom
  REAL(8) :: delta_chi2_1                                                                    ! Chi2 interval in which "good models" are included at 1 sigma (68.3%)
  REAL(8) :: delta_chi2_2                                                                    ! Chi2 interval in which "good models" are included at 2 sigma (95.5%)
  REAL(8) :: delta_chi2_3                                                                    ! Chi2 interval in which "good models" are included at 3 sigma (99.7%)
  REAL(8) :: Chi2_min = 1.d21                                                                ! Lowest Chi2 found
  REAL(8) :: Chi2_best = 1.d21                                                               ! Best Chi2 found (used in post processing to redefine acceptable interval of Chi2)
!  REAL(8) :: Chi2_old = 1.d21                                                                ! Old Chi2 used for MCMC
!  REAL(8) :: Chi2_new = 1.d21                                                                ! New Chi2 used for MCMC  
  INTEGER :: good_model=0, Nb_good_models=0, Nb_models=0                                     ! Number of models and good models (defined above)
  INTEGER :: Nb_good_models_to_pp=0, Nb_good_models_pp=0                                     ! Number of good models to post-process and number that have been post-processed so far


  ! ------- Likelihood ------- !
  ! -------------------------- !

  REAL(8), DIMENSION(0:N_Constraints) :: lnL     = 0.d0                           ! Log of the likelihood function
  REAL(8), DIMENSION(0:N_Constraints) :: lnL_max_test = 0.d0                      ! Log of the likelihood function
  REAL(8), DIMENSION(0:N_Constraints) :: lnL_empty_test = 0.d0                    ! Log of the likelihood function
  REAL(8), DIMENSION(0:N_Constraints) :: lnL_old = 0.d0                           ! Log of the likelihood function (used in MCMC)
  REAL(8)                             :: epsilon = 1.d-3                          ! Small addition to avoid empty bins in log likelihood
  REAL(8)                             :: lnL_max = -999.d0                        ! Max value found for lnL
  REAL(8), DIMENSION(0:N_Constraints) :: lnL_weight = 1.d0                        ! Weight given to each constraint in lnL calculation


  
  REAL(8) :: t
  INTEGER :: ii, jj, kk                                                                          ! Indices for loops
  REAL(8) :: pseudo_collapse_rate   = 0.d0
  REAL(8) :: collapse_rate_from_SFR = 5.66786d8                               ! From SFR Vangioni+15 (based on Springel Hernquist)
    


  ! ----- Variables for the GRB draws ----- !

  ! Intrinsic
  REAL(8) :: L,z,D_L,Ep,alpha,beta              ! variables used for each GRB draw
  REAL(8) :: ktild                              ! Normalization of spectrum
  REAL(8) :: t_star                             ! Luminosity BPL_evol draw
  
  ! Assumptions
  REAL(8) :: z0,zmax,a,b,zm,R0
  REAL(8) :: Ep0, sigmaLog
  
  ! Observed
  REAL(8) :: softness
  REAL(8) :: Epobs
  
  
   

  ! ------ MCMC stuff ----- !

  INTEGER :: N_MCMC_iter = 25000                      ! Number of iterations to run the MCMC chain for
  INTEGER :: MCMC_run_nb                              ! Current run number
  INTEGER :: N_walkers = 3                            ! Number of walkers to try
  REAL(8) :: a_ratio = 1.d0                           ! Ratio used in the Metropolis-Hastings algorithm 
  REAL(8) :: accepted = 0.d0                          ! Number of accepted jumps
  REAL(8) :: acceptance_ratio = 0.d0                  ! Ratio of accepted jumps to total jumps
  INTEGER :: steps_since_reset = 0                    ! Number of steps since reset
  INTEGER :: rejected = 0                             ! Counts how many times a step was rejected to avoid getting stuck somewhere
  INTEGER :: rejected_limit = 1000                    ! Maximum number of steps rejected after which the chain is reset
  INTEGER :: accepted_rec = 0                         ! Acceptance recorded (saved). 0 if step was rejected, 1 if step was accepted (used for saving and post processing purposes)
  INTEGER :: MCMC_step_check = 250                    ! Number of steps after which the codes checks if the step size is appropriate
  INTEGER :: temperature_check = 100                  ! Number of steps after which the codes updates the temperature for the annealing algorithm.
  REAL(8) :: tau       = 0.d0                         ! Temperature used in the annealing algorithm
  REAL(8) :: tau_0     = 1.d3                         ! Starting temperature
  REAL(8) :: tau_lim   = 1.d0                         ! Temperature limit towards which tau converges 
  REAL(8) :: kappa_tau = 0.d0                         ! Cool down speed
  
  ! --------------------------------- !
  ! ------- End of declaration ------ !
  ! --------------------------------- !


  ! ----- Initialization ----- !

  CALL MPI_INIT(code)                                     ! MPI environment starts here
  !CALL PTIM_start(user_label)                             ! times stuff
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)          ! gives the rank of the processor
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)      ! gives the number of processors being used
  WRITE(str_rank,'(I1)') rank                             ! gives a string version of the processor rank
  
  IF(rank == master_proc) SHOW_READ = .TRUE.         ! Only master proc has verbosity (otherwise unreadable)
  CALL InitPhysics(.FALSE.)                          ! Initializes Omega_M, Omega_Lambda etc... and filters
  CALL InitCosmology(SHOW_READ)                      ! Creates precise tables for z, dVdz, D_L and Vz
  CALL Generate_Names()                              ! Creates TabSample_name
  CALL ReadInitFile(InitFile)                        ! Reads the input file and sets the models according to the user's choice
 
  IF(run_mode > 0) THEN
     CALL ReadParamSearchFile()                         ! Reads the input file for parameter search mode
  ELSE
     CALL Calculate_dof()
  END IF
  CALL Generate_Paths()                              ! Creates the paths for the save files
  CALL Init_RNG()                                    ! Initializes random number generator  
  CALL InitBoundaries()                              ! Creates various boundaries for the samples 
  IF(run_mode < 2) CALL Reset_eBAT6_output_files()
  IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " Initializion OK"
  ! CALL TestKiss()

  ! --- Observational constraints --- !

  CALL Prepare_Constraints()                        ! Generate the tables with observational data
  IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " Prepare_Constraints OK"
  delta_chi2_1 = Calc_DeltaChi2(0.6827d0, REAL(dof,8)) 
  delta_chi2_2 = Calc_DeltaChi2(0.9545d0, REAL(dof,8)) 
  delta_chi2_3 = Calc_DeltaChi2(0.9973d0, REAL(dof,8)) 

  IF ((reprise .EQV. .TRUE.) .OR. (run_mode == post_process)) THEN
     ! Do nothing 
  ELSE
     IF(rank == master_proc) THEN
        CALL WRITE_INFO()
     END IF
  END IF
  
  ! Reprise mode 
  CALL Reprise_mode()
  IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " Reprise_mode OK"

 
  IF(run_mode == param_search) THEN
     
     DO ii=1, 100
        good_model = 0
       
        ! Save random generator indexes for reproductibility
        CALL Save_KISS()
        IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " Save_KISS OK"

        
        ! Random draw of parameters
        IF(rank == master_proc) THEN
           !CALL Draw_model_parameters()
           !TabParam_Lum(Param_Lum_Lbreak) = 10.d0**(50.0d0 + ii/100.d0 * 3.d0)
           !TabParam_Ep(Param_Ep_alpha_amati) = 0.0d0 + ii/100.d0 * 1.d0
           TabParam_z(Param_z_zeta) = 0.25d0 + ii/100.d0 * .25d0
           
        END IF
        IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " Draw_model_parameters OK"
        CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)

        CALL MonteCarlo(hist_flag)
        IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " MonteCarlo  OK"

      
        !IF(Chi2(0) < Chi2_min) THEN
!!$        IF(lnL(0) > lnL_max) THEN
!!$           ! Rerun model with correct seeds and save histograms
!!$           iKiss = TabSave_ijk_KISS(1,rank)
!!$           jKiss = TabSave_ijk_KISS(2,rank)
!!$           kKiss = TabSave_ijk_KISS(3,rank)
!!$           IF(rank == master_proc) THEN
!!$              CALL Draw_model_parameters()
!!$           END IF
!!$           CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
!!$           CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
!!$           CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
!!$           CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
!!$           hist_flag = 2
!!$           CALL Reset_eBAT6_output_files()
!!$           CALL MonteCarlo(hist_flag)
!!$           hist_flag = 0
!!$           IF(rank == master_proc) THEN
!!$              !WRITE(*,*) "New best model :", Chi2(0)," (compared to :",Chi2_min,"). Histograms saved."
!!$              WRITE(*,*) "New best model :", lnL(0)," (compared to :",lnL_max,"). Histograms saved."
!!$           END IF
!!$           !Chi2_min = Chi2(0)
!!$           lnL_max = lnL(0)
!!$           IF(rank == master_proc) CALL Update_best_parameters()
!!$     END IF
        
!!$        IF(Chi2(0) <= dof + delta_chi2_3) THEN ! 3 sigma interval for now 
!!$           good_model = 1
!!$           Nb_good_models = Nb_good_models + 1
!!$        END IF
!!$       
        Nb_models = Nb_models + 1
 
!!$        ! Write in reprise
!!$        OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted', POSITION='append')
!!$        IF(rank == master_proc) WRITE(43) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
!!$             & Chi2, lnL, k_Kommers, k_Stern, k_Preece, k_EpGBM, Frac_XRFHETE2, dof, lnL_max, Nb_good_models, Nb_models
!!$        CLOSE(43) ! Close reprise file
        OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted', POSITION='append')                     
        IF(rank == master_proc) WRITE(43) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &           
             & Chi2, lnL, k_Kommers, k_Stern, k_Preece, k_EpGBM, k_eBAT6, dof, accepted_rec, tau, Step_Lum, Step_z, Step_Ep,&
             & TabHistLogL_master, TabHistz_master, TabHistLogEp_master, TabHistLogEpobs_master, TabHistLogP_master, &
             & TabHistKomm_DRDP, TabHistPreece_Ep_master, TabHistStern_P23_master, TabHistEpGBM_Epobs_master, TabHisteBAT6_z_master
        CLOSE(43) ! Close reprise file

        IF(Nb_models >= 25000) THEN
           IF(rank == master_proc)   WRITE(*,*) "Reached 25000 models"
           EXIT
        END IF
        
        IF(rank == master_proc)THEN  
           IF(MOD(Nb_models, 5) == 0) WRITE(*,'(I6x,I5,A,F5.1)') Nb_models, Nb_good_models, ' Best lnL : ',lnL_max
        END IF
       
     END DO

    
    

  ! Markov Chain Monte Carlo
  ELSE IF(run_mode == MCMC) THEN
     IF (rank == master_proc) WRITE(*,*) " You chose MCMC mode"

     DO kk = 1, N_walkers
        Tabparam_Lum_old = TabParam_Lum
        TabParam_Ep_old  = TabParam_Ep
        TabParam_z_old   = TabParam_z

        ! Save random generator indexes for reproductibility
        CALL Save_KISS()

        !CALL Reset_Markov_Chain()
        CALL Reset_MCMC_Step_Size()
        CALL Reset_temperature()
        CALL MonteCarlo(hist_flag)                 ! Calculate Chi2
        rejected = 0                               ! Reset consecutive rejected jumps
        accepted_rec = 1                           ! By design this jump is accepted
        accepted = 0.d0                            ! Reset the accepted rate
        acceptance_ratio = 0.d0                    ! Reset ratio of accepted jumps
        steps_since_reset = 0                      ! Reset steps count
        MCMC_run_nb = MCMC_run_nb + 1
        
        kappa_tau = Calc_kappa_tau(tau_0, 1.d-3) ! Use to calculate the appropriate cooling rate for the annealing algorithm
        starting_i = Calc_starting_i()
        DO ii = starting_i, N_MCMC_iter
           ! Save random generator indexes for reproductibility
           CALL Save_KISS()


           IF (MOD(ii, temperature_check) == 0) CALL Update_temperature(kappa_tau, tau_lim)
           
           IF (tau < 1.001d0) THEN ! Don't check acceptance until chain has cooled
              rejected = 0
              accepted = 0.d0
              acceptance_ratio = 0.d0
           ELSE
              CALL Reset_MCMC_Step_Size() ! TEMP
              Step_Lum = Step_Lum * (1.d0 + 5.d0 * tau/tau_0)
              Step_Ep  = Step_Ep *(1.d0 + 5.d0 * tau/tau_0)
              Step_z   = Step_z *(1.d0 + 5.d0 * tau/tau_0)
           END IF
           
           !IF ( (MOD(ii, MCMC_step_check) == 0) .AND. (tau < 1.001d0) )CALL Update_MCMC_Step_Size(0.5d0, 2.d0)
          
         
           IF (rejected >= rejected_limit) THEN ! If after a certain amount of jumps the chain hasn't moved, manually reset it (avoids getting stuck in areas of low likelihood)
              CALL Reset_Markov_Chain()                  ! Randomly assign a new position in parameter space
           END IF

           steps_since_reset = steps_since_reset + 1
           
           ! Save current state 
           TabParam_Lum_old = TabParam_Lum
           TabParam_z_old   = TabParam_z
           TabParam_Ep_old  = TabParam_Ep
           Chi2_old = Chi2
           lnL_old  = lnL
           
           ! Draw Markov jump in parameter space
           IF(rank == master_proc) THEN
              CALL Draw_Markov_Jump()
           END IF
           CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           
           CALL MonteCarlo(hist_flag)

           MCMC_run_nb = MCMC_run_nb + 1
           IF((rank == master_proc) .AND. (verbose >= 1)) WRITE(*,*) " run : ",MCMC_run_nb," MonteCarlo OK"
          

           IF (lnL(0) >= lnL_old(0)) THEN
              accepted = accepted + 1.d0
              rejected = 0              
           ELSE
              a_ratio = (lnL(0)- lnL_old(0)) / tau
              IF(a_ratio <= -50.d0) THEN               ! because exp(-50) ~ 1e-22, essentially zero
                 a_ratio = 0.d0
              ELSE
                 a_ratio = EXP(a_ratio)        
              END IF
              IF(rank == master_proc) t = uniform()
              CALL MPI_BCAST(t, 1, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
              IF(t <= a_ratio ) THEN
                 accepted = accepted + 1.d0
                 rejected = 0
              ELSE ! if jump not accepted
                 rejected = rejected + 1
                 CALL Update_Markov_Chain()    ! If the Markov jump is not accepted, revert to previous state                                                                                        
              END IF
              
           END IF ! end of better or worse likelihood
           acceptance_ratio = accepted / REAL(steps_since_reset,8)
           
           CALL Write_reprise_MCMC()
           
        END DO ! End MCMC iterations
     END DO ! end walkers loop 
     
  ELSE IF (run_mode == post_process) THEN
     ! This mode is for recomputing the good models and generating error bars on the observables

     CALL Define_Nb_lines()
     Chi2_best = Chi2_min
     IF(rank == master_proc) THEN
        WRITE(*,*) "Chi2_min found = ", Chi2_best
        OPEN(UNIT=46, FILE=TRIM(path)//'reprise.dat', FORM='unformatted')
     END IF
     Nb_good_models_pp = 0
     pp_flag = 0
     DO WHILE(pp_flag == 0)
        IF(rank == master_proc) THEN
           DO
              READ(46, IOSTAT=pp_flag) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
                   & Chi2, k_Kommers, k_Stern, k_Preece, k_EpGBM, Frac_XRFHETE2, dof, Chi2_min, Nb_good_models, Nb_models
              jj = jj + 1
              IF(ABS(Chi2(0)) <= Chi2_best + delta_chi2_3) THEN
                 Nb_good_models_pp = Nb_good_models_pp + 1
                 WRITE(*,'(A,I5,A,F4.0,A)') "Nb_good_models_pp : ", Nb_good_models_pp, " (",100.d0*REAL(jj,8)/REAL(Nb_lines,8),"% of file)"
                 !WRITE(*,*) "Good Chi2 = ", Chi2(0), Chi2(Constraint_EpGBM), Chi2(Constraint_Stern)
                 EXIT
              END IF
           END DO
        END IF

        CALL MPI_BCAST(pp_flag, 1, MPI_INT, master_proc, MPI_COMM_WORLD, code)        
        
        !WRITE(*,*) "Before BCAST ; rank : ", rank," ijk Kiss : ", TabSave_ijk_KISS
        CALL MPI_BCAST(TabSave_ijk_KISS, 3*nb_procs, MPI_INT, master_proc, MPI_COMM_WORLD, code)
        !WRITE(*,*) "After BCAST ; rank : ", rank," ijk Kiss : ", TabSave_ijk_KISS
        
        iKiss = TabSave_ijk_KISS(1,rank)
        jKiss = TabSave_ijk_KISS(2,rank)
        kKiss = TabSave_ijk_KISS(3,rank)
  
        !WRITE(*,*) "Before draw ; rank : ", rank," Lmin, Lmax, slope : ", TabParam_Lum(Param_Lum_Lmin),  TabParam_Lum(Param_Lum_Lmax),  TabParam_Lum(Param_Lum_slope)
        !WRITE(*,*) "Before draw ; rank : ", rank," Tabparamz : ",TabParam_z   
        IF(rank == master_proc) THEN
           CALL Draw_model_parameters()
        END IF      
        
        CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        !WRITE(*,*) "After draw ; rank : ", rank," Tabparamz : ",TabParam_z   
        !WRITE(*,*) "After  draw ; rank : ", rank," Lmin, Lmax, slope : ", TabParam_Lum(Param_Lum_Lmin),  TabParam_Lum(Param_Lum_Lmax),  TabParam_Lum(Param_Lum_slope)
        
        CALL MonteCarlo(hist_flag)

        !WRITE(*,*) "PP Chi2 = ", Chi2(0), Chi2(Constraint_EpGBM), Chi2(Constraint_Stern)
        OPEN(UNIT=43, FILE=TRIM(path)//'reprise_pp.dat', FORM='unformatted', POSITION='append')
        IF(rank == master_proc) WRITE(43) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
             & Chi2, k_Kommers, k_Stern, k_Preece, k_EpGBM, Frac_XRFHETE2, dof, chi2_min, Nb_good_models, Nb_models
        CLOSE(43) ! Close post-process reprise file
     END DO

     
     IF(rank == master_proc) THEN
        CLOSE(46) ! Close reprise file
     END IF

  ELSE ! IF(run_mode == one_run)

  IF (rank == master_proc) WRITE(*,*) " You chose one run mode"   

  IF (rank == master_proc)  CALL InitECLAIRs(nsigmasECLAIRs)
  
  CALL MPI_BCAST(bkgECLAIRsB1,  1,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  Threshold(Sample_SVOM) = bkgECLAIRsB1
  CALL MPI_BCAST(omegaECLAIRs,  1,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABeffECLAIRsE,  NeffECLAIRs,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABeffECLAIRsA,  NeffECLAIRs,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABbkgECLAIRsE,  NbkgECLAIRs,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABbkgECLAIRsB,  NbkgECLAIRs,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABomegaECLAIRs,  2*NoffECLAIRs+1,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_BCAST(TABoffECLAIRs  ,  2*NoffECLAIRs+1,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
 

  CALL Save_KISS()

  CALL MonteCarlo(hist_flag)
  
END IF

IF(rank == master_proc) WRITE(*,*) "[------------------- END OF CODE ------------------]"
CALL MPI_FINALIZE(code)

CONTAINS  
  



  SUBROUTINE MonteCarlo(hist_flag)
    INTEGER, INTENT(in) :: hist_flag
    ! Main subroutine of the code, is organized as follows : 
    !
    !    1. Select models (Lum, z, Ep...)
    !
    !    (1.bis) Prepare histograms (x-axis, with the limits given from the models, reset histograms)
    !
    !    2. Generate GRBs (random draw) 
    !        a. Properties
    !           - Intrinsic (Lum, z, Ep...)
    !           - Observed  (Peak flux for various instruments, Peak energy observed...)
    !        b. Detection by instrument (BATSE23, SWIFT...)
    !        c. Preparation of indices (Lum, z, Ep, Peak flux...)
    !        d. Filling of histograms per sample (Intrinsic, BATSE23, SWIFT...)
    !
    !    3. Normalizations (GRB rate, k_stern, etc...)
    !
    !    4. Chi2 Calculation (using Kommers, Preece, Stern...)
    !
    !    (5.) Save histograms (Lum, z, Ep...)
    !
    !    (6.) Save constraints (Kommers, Preece, Stern...)
    !

  
    INTEGER, PARAMETER                              :: N_saved_Prop = 6                       ! Saving Log(L), z, Ep, alpha, beta
    INTEGER                                         :: i_Prop                                 ! Index for loops
    CHARACTER(*), PARAMETER                         :: GRB_PropFile = "GRB_Properties"        ! Filename for saving the properties
    REAL(8), DIMENSION(0:N_Samples, 0:N_saved_Prop) :: GRB_Prop                               ! Table that saves the properties
    REAL(8), DIMENSION(0:N_Samples, 0:N_saved_Prop) :: GRB_Prop_master                        ! Table that saves the properties for master proc
    REAL(8), DIMENSION(0:N_Samples, 0:N_saved_Prop) :: average_Prop                           ! Table for the average of the properties
    REAL(8), DIMENSION(0:N_Samples, 0:N_saved_Prop) :: average_Prop_master                    ! Table for the average of the properties for master proc

    REAL(8) :: t, x1, x2
    INTEGER :: imin,imax,i,j,M                                                                ! variables for loops and dichotomy
   
    IF(rank == master_proc) THEN
       IF(run_mode == one_run) WRITE(*,*) " "
       IF(run_mode == one_run) WRITE(*,'(A)') " ====================================== Monte Carlo ================================================= "
       IF(run_mode == one_run) WRITE(*,*) " "
       IF(run_mode == one_run) WRITE(*, '(A)', advance='no') "[     (%) =  0"
    END IF
    
    ! --------------------------------- !
    ! -------- 1. Select model -------- !
    ! --------------------------------- !
    
    CALL Set_Fixed_Parameters()                ! Calculate certain quantities once for the whole population (t_star, ktild, L, z if fixed...)
    CALL Prepare_Redshift()                    ! Prepare distribution function for redshift draws
    CALL Prepare_Luminosity()                  ! Prepare distribution function for Luminosity draws  
    CALL Prepare_Spec()                        ! Prepare distribution function for alpha and beta draws (currently only from GBM catalog)

    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo Select model  OK"
    
    ! ---------------------------------- !
    ! --- (1.bis) Prepare histograms --- !
    ! ---------------------------------- !
    
    IF(hist_flag >= 1) THEN
       CALL Prepare_Histogram_Prop(TabHistlim_inf, TabHistlim_sup, n_call)        ! Create histogram limits and resets them
    END IF
    
    IF(hist_flag >= 1) THEN
       CALL Reset_Histogram_Samples()                                             ! Resets histograms not used in Chi2 calculations but still included
    END IF
  
    CALL Reset_Histogram_Chi2()                                                   ! Resets histograms used in Chi2 calculation
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo Prepare histograms  OK"

    CALL Reset_GRB_seeds()                                                        ! Reset GRB seeds to insure Chi2 variations aren't due to different realizations of the intrinsic population
    
    ! -------------------------------------- !
    ! ---------- 2. Generate GRBs ---------- !
    ! -------------------------------------- !
    
    IF(run_mode == one_run) THEN
       IF(Save_all_GRB .EQV. .TRUE.) THEN
          IF(rank == master_proc) OPEN(UNIT=66, FILE=TRIM(path)//GRB_PropFile//".dat", FORM='unformatted')                ! Save all GRB properties, be careful the file is heavy (several GB)
       END IF
       GRB_Prop = 0.d0                                    ! Reset GRB properties
       average_Prop = 0.d0                                ! Reset GRB properties
    END IF
    NaNtest_hist = .FALSE.                             ! Reset all NaN tests to false
    
    ! ----------------------------------------------------------------------------------- MAIN LOOP ----------------------------------------------------------------------------------------- !
    
    
    DO i=1, Nb_GRB/nb_procs
      
       NaNtest_prop = .FALSE. 
       
       ! ---------- a. Properties ---------- !
       ! ----------------------------------- !

       ! ----------- Intrinsic properties ----------- !
       CALL Draw_Redshift(z)
!       z = 1.0d0
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo z  OK : ", z
       CALL Draw_Luminosity(L)
!       L = 1.0d52
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo Lum  OK : ", L
       CALL Draw_alpha_beta(alpha, beta)
!       alpha = 1.0d0
!       beta = 2.5d0
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo alpha, beta  OK : ", alpha, beta
       CALL Draw_Ep(Ep)
!       Ep = 600.d0
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo Ep  OK  : ", Ep
       
       ! ------------ Observed properties ----------- !      
       Epobs = Calc_Epobs(Ep,z)
      
       CALL Calculate_Peakflux_Instrument()
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo Calculate_Peakflux_Instrument  OK"
       CALL Calculate_Detection_Probability()
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo Calculate_Detection_Probability  OK"
       
       ! --------------- Histograms --------------- !
       ! ------------------------------------------ !
       
       IF(hist_flag >= 1) THEN
          CALL Fill_Histogram_Prop()                       ! Find the bins and add the probability of detection to the histograms of L, z, Ep etc...
       END IF

       CALL Fill_Histogram_Chi2()                          ! Find the bins and add the probability of detection to the histograms to use in Chi2 calculation

       IF(hist_flag >= 1) CALL Fill_Samples()                ! Find the bins and add the probability of detection to the histograms included but not used in Chi2 calculation
       
       IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo filling of Histograms  OK"
       
       IF(run_mode == one_run) THEN                        ! Calculate the GRB population properties 
       
          ! --- GRB properties --- ! 
          DO i_Sample=0, N_Samples
             GRB_Prop(i_Sample,          0) = GRB_Prop(i_Sample,          0) + Prob_det(i_Sample)
             GRB_Prop(i_Sample,     Prop_L) = GRB_Prop(i_Sample,     Prop_L) + L     * Prob_det(i_Sample)
             GRB_Prop(i_Sample,     Prop_z) = GRB_Prop(i_Sample,     Prop_z) + z     * Prob_det(i_Sample)
             GRB_Prop(i_Sample,    Prop_Ep) = GRB_Prop(i_Sample,    Prop_Ep) + Epobs * Prob_det(i_Sample)
             GRB_Prop(i_Sample, Prop_alpha) = GRB_Prop(i_Sample, Prop_alpha) + alpha * Prob_det(i_Sample)
             GRB_Prop(i_Sample,  Prop_beta) = GRB_Prop(i_Sample,  Prop_beta) + beta  * Prob_det(i_Sample)
             
          END DO
        
          ! Only record Intrinsic sample
          !      1                  2                3                   4                       5                        6
          ! [GRB number]   [Luminosity erg/s]    [Redshift]   [Observed Peak energy]   [spectral slope alpha]   [spectral slope beta]
          ! i4, 7f8, 5f8, 11f8
          
          IF(Save_all_GRB .EQV. .TRUE.)  WRITE(66) i, L, z, D_L, Epobs, alpha, beta, ktild, Peakflux_Instrument, Prob_det
          
         
          IF (Mode_test) THEN
             IF(rank == master_proc) THEN
                WRITE(*,*) " "
                WRITE(*,*) "GRB_nb =", i
                WRITE(*,*) "z =", z
                WRITE(*,*) "D_L =", D_L*Mpc
                WRITE(*,'(A,1ES12.5)') " L =", L
                WRITE(*,*) "alpha =", alpha
                WRITE(*,*) "beta =",beta
                WRITE(*,*) "Ep =",Ep
                WRITE(*,*) "IBand =", 1.d0/ktild
                WRITE(*,*) "ktild =", ktild
                DO i_Instrument=1, N_Instruments
                   WRITE(*,*) "Instrument ",TRIM(TabInstrument_name(i_Instrument)), " : Peakflux =",Peakflux_Instrument(i_Instrument)
                END DO
                WRITE(*,*) "softness =",softness
                DO i_Sample=1, N_Samples
                   WRITE(*,*) "Prob_det "//TRIM(TabSample_name(i_Sample))//" =", Prob_det(i_Sample)
                END DO
                WRITE(*,*) "Epobs =", Epobs
                READ(*,*)
             END IF
             IF(verbose==1 .OR. verbose==2) READ(*,*)
          END IF
         
 
          ! --- Nan testing --- !
          IF (NaNtest_prop) THEN
             WRITE(*,*) "[       z         D_L        Lum        alpha         beta         ktild        Ep         Epobs    ]"
             WRITE(*,*) "[  ", z , D_L, L, alpha, beta, ktild, Ep, Epobs, "   ]"
             WRITE(*,*) "[                                      Peak Flux Instrument                                         ]"
             WRITE(*,*) "[  ", Peakflux_Instrument,"  ]"
             WRITE(*,*) "[                                        Peak Flux Sample                                           ]"
             WRITE(*,*) "[  ", Peakflux,"  ]"
             WRITE(*,*) "[                                       Prob det Instrument                                         ]"
             WRITE(*,*) "[  ", Prob_det,"  ]"  
          END IF
                 
                
          ! Track progress
          IF(rank == master_proc) THEN
             IF( MOD(100.d0 * REAL(i,8)*REAL(nb_procs,8)/REAL(Nb_GRB,8), 5.d0) == 0 ) THEN
                WRITE(*,'(I3,1x)', advance='no') INT(100.d0 * REAL(i,8)*REAL(nb_procs,8)/REAL(Nb_GRB))
             END IF
          END IF
       END IF ! End of one_run mode
       
    END DO ! end of main loop

    IF(run_mode == one_run) THEN
       IF((Save_all_GRB .EQV. .TRUE.) .AND. (rank == master_proc)) CLOSE(66)     
       IF(rank == master_proc) WRITE(*, '(A)') "      ]"
    END IF
   ! --------------------------------- !
   ! ----- End of : Generate GRB ----- !
   ! --------------------------------- !
   
   ! ----------------------------------------------------------------------------------- END OF MAIN LOOP ----------------------------------------------------------------------------------- !
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Main Loop ended OK"
    ! WRITE(*,*) " in Monte Carlo Main Loop ended OK for proc : ", rank
    
   ! -------- Combine histograms of each procs -------- !
    IF(Sample_Included(Sample_Kommers)) CALL MPI_REDUCE(TabHistKomm_P23,    TabHistKomm_P23_master,    N_Komm,   MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Kommers reduce OK"
    IF(Sample_Included(Sample_Stern))   CALL MPI_REDUCE(TabHistStern_P23,   TabHistStern_P23_master,   N_Stern,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Stern reduce OK"
    IF(Sample_Included(Sample_Preece))  CALL MPI_REDUCE(TabHistPreece_Ep,   TabHistPreece_Ep_master,   N_Preece, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Preece reduce OK"
    IF(Sample_Included(Sample_EpGBM))   CALL MPI_REDUCE(TabHistEpGBM_Epobs, TabHistEpGBM_Epobs_master, N_EpGBM,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo EpGBM reduce OK"
    IF(Sample_Included(Sample_eBAT6))   CALL MPI_REDUCE(TabHisteBAT6_z,     TabHisteBAT6_z_master,     N_eBAT6,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo eBAT6 reduce OK"
    
    IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Reduced constraints OK"
    IF(hist_flag >= 1) THEN
       DO i_Sample = 0, N_Samples
          IF(Sample_Included(i_Sample)) THEN
             CALL MPI_REDUCE(    TabHistLogL(i_Sample,:),     TabHistLogL_master(i_Sample,:), N_L,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             CALL MPI_REDUCE(       TabHistz(i_Sample,:),        TabHistz_master(i_Sample,:), N_z,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             CALL MPI_REDUCE(   TabHistLogEp(i_Sample,:),    TabHistLogEp_master(i_Sample,:), N_Ep, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             CALL MPI_REDUCE(TabHistLogEpobs(i_Sample,:), TabHistLogEpobs_master(i_Sample,:), N_Ep, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             CALL MPI_REDUCE(   TabHistalpha(i_Sample,:),    TabHistalpha_master(i_Sample,:), N_spec_a,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             CALL MPI_REDUCE(    TabHistbeta(i_Sample,:),     TabHistbeta_master(i_Sample,:), N_spec_b,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             IF (i_Sample >= 1) THEN
                CALL MPI_REDUCE( TabHistLogP(i_Sample,:),     TabHistLogP_master(i_Sample,:), N_P,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
             END IF
          END IF
       END DO
       ! eBAT6 Ep-L plane
       CALL MPI_REDUCE(TabHisteBAT6_EpL, TabHisteBAT6_EpL_master, N_eBAT6_EpL**2,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
       IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in Monte Carlo Reduced prop histograms OK"
   END IF
   !WRITE(*,*) 'rank : ',rank,' tot, ind : ', GRB_Prop_master(0, 0), GRB_Prop(0, 0)
   IF(run_mode == one_run) THEN
      DO i_Sample = 0, N_Samples
         DO i_Prop=0, N_saved_Prop
            CALL MPI_REDUCE(    GRB_Prop(i_Sample,i_Prop),     GRB_Prop_master(i_Sample,i_Prop), 1, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD,code)
            CALL MPI_REDUCE(average_Prop(i_Sample,i_Prop), average_Prop_master(i_Sample,i_Prop), 1, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD,code)
         END DO
      END DO
   END IF
  
   
   ! Master does calculations
   IF(rank == master_proc) THEN
      
   ! --- Print average quantities --- !
      IF(run_mode == one_run) THEN
         WRITE(*,*) " "
         WRITE(*,'(A)') " ======================================= Average quantities ========================================="
         WRITE(*,'(A)') "[  Sample              Nb_GRB           Luminosity   Redshift       Ep        alpha       beta      ]"
         WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
         
         DO i_Sample=0, N_Samples
            IF(Sample_Included(i_Sample)) THEN
               average_Prop_master(i_Sample, 0) = GRB_Prop_master(i_Sample, 0)
               IF ( average_Prop_master(i_Sample, 0) .NE. 0.d0 ) THEN
                  !WRITE(*,*) 'rank : ',rank,' L, tot_av : ', GRB_Prop_master(i_Sample,     Prop_L), average_Prop_master(i_Sample, 0)
                  average_Prop_master(i_Sample,     Prop_L) = GRB_Prop_master(i_Sample,     Prop_L) / average_Prop_master(i_Sample, 0) 
                  average_Prop_master(i_Sample,     Prop_z) = GRB_Prop_master(i_Sample,     Prop_z) / average_Prop_master(i_Sample, 0) 
                  average_Prop_master(i_Sample,    Prop_Ep) = GRB_Prop_master(i_Sample,    Prop_Ep) / average_Prop_master(i_Sample, 0) 
                  average_Prop_master(i_Sample, Prop_alpha) = GRB_Prop_master(i_Sample, Prop_alpha) / average_Prop_master(i_Sample, 0) 
                  average_Prop_master(i_Sample,  Prop_beta) = GRB_Prop_master(i_Sample,  Prop_beta) / average_Prop_master(i_Sample, 0)
               ELSE 
                  average_Prop_master(i_Sample, :) = 0.d0
               END IF
               IF(run_mode == one_run)  WRITE(*,'(A,A15,A,ES12.5,F6.1,A,5ES12.5,A)') "[  ",TabSample_name(i_Sample),":", &
                    & average_Prop_master(i_Sample, 0), 100.d0 * average_Prop_master(i_Sample, 0)/average_Prop_master(Sample_Intrinsic, 0), "%", &
                    & average_Prop_master(i_Sample,     Prop_L), &
                    & average_Prop_master(i_Sample,     Prop_z), &
                    & average_Prop_master(i_Sample,    Prop_Ep), &
                    & average_Prop_master(i_Sample, Prop_alpha), &
                    & average_Prop_master(i_Sample,  Prop_beta), &
                    & "  ]"
            END IF
         END DO
         WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      END IF ! end one_run mode printing

      ! ---------- 3. Normalizations ---------- !
      ! --------------------------------------- !

      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                                                                                   ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                    Normalization coefficients                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      CALL Normalize_Model()
      IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo Normalize_Model  OK"

      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                                                                                   ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                          Log Likelihood                                           ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      CALL Calculate_unnormalized_Likelihood()
      IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo Calculate_unnormed_lnL  OK"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                            lnL = ", lnL(0), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                        max lnL = ", lnL_max_test(0), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                      empty lnL = ", lnL_empty_test(0), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      ! ----- 4. Chi squared calculation ----- !
      ! -------------------------------------- !
      
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                                                                                   ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                          Chi squared                                              ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      CALL Calculate_Chi2()
      IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo Calculate_Chi2  OK"
      !IF(rank == master_proc)  WRITE(*,*) " new Chi2 = ", Chi2(0), Chi2(Constraint_EpGBM), Chi2(Constraint_Stern)
      
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                           Chi2 = ", Chi2(0), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                           Delta Chi2 (3 sigma) = ", dof + delta_chi2_3, "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,I4.0,A)')   "[                             Degrees of freedom =   ", dof, "                                           ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                   Reduced Chi2 = ", Chi2(0)/REAL(dof,8), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      
      ! ----- (5.) Save histograms ----- !
      ! -------------------------------- !
      
      IF(hist_flag == 2) CALL Save_Histograms()
      IF(hist_flag == 2) CALL Save_eBAT6()
      
      ! ----- (6.) Save constraints ----- !
      ! --------------------------------- !
      
      IF(hist_flag == 2)  CALL Save_Constraints()
      
      IF(run_mode == post_process) CALL Post_process_Constraints()
      
      IF(run_mode == one_run)  WRITE(*,*) " "
      IF(run_mode == one_run)  WRITE(*,'(A)') " ==================================== Monte Carlo done ============================================== "
      IF(run_mode == one_run)  WRITE(*,*) " "
      
      
   END IF ! Ends master rank computation
   
   ! Broadcast chi2
   CALL MPI_BCAST(Chi2, N_Constraints, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   ! Broadcast lnL
   CALL MPI_BCAST(lnL, N_Constraints, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   

 END SUBROUTINE MonteCarlo
 

 
 SUBROUTINE ReadInitFile(Name)
   INTEGER :: s_L,s_z,s_Spec,s_Ep,s_P,s_Constraint,i,incl,temp
   REAL(8) :: temp_real
   ! Initializes the parameters for the code
   CHARACTER(Len=*), INTENT(in) :: Name

   IF(SHOW_READ) WRITE(*,'(A)')     " "
   IF(SHOW_READ) WRITE(*,'(A)')     " ========================= ReadInitFile ========================= "
   IF(SHOW_READ) WRITE(*,'(A)')     " "
   IF(SHOW_READ) WRITE(*,'(A64A)')   "[           Reading parameters in : "//TRIM(Name)//"        "," ]"
   IF(SHOW_READ) WRITE(*,'(A)')     "[ -------------------------------------------------------------- ]"
   OPEN(UNIT=50, FILE=TRIM(Name))

   IF(SHOW_READ) WRITE(*,'(A)')     "[                                                                ]"
  
   IF(SHOW_READ) WRITE(*,'(A)')     "[                     Input model parameters                     ]"
   IF(SHOW_READ) WRITE(*,'(A)')     "[ -------------------------------------------------------------- ]"   
  
   ! Path

   CALL ReadLine(50, path)
   path = '../Model_outputs/' // path
   IF(SHOW_READ) WRITE(*,'(A,A47,A)')           "[  Output path : ",ADJUSTL(path)," ]"
   
   ! Reprise mode

   CALL ReadInteger(50, temp)
   IF (temp == 1) THEN
      reprise = .TRUE. 
      IF(SHOW_READ) WRITE(*,'(A)') "[                       reprise mode : TRUE                      ]"
   ELSE
      IF(SHOW_READ) WRITE(*,'(A)') "[                       reprise mode : FALSE                     ]"
   END IF
   
   ! RNG
   CALL ReadInteger(50, RNG)
   IF (RNG == MT19937) THEN 
      IF(SHOW_READ) WRITE(*,'(A)') "[                         RNG chosen : MT                        ]"
   ELSE IF (RNG == Kiss_rng) THEN
      IF(SHOW_READ) WRITE(*,'(A)') "[                         RNG chosen : KISS                      ]"
   ELSE
      IF(SHOW_READ) WRITE(*,'(A)') "[                         RNG chosen : INVALID                   ]"
   END IF
   
   ! Run mode
   CALL ReadInteger(50, run_mode)
   IF (run_mode == one_run) THEN
      hist_flag = 2
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : one run                   ]"
   ELSE IF (run_mode == param_search) THEN
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : parameter search          ]"
   ELSE IF (run_mode == post_process) THEN
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : post processing           ]"
   ELSE IF (run_mode == MCMC) THEN
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : MCMC                      ]"
     ! hist_flag = 2
   ELSE
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : INVALID                   ]"
      STOP "INVALID RUN MODE"
   END IF


   ! Number of GRBs 
   
   CALL ReadReal(50, temp_real) 
   Nb_GRB = INT(temp_real)
   IF(SHOW_READ) WRITE(*,'(A,1ES12.5,A)') "[           Number of GRBs simulated : ", REAL(Nb_GRB,8), "              ]"
   IF(SHOW_READ) WRITE(*,'(A)')           "[                                                                ]"
  
   ! Luminosity Function 
   
   CALL ReadInteger(50, Model_Lum)
   
   SELECT CASE(Model_Lum)
   CASE(Model_LumFix) ! Fixed Luminosity
      CALL ReadReal(50, TabParam_Lum(Param_Lum_L0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Fix ->", "L0 =", TabParam_Lum(Param_Lum_L0),"  ]"
      
   CASE(Model_LumPL) ! Power Law distribution
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmin))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmin =",TabParam_Lum(Param_Lum_Lmin),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmax =",TabParam_Lum(Param_Lum_Lmax),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slope))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "slope =",TabParam_Lum(Param_Lum_slope),"  ]"
      
   CASE(Model_LumBPL_evol) ! Broken Power Law distribution
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmin))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmin =",TabParam_Lum(Param_Lum_Lmin),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmax =",TabParam_Lum(Param_Lum_Lmax),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lbreak))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lbreak =",TabParam_Lum(Param_Lum_Lbreak),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slopeL))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeL =",TabParam_Lum(Param_Lum_slopeL),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slopeH))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeH =",TabParam_Lum(Param_Lum_slopeH),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_k_evol))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "k_evol =",TabParam_Lum(Param_Lum_k_evol),"  ]"

   CASE(Model_LumPL_evol) ! Power Law distribution that evolves
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmin))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "Lmin =",TabParam_Lum(Param_Lum_Lmin),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "Lmax =",TabParam_Lum(Param_Lum_Lmax),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slope))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "slope =",TabParam_Lum(Param_Lum_slope),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_k_evol))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "k_evol =",TabParam_Lum(Param_Lum_k_evol),"  ]"
   CASE(Model_LumSch) ! Schechter distribution
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmin))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lmin =",TabParam_Lum(Param_Lum_Lmin),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lbreak))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lbreak =",TabParam_Lum(Param_Lum_Lbreak),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slope))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "slope =",TabParam_Lum(Param_Lum_slope),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_k_evol))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "k_evol =",TabParam_Lum(Param_Lum_k_evol),"  ]"
      
   CASE DEFAULT
      STOP "Error : Luminosity Function unknown (check .init file)" 
   END SELECT
   
   ! Redshift Distribution 
   
   CALL ReadInteger(50, Model_z)
   SELECT CASE(Model_z)
   CASE(Model_zFix) ! Fixed z
      CALL ReadReal(50, TabParam_z(Param_z_z0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Fix ->", "z0 =", TabParam_z(Param_z_z0),"  ]"
   CASE(Model_zUniform) ! Uniform z up to zmax
      CALL ReadReal(50, TabParam_z(Param_z_zmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Uniform ->", "zmax =", TabParam_z(Param_z_zmax),"  ]"
   CASE(Model_zSH) ! Springel-Hernquist distribution
      CALL ReadReal(50, TabParam_z(Param_z_zmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zmax =", TabParam_z(Param_z_zmax),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_zm))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zm =", TabParam_z(Param_z_zm),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "b =", TabParam_z(Param_z_b),"  ]"
   CASE(Model_zDaigne)
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Daigne 2006 ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Daigne 2006 ->", "b =", TabParam_z(Param_z_b),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_c))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Daigne 2006 ->", "c =", TabParam_z(Param_z_c),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_d))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Daigne 2006 ->", "d =", TabParam_z(Param_z_d),"  ]"
   CASE(Model_z_evol) ! Springel-Hernquist distribution
      CALL ReadReal(50, TabParam_z(Param_z_zmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "zmax =", TabParam_z(Param_z_zmax),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_zm))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "zm =", TabParam_z(Param_z_zm),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "b =", TabParam_z(Param_z_b),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_zeta))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "zeta =", TabParam_z(Param_z_zeta),"  ]"
   CASE(Model_zLi)
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Li 2008 ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Li 2008 ->", "b =", TabParam_z(Param_z_b),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_c))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Li 2008 ->", "c =", TabParam_z(Param_z_c),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_d))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Li 2008 ->", "d =", TabParam_z(Param_z_d),"  ]"
   CASE(Model_zPesc)
      STOP "Pescalli redshift distribution not implemented yet"
   CASE(Model_zBPL)
      CALL ReadReal(50, TabParam_z(Param_z_zmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Power Law ->", "zmax =", TabParam_z(Param_z_zmax),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_zm))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Power Law ->", "zm =", TabParam_z(Param_z_zm),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Power Law ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Power Law ->", "b =", TabParam_z(Param_z_b),"  ]"
   CASE(Model_zBExp)
      CALL ReadReal(50, TabParam_z(Param_z_zmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "zmax =", TabParam_z(Param_z_zmax),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_zm))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "zm =", TabParam_z(Param_z_zm),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_a))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "a =", TabParam_z(Param_z_a),"  ]"
      CALL ReadReal(50, TabParam_z(Param_z_b))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "b =", TabParam_z(Param_z_b),"  ]"
   CASE DEFAULT
      STOP "Error : Redshift Distribution unknown (check .init file)"
   END SELECT
   
   
   ! Spectral Model (SP)
   
   CALL ReadInteger(50, Model_Spec)
   SELECT CASE(Model_Spec)
   CASE(Model_SpecBPLFix) ! Fixed Broken Power Law
      CALL ReadReal(50, TabParam_Spec(Param_spec_alpha))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "BPL Fix ->", "alpha =", TabParam_Spec(Param_spec_alpha),"  ]"
      CALL ReadReal(50, TabParam_Spec(Param_spec_beta))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "BPL Fix ->", "beta =", TabParam_Spec(Param_spec_beta),"  ]"
   CASE (Model_SpecBPLK) ! Kaneko Broken Power Law
      !IF(SHOW_READ) WRITE(*,Format_RIF) "Model_Spec :", "BPL Kaneko ->", "alpha =", TabParam_Spec(Param_spec_alpha)
      STOP "Model_Spec : Kaneko BPL Function not coded yet"
   CASE(Model_SpecBandFix) ! Fixed Band distribution
      CALL ReadReal(50, TabParam_Spec(Param_spec_alpha))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band Fix ->", "alpha =", TabParam_Spec(Param_spec_alpha),"  ]"
      CALL ReadReal(50, TabParam_Spec(Param_spec_beta))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band Fix ->", "beta =", TabParam_Spec(Param_spec_beta),"  ]"
   CASE(Model_SpecBandK) ! Kaneko Band distribution
      STOP "Model_Spec : Kaneko Band Function not coded yet"
   CASE(Model_SpecBandD) ! Daigne Band distribution
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band Daigne ->", "alpha =", 1.,"  ]"
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band Daigne ->", "sigma_a =", 0.5,"  ]"
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band Daigne ->", "beta >=", 2.,"  ]"
   CASE(Model_SpecBandGBM) ! Band distribution from GBM catalog
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band GBM ->", "alpha ~", .6,"  ]"
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Spec :", "Band GBM ->", "beta ~", 4.2,"  ]"
   CASE DEFAULT
      STOP "Error : Spectral Model unknown (check .init file)"
   END SELECT
   
   
   ! Ep Model
   
   CALL ReadInteger(50, Model_Ep)
   SELECT CASE(Model_Ep)
   CASE(Model_EpFix) ! Fixed Ep
      CALL ReadReal(50, TabParam_Ep(Param_Ep_Ep0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Fix ->", "Ep0 =", TabParam_Ep(Param_Ep_Ep0),"  ]"
   CASE(Model_EpLogNormal) ! Log Normal Ep distribution
      CALL ReadReal(50, TabParam_Ep(Param_Ep_Ep0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Log Normal ->", "Ep0 =", TabParam_Ep(Param_Ep_Ep0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_sigmaLog))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Log Normal ->", "sigmaLog =", TabParam_Ep(Param_Ep_sigmaLog),"  ]"
   CASE(Model_EpY) ! Yonetoku 
      STOP "Model_Ep : Yonetoku not coded yet"     
   CASE(Model_EpAmati)
      CALL ReadReal(50, TabParam_Ep(Param_Ep_Ep0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "Ep0 =", TabParam_Ep(Param_Ep_Ep0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_sigmaLog))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "sigmaLog =", TabParam_Ep(Param_Ep_sigmaLog),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_L0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "L0 =", TabParam_Ep(Param_Ep_L0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_alpha_amati))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "alpha_amati =", TabParam_Ep(Param_Ep_alpha_amati),"  ]"
   CASE DEFAULT
      STOP "Error : Ep Model unknown (check .init file)"
   END SELECT
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   

   ! Include Sample
   IF(SHOW_READ) WRITE(*,'(A)') "[                                                                ]"
   IF(SHOW_READ) WRITE(*,'(A)') "[                           Samples                              ]"
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"
   Sample_Included(0) = .TRUE.  ! Always include intrinsic sample

   DO i_Sample = 1, N_Samples
      CALL ReadInteger(50, incl)
      IF(incl == 1) THEN
         Sample_Included(i_Sample) = .TRUE.
         IF(SHOW_READ) WRITE(*,'(A,A20,A)') '[                 Including sample : ', TabSample_name(i_Sample),"        ]"                 
         SELECT CASE(i_Sample)
         CASE(Sample_Kommers)
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Sample_Preece)
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Sample_Stern)
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Sample_SWIFTweak)
            Instrument_Included(Instrument_BAT  ) = .TRUE.
         CASE(Sample_SWIFT)
            Instrument_Included(Instrument_BAT  ) = .TRUE.
         CASE(Sample_SWIFTbright)
            Instrument_Included(Instrument_BAT  ) = .TRUE.
         CASE(Sample_HETE2)
            Instrument_Included(Instrument_FREGATE) = .TRUE.
            Instrument_Included(Instrument_WXM) = .TRUE.
         CASE(Sample_eBAT6)
            Instrument_Included(Instrument_BAT  ) = .TRUE.
         CASE(Sample_EpGBM)
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Sample_SVOM)
            Instrument_Included(Instrument_ECLAIRs) = .TRUE.
         END SELECT
      END IF
   END DO

   ! Include Constraint
  
   DO i_Constraint=1, N_Constraints
      CALL ReadInteger(50, incl)
      IF(incl == 1) THEN
         Constraint_Included(i_Constraint) = .TRUE.
         IF(SHOW_READ) WRITE(*,'(A,A20,A)') '[             Including constraint : ', TabConstraint_name(i_Constraint),"        ]"
         SELECT CASE(i_Constraint)
         CASE(Constraint_Kommers)
            Sample_Included(Sample_Kommers) = .TRUE.
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Constraint_Preece)
            Sample_Included(Sample_Preece) = .TRUE.
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Constraint_Stern)
            Sample_Included(Sample_Stern) = .TRUE.
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Constraint_HETE2)
            Sample_Included(Sample_HETE2) = .TRUE.
            Instrument_Included(Instrument_FREGATE) = .TRUE.
            Instrument_Included(Instrument_WXM) = .TRUE.
         CASE(Constraint_EpGBM)
            Sample_Included(Sample_EpGBM) = .TRUE.
            Instrument_Included(Instrument_BATSE) = .TRUE.
         CASE(Constraint_eBAT6)
            Sample_Included(Sample_eBAT6) = .TRUE.
            Instrument_Included(Instrument_BAT) = .TRUE.
         END SELECT 
      END IF
   END DO
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   

   ! Save histograms

   CALL ReadInteger(50, s_L)
   IF(s_L == 1) THEN 
      Lsave = .TRUE.
      IF(SHOW_READ) WRITE(*,'(A,A)') '[            Saving Luminosity as :','   luminosity_[sample].dat    ]'
   END IF
   CALL ReadInteger(50, s_z)
   IF(s_z == 1) THEN 
      zsave = .TRUE.
      IF(SHOW_READ) WRITE(*,'(A,A)') '[              Saving Redshift as :','     redshift_[sample].dat    ]'
   END IF
   CALL ReadInteger(50, s_Spec)
   IF(s_Spec == 1) THEN 
      Specsave = .TRUE.
      IF(SHOW_READ) WRITE(*,'(A,A)') '[            Saving Model Spec as :','    modelspec_[sample].dat    ]'
   END IF
   CALL ReadInteger(50, s_Ep)
   IF(s_Ep == 1) THEN 
      Epsave = .TRUE.
      IF(SHOW_READ) WRITE(*,'(A,A)') '[           Saving Peak Energy as :','   peakenergy_[sample].dat    ]'
   END IF
   CALL ReadInteger(50, s_P)
   IF(s_P == 1) THEN 
      Psave = .TRUE.
      IF(SHOW_READ) WRITE(*,'(A,A)') '[             Saving Peak Flux as :','     peakflux_[sample].dat    ]'
   END IF

   ! Save constraints
   DO i = 1, N_Constraints
      CALL ReadInteger(50, s_Constraint) 
      IF(s_Constraint == 1) THEN
         Constraint_save(i) = .TRUE.
         IF(SHOW_READ) WRITE(*,'(A,A34,A31)') '[','Saving '// TRIM(TabConstraint_name(i)) // &
              ' constraint as :',' ' // TRIM(TabConstraint_name(i)) // '_constraint[err].dat  ]'
      END IF
   END DO
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') " ====================== ReadInitFile done ======================== "
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') " "
   CLOSE(50)

 END SUBROUTINE ReadInitFile


 SUBROUTINE ReadParamSearchFile()
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') " ==================== Parameter Search File ====================== "
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   


   ! Luminosity
   OPEN(UNIT=51, FILE=TRIM(param_search_path)//'lum_param_search.init')
   CALL ReadInteger(51, lum_explore)

   IF(lum_explore == 1) THEN
      IF(SHOW_READ) WRITE(*,'(A)') " "
      IF(SHOW_READ) WRITE(*,'(A)') "[                           Luminosity                           ]"   
      IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
      SELECT CASE(Model_Lum)
      CASE(Model_LumFix) ! Fixed Luminosity
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_L0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Fix ->", "L0_min =", TabParam_Lum_min(Param_Lum_L0),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_L0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Fix ->", "L0_max =", TabParam_Lum_max(Param_Lum_L0),"  ]"
         
      CASE(Model_LumPL) ! Power Law distribution
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmin_min =",TabParam_Lum_min(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmin_max =",TabParam_Lum_max(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmax_min =",TabParam_Lum_min(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "Lmax_max =",TabParam_Lum_max(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "slope_min =",TabParam_Lum_min(Param_Lum_slope),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Power Law ->", "slope_max =",TabParam_Lum_max(Param_Lum_slope),"  ]"

          
      CASE(Model_LumBPL_evol) ! Evolving BPL distribution
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmin_min =",TabParam_Lum_min(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmin_max =",TabParam_Lum_max(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmax_min =",TabParam_Lum_min(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lmax_max =",TabParam_Lum_max(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lbreak_min =",TabParam_Lum_min(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "Lbreak_max =",TabParam_Lum_max(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slopeL))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeL_min =",TabParam_Lum_min(Param_Lum_slopeL),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slopeL))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeL_max =",TabParam_Lum_max(Param_Lum_slopeL),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slopeH))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeH_min =",TabParam_Lum_min(Param_Lum_slopeH),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slopeH))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "slopeH_max =",TabParam_Lum_max(Param_Lum_slopeH),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "k_evol_min =",TabParam_Lum_min(Param_Lum_k_evol),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving BPL ->", "k_evol_max =",TabParam_Lum_max(Param_Lum_k_evol),"  ]"
         
      CASE(Model_LumPL_evol) ! Evolving power Law distribution
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "Lmin_min =",TabParam_Lum_min(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "Lmin_max =",TabParam_Lum_max(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "Lmax_min =",TabParam_Lum_min(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "Lmax_max =",TabParam_Lum_max(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "slope_min =",TabParam_Lum_min(Param_Lum_slope),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "slope_max =",TabParam_Lum_max(Param_Lum_slope),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power law ->", "k_evol_min =",TabParam_Lum_min(Param_Lum_k_evol),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Evolving power Law ->", "k_evol_max =",TabParam_Lum_max(Param_Lum_k_evol),"  ]"

      CASE(Model_LumSch) ! Schechter function
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lmin_min =",TabParam_Lum_min(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lmin_max =",TabParam_Lum_max(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lbreak_min =",TabParam_Lum_min(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "Lbreak_max =",TabParam_Lum_max(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "slope_min =",TabParam_Lum_min(Param_Lum_slope),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slope))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "slope_max =",TabParam_Lum_max(Param_Lum_slope),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "k_evol_min =",TabParam_Lum_min(Param_Lum_k_evol),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_k_evol))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Schechter ->", "k_evol_max =",TabParam_Lum_max(Param_Lum_k_evol),"  ]"

   END SELECT
   END IF
         
   CLOSE(51)

   ! Redshift
   OPEN(UNIT=52, FILE=TRIM(param_search_path)//'redshift_param_search.init')
   CALL ReadInteger(52, z_explore)

   IF(z_explore == 1) THEN
      IF(SHOW_READ) WRITE(*,'(A)') " "
      IF(SHOW_READ) WRITE(*,'(A)') "[                            Redshift                            ]"   
      IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
      SELECT CASE(Model_z)
      CASE(Model_zFix) ! Fixed redshift
         CALL ReadReal(52, TabParam_z_min(Param_z_z0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Fix ->", "z0_min =", TabParam_z_min(Param_z_z0),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_z0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Fix ->", "z0_max =", TabParam_z_max(Param_z_z0),"  ]"
      CASE(Model_zSH)
         CALL ReadReal(52, TabParam_z_min(Param_z_zmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zmax_min =", TabParam_z_min(Param_z_zmax),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_zmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zmax_max =", TabParam_z_max(Param_z_zmax),"  ]"
         CALL ReadReal(52, TabParam_z_min(Param_z_zm))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zm_min =", TabParam_z_min(Param_z_zm),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_zm))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "zm_max =", TabParam_z_max(Param_z_zm),"  ]"
         CALL ReadReal(52, TabParam_z_min(Param_z_a))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "a_min =", TabParam_z_min(Param_z_a),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_a))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "a_max =", TabParam_z_max(Param_z_a),"  ]"
         CALL ReadReal(52, TabParam_z_min(Param_z_b))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "b_min =", TabParam_z_min(Param_z_b),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_b))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Springel-Hernquist ->", "b_max =", TabParam_z_max(Param_z_b),"  ]"
      CASE(Model_z_evol)
         CALL ReadReal(52, TabParam_z_min(Param_z_zeta))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "zeta_min =", TabParam_z_min(Param_z_zeta),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_zeta))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "z evolution (SH) ->", "zeta_max =", TabParam_z_max(Param_z_zeta),"  ]"
      CASE(Model_zBExp)
         CALL ReadReal(52, TabParam_z_min(Param_z_zm))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "zm_min =", TabParam_z_min(Param_z_zm),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_zm))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "zm_max =", TabParam_z_max(Param_z_zm),"  ]"
         CALL ReadReal(52, TabParam_z_min(Param_z_a))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "a_min =", TabParam_z_min(Param_z_a),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_a))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "a_max =", TabParam_z_max(Param_z_a),"  ]"
         CALL ReadReal(52, TabParam_z_min(Param_z_b))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "b_min =", TabParam_z_min(Param_z_b),"  ]"
         CALL ReadReal(52, TabParam_z_max(Param_z_b))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_z :", "Broken Exponential ->", "b_max =", TabParam_z_max(Param_z_b),"  ]"
      CASE DEFAULT
         IF(SHOW_READ) WRITE(*,'(A)') "[           No parameter search for this redshift model          ]"
         IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
      END SELECT
   END IF
   CLOSE(52)
   
   ! Peak energy
   OPEN(UNIT=53, FILE=TRIM(param_search_path)//'Ep_param_search.init')
   CALL ReadInteger(53, Ep_explore)

   IF(Ep_explore == 1) THEN
      IF(SHOW_READ) WRITE(*,'(A)') " "
      IF(SHOW_READ) WRITE(*,'(A)') "[                          Peak Energy                           ]"   
      IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
      SELECT CASE(Model_Ep)
      CASE(Model_EpFix) ! Fixed Ep
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Fix ->", "Ep0_min =", TabParam_Ep_min(Param_Ep_Ep0),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Fix ->", "Ep0_max =", TabParam_Ep_max(Param_Ep_Ep0),"  ]"
      CASE(Model_EpLogNormal) ! Lognormal Ep
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "LogNormal ->", "Ep0_min =", TabParam_Ep_min(Param_Ep_Ep0),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "LogNormal ->", "Ep0_max =", TabParam_Ep_max(Param_Ep_Ep0),"  ]"
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_sigmaLog))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "LogNormal ->", "sigmaLog_min =", TabParam_Ep_min(Param_Ep_sigmaLog),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_sigmaLog))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "LogNormal ->", "sigmaLog_max =", TabParam_Ep_max(Param_Ep_sigmaLog),"  ]"
      CASE(Model_EpAmati) 
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "Ep0_min =", TabParam_Ep_min(Param_Ep_Ep0),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_Ep0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "Ep0_max =", TabParam_Ep_max(Param_Ep_Ep0),"  ]"
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_sigmaLog))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "sigmaLog_min =", TabParam_Ep_min(Param_Ep_sigmaLog),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_sigmaLog))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "sigmaLog_max =", TabParam_Ep_max(Param_Ep_sigmaLog),"  ]"
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_L0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " L0_min =", TabParam_Ep_min(Param_Ep_L0),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_L0))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " L0_max =", TabParam_Ep_max(Param_Ep_L0),"  ]"
         CALL ReadReal(53, TabParam_Ep_min(Param_Ep_alpha_amati))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " alpha_min =", TabParam_Ep_min(Param_Ep_alpha_amati),"  ]"
         CALL ReadReal(53, TabParam_Ep_max(Param_Ep_alpha_amati))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " alpha_max =", TabParam_Ep_max(Param_Ep_alpha_amati),"  ]"
      CASE DEFAULT
         STOP "Error : Ep Model unknown (check .init file)"
      END SELECT
   END IF
      CLOSE(53)
        
   CALL Calculate_dof()
    
   IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"
   IF(SHOW_READ) WRITE(*,'(A,I2,A)') "[         (data to fit - degrees of freedom) : ",dof,"                ]"   
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') " ================== Parameter Search File done =================== "
   IF(SHOW_READ) WRITE(*,'(A)') " "
   IF(SHOW_READ) WRITE(*,'(A)') " "
   
 END SUBROUTINE ReadParamSearchFile

 SUBROUTINE Calculate_dof()
   dof = 0
   IF(lum_explore == 1) THEN
      SELECT CASE(Model_Lum)
      CASE(Model_LumFix)
         dof = dof - 1
      CASE(Model_LumPL)
         IF(TabParam_Lum_min(Param_Lum_Lmin)  /= TabParam_Lum_max(Param_Lum_Lmin))  dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lmax)  /= TabParam_Lum_max(Param_Lum_Lmax))  dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slope) /= TabParam_Lum_max(Param_Lum_slope)) dof = dof - 1
      CASE(Model_LumBPL_evol)
         IF(TabParam_Lum_min(Param_Lum_Lmin)   /= TabParam_Lum_max(Param_Lum_Lmin))   dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lmax)   /= TabParam_Lum_max(Param_Lum_Lmax))   dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lbreak) /= TabParam_Lum_max(Param_Lum_Lbreak)) dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slopeL) /= TabParam_Lum_max(Param_Lum_slopeL)) dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slopeH) /= TabParam_Lum_max(Param_Lum_slopeH)) dof = dof - 1
      CASE(Model_LumPL_evol)
         IF(TabParam_Lum_min(Param_Lum_Lmin)    /= TabParam_Lum_max(Param_Lum_Lmin))    dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lmax)    /= TabParam_Lum_max(Param_Lum_Lmax))    dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slope)   /= TabParam_Lum_max(Param_Lum_slope))   dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_k_evol) /= TabParam_Lum_max(Param_Lum_k_evol)) dof = dof - 1
      END SELECT
   END IF
   
   IF(z_explore == 1) THEN
      SELECT CASE(Model_z)
      CASE(Model_zFix)
         dof = dof - 1
      CASE(Model_zSH)
         IF(TabParam_z_min(Param_z_zmax) /= TabParam_z_max(Param_z_zmax)) dof = dof - 1
         IF(TabParam_z_min(Param_z_zm)   /= TabParam_z_max(Param_z_zm))   dof = dof - 1
         IF(TabParam_z_min(Param_z_a)    /= TabParam_z_max(Param_z_a))    dof = dof - 1
         IF(TabParam_z_min(Param_z_b)    /= TabParam_z_max(Param_z_b))    dof = dof - 1
      END SELECT
   END IF
   
   IF(Ep_explore == 1) THEN
      SELECT CASE(Model_Ep)
      CASE(Model_EpFix)
         dof = dof - 1
      CASE(Model_EpLogNormal)
         IF(TabParam_Ep_min(Param_Ep_Ep0)       /= TabParam_Ep_max(Param_Ep_Ep0))       dof = dof - 1
         IF(TabParam_Ep_min(Param_Ep_sigmaLog)  /= TabParam_Ep_max(Param_Ep_sigmaLog))  dof = dof - 1
      CASE(Model_EpAmati)
         IF(TabParam_Ep_min(Param_Ep_Ep0)         /= TabParam_Ep_max(Param_Ep_Ep0))         dof = dof - 1
         IF(TabParam_Ep_min(Param_Ep_sigmaLog)    /= TabParam_Ep_max(Param_Ep_sigmaLog))    dof = dof - 1
         IF(TabParam_Ep_min(Param_Ep_L0)          /= TabParam_Ep_max(Param_Ep_L0))          dof = dof - 1
         IF(TabParam_Ep_min(Param_Ep_alpha_amati) /= TabParam_Ep_max(Param_Ep_alpha_amati)) dof = dof - 1
      END SELECT
   END IF
   
   DO i_Constraint = 1, N_Constraints
      IF(Constraint_included(i_Constraint)) THEN
         IF(i_Constraint == Constraint_Kommers)  dof = dof + N_Komm   - 1  ! -1 for the normalization
         IF(i_Constraint == Constraint_Preece )  dof = dof + N_Preece - 1 
         IF(i_Constraint == Constraint_Stern  )  dof = dof + N_Stern  - 1 
         IF(i_Constraint == Constraint_HETE2) dof = dof + 1
         IF(i_Constraint == Constraint_EpGBM)    dof = dof + N_EpGBM  - 1
         !IF(i_Constraint == Constraint_eBAT6)    dof = dof + N_eBAT6  - 1    
      END IF
   END DO

 END SUBROUTINE Calculate_dof
 
 SUBROUTINE Set_Fixed_Parameters()
   INTEGER              :: imin, imax, j
   ! --- Redshift --- !
   
   IF(Model_z == Model_zFix) THEN
      z = TabParam_z(Param_z_z0)  
      imin = 1
      imax = INT(SIZE(Tabprecisez)-1)
      IF(z <= 0.d0) WRITE(*,*) "WARNING in MonteCarlo : z <= 0 not allowed, results will be absurd"
      DO
         j = (imin+imax)/2
         IF (z >= Tabprecisez(j)) THEN
            imin = j+1
         ELSE IF (z < Tabprecisez(j-1)) THEN
            imax = j
         ELSE
            EXIT
         END IF
      END DO
      
      D_L = TabD_L(j-1) + (TabD_L(j)-TabD_L(j-1)) * (z-Tabprecisez(j-1)) / (Tabprecisez(j)-Tabprecisez(j-1)) 
   END IF
   
   
   ! -- Luminosity [erg/s] -- !
   
   IF(Model_Lum == Model_LumFix) L = TabParam_Lum(Param_Lum_L0)
   IF(Model_Lum == Model_LumBPL_evol) t_star = Calc_t_star(TabParam_Lum)
   
   
   ! --- Spectral Model --- !
   
   IF(Model_Spec == Model_SpecBPLFix .OR. Model_Spec == Model_SpecBandFix) THEN
      alpha = TabParam_Spec(Param_spec_alpha)
      beta  = TabParam_Spec(Param_spec_beta)
      ktild = Calc_ktild(alpha, beta, Model_Spec)
   END IF
   IF(verbose == 2) PRINT*, "ktild =", ktild
   
   
   ! --- Peak Energy [keV] --- !
   
   IF(Model_Ep == Model_EpFix) Ep = TabParam_Ep(Param_Ep_Ep0)
   
   
 END SUBROUTINE Set_Fixed_Parameters

 SUBROUTINE Draw_model_parameters()
   IF(Lum_Explore == 1) THEN
      SELECT CASE(Model_lum)
      CASE(Model_LumPL)
         IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmin) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmin))-LOG10(TabParam_Lum_min(Param_Lum_Lmin))) + LOG10(TabParam_Lum_min(Param_Lum_Lmin)))
         ELSE
            TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmax) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmax))-LOG10(TabParam_Lum_min(Param_Lum_Lmax))) + LOG10(TabParam_Lum_min(Param_Lum_Lmax)))
         ELSE
            TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN   
            t = uniform()
            TabParam_Lum(Param_Lum_slope) = t*(TabParam_Lum_max(Param_Lum_slope)-TabParam_Lum_min(Param_Lum_slope)) + TabParam_Lum_min(Param_Lum_slope)
         ELSE
            TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
         END IF
         
      CASE(Model_LumBPL_evol)
         IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmin) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmin))-LOG10(TabParam_Lum_min(Param_Lum_Lmin))) + LOG10(TabParam_Lum_min(Param_Lum_Lmin)))
         ELSE
            TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmax) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmax))-LOG10(TabParam_Lum_min(Param_Lum_Lmax))) + LOG10(TabParam_Lum_min(Param_Lum_Lmax)))
         ELSE
            TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lbreak) /= TabParam_Lum_min(Param_Lum_Lbreak)) THEN          
            t = uniform()
            TabParam_Lum(Param_Lum_Lbreak) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lbreak))-LOG10(TabParam_Lum_min(Param_Lum_Lbreak))) + LOG10(TabParam_Lum_min(Param_Lum_Lbreak)))
         ELSE
            TabParam_Lum(Param_Lum_Lbreak) = TabParam_Lum_min(Param_Lum_Lbreak)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slopeL) /= TabParam_Lum_min(Param_Lum_slopeL)) THEN          
            t = uniform()
            TabParam_Lum(Param_Lum_slopeL) = t*(TabParam_Lum_max(Param_Lum_slopeL)-TabParam_Lum_min(Param_Lum_slopeL)) + TabParam_Lum_min(Param_Lum_slopeL)
         ELSE
            TabParam_Lum(Param_Lum_slopeL) = TabParam_Lum_min(Param_Lum_slopeL)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slopeH) /= TabParam_Lum_min(Param_Lum_slopeH)) THEN          
            t = uniform()
            TabParam_Lum(Param_Lum_slopeH) = t*(TabParam_Lum_max(Param_Lum_slopeH)-TabParam_Lum_min(Param_Lum_slopeH)) + TabParam_Lum_min(Param_Lum_slopeH)
         ELSE
            TabParam_Lum(Param_Lum_slopeH) = TabParam_Lum_min(Param_Lum_slopeH)
         END IF
         
      CASE(Model_LumPL_evol)
         IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmin) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmin))-LOG10(TabParam_Lum_min(Param_Lum_Lmin))) + LOG10(TabParam_Lum_min(Param_Lum_Lmin)))
         ELSE
            TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmax) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmax))-LOG10(TabParam_Lum_min(Param_Lum_Lmax))) + LOG10(TabParam_Lum_min(Param_Lum_Lmax)))
         ELSE
            TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN   
            t = uniform()
            TabParam_Lum(Param_Lum_slope) = t*(TabParam_Lum_max(Param_Lum_slope)-TabParam_Lum_min(Param_Lum_slope)) + TabParam_Lum_min(Param_Lum_slope)
         ELSE
            TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
         END IF
         IF(TabParam_Lum_max(Param_Lum_k_evol) /= TabParam_Lum_min(Param_Lum_k_evol)) THEN   
            t = uniform()
            TabParam_Lum(Param_Lum_k_evol) = t*(TabParam_Lum_max(Param_Lum_k_evol)-TabParam_Lum_min(Param_Lum_k_evol)) + TabParam_Lum_min(Param_Lum_k_evol)
         ELSE
            TabParam_Lum(Param_Lum_k_evol) = TabParam_Lum_min(Param_Lum_k_evol)
         END IF
      CASE(Model_LumSch)
         IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lmin) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lmin))-LOG10(TabParam_Lum_min(Param_Lum_Lmin))) + LOG10(TabParam_Lum_min(Param_Lum_Lmin)))
         ELSE
            TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lbreak) /= TabParam_Lum_min(Param_Lum_Lbreak)) THEN
            t = uniform()
            TabParam_Lum(Param_Lum_Lbreak) = 10.d0**(t*(LOG10(TabParam_Lum_max(Param_Lum_Lbreak))-LOG10(TabParam_Lum_min(Param_Lum_Lbreak))) + LOG10(TabParam_Lum_min(Param_Lum_Lbreak)))
         ELSE
            TabParam_Lum(Param_Lum_Lbreak) = TabParam_Lum_min(Param_Lum_Lbreak)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN   
            t = uniform()
            TabParam_Lum(Param_Lum_slope) = t*(TabParam_Lum_max(Param_Lum_slope)-TabParam_Lum_min(Param_Lum_slope)) + TabParam_Lum_min(Param_Lum_slope)
         ELSE
            TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
         END IF
         IF(TabParam_Lum_max(Param_Lum_k_evol) /= TabParam_Lum_min(Param_Lum_k_evol)) THEN   
            t = uniform()
            TabParam_Lum(Param_Lum_k_evol) = t*(TabParam_Lum_max(Param_Lum_k_evol)-TabParam_Lum_min(Param_Lum_k_evol)) + TabParam_Lum_min(Param_Lum_k_evol)
         ELSE
            TabParam_Lum(Param_Lum_k_evol) = TabParam_Lum_min(Param_Lum_k_evol)
         END IF
      END SELECT
   END IF

   IF(Ep_explore == 1) THEN
      SELECT CASE(Model_Ep)
      CASE(Model_EpFix)
         IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_Ep0) = 10.d0**(t*(LOG10(TabParam_Ep_max(Param_Ep_Ep0))-LOG10(TabParam_Ep_min(Param_Ep_Ep0))) + LOG10(TabParam_Ep_min(Param_Ep_Ep0)))
         ELSE
            TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
         END IF
      CASE(Model_EpLogNormal)
         IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_Ep0) = 10.d0**(t*(LOG10(TabParam_Ep_max(Param_Ep_Ep0))-LOG10(TabParam_Ep_min(Param_Ep_Ep0))) + LOG10(TabParam_Ep_min(Param_Ep_Ep0)))
         ELSE
            TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
         END IF
         IF(TabParam_Ep_max(Param_Ep_sigmaLog) /= TabParam_Ep_min(Param_Ep_sigmaLog)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_sigmaLog) = t*(TabParam_Ep_max(Param_Ep_sigmaLog)-TabParam_Ep_min(Param_Ep_sigmaLog)) + TabParam_Ep_min(Param_Ep_sigmaLog)
         ELSE
            TabParam_Ep(Param_Ep_sigmaLog) = TabParam_Ep_min(Param_Ep_sigmaLog)
         END IF
      CASE(Model_EpAmati)
         IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_Ep0) = 10.d0**(t*(LOG10(TabParam_Ep_max(Param_Ep_Ep0))-LOG10(TabParam_Ep_min(Param_Ep_Ep0))) + LOG10(TabParam_Ep_min(Param_Ep_Ep0)))
         ELSE
            TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
         END IF
         IF(TabParam_Ep_max(Param_Ep_sigmaLog) /= TabParam_Ep_min(Param_Ep_sigmaLog)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_sigmaLog) = t*(TabParam_Ep_max(Param_Ep_sigmaLog)-TabParam_Ep_min(Param_Ep_sigmaLog)) + TabParam_Ep_min(Param_Ep_sigmaLog)
         ELSE
            TabParam_Ep(Param_Ep_sigmaLog) = TabParam_Ep_min(Param_Ep_sigmaLog)
         END IF
         IF(TabParam_Ep_max(Param_Ep_L0) /= TabParam_Ep_min(Param_Ep_L0)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_L0) = 10.d0**(t*(LOG10(TabParam_Ep_max(Param_Ep_L0))-LOG10(TabParam_Ep_min(Param_Ep_L0))) + LOG10(TabParam_Ep_min(Param_Ep_L0)))
         ELSE
            TabParam_Ep(Param_Ep_L0) = TabParam_Ep_min(Param_Ep_L0)
         END IF
         IF(TabParam_Ep_max(Param_Ep_alpha_amati) /= TabParam_Ep_min(Param_Ep_alpha_amati)) THEN          
            t = uniform()
            TabParam_Ep(Param_Ep_alpha_amati) = t*(TabParam_Ep_max(Param_Ep_alpha_amati)-TabParam_Ep_min(Param_Ep_alpha_amati)) + TabParam_Ep_min(Param_Ep_alpha_amati)
         ELSE
            TabParam_Ep(Param_Ep_alpha_amati) = TabParam_Ep_min(Param_Ep_alpha_amati)
         END IF
      END SELECT
   END IF
   
   IF(z_explore == 1) THEN
      SELECT CASE(Model_z)
      CASE(Model_z_evol)
         IF(TabParam_z_max(Param_z_zeta) /= TabParam_z_min(Param_z_zeta)) THEN          
            t = uniform()
            TabParam_z(Param_z_zeta) = t*(TabParam_z_max(Param_z_zeta)-TabParam_z_min(Param_z_zeta)) + TabParam_z_min(Param_z_zeta)
         ELSE
            TabParam_z(Param_z_zeta) = TabParam_z_min(Param_z_zeta)
         END IF
      CASE DEFAULT
         STOP "Error : no exploration implemented with this redshift distribution (check .init file)" 
      END SELECT
   END IF
   
 END SUBROUTINE Draw_model_parameters

 SUBROUTINE Draw_Markov_Jump()
   REAL(8) :: temp = 0.d0
   INTEGER(4) :: safeguard = 0
   
   IF(Lum_Explore == 1) THEN
      SELECT CASE(Model_lum)
       
      CASE(Model_LumPL)
         IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
            safeguard = 0
            DO
               temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmin)), Step_Lum(Param_Lum_Lmin)))
               safeguard = safeguard + 1
               IF(safeguard > 10000) STOP 'Looping too many times for Lmin in Draw_Markov_Jump()'
               IF( (temp >= TabParam_Lum_min(Param_Lum_Lmin)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmin))) THEN
                  TabParam_Lum(Param_Lum_Lmin) = temp
                  EXIT
               END IF
            END DO
        ELSE
            TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
         END IF
         IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
            safeguard = 0
            DO
               temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmax)), Step_Lum(Param_Lum_Lmax)))
               safeguard = safeguard + 1
               IF(safeguard > 10000) STOP 'Looping too many times for Lmax in Draw_Markov_Jump()'
               IF( (temp >= TabParam_Lum_min(Param_Lum_Lmax)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmax))) THEN
                  TabParam_Lum(Param_Lum_Lmax) = temp
                  EXIT
               END IF
            END DO
         ELSE
            TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN
            safeguard = 0
            DO
               temp = gaussian(TabParam_Lum(Param_Lum_slope), Step_Lum(Param_Lum_slope))
               safeguard = safeguard + 1
               IF(safeguard > 10000) STOP 'Looping too many times for slope in Draw_Markov_Jump()'
               IF( (temp >= TabParam_Lum_min(Param_Lum_slope)) .AND. (temp <= TabParam_Lum_max(Param_Lum_slope))) THEN
                  TabParam_Lum(Param_Lum_slope) = temp
                  EXIT
               END IF
            END DO
            ELSE
            TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
         END IF
         
         CASE(Model_LumBPL_evol)
            IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmin)), Step_Lum(Param_Lum_Lmin)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lmin in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lmin)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmin))) THEN
                     TabParam_Lum(Param_Lum_Lmin) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
            END IF
            IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmax)), Step_Lum(Param_Lum_Lmax)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lmax in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lmax)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmax))) THEN
                     TabParam_Lum(Param_Lum_Lmax) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
            END IF
            IF(TabParam_Lum_max(Param_Lum_Lbreak) /= TabParam_Lum_min(Param_Lum_Lbreak)) THEN          
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lbreak)), Step_Lum(Param_Lum_Lbreak)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lbreak in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lbreak)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lbreak))) THEN
                     TabParam_Lum(Param_Lum_Lbreak) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lbreak) = TabParam_Lum_min(Param_Lum_Lbreak)
            END IF
            IF(TabParam_Lum_max(Param_Lum_slopeL) /= TabParam_Lum_min(Param_Lum_slopeL)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Lum(Param_Lum_slopeL), Step_Lum(Param_Lum_slopeL))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for slopeL in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_slopeL)) .AND. (temp <= TabParam_Lum_max(Param_Lum_slopeL))) THEN
                     TabParam_Lum(Param_Lum_slopeL) = temp
                     EXIT
                  END IF
            END DO
            ELSE
               TabParam_Lum(Param_Lum_slopeL) = TabParam_Lum_min(Param_Lum_slopeL)
            END IF
            IF(TabParam_Lum_max(Param_Lum_slopeH) /= TabParam_Lum_min(Param_Lum_slopeH)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Lum(Param_Lum_slopeH), Step_Lum(Param_Lum_slopeH))
                   safeguard = safeguard + 1
                   IF(safeguard > 10000) STOP 'Looping too many times for slopeH in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_slopeH)) .AND. (temp <= TabParam_Lum_max(Param_Lum_slopeH))) THEN
                     TabParam_Lum(Param_Lum_slopeH) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_slopeH) = TabParam_Lum_min(Param_Lum_slopeH)
            END IF

         CASE(Model_LumPL_evol)
            IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmin)), Step_Lum(Param_Lum_Lmin)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lmin in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lmin)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmin))) THEN
                     TabParam_Lum(Param_Lum_Lmin) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
            END IF
            IF(TabParam_Lum_max(Param_Lum_Lmax) /= TabParam_Lum_min(Param_Lum_Lmax)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmax)), Step_Lum(Param_Lum_Lmax)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lmax in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lmax)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmax))) THEN
                     TabParam_Lum(Param_Lum_Lmax) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lmax) = TabParam_Lum_min(Param_Lum_Lmax)
            END IF
            IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Lum(Param_Lum_slope), Step_Lum(Param_Lum_slope))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for slope in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_slope)) .AND. (temp <= TabParam_Lum_max(Param_Lum_slope))) THEN
                     TabParam_Lum(Param_Lum_slope) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
            END IF
            IF(TabParam_Lum_max(Param_Lum_k_evol) /= TabParam_Lum_min(Param_Lum_k_evol)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Lum(Param_Lum_k_evol), Step_Lum(Param_Lum_k_evol))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for k_evol in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_k_evol)) .AND. (temp <= TabParam_Lum_max(Param_Lum_k_evol))) THEN
                     TabParam_Lum(Param_Lum_k_evol) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_k_evol) = TabParam_Lum_min(Param_Lum_k_evol)
            END IF
            
         CASE(Model_LumSch)
            IF(TabParam_Lum_max(Param_Lum_Lmin) /= TabParam_Lum_min(Param_Lum_Lmin)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lmin)), Step_Lum(Param_Lum_Lmin)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Lmin in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_Lmin)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lmin))) THEN
                     TabParam_Lum(Param_Lum_Lmin) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_Lmin) = TabParam_Lum_min(Param_Lum_Lmin)
            END IF
         IF(TabParam_Lum_max(Param_Lum_Lbreak) /= TabParam_Lum_min(Param_Lum_Lbreak)) THEN
            safeguard = 0
            DO
               temp = 10.d0**(gaussian(LOG10(TabParam_Lum(Param_Lum_Lbreak)), Step_Lum(Param_Lum_Lbreak)))
               safeguard = safeguard + 1
               IF(safeguard > 10000) STOP 'Looping too many times for Lbreak in Draw_Markov_Jump()'
               IF( (temp >= TabParam_Lum_min(Param_Lum_Lbreak)) .AND. (temp <= TabParam_Lum_max(Param_Lum_Lbreak))) THEN
                  TabParam_Lum(Param_Lum_Lbreak) = temp
                  EXIT
               END IF
            END DO
         ELSE
            TabParam_Lum(Param_Lum_Lbreak) = TabParam_Lum_min(Param_Lum_Lbreak)
         END IF
         IF(TabParam_Lum_max(Param_Lum_slope) /= TabParam_Lum_min(Param_Lum_slope)) THEN
            safeguard = 0
            DO
               temp = gaussian(TabParam_Lum(Param_Lum_slope), Step_Lum(Param_Lum_slope))
               safeguard = safeguard + 1
               IF(safeguard > 10000) STOP 'Looping too many times for slope in Draw_Markov_Jump()'
               IF( (temp >= TabParam_Lum_min(Param_Lum_slope)) .AND. (temp <= TabParam_Lum_max(Param_Lum_slope))) THEN
                  TabParam_Lum(Param_Lum_slope) = temp
                  EXIT
               END IF
            END DO
         ELSE
            TabParam_Lum(Param_Lum_slope) = TabParam_Lum_min(Param_Lum_slope)
         END IF
         IF(TabParam_Lum_max(Param_Lum_k_evol) /= TabParam_Lum_min(Param_Lum_k_evol)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Lum(Param_Lum_k_evol), Step_Lum(Param_Lum_k_evol))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for k_evol in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Lum_min(Param_Lum_k_evol)) .AND. (temp <= TabParam_Lum_max(Param_Lum_k_evol))) THEN
                     TabParam_Lum(Param_Lum_k_evol) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Lum(Param_Lum_k_evol) = TabParam_Lum_min(Param_Lum_k_evol)
            END IF
         END SELECT
      END IF
      
      IF(Ep_explore == 1) THEN
         SELECT CASE(Model_Ep)
         CASE(Model_EpFix)
            IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Ep(Param_Ep_Ep0)), Step_Ep(Param_Ep_Ep0)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Ep0 in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_Ep0)) .AND. (temp <= TabParam_Ep_max(Param_Ep_Ep0))) THEN
                     TabParam_Ep(Param_Ep_Ep0) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
            END IF
         CASE(Model_EpLogNormal)
            IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Ep(Param_Ep_Ep0)), Step_Ep(Param_Ep_Ep0)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Ep0 in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_Ep0)) .AND. (temp <= TabParam_Ep_max(Param_Ep_Ep0))) THEN
                     TabParam_Ep(Param_Ep_Ep0) = temp
                     EXIT
                  END IF
               END DO
               ELSE
               TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
            END IF
            IF(TabParam_Ep_max(Param_Ep_sigmaLog) /= TabParam_Ep_min(Param_Ep_sigmaLog)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Ep(Param_Ep_sigmaLog), Step_Ep(Param_Ep_sigmaLog))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for sigmaLog in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_sigmaLog)) .AND. (temp <= TabParam_Ep_max(Param_Ep_sigmaLog))) THEN
                     TabParam_Ep(Param_Ep_sigmaLog) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Ep(Param_Ep_sigmaLog) = TabParam_Ep_min(Param_Ep_sigmaLog)
            END IF
         CASE(Model_EpAmati)
            IF(TabParam_Ep_max(Param_Ep_Ep0) /= TabParam_Ep_min(Param_Ep_Ep0)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Ep(Param_Ep_Ep0)), Step_Ep(Param_Ep_Ep0)))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for Ep0 in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_Ep0)) .AND. (temp <= TabParam_Ep_max(Param_Ep_Ep0))) THEN
                     TabParam_Ep(Param_Ep_Ep0) = temp
                     EXIT
                  END IF
               END DO
            ELSE
                  TabParam_Ep(Param_Ep_Ep0) = TabParam_Ep_min(Param_Ep_Ep0)
            END IF
            IF(TabParam_Ep_max(Param_Ep_sigmaLog) /= TabParam_Ep_min(Param_Ep_sigmaLog)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Ep(Param_Ep_sigmaLog), Step_Ep(Param_Ep_sigmaLog))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for sigmaLog in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_sigmaLog)) .AND. (temp <= TabParam_Ep_max(Param_Ep_sigmaLog))) THEN
                     TabParam_Ep(Param_Ep_sigmaLog) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Ep(Param_Ep_sigmaLog) = TabParam_Ep_min(Param_Ep_sigmaLog)
            END IF
            IF(TabParam_Ep_max(Param_Ep_L0) /= TabParam_Ep_min(Param_Ep_L0)) THEN
               safeguard = 0
               DO
                  temp = 10.d0**(gaussian(LOG10(TabParam_Ep(Param_Ep_L0)), Step_Ep(Param_Ep_L0)))
                   safeguard = safeguard + 1
                   IF(safeguard > 10000) STOP 'Looping too many times for L0 in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_L0)) .AND. (temp <= TabParam_Ep_max(Param_Ep_L0))) THEN
                     TabParam_Ep(Param_Ep_L0) = temp
                     EXIT
                  END IF
               END DO
               ELSE
               TabParam_Ep(Param_Ep_L0) = TabParam_Ep_min(Param_Ep_L0)
            END IF
            IF(TabParam_Ep_max(Param_Ep_alpha_amati) /= TabParam_Ep_min(Param_Ep_alpha_amati)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_Ep(Param_Ep_alpha_amati), Step_Ep(Param_Ep_alpha_amati))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for alpha_amati in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_Ep_min(Param_Ep_alpha_amati)) .AND. (temp <= TabParam_Ep_max(Param_Ep_alpha_amati))) THEN
                     TabParam_Ep(Param_Ep_alpha_amati) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_Ep(Param_Ep_alpha_amati) = TabParam_Ep_min(Param_Ep_alpha_amati)
            END IF
         END SELECT
      END IF

      IF(z_explore == 1) THEN
         SELECT CASE(Model_z)
         CASE(Model_z_evol)
            IF(TabParam_z_max(Param_z_zeta) /= TabParam_z_min(Param_z_zeta)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_z(Param_z_zeta), Step_z(Param_z_zeta))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for zeta in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_z_min(Param_z_zeta)) .AND. (temp <= TabParam_z_max(Param_z_zeta))) THEN
                     TabParam_z(Param_z_zeta) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_z(Param_z_zeta) = TabParam_z_min(Param_z_zeta)
            END IF
         CASE(Model_zBExp)
             IF(TabParam_z_max(Param_z_zm) /= TabParam_z_min(Param_z_zm)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_z(Param_z_zm), Step_z(Param_z_zm))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for zm in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_z_min(Param_z_zm)) .AND. (temp <= TabParam_z_max(Param_z_zm))) THEN
                     TabParam_z(Param_z_zm) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_z(Param_z_zm) = TabParam_z_min(Param_z_zm)
            END IF
            
            IF(TabParam_z_max(Param_z_a) /= TabParam_z_min(Param_z_a)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_z(Param_z_a), Step_z(Param_z_a))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for a in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_z_min(Param_z_a)) .AND. (temp <= TabParam_z_max(Param_z_a))) THEN
                     TabParam_z(Param_z_a) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_z(Param_z_a) = TabParam_z_min(Param_z_a)
            END IF
            
            IF(TabParam_z_max(Param_z_b) /= TabParam_z_min(Param_z_b)) THEN
               safeguard = 0
               DO
                  temp = gaussian(TabParam_z(Param_z_b), Step_z(Param_z_b))
                  safeguard = safeguard + 1
                  IF(safeguard > 10000) STOP 'Looping too many times for b in Draw_Markov_Jump()'
                  IF( (temp >= TabParam_z_min(Param_z_b)) .AND. (temp <= TabParam_z_max(Param_z_b))) THEN
                     TabParam_z(Param_z_b) = temp
                     EXIT
                  END IF
               END DO
            ELSE
               TabParam_z(Param_z_b) = TabParam_z_min(Param_z_b)
            END IF
         END SELECT
      END IF
    
 END SUBROUTINE Draw_Markov_Jump


 SUBROUTINE Reset_Markov_Chain()

   CALL Reset_MCMC_Step_Size()
   CALL Reset_temperature()
   IF(rank == master_proc) CALL Draw_model_Parameters()               ! Move to random place in parameter space
   CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   CALL MonteCarlo(hist_flag)                 ! Calculate Chi2
   rejected = 0                               ! Reset consecutive rejected jumps
   accepted_rec = 1                           ! By design this jump is accepted
   accepted = 0.d0                            ! Reset the accepted rate
   acceptance_ratio = 0.d0                    ! Reset ratio of accepted jumps
   steps_since_reset = 0                      ! Reset steps count
   MCMC_run_nb = MCMC_run_nb + 1
   IF((rank == master_proc) .AND. (verbose >= 1)) WRITE(*,*) "[RESET CHAIN] run : ",MCMC_run_nb," MonteCarlo OK"
   !IF(rank == master_proc) WRITE(*,*) "[RESET CHAIN] run : ",MCMC_run_nb," Chi2 : ", Chi2(0), " Lmin, slope : ",TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_slope)
   IF(rank == master_proc) WRITE(*,'(A,I6,A,F8.2,A,F7.4x,F6.4,A,2F6.3)') "[RESET CHAIN] run : ",MCMC_run_nb," Chi2 : ", Chi2(0), " Lmin, slope : ",LOG10(TabParam_Lum(Param_Lum_Lmin)),TabParam_Lum(Param_Lum_slope), " steps : Lmin, slope ", Step_Lum(Param_Lum_Lmin), Step_Lum(Param_Lum_slope)
    
   CALL Write_reprise_MCMC()
   
 END SUBROUTINE Reset_Markov_Chain

 
 SUBROUTINE Update_Markov_Chain()
   ! If the Markov jump is not accepted, revert to previous state
   IF (lum_explore == 1) TabParam_Lum = TabParam_Lum_old
   IF (Ep_explore == 1)  TabParam_Ep = TabParam_Ep_old
   IF (z_explore == 1 )  TabParam_z = TabParam_z_old
   Chi2 = Chi2_old
   lnL = lnL_old
   
 END SUBROUTINE Update_Markov_Chain

 SUBROUTINE Reset_temperature()
   tau = 1.d3
 END SUBROUTINE Reset_temperature

 SUBROUTINE Reset_MCMC_Step_Size()

   Step_Lum = 0.d0
   Step_Ep  = 0.d0
   Step_z   = 0.d0
   
   Step_Lum(Param_Lum_Lmin) = 0.02d0
   Step_Lum(Param_Lum_Lmax) = 0.025d0
   Step_Lum(Param_Lum_Lbreak) = 0.025d0
   Step_Lum(Param_Lum_slope) = 0.02d0
   Step_Lum(Param_Lum_slopeL) = 0.05d0
   Step_Lum(Param_Lum_slopeH) = 0.02d0
   Step_Lum(Param_Lum_k_evol) = 0.04d0
   
   Step_Ep(Param_Ep_Ep0) = 0.0075d0
   Step_Ep(Param_Ep_sigmaLog) = 0.01d0
   Step_Ep(Param_Ep_alpha_amati) = 0.01d0

   Step_z(Param_z_zeta) = 0.01d0
   Step_z(Param_z_zm)   = 0.03d0
   Step_z(Param_z_a) = 0.03d0
   Step_z(Param_z_b) = 0.01d0
   
 END SUBROUTINE Reset_MCMC_Step_Size

 SUBROUTINE Update_MCMC_Step_Size(downscale, upscale)
   REAL(8), intent(in) :: downscale, upscale

   IF (acceptance_ratio > 0.234d0) THEN
      Step_Lum = Step_Lum * upscale
      
      Step_Ep = Step_Ep * upscale
      
      Step_z = Step_z * upscale
      
   ELSE IF((acceptance_ratio <= 0.234d0) .AND. (acceptance_ratio > 0.d0)) THEN 
      Step_Lum = Step_Lum * downscale
      
      Step_Ep = Step_Ep * downscale
      Step_z = Step_z * downscale
   END IF
      
 END SUBROUTINE Update_MCMC_Step_Size

 SUBROUTINE Update_temperature(kappa_tau, tau_lim)
   REAL(8), intent(in) :: kappa_tau
   REAL(8), intent(in) :: tau_lim
   
   tau = kappa_tau * (tau - tau_lim) + tau_lim
   
 END SUBROUTINE Update_temperature
 
 SUBROUTINE Write_reprise_MCMC()
   ! Write in reprise                                                                                            
   OPEN(UNIT=43, FILE=TRIM(path)//'reprise_MCMC.dat', FORM='unformatted', POSITION='append')                     
   IF(rank == master_proc) WRITE(43) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &           
        & Chi2, lnL, k_Kommers, k_Stern, k_Preece, k_EpGBM, k_eBAT6, dof, accepted_rec, tau, Step_Lum, Step_z, Step_Ep,&
        & TabHistLogL_master, TabHistz_master, TabHistLogEp_master, TabHistLogEpobs_master, TabHistLogP_master, &
        & TabHistKomm_DRDP, TabHistPreece_Ep_master, TabHistStern_P23_master, TabHistEpGBM_Epobs_master, TabHisteBAT6_z_master
   CLOSE(43) ! Close reprise file
 END SUBROUTINE Write_reprise_MCMC

 
 SUBROUTINE Update_best_parameters()
   INTEGER :: i
   
   OPEN(UNIT=876, FILE=TRIM(path)//'best_chi2.txt')
   WRITE(876,'(A)') "# This file contains the information of the best model"
   WRITE(876,'(A)') "# "
   WRITE(876,'(A,ES12.5)') "# Best Chi2 : ", Chi2_min
   WRITE(876,'(A,I3)') "# degrees of freedom : ", dof
   WRITE(876,'(A,F6.3)') "# 3 sigma interval : ", delta_chi2_3
   WRITE(876,'(A)') "# "
   WRITE(876,'(A,I2,A)') "# seeds (KISS i,j,k) for ",nb_procs," cores : "
   DO i=0, nb_procs-1
      WRITE(876, *) TabSave_ijk_KISS(1,i),TabSave_ijk_KISS(2,i),TabSave_ijk_KISS(3,i)
   END DO
   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Constraints used : "
   DO i_Constraint=1, N_Constraints
      IF(Constraint_Included(i_Constraint)) WRITE(876,'(a)') "# - "//TRIM(TabConstraint_name(i_Constraint))
   END DO
   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Parameters : "
   WRITE(876,'(A)') "# "

   WRITE(876,'(A)') "# Luminosity function : (0) Fix (1 arg : L0)"
   WRITE(876,'(A)') "# Luminosity function : (1) Power Law (3 args : Lmin, Lmax, slope)"
   WRITE(876,'(A)') "# Luminosity function : (2) Evolving BPL (6 args : Lmin, Lmax, Lbreak, slopeL, slopeH, k_evol)"
   WRITE(876,'(A)') "# Luminosity function : (3) Evol. Power Law (3 args : Lmin, Lmax, slope, k_evol)"
   SELECT CASE(Model_Lum)
   CASE(Model_LumFix)
      TabBestParam_Lum(Param_Lum_L0) = TabParam_Lum(Param_Lum_L0)
      WRITE(876,'(I1)') Model_LumFix
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_L0)
   CASE(Model_LumPL)
      TabBestParam_Lum(Param_Lum_Lmin) = TabParam_Lum(Param_Lum_Lmin)
      TabBestParam_Lum(Param_Lum_Lmax) = TabParam_Lum(Param_Lum_Lmax)
      TabBestParam_Lum(Param_Lum_slope) = TabParam_Lum(Param_Lum_slope)
      WRITE(876,'(I1)') Model_LumPL
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmin)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmax)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slope)
   CASE(Model_LumBPL_evol)
      TabBestParam_Lum(Param_Lum_Lmin) = TabParam_Lum(Param_Lum_Lmin)
      TabBestParam_Lum(Param_Lum_Lmax) = TabParam_Lum(Param_Lum_Lmax)
      TabBestParam_Lum(Param_Lum_Lbreak) = TabParam_Lum(Param_Lum_Lbreak)
      TabBestParam_Lum(Param_Lum_slopeL) = TabParam_Lum(Param_Lum_slopeL)
      TabBestParam_Lum(Param_Lum_slopeH) = TabParam_Lum(Param_Lum_slopeH)
      TabBestParam_Lum(Param_Lum_k_evol) = TabParam_Lum(Param_Lum_k_evol)
      WRITE(876,'(I1)') Model_LumBPL_evol
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmin)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmax)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lbreak)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slopeL)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slopeH)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_k_evol)
   CASE(Model_LumPL_evol)
      TabBestParam_Lum(Param_Lum_Lmin)   = TabParam_Lum(Param_Lum_Lmin)
      TabBestParam_Lum(Param_Lum_Lmax)   = TabParam_Lum(Param_Lum_Lmax)
      TabBestParam_Lum(Param_Lum_slope)  = TabParam_Lum(Param_Lum_slope)
      TabBestParam_Lum(Param_Lum_k_evol) = TabParam_Lum(Param_Lum_k_evol)
      WRITE(876,'(I1)') Model_LumPL_evol
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmin)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmax)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slope)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_k_evol)
   END SELECT
   
   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Peak energy : (0) Fix (1 arg : Ep0)"
   WRITE(876,'(A)') "# Peak energy : (1) LogNormal  (2 args : Ep0[keV], sigmaLog)"
   WRITE(876,'(A)') "# Peak energy : (3) Amati-like (4 args : Ep0[keV], sigmaLog, L0[erg/s], alpha_amati)"
   ! DONT FORGET TO ADD THIS LATER (the other possible models)
   SELECT CASE(Model_Ep)
   CASE(Model_EpFix)
      TabBestParam_Ep(Param_Ep_Ep0) = TabParam_Ep(Param_Ep_Ep0)
      WRITE(876,'(I1)') Model_EpFix
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_Ep0)
   CASE(Model_EpLogNormal)
      TabBestParam_Ep(Param_Ep_Ep0) = TabParam_Ep(Param_Ep_Ep0)
      TabBestParam_Ep(Param_Ep_sigmaLog) = TabParam_Ep(Param_Ep_sigmaLog)
      WRITE(876,'(I1)') Model_EpLogNormal
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_Ep0)
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_sigmaLog)
   CASE(Model_EpAmati)
      TabBestParam_Ep(Param_Ep_Ep0)         = TabParam_Ep(Param_Ep_Ep0)
      TabBestParam_Ep(Param_Ep_sigmaLog)    = TabParam_Ep(Param_Ep_sigmaLog)
      TabBestParam_Ep(Param_Ep_L0)          = TabParam_Ep(Param_Ep_L0)
      TabBestParam_Ep(Param_Ep_alpha_amati) = TabParam_Ep(Param_Ep_alpha_amati)    
      WRITE(876,'(I1)') Model_EpAmati
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_Ep0)
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_sigmaLog)
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_L0)
      WRITE(876,'(ES12.5)') TabBestParam_Ep(Param_Ep_alpha_amati)
     
   END SELECT
 
   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Redshift : (0) z evolution (1 arg : zeta)"
   SELECT CASE(Model_z)
   CASE(Model_z_evol)
      TabBestParam_z(Param_z_zeta) = TabParam_z(Param_z_zeta)
      WRITE(876,'(I1)') Model_z_evol
      WRITE(876,'(ES12.5)') TabBestParam_z(Param_z_zeta)
   END SELECT

   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Normalizations : "
   WRITE(876,*) " Number of GRBs : ", Nb_GRB
   WRITE(876,*) " Kommers : ", k_Kommers
   WRITE(876,*) " Preece : ", k_Preece
   WRITE(876,*) " Stern : ", k_Stern
   WRITE(876,*) " EpGBM : ", k_EpGBM
   WRITE(876,*) " Pseudo collapse rate : ", pseudo_collapse_rate
    
   CLOSE(876)

 END SUBROUTINE Update_best_parameters

 SUBROUTINE Draw_Redshift(z)
   ! --- Redshift ---!
   REAL(8), INTENT(out) :: z
   INTEGER              :: imin, imax, j

   IF(Model_z /= Model_zFix) THEN
      DO
         t = uniform_GRB() ! Different RNG for reproductibility
         IF((t > 0.d0) .AND. (t < 1.0d0))  EXIT
      END DO
      imin = 1
      imax = izmax
      
      DO
         j = (imin+imax)/2
         IF (t >= TabFctDistrz(j)) THEN
            imin = j+1
         ELSE IF (t < TabFctDistrz(j-1)) THEN
            imax = j
         ELSE
            EXIT
         END IF
      END DO
      
      z   = TabPrecisez(j-1) + (TabPrecisez(j) - TabPrecisez(j-1)) * (t-TabFctDistrz(j-1)) / (TabFctDistrz(j)-TabFctDistrz(j-1))
      D_L = TabD_L(j-1)      + (TabD_L(j)      - TabD_L(j-1)     ) * (t-TabFctDistrz(j-1)) / (TabFctDistrz(j)-TabFctDistrz(j-1))                       
   END IF
   
   IF (verbose==2) PRINT*, "z =", z
   IF (verbose==2) PRINT*, "D_L = ", D_L, "[Mpc]"
   IF (ISNAN(z)) NaNtest_prop = .TRUE.
   
 END SUBROUTINE Draw_Redshift
 
 SUBROUTINE Draw_Luminosity(L)
   ! -- Luminosity [erg/s] -- !
   REAL(8), INTENT(out) :: L
   INTEGER              :: imin, imax, j
   
   SELECT CASE(Model_Lum)
   CASE(Model_LumFix)
      
   CASE(Model_LumPL) 
      t = uniform_GRB()
      IF(TabParam_Lum(Param_Lum_slope) == 1.d0) THEN
         L = TabParam_Lum(Param_Lum_Lmin) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**t
      ELSE
         L = TabParam_Lum(Param_Lum_Lmin) * &
              & ( 1.d0 - t*( 1.d0 - (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**(1.d0-TabParam_Lum(Param_Lum_slope)) )&
              & )**(1.d0/ (1.d0-TabParam_Lum(Param_Lum_slope)) )
      END IF
      
   CASE(Model_LumBPL_evol)
      t = uniform_GRB()
      IF(TabParam_Lum(Param_Lum_slopeL) == 1) THEN
         IF(TabParam_Lum(Param_Lum_slopeH) == 1) THEN
            L = TabParam_Lum(Param_Lum_Lmin)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**t
         ELSE
            IF(t < t_star) THEN 
               L = TabParam_Lum(Param_Lum_Lmin)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))**(t/t_star)
            ELSE
               L = TabParam_Lum(Param_Lum_Lbreak)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * ( 1.d0 - (t-t_star) * (1.d0 -&
                    & ( TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak) )**(1.d0-TabParam_Lum(Param_Lum_slopeH)) )  / (1.d0-t_star) )&
                    & **( 1.d0/(1.d0 - TabParam_Lum(Param_Lum_slopeH)) )
            END IF
         END IF
      ELSE
         IF(TabParam_Lum(Param_Lum_slopeH) == 1) THEN
            t = uniform_GRB()
            IF(t < t_star) THEN
               L = TabParam_Lum(Param_Lum_Lmin)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (1.d0 - t/t_star * (1.d0 - (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))&
                    & **(1.d0 -TabParam_Lum(Param_Lum_slopeL))) )**(1.d0/(1.d0-TabParam_Lum(Param_Lum_slopeL)))
            ELSE 
               L = TabParam_Lum(Param_Lum_Lbreak)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak))&
                    & **((t-t_star)/(1.d0-t_star))
            END IF
         ELSE
            t = uniform_GRB()
            IF(t < t_star) THEN
               L = TabParam_Lum(Param_Lum_Lmin)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (1.d0 - t/t_star * (1.d0 - (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))&
                    & **(1.d0 -TabParam_Lum(Param_Lum_slopeL))) )**(1.d0/(1.d0-TabParam_Lum(Param_Lum_slopeL)))
            ELSE
               L = TabParam_Lum(Param_Lum_Lbreak)* (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * ( 1.d0 - (t-t_star) * (1.d0 -&
                    & ( TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak) )**(1.d0-TabParam_Lum(Param_Lum_slopeH)) )  / (1.d0-t_star) )&
                    & **( 1.d0/(1.d0 - TabParam_Lum(Param_Lum_slopeH)) )
            END IF
         END IF
      END IF

   CASE(Model_LumPL_evol) 
      t = uniform_GRB()
      IF(TabParam_Lum(Param_Lum_slope) == 1.d0) THEN
         L = TabParam_Lum(Param_Lum_Lmin) * (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**t
      ELSE
         L = TabParam_Lum(Param_Lum_Lmin) * (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * &
              & ( 1.d0 - t*( 1.d0 - (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**(1.d0-TabParam_Lum(Param_Lum_slope)) )&
              & )**(1.d0/ (1.d0-TabParam_Lum(Param_Lum_slope)) )
      END IF
   CASE(Model_LumSch)    
      DO
         t = uniform_GRB() ! Different RNG for reproductibility
         IF((t > 0.d0) .AND. (t < 1.0d0))  EXIT
      END DO
      imin = 1
      imax = M_L
      
      DO
         j = (imin+imax)/2
         IF (t >= TabFctDistrLum(j)) THEN
            imin = j+1
         ELSE IF (t < TabFctDistrLum(j-1)) THEN
            imax = j
         ELSE
            EXIT
         END IF
      END DO

      L = TabLogL_CDF(j-1) + (TabLogL_CDF(j) - TabLogL_CDF(j-1)) * (t-TabFctDistrLum(j-1)) / (TabFctDistrLum(j)-TabFctDistrLum(j-1))
      L = (1.d0+z)**TabParam_Lum(Param_Lum_k_evol) * 10.d0**(L)
      
   CASE DEFAULT
      STOP "Monte Carlo (2) : ERROR Lum, not implemented yet"
      
   END SELECT
   IF (verbose==2) PRINT*, "L =", L, "[erg/s]"
   IF (ISNAN(L)) NaNtest_prop = .TRUE.
   
 END SUBROUTINE Draw_Luminosity

 SUBROUTINE Draw_alpha_beta(alpha, beta)
   ! --- Spectral Model --- ! 
   REAL(8), INTENT(out) :: alpha, beta
   INTEGER              :: imin, imax, j
   
   IF(Model_Spec == Model_SpecBandD) THEN
      
      DO 
         alpha = gaussian_GRB(1.d0, 0.5d0)
         IF (alpha < 2.0d0) EXIT
      END DO
      
      DO 
         t = uniform_GRB()
         IF( t <= 0.15625d0) THEN
            beta =  2.d0 + SQRT(t/2.5d0)
         ELSE
            beta = (9.d0 - SQRT(13.5d0*(1.d0-t)) )/2.5d0
         END IF
         IF (beta > 2.0d0) EXIT
      END DO
      
   ELSE IF (Model_Spec == Model_SpecBandGBM) THEN
      DO
         t = uniform_GRB()
         IF((t > 0.d0) .AND. (t < 1.0d0))  EXIT
      END DO
      imin = 1
      imax = N_GBM_alpha
      DO
         j = (imin+imax)/2
         IF (t >= TabFctDistrGBM_alpha(j)) THEN
            imin = j+1
         ELSE IF (t < TabFctDistrGBM_alpha(j-1)) THEN
            imax = j
         ELSE
            EXIT
         END IF
      END DO
 
      alpha = TabGBM_alpha(j-1) + (TabGBM_alpha(j) - TabGBM_alpha(j-1)) * (t-TabFctDistrGBM_alpha(j-1)) / (TabFctDistrGBM_alpha(j)-TabFctDistrGBM_alpha(j-1))
      alpha = -alpha ! (because convention is opposite in GBM)
      IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, alpha   OK : ", alpha
      DO
         t = uniform_GRB()
         IF((t > 0.d0) .AND. (t < 1.0d0)) EXIT
      END DO
     IF((rank == master_proc) .AND. (verbose == 2))  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, t for beta   OK : ", t
      imin = 1
      imax = N_GBM_beta
      DO
         j = (imin+imax)/2
         IF (t >= TabFctDistrGBM_beta(j)) THEN
            imin = j+1
         ELSE IF (t < TabFctDistrGBM_beta(j-1)) THEN
            imax = j
         ELSE
            EXIT
         END IF
      END DO

      beta = TabGBM_beta(j-1) + (TabGBM_beta(j) - TabGBM_beta(j-1)) * (t-TabFctDistrGBM_beta(j-1)) / (TabFctDistrGBM_beta(j)-TabFctDistrGBM_beta(j-1))
      beta = -beta ! (because convention is opposite in GBM)
       IF((rank == master_proc) .AND. (verbose >= 1))  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, beta   OK : ", beta
   END IF
   ktild = Calc_ktild(alpha, beta, Model_Spec) 
   IF (ISNAN(ktild)) NaNtest_prop = .TRUE.
   IF (ISNAN(alpha)) NaNtest_prop = .TRUE.
   IF (ISNAN(beta)) NaNtest_prop = .TRUE.
   IF(ISNAN(ktild)) WRITE(*,'(A, 3ES12.5)')  "param(ktild, alpha, beta-alpha)  :  ", ktild, alpha, (beta-alpha) 
   IF(verbose == 2) THEN
      PRINT*, 'alpha = ', alpha
      PRINT*, 'beta = ', beta        
   END IF
   
   
 END SUBROUTINE Draw_alpha_beta
 
 SUBROUTINE Draw_Ep(Ep)
   ! --- Peak Energy [keV] --- !
   REAL(8), INTENT(out) :: Ep
   
   SELECT CASE(Model_Ep)
   CASE(Model_EpLogNormal)
      DO
         Ep = 10**( gaussian_GRB(LOG10(TabParam_Ep(Param_Ep_Ep0)), TabParam_Ep(Param_Ep_sigmaLog)) )
         IF(Ep > 0.d0) EXIT
      END DO
   CASE(Model_EpAmati)
          DO
         t  = gaussian_GRB(0.d0, TabParam_Ep(Param_Ep_sigmaLog))
         !IF(rank == master_proc)  WRITE(*,*) " in Draw_Ep  t : ", t
         !IF(rank == master_proc)  WRITE(*,*) " in Draw_Ep  L : ", L
         Ep = TabParam_Ep(Param_Ep_Ep0) * (L/TabParam_Ep(Param_Ep_L0))**TabParam_Ep(Param_Ep_alpha_amati) * &
              & 10.d0**( SQRT(1.d0 + TabParam_Ep(Param_Ep_alpha_amati)**2) * t )
         IF(Ep > 0.d0) EXIT
      END DO
      !IF(rank == master_proc)  WRITE(*,*) " in Draw_Ep  Ep0 : ", TabParam_Ep(Param_Ep_Ep0)
   END SELECT
   IF (ISNAN(Ep)) NaNtest_prop = .TRUE.
   IF(verbose==2) PRINT*, "Ep : ", Ep, "[keV]"
   
 END SUBROUTINE Draw_Ep
 
 
 SUBROUTINE InitBoundaries()
   ! Creates the limits for the histograms
   ! Fills TabEmin and TabEmax with the values for each instrument
   ! Creates the detection threshold for each sample
   
   ! --- Histogram limits --- !
   ! ------------------------ !
   
   ! --- Luminosity [erg/s] --- !
   TabHistlim_inf(Prop_LogL) = 46.d0
   TabHistlim_sup(Prop_LogL) = 56.d0
   
   ! --- Redshift --- !
   TabHistlim_inf(Prop_z) = 0.d0
   TabHistlim_sup(Prop_z) = z_maximum

   ! --- Peak energy [keV] --- !
   TabHistlim_inf(Prop_Ep) = 1.d-1
   TabHistlim_sup(Prop_Ep) = 1.d4

   ! --- Peak flux [ph/cm^2/s] --- !
   TabHistlim_inf(Prop_LogP) = -4.d0
   TabHistlim_sup(Prop_LogP) = 4.d0

   ! --- Spec alpha --- !
   TabHistlim_inf(Prop_alpha) = -1.d0
   TabHistlim_sup(Prop_alpha) = 1.9d0

   ! --- Spec beta --- !
   TabHistlim_inf(Prop_beta) = 2.d0
   TabHistlim_sup(Prop_beta) = 22.d0

   
   ! --- Filling TabEmin and TabEmax --- !
   ! ----------------------------------- !
      
   ! BATSE
   TabEmin(Instrument_BATSE) = 50.d0
   TabEmax(Instrument_BATSE) = 300.d0

   ! SWIFT
   TabEmin(Instrument_BAT  ) = 15.d0
   TabEmax(Instrument_BAT  ) = 150.d0

   ! FREGATE
   TabEmin(Instrument_FREGATE) = 30.d0
   TabEmax(Instrument_FREGATE) = 400.d0

   ! WXM
   TabEmin(Instrument_WXM) = 2.d0
   TabEmax(Instrument_WXM) = 10.d0

   ! ECLAIRs
   TabEmin(Instrument_ECLAIRs) = 4.d0
   TabEmax(Instrument_ECLAIRs) = 120.d0


   ! --- Detection Threshold [ph/cm2/s] --- !
   ! -------------------------------------- !

   Threshold(Sample_Intrinsic)   = 0.d0
   Threshold(Sample_Kommers)     = 0.d0   ! 50-300 keV
   Threshold(Sample_Preece)      = 5.d0   ! 50-300 keV
   Threshold(Sample_Stern)       = 0.066825d0 ! 50-300 keV
   Threshold(Sample_SWIFTweak)   = 0.01d0 ! 15-150 keV
   Threshold(Sample_SWIFT)       = 0.2d0  ! 15-150 keV
   Threshold(Sample_SWIFTbright) = 1.d0   ! 15-150 keV
   Threshold(Sample_HETE2)       = 1.d0   ! 2-10 keV or 30-400 keV
   Threshold(Sample_eBAT6)       = 2.6d0  ! 15-150 keV
   Threshold(Sample_EpGBM)       = 0.9d0  ! 50-300 keV
   Threshold(Sample_SVOM)        = bkgECLAIRsB1  ! 4-150 keV

   
 END SUBROUTINE InitBoundaries
 
 SUBROUTINE Init_RNG()
   INTEGER :: i,j
   CHARACTER :: skip
   
   IF(RNG == Kiss_rng) THEN
      ! This is the parallel implementation
      ! These seeds were generated to garantee 10**12 draws from each core before overlapping
      ! 1st core
      TabInit_ijk_Kiss(1,0) =123456789 
      TabInit_ijk_Kiss(2,0) =362436069 
      TabInit_ijk_Kiss(3,0) =521288629 
      ! 2nd core
      TabInit_ijk_Kiss(1,1) =293731605 
      TabInit_ijk_Kiss(2,1) =-985101706 
      TabInit_ijk_Kiss(3,1) =909351951 
      ! 3rd core
      TabInit_ijk_Kiss(1,2) =-408408811 
      TabInit_ijk_Kiss(2,2) =1991222469 
      TabInit_ijk_Kiss(3,2) =1423132033 
      ! 4th core
      TabInit_ijk_Kiss(1,3) =-1982964459 
      TabInit_ijk_Kiss(2,3) =-1145452629 
      TabInit_ijk_Kiss(3,3) =946873952  
      ! 5th core
      TabInit_ijk_Kiss(1,4) =-134968043 
      TabInit_ijk_Kiss(2,4) =27302573  
      TabInit_ijk_Kiss(3,4) =1817789790 
      ! 6th core
      TabInit_ijk_Kiss(1,5) =840613141  
      TabInit_ijk_Kiss(2,5) =-1675529275 
      TabInit_ijk_Kiss(3,5) =1643740297  
      ! 7th core
      TabInit_ijk_Kiss(1,6) =943779093  
      TabInit_ijk_Kiss(2,6) =1835649319  
      TabInit_ijk_Kiss(3,6) =1356939134  
      ! 8th core
      TabInit_ijk_Kiss(1,7) =174529813  
      TabInit_ijk_Kiss(2,7) =1388453552  
      TabInit_ijk_Kiss(3,7) =814083514



!!$      ! TEMPPPP
!!$        ! 1st core
!!$      TabInit_ijk_Kiss(1,0) =491832654
!!$      TabInit_ijk_Kiss(2,0) =918875735
!!$      TabInit_ijk_Kiss(3,0) =-98911476 
!!$      ! 2nd core
!!$      TabInit_ijk_Kiss(1,1) =356972844
!!$      TabInit_ijk_Kiss(2,1) =-95462569 
!!$      TabInit_ijk_Kiss(3,1) =847552878 
!!$      ! 3rd core
!!$      TabInit_ijk_Kiss(1,2) =-36149117 
!!$      TabInit_ijk_Kiss(2,2) =266254988 
!!$      TabInit_ijk_Kiss(3,2) =589704982
!!$      ! 4th core
!!$      TabInit_ijk_Kiss(1,3) =-013698126 
!!$      TabInit_ijk_Kiss(2,3) =-128796635 
!!$      TabInit_ijk_Kiss(3,3) =462139459  


      
      IF(reprise .EQV. .TRUE.) THEN
         OPEN(UNIT=74, FILE=TRIM(path)//'save_KISS.dat')
         READ(74,*) skip ! skip line
         DO i=0, nb_procs-1
            READ(74, *) TabInit_ijk_Kiss(1,i), TabInit_ijk_Kiss(2,i), TabInit_ijk_Kiss(3,i)
         END DO
         CLOSE(74)
      END IF
     
      
      ! Set the different seeds 
      iKiss = TabInit_ijk_Kiss(1,rank)
      jKiss = TabInit_ijk_Kiss(2,rank)
      kKiss = TabInit_ijk_Kiss(3,rank)


   ELSE IF (RNG == MT19937) THEN
      TabInit_seed_MT(0) = 85105
      TabInit_seed_MT(1) = 12345
      TabInit_seed_MT(2) = 71123
      TabInit_seed_MT(3) = 2153
      TabInit_seed_MT(4) = 2862
      TabInit_seed_MT(5) = 35836
      TabInit_seed_MT(6) = 58214
      TabInit_seed_MT(7) = 457622
      
      ! Set the different seeds
      MT_seed = TabInit_seed_MT(rank)
      CALL sgrnd(MT_seed)

   END IF

         
 END SUBROUTINE Init_RNG

 SUBROUTINE Reset_GRB_seeds()
   ! Subroutine to reset the seeds for the GRB property drawings to ensure a smooth chi2 surface
   iKiss_GRB = TabInit_ijk_Kiss(1,rank)
   jKiss_GRB = TabInit_ijk_Kiss(2,rank)
   kKiss_GRB = TabInit_ijk_Kiss(3,rank)
   
 END SUBROUTINE Reset_GRB_seeds
 
 SUBROUTINE TestKiss()
  INTEGER :: i
  REAL(8) :: x
  
  DO i=1, 1000000000
     x = uniform()
     IF (MODULO(i,100000000)==0) WRITE(*,'(4(A,I12))') "i=",i,"  Kiss : i=",iKiss," j=",jKiss," k=",kKiss
  END DO
END SUBROUTINE TestKiss

SUBROUTINE Save_KISS()
  INTEGER :: i
  TabSave_ijk_KISS = 0
  TabSave_ijk_KISS(1,rank) = iKiss
  TabSave_ijk_KISS(2,rank) = jKiss
  TabSave_ijk_KISS(3,rank) = kKiss
  CALL MPI_REDUCE(TabSave_ijk_KISS(1,:), TabSave_ijk_KISS(1,:), nb_procs, MPI_INT, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_REDUCE(TabSave_ijk_KISS(2,:), TabSave_ijk_KISS(2,:), nb_procs, MPI_INT, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
  CALL MPI_REDUCE(TabSave_ijk_KISS(3,:), TabSave_ijk_KISS(3,:), nb_procs, MPI_INT, MPI_SUM, master_proc, MPI_COMM_WORLD, code)

  IF(rank == master_proc) THEN
     OPEN(UNIT=325, FILE=TRIM(path)//"save_KISS.dat")
     WRITE(325, '(A,I2,A)') '# Last saved KISS seeds i,j,k for ',nb_procs,' cores : ' 
     DO i=0, nb_procs-1
        WRITE(325, *) TabSave_ijk_KISS(1,i),TabSave_ijk_KISS(2,i),TabSave_ijk_KISS(3,i)
     END DO
     CLOSE(325)
  END IF
END SUBROUTINE Save_KISS

 ! -------------------------- !
 ! ----- Functions used ----- !
 ! -------------------------- !

FUNCTION Calc_DeltaChi2(p,nu)
  ! Frederic provided this, might need to ask him to understand it better
  REAL(8), INTENT(in) :: p,nu
  REAL(8)             :: Calc_DeltaChi2
  REAL(8) :: xmin,x,xmax,ymin,ymax,y
  INTEGER :: i
  
  xmin = 1.d0
  ymin = GammaQ((nu/2.d0),(xmin/2.d0))-1.d0+p
  xmax = 1000.d0
  ymax = GammaQ((nu/2.d0),(xmax/2.d0))-1.d0+p

  DO i=1, 100
    x = (xmin+xmax)/2.d0
    y = GammaQ((nu/2.d0),(x/2.d0))-1.d0+p
    IF (y*ymax>0.d0) THEN
      xmax = x
      ymax = y
    ELSE
      xmin = x
      ymin = y
    END IF
  END DO
  Calc_DeltaChi2 = x

END FUNCTION Calc_DeltaChi2

 ! Uniform variable on [0,1]
 FUNCTION uniform()
   REAL(8) :: uniform
   INTEGER :: j
   INTEGER, PARAMETER :: Maxi=2147483647
   IF(RNG == Kiss_rng) THEN
      j=Kiss()
      uniform=(REAL(j,8)+REAL(Maxi,8)+1.d0)/(2.d0*REAL(Maxi,8)+1.d0)
   ELSE IF( RNG == MT19937) THEN
   uniform = grnd()   ! This is the MT19937 RNG 
   ELSE
      CALL RANDOM_NUMBER(uniform)
   END IF
 END FUNCTION uniform

! Uniform variable on [0,1]
 FUNCTION uniform_GRB()
   ! Different RNG for GRB sampling 
   REAL(8) :: uniform_GRB
   INTEGER :: j
   INTEGER, PARAMETER :: Maxi=2147483647
  
   j=Kiss_GRB()
   uniform_GRB=(REAL(j,8)+REAL(Maxi,8)+1.d0)/(2.d0*REAL(Maxi,8)+1.d0)
  
 END FUNCTION uniform_GRB


 FUNCTION Kiss()
   INTEGER :: Kiss
   iKiss=69069*iKiss+32606797
   jKiss=Melange(Melange(jKiss,17),-15)
   kKiss=Melange(IAND(Melange(kKiss,18),2147483647),-13)
   Kiss=iKiss+jKiss+kKiss
 END FUNCTION Kiss
 
 FUNCTION Kiss_GRB()
   INTEGER :: Kiss_GRB
   iKiss_GRB=69069*iKiss_GRB+32606797
   jKiss_GRB=Melange(Melange(jKiss_GRB,17),-15)
   kKiss_GRB=Melange(IAND(Melange(kKiss_GRB,18),2147483647),-13)
   Kiss_GRB=iKiss_GRB+jKiss_GRB+kKiss_GRB
 END FUNCTION Kiss_GRB

 
 FUNCTION Melange(k,n)
   INTEGER, INTENT(In) :: k,n
   INTEGER             :: Melange
   Melange=IEOR(k,ISHFT(k,n))
 END FUNCTION Melange



 REAL(8) FUNCTION gaussian(mu, sigma)
   ! Random draw with a gaussian distribution
   REAL(8), INTENT(in) :: mu, sigma
   REAL(8)             :: t, theta

   DO
      t = uniform()
      IF(t>0) EXIT
   END DO
   theta = uniform()
   theta    = theta * 2*Pi
   gaussian = mu + SQRT(2.d0) * sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   ! sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   ! for comparison with Daigne 2006 remove sqrt(2)
 END FUNCTION gaussian
 
 REAL(8) FUNCTION gaussian_GRB(mu, sigma)
   ! Random draw with a gaussian distribution
   REAL(8), INTENT(in) :: mu, sigma
   REAL(8)             :: t, theta

   DO
      t = uniform_GRB()
      IF(t>0) EXIT
   END DO
   theta = uniform_GRB()
   theta    = theta * 2*Pi
   gaussian_GRB = mu + SQRT(2.d0) * sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   ! sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   ! for comparison with Daigne 2006 remove sqrt(2)
 END FUNCTION gaussian_GRB


 REAL(8) FUNCTION Rz(z, Model_z, Tabparam_z)
   ! GRB rate as a function of redshift [GRB/yr/Mpc^3]
   ! b   : slope at low z
   ! b-a : slope at high z
   
   REAL(8), INTENT(in)                  :: z
   INTEGER, INTENT(in)                  :: Model_z
   REAL(8), DIMENSION(1:10), INTENT(in) :: Tabparam_z
   REAL(8)                              :: a, b, c, d, zm, zeta
   
   SELECT CASE(Model_z)
   CASE(Model_zFix)
      Rz = 1.d0
   CASE(Model_zUniform)
      Rz = 1.d0
   CASE(Model_zSH)
      zm = TabParam_z(Param_z_zm)
      a  = TabParam_z(Param_z_a)
      b  = TabParam_z(Param_z_b)
      Rz = 0.00132d0 * a * EXP(b*(z-zm)) / ( a-b + b*EXP(a*(z-zm)) )
   CASE(Model_zDaigne)
      a = TabParam_z(Param_z_a)
      b = TabParam_z(Param_z_b)
      c = TabParam_z(Param_z_c)
      d = TabParam_z(Param_z_d)
      Rz = 0.0122d0 * a * EXP(b*z) / ( EXP(c*z) + d )
   CASE(Model_z_evol)
      a  = TabParam_z(Param_z_a)
      b  = TabParam_z(Param_z_b)
      zm   = TabParam_z(Param_z_zm)
      zeta = TabParam_z(Param_z_zeta)
      ! normalization comes from IMF Salpeter
      Rz = 0.00132d0 * z_evolution(z, zm, zeta) * a * EXP(b*(z-zm)) / ( a-b + b*EXP(a*(z-zm)) )
   CASE(Model_zLi)
      a = TabParam_z(Param_z_a)
      b = TabParam_z(Param_z_b)
      c = TabParam_z(Param_z_c)
      d = TabParam_z(Param_z_d)
      ! Normalization from Salpeter IMF
      Rz = 0.007422d0 * (a + b*z) / ( 1.d0 + (z/c)**d )
   CASE(Model_zBPL)
      a  = TabParam_z(Param_z_a)
      b  = TabParam_z(Param_z_b)
      zm = TabParam_z(Param_z_zm)
      ! No normalization (don't need for pdf)
      IF(z <= zm) THEN
         Rz = (1.d0 + z)**a
      ELSE
         Rz = (1.d0 + z)**b * (1.d0 + zm)**(a-b)
      END IF
   CASE(Model_zBExp)
      a  = TabParam_z(Param_z_a)
      b  = TabParam_z(Param_z_b)
      zm = TabParam_z(Param_z_zm)
      ! Normalization from fit to SH with Salpeter IMF  
      IF(z <= zm) THEN
         Rz =0.00033313d0 * EXP(a*z)
      ELSE
         Rz =0.00033313d0 * EXP(b*z) * EXP((a-b)*zm)
      END IF
   END SELECT
   
 END FUNCTION Rz

 REAL(8) FUNCTION z_evolution(z, zm, zeta)
   REAL(8), INTENT(in) :: z, zm, zeta
   ! Function that returns the redshift dependence of the "efficiency" of collapses to form GRBs
   !IF(z <= zm) THEN
   !   z_evolution = 1.d0
   !ELSE
      !z_evolution = ( (1.d0+z)/(1.d0+zm) )**zeta
   z_evolution = EXP(zeta*z)
   !END IF
 END FUNCTION z_evolution

 REAL(8) FUNCTION Calc_Epobs(Ep,z)
   REAL(8), INTENT(in) :: Ep, z
   ! --- Peak Energy (obs) [keV] --- !
   Calc_Epobs = Ep / (1.d0 + z)
   IF (ISNAN(Epobs)) NaNtest_prop = .TRUE.
   IF(verbose==2) PRINT*, "Epobs : ", Calc_Epobs, "[keV]"
   
 END FUNCTION Calc_Epobs


 REAL(8) FUNCTION Calc_t_star(TabParam_Lum)
   REAL(8), DIMENSION(1:10), INTENT(in) :: TabParam_Lum
   REAL(8)                              :: slopeL, slopeH, Lbreak, Lmin, Lmax

   Lmin   = TabParam_Lum(Param_Lum_Lmin)
   Lmax   = TabParam_Lum(Param_Lum_Lmax)
   Lbreak = TabParam_Lum(Param_Lum_Lbreak)
   slopeL = TabParam_Lum(Param_Lum_slopeL)
   slopeH = TabParam_Lum(Param_Lum_slopeH)

   
   IF(slopeL == 1) THEN
      IF(slopeH == 1) THEN
         STOP "Problem in Calc_t_star : slopeL = 1 and slopeH = 1"
      ELSE
         Calc_t_star = LOG(Lbreak/Lmin) / ( LOG(Lbreak/Lmin) + (slopeH-1.d0)**(-1) * ( 1.d0 - (Lmax/Lbreak)**(1.d0-slopeH) ) ) 
      END IF
   ELSE
      IF(slopeH == 1) THEN
         Calc_t_star = ( Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL) ) / ( (Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL)) + (1.d0-slopeL)*Lbreak**(1.d0-slopeL)*LOG(Lmax/Lbreak) )
 
      ELSE
         Calc_t_star = ( Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL) ) / &
              & ( (Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL)) + Lbreak**(slopeH-slopeL) * (1.d0-slopeL)/(1.d0-slopeH) * (Lmax**(1.d0-slopeH) - Lbreak**(1.d0-slopeH)) )
      END IF
    END IF
  END FUNCTION Calc_t_star

 REAL(8) FUNCTION Calc_Lum(L, Model_Lum, TabParam_Lum)
   ! Calc_Lum returns the probability of getting the luminosity L with a given model
   ! Luminosity in [erg/s]
   REAL(8), INTENT(in)                  :: L
   INTEGER, INTENT(in)                  :: Model_Lum
   REAL(8), DIMENSION(1:10), INTENT(in) :: TabParam_Lum
   REAL(8)                              :: slope, Lmin, Lmax
   REAL(8)                              :: slopeL, slopeH, Lbreak, t_star

   SELECT CASE(Model_Lum)

   CASE(Model_LumFix)
      Calc_Lum = 0.d0

   CASE(Model_LumPL)
      Lmin  = TabParam_Lum(Param_Lum_Lmin)
      Lmax  = TabParam_Lum(Param_Lum_Lmax)
      slope = TabParam_Lum(Param_Lum_slope)

      IF(L < Lmin .OR. L > Lmax) THEN
         Calc_Lum = 0.d0
      ELSE
         IF(slope == 1.d0) THEN
            Calc_Lum = 1.d0/(LOG(Lmax/Lmin)*L)
         ELSE 
            Calc_Lum = L**(-slope) * (slope-1.d0) / ( Lmin**(1.d0-slope) - Lmax**(1.d0-slope) )
         END IF
      END IF

   CASE(Model_LumBPL_evol)
      Lmin   = TabParam_Lum(Param_Lum_Lmin)
      Lmax   = TabParam_Lum(Param_Lum_Lmax)
      Lbreak = TabParam_Lum(Param_Lum_Lbreak)
      slopeL = TabParam_Lum(Param_Lum_slopeL)
      slopeH = TabParam_Lum(Param_Lum_slopeH)

      IF(L <= Lmin .OR. L >= Lmax) THEN
         Calc_Lum = 0.d0
      ELSE
         IF (slopeL == 1.d0 .AND. slopeH == 1.d0) THEN
            Calc_Lum = 1.d0/(LOG(Lmax/Lmin)*L)
            
         ELSE IF (slopeL == 1.d0 .AND. slopeH /= 1.d0) THEN
            
            t_star = LOG(Lbreak/Lmin) / ( LOG(Lbreak/Lmin) + (slopeH-1.d0)**(-1) * ( 1.d0 - (Lmax/Lbreak)**(1.d0-slopeH) ) )
            IF (L <= Lbreak) THEN
               Calc_Lum = t_star / ( L* LOG(Lbreak/Lmin) )
            ELSE IF (L > Lbreak) THEN
               Calc_Lum = (1.d0 - t_star) * (slopeH - 1.d0) * L**(-slopeH) / ( Lbreak**(1.d0 - slopeH) - Lmax**(1.d0 - slopeH) )
            END IF
            
         ELSE IF (slopeL /= 1.d0 .AND. slopeH == 1.d0) THEN
            
            t_star = ( L**(1.d0-slopeL) - Lmin**(1.d0-slopeL) ) / ( (Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL)) + (1.d0-slopeL)*Lbreak**(1.d0-slopeL)*LOG(Lmax/Lbreak) )
            IF (L <= Lbreak) THEN
               Calc_Lum = t_star * (slopeH - 1.d0) * L**(-slopeH) / ( Lbreak**(1.d0 - slopeH) - Lmax**(1.d0 - slopeH) )
            ELSE IF (L > Lbreak) THEN
               Calc_Lum = (1.d0 - t_star) / ( L* LOG(Lbreak/Lmin) )
            END IF
            
         ELSE IF (slopeL /= 1.d0 .AND. slopeH /= 1.d0) THEN
            
            t_star = ( L**(1.d0-slopeL) - Lmin**(1.d0-slopeL) ) / &
                 & ( (Lbreak**(1.d0-slopeL) - Lmin**(1.d0-slopeL)) + Lbreak**(slopeH-slopeL) * (1.d0-slopeL)/(1.d0-slopeH) * (Lmax**(1.d0-slopeH) - Lbreak**(1.d0-slopeH)) ) 
            IF (L <= Lbreak) THEN
               Calc_Lum = t_star * (slopeH - 1.d0) * L**(-slopeH) / ( Lbreak**(1.d0 - slopeH) - Lmax**(1.d0 - slopeH) )
            ELSE IF (L > Lbreak) THEN
               Calc_Lum = (1.d0 - t_star) * (slopeH - 1.d0) * L**(-slopeH) / ( Lbreak**(1.d0 - slopeH) - Lmax**(1.d0 - slopeH) )
            END IF
            
         ELSE 
            STOP "Problem in Calc_Lum"
         END IF

      END IF
   CASE(Model_LumSch)
      slope = TabParam_Lum(Param_Lum_slope)
      Lbreak = TabParam_Lum(Param_Lum_Lbreak)
      Calc_Lum = (L / Lbreak )**(-slope) * EXP(-L/Lbreak)
   CASE DEFAULT
      STOP "Calc_Lum : error, model not defined."
      
   END SELECT
  
 END FUNCTION Calc_Lum
 
 REAL(8) FUNCTION Calc_ktild(alpha, beta, Model_Spec)
   ! Calculates the normalization factor for the spectral shape
   ! Uses the Gamma functions for optimized speed. See stats.f90 for details
   REAL(8), INTENT(in) :: alpha, beta 
   INTEGER, INTENT(in) :: Model_Spec
   REAL(8)             :: x_c, xBx, Gamma, gln

   IF(Model_Spec == Model_SpecBPLFix) THEN
      Calc_ktild = (2.d0-alpha) * (beta-2.d0) / (beta-alpha)
   ELSE
      x_c     = (beta - alpha) / (2.d0 - alpha)
      CALL GammaSeries(Gamma, 2.d0-alpha, beta-alpha, gln)
      xBx = Gamma * EXP(gln) /( (2.d0 - alpha)**(2.d0-alpha) )
      IF(verbose == 2) PRINT*, "ktild integral up to x_c = ", xBx
      xBx = xBx +  x_c**(2.d0-alpha) * EXP(alpha-beta) / (beta - 2.d0)
      Calc_ktild = 1.d0/xBx
   END IF

 END FUNCTION Calc_ktild
 
 SUBROUTINE Calculate_Peakflux_Instrument()
   ! --- Number of Photons between Emin and Emax [ph/cm^2/s] --- !

   DO i_Instrument=1, N_Instruments
      IF(Instrument_Included(i_Instrument)) THEN
         Peakflux_Instrument(i_Instrument) = Nb_ph(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(i_Instrument), TabEmax(i_Instrument))
         IF (verbose==2) THEN 
            WRITE(*,'(A,A20,ES12.5,A,F3.0,A,F4.0,A)')&
                 & " Instrument : " , TRIM(TabInstrument_name(i_Instrument)) // " Peakflux =",Peakflux_Instrument(i_Instrument),&
                 & " [ph/cm2/s between ",TabEmin(i_Instrument)," and ", TabEmax(i_Instrument), " keV]"
         END IF
         IF (ISNAN(Peakflux_Instrument(i_Instrument))) NaNtest_prop = .TRUE.
      END IF
   END DO
   
 END SUBROUTINE Calculate_Peakflux_Instrument

 SUBROUTINE Calculate_Detection_Probability()
   REAL(8) :: cts_ECLAIRs, Prob_det_ij=0.d0
   INTEGER :: i,j
   
   ! ---------------- Detection Probability [0,1] ---------------- !
   
   DO i_Sample=0, N_Samples
      IF(i_Sample >= 1) THEN
         IF(Sample_Included(i_Sample)) THEN      
            SELECT CASE(i_Sample)
            CASE(Sample_Kommers)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               ! BATSE23 (Kommers et al. 2000, eq (4) ) :
               !PRINT*,Peakflux_Instrument(Instrument_BATSE), ERF(-4.801d0 + 29.868d0 *1.d0 )
               IF(Peakflux_Instrument(Instrument_BATSE) <= 1.d0) THEN ! Had to add this because above 1.d0, ERF function gets Floating Point Exception : erroneous arithmetic operation.
                  Prob_det(Sample_Kommers) = 0.5d0 * ( 1.d0 + ERF(-4.801d0 + 29.868d0 * Peakflux_Instrument(Instrument_BATSE) ))
               ELSE
                  Prob_det(Sample_Kommers) = 1.d0
               END IF
               IF (ISNAN(Prob_det(Sample_Kommers))) NaNtest_prop = .TRUE.
               
            CASE(Sample_Preece)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               IF(Peakflux_Instrument(Instrument_BATSE) >= Threshold(Sample_Preece)) THEN
                  Prob_det(Sample_Preece) = 1.d0
               ELSE
                  Prob_det(Sample_Preece) = 0.d0
               END IF

            CASE(Sample_Stern)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               IF(Peakflux_Instrument(Instrument_BATSE) >= Threshold(Sample_Stern)) THEN
                  Prob_det(Sample_Stern) = 1.d0
               ELSE
                  Prob_det(Sample_Stern) = 0.d0
               END IF
               
            CASE(Sample_SWIFTweak)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BAT  )
               IF(Peakflux_Instrument(Instrument_BAT  ) >= Threshold(Sample_Swiftweak)) THEN
                  Prob_det(Sample_Swiftweak) = 1.d0
               ELSE
                  Prob_det(Sample_Swiftweak) = 0.d0
               END IF
               
            CASE(Sample_SWIFT)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BAT  )
               IF(Peakflux_Instrument(Instrument_BAT  ) >= Threshold(Sample_Swift)) THEN
                  Prob_det(Sample_Swift) = 1.d0
               ELSE
                  Prob_det(Sample_Swift) = 0.d0
               END IF
               
            CASE(Sample_SWIFTbright) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BAT  )
               IF(Peakflux_Instrument(Instrument_BAT  ) >= Threshold(Sample_Swiftbright)) THEN
                  Prob_det(Sample_Swiftbright) = 1.d0
               ELSE
                  Prob_det(Sample_Swiftbright) = 0.d0
               END IF
               
            CASE(Sample_HETE2)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_FREGATE)
               IF(Peakflux_Instrument(Instrument_WXM) >= Threshold(Sample_HETE2) .OR. Peakflux_Instrument(Instrument_FREGATE) >= Threshold(Sample_HETE2) ) THEN
                  Prob_det(Sample_HETE2) = 1.d0
               ELSE
                  Prob_det(Sample_HETE2) = 0.d0
               END IF
            CASE(Sample_eBAT6) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BAT  )
               IF(Peakflux_Instrument(Instrument_BAT  ) >= Threshold(Sample_eBAT6)) THEN
                  Prob_det(Sample_eBAT6) = 1.d0
               ELSE
                  Prob_det(Sample_eBAT6) = 0.d0
               END IF
            CASE(Sample_EpGBM) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               IF(Peakflux_Instrument(Instrument_BATSE) >= Threshold(Sample_EpGBM)) THEN
                  Prob_det(Sample_EpGBM) = 1.d0
               ELSE
                  Prob_det(Sample_EpGBM) = 0.d0
               END IF
            CASE(Sample_SVOM)
               cts_ECLAIRs = Nb_cts_ECLAIRs(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(Instrument_ECLAIRs), TabEmax(Instrument_ECLAIRs)) ! on-axis
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_ECLAIRs)
               Prob_det(Sample_SVOM) = 0.D0
               DO i=-NoffECLAIRs, NoffECLAIRs
                  DO j=-NoffECLAIRs, NoffECLAIRs
                     IF(cts_ECLAIRs * TABoffECLAIRs(i,j) >= nsigmasECLAIRs * SQRT(Threshold(Sample_SVOM) / Delta_t_pflx) ) THEN
                        Prob_det_ij = 1.d0
                     ELSE
                        Prob_det_ij = 0.d0
                     END IF
                     Prob_det(Sample_SVOM) = Prob_det(Sample_SVOM) + TABomegaECLAIRs(i,j) * Prob_det_ij 
                  END DO
               END DO
               Prob_det(Sample_SVOM) = Prob_det(Sample_SVOM) / omegaECLAIRs

               IF (verbose == 2) PRINT*, "cts ECLAIRS =",cts_ECLAIRs
               IF (verbose == 2) PRINT*, "detec limit =", nsigmasECLAIRs * SQRT(Threshold(Sample_SVOM) / Delta_t_pflx)
               !IF(cts_ECLAIRs >= nsigmasECLAIRs * SQRT(Threshold(Sample_SVOM) / Delta_t_pflx) ) THEN
               !   Prob_det(Sample_SVOM) = 1.d0
               !ELSE
               !   Prob_det(Sample_SVOM) = 0.d0
               !END IF  
            CASE DEFAULT
               STOP "Error, undefined sample."
               
            END SELECT
            
            IF (ISNAN(Peakflux(i_Sample))) NaNtest_prop = .TRUE.
            
            IF (verbose==2) PRINT*, "Prob_det of ",TRIM(TabSample_name(i_Sample)) ," = ", Prob_det(i_Sample)       
         END IF
      ELSE
         Prob_det(i_Sample) = 1.d0
      END IF
   END DO
    
 END SUBROUTINE Calculate_Detection_Probability




 REAL(8) FUNCTION Nb_ph(L, z, Ep,  D_L, alpha, beta, ktild, Emin, Emax) 
   ! Returns the Number of photons [ph/cm^2/s] between Emin and Emax 
   ! L         [erg/s]
   ! Ep        [keV] (source frame)
   ! D_L       [Mpc]
   ! Emin,Emax [keV] (observer frame)
   REAL(8), INTENT(in) :: L, z, Ep, D_L, alpha, beta, ktild, Emin, Emax
   REAL(8)             :: xmin,xmax,x1, x2, B1, B2, B, dx, dlogx
   INTEGER             :: M2,j
   
   xmin = (1.d0+z)*Emin/Ep
   xmax = (1.d0+z)*Emax/Ep
   dlogx = 1.d-2
   M2 = MAX(1, INT(LOG10(xmax/xmin)/dlogx))
   IF(verbose == 2) PRINT*, " Number of steps for Nb_ph integration : ", M2
   x1 = xmin 
   B1 = Btild(x1, ktild, Model_Spec, alpha, beta)
   B = 0.d0
   DO j=1,M2
      x2 = 10.d0**(LOG10(xmin)+LOG10(xmax/xmin)*REAL(j,8)/REAL(M2,8))
      B2 = Btild(x2, ktild, Model_Spec, alpha, beta)
      B  = B + 0.5d0*(x2-x1)*(B2+B1)
      x1 = x2
      B1 = B2
   END DO
   IF(verbose == 2) PRINT*, " Btild integrated value : ", B
   Nb_ph = (1.d0+z) * L /(4.d0*Pi*(Ep*keV)*(D_L*Mpc)**2) * B
 END FUNCTION Nb_ph

 REAL(8) FUNCTION Nb_cts_ECLAIRs(L, z, Ep,  D_L, alpha, beta, ktild, Emin, Emax) 
   ! Returns the Number of photons [ph/cm^2/s] between Emin and Emax 
   ! L         [erg/s]
   ! Ep        [keV] (source frame)
   ! D_L       [Mpc]
   ! Emin,Emax [keV] (observer frame)
   REAL(8), INTENT(in) :: L, z, Ep, D_L, alpha, beta, ktild, Emin, Emax
   REAL(8)             :: xmin,xmax,x1, x2, B1, B2, Aeff1, Aeff2, B, dx, dlogx
   INTEGER             :: M2,j,i,ieff
   
   xmin = (1.d0+z)*Emin/Ep
   xmax = (1.d0+z)*Emax/Ep
   dlogx = 1.d-2
   M2 = MAX(1, INT(LOG10(xmax/xmin)/dlogx))
   IF(verbose == 2) PRINT*, " Number of steps for Nb_ph integration : ", M2
   x1 = xmin 
   B1 = Btild(x1, ktild, Model_Spec, alpha, beta)
   i = 1
   DO WHILE(TABeffECLAIRsE(i)<= Emin)
      ieff = i
      i = i + 1
      IF(i > NeffECLAIRs) THEN
         PRINT*,"Warning in Nb_ph_ECLAIRs: Emin=",Emin," larger than TABeffECLAIRsE(NeffECLAIRs)=",TABeffECLAIRsE(NeffECLAIRs)
         EXIT
      END IF
   END DO
   Aeff1 = TABeffECLAIRsA(ieff)
   B = 0.d0
   DO j=1,M2
      x2 = 10.d0**(LOG10(xmin)+LOG10(xmax/xmin)*REAL(j,8)/REAL(M2,8))
      B2 = Btild(x2, ktild, Model_Spec, alpha, beta)
      i = 1
      DO WHILE(TABeffECLAIRsE(i)<= Ep*x2/(1.d0+z) )
         ieff = i
         i = i + 1
         IF(i > NeffECLAIRs) THEN
            PRINT*,"Warning in Nb_ph_ECLAIRs: E2=",Ep*x2/(1.d0+z)," larger than TABeffECLAIRsE(NeffECLAIRs)=",TABeffECLAIRsE(NeffECLAIRs)
            EXIT         
         END IF
      END DO
      Aeff2 = TABeffECLAIRsA(ieff)
      B  = B + 0.5d0*(x2-x1)*(Aeff2*B2 + Aeff1*B1)
      x1 = x2
      B1 = B2
      Aeff1 = Aeff2
   END DO
   IF(verbose == 2) PRINT*, " Btild integrated value : ", B
   Nb_cts_ECLAIRs = (1.d0+z) * L /(4.d0*Pi*(Ep*keV)*(D_L*Mpc)**2) * B
 END FUNCTION Nb_cts_ECLAIRs


 REAL(8) FUNCTION F_12(L, z, Ep,  D_L, alpha, beta, ktild, Emin, Emax)
   ! Returns the Energy flux [erg/cm^2/s] between Emin and Emax 
   ! L         [erg/s]
   ! Ep        [keV] (source frame)
   ! D_L       [Mpc]
   ! Emin,Emax [keV] (observer frame) 
   REAL(8), INTENT(in) :: L, z, Ep, D_L, alpha, beta, ktild, Emin, Emax
   REAL(8)             :: xmin,xmax,x1, x2, B1, B2, B, dx, dlogx
   INTEGER             :: M2,j
   
   xmin = (1.d0+z)*Emin/Ep
   xmax = (1.d0+z)*Emax/Ep
   dlogx = 0.01d0
   M2 = MAX(1, INT(LOG10(xmax/xmin)/dlogx))
   x1 = xmin 
   B1 = xBtild(x1, ktild, Model_Spec, alpha, beta)
   B = 0.d0
   DO j=1,M2
      x2 = 10.d0**(LOG10(xmin)+LOG10(xmax/xmin)*REAL(j,8)/REAL(M2,8))
      B2 = xBtild(x2, ktild, Model_Spec, alpha, beta)
      B  = B + 0.5d0*(x2-x1)*(B2+B1)
      x1 = x2
      B1 = B2
   END DO
   IF(verbose==2) PRINT*, "B =", B
   F_12 = L /(4.d0*Pi*(D_L*Mpc)**2) * B
 END FUNCTION F_12
 
 
 REAL(8) FUNCTION Btild(x, ktild, Model_Spec, alpha, beta) ! No units
   ! Returns the unitless spectral shape. Used by Nb_ph  
   REAL(8), INTENT(in) :: x, ktild
   INTEGER, INTENT(in) :: Model_Spec
   REAL(8)             :: alpha, beta, x_c
   INTEGER             :: j, M
   
   IF(Model_Spec == Model_SpecBPLFix) THEN
      IF(x<1.d0) THEN
         Btild = ktild * x**(-alpha)
      ELSE IF(x>1.d0) THEN
         Btild = ktild * x**(-beta)
      ELSE
         Btild = ktild
         PRINT*, "In Btild, E = Ep"
      END IF
   ELSE IF (Model_Spec == Model_SpecBandFix .OR. Model_Spec == Model_SpecBandK .OR. Model_Spec == Model_SpecBandD .OR. Model_Spec == Model_SpecBandGBM) THEN
      x_c     = (beta - alpha) / (2.d0 - alpha)
      IF(x < x_c) THEN
         Btild = ktild * x**(-alpha)   * EXP( (alpha - 2.d0) * x )
      ELSE IF (x > x_c) THEN
         Btild = ktild * x**(-beta )   * EXP(alpha - beta) * x_c**(beta - alpha)
      ELSE
         Btild = ktild * x_c**(-alpha) * EXP(alpha - beta)
         PRINT*, "In Btild, x = x_c"
      END IF
   END IF
 END FUNCTION Btild

 ! Returns the unitless spectral shape. Used by F_12 (energy flux) 
 REAL(8) FUNCTION xBtild(x, ktild, Model_Spec, alpha, beta) ! No units
   REAL(8), INTENT(in) :: x, ktild
   INTEGER, INTENT(in) :: Model_Spec
   REAL(8)             :: alpha, beta, x_c
   INTEGER             :: j, M
   
   IF(Model_Spec == Model_SpecBPLFix) THEN
      IF(x<1.d0) THEN
         xBtild = ktild * x**(1.d0-alpha)
      ELSE IF(x>1.d0) THEN
         xBtild = ktild * x**(1.d0-beta)
      ELSE
         xBtild = ktild
         PRINT*, "In xBtild, E = Ep"
      END IF
   ELSE IF (Model_Spec == Model_SpecBandFix .OR. Model_Spec == Model_SpecBandK .OR. Model_Spec == Model_SpecBandD .OR. Model_Spec == Model_SpecBandGBM) THEN
      x_c     = (beta - alpha) / (2.d0 - alpha)
      IF(x < x_c) THEN
         xBtild = ktild * x**(1.d0-alpha)   * EXP( (alpha - 2.d0) * x )
      ELSE IF (x > x_c) THEN
         xBtild = ktild * x**(1.d0-beta )   * EXP(alpha - beta) * x_c**(beta - alpha)
      ELSE
         xBtild = ktild * x_c**(1.d0-alpha) * EXP(alpha - beta)
         PRINT*, "In xBtild, x = x_c"
      END IF
   END IF
   
 END FUNCTION xBtild
 
 ! -----------------------
 ! End of : Functions used
 
 SUBROUTINE Generate_Names()
   ! Generates TabSample_name, TabConstraint_name
   
   TabSample_name(Sample_Intrinsic)   = "intr_sample"
   TabSample_name(Sample_Kommers)     = "Kommers_sample"
   TabSample_name(Sample_Preece)      = "Preece_sample"
   TabSample_name(Sample_Stern)       = "Stern_sample"
   TabSample_name(Sample_SWIFTweak)   = "SWIFTweak_sample"
   TabSample_name(Sample_SWIFT)       = "SWIFT_sample"
   TabSample_name(Sample_SWIFTbright) = "SWIFTbright_sample"
   TabSample_name(Sample_HETE2)       = "HETE2_sample"
   TabSample_name(Sample_eBAT6)       = "eBAT6_sample"
   TabSample_name(Sample_EpGBM)       = "GBM_sample"
   TabSample_name(Sample_SVOM)        = "SVOM_sample"
   
   TabInstrument_name(Instrument_BATSE)   = "BATSE"
   TabInstrument_name(Instrument_BAT  )   = "BAT"
   TabInstrument_name(Instrument_FREGATE) = "FREGATE"
   TabInstrument_name(Instrument_WXM)     = "WXM"
   TabInstrument_name(Instrument_ECLAIRs) = "ECLAIRs"

   TabConstraint_name(Constraint_Kommers)  = "Kommers"
   TabConstraint_name(Constraint_Preece)   = "Preece"
   TabConstraint_name(Constraint_Stern)    = "Stern"
   TabConstraint_name(Constraint_HETE2)    = "XRFHETE2"
   TabConstraint_name(Constraint_EpGBM)    = "EpGBM"
   TabConstraint_name(Constraint_eBAT6)    = "eBAT6"
   
 END SUBROUTINE Generate_Names

 SUBROUTINE Generate_Paths()
   INTEGER :: i
   DO i=0,N_Samples
      LFile(i)      = TRIM(path) // "luminosity_"  // TRIM(TabSample_name(i))
      zFile(i)      = TRIM(path) // "redshift_"    // TRIM(TabSample_name(i))
      EpFile(i)     = TRIM(path) // "peakenergy_"  // TRIM(TabSample_name(i))
      PFile(i)      = TRIM(path) // "peakflux_"    // TRIM(TabSample_name(i))
      SpecFile_a(i) = TRIM(path) // "spec_alpha_"  // TRIM(TabSample_name(i))
      SpecFile_b(i) = TRIM(path) // "spec_beta_"   // TRIM(TabSample_name(i))

      zFile_cumul(i)      = TRIM(path) // "cumul_redshift_"    // TRIM(TabSample_name(i))
      
      LErrorFile(i)      = TRIM(path) // "luminosity_err_" // TRIM(TabSample_name(i))
      zErrorFile(i)      = TRIM(path) // "redshift_err_"   // TRIM(TabSample_name(i))
      EpErrorFile(i)     = TRIM(path) // "peakenergy_err_" // TRIM(TabSample_name(i))
      PErrorFile(i)      = TRIM(path) // "peakflux_err_"   // TRIM(TabSample_name(i))
      SpecErrorFile_a(i) = TRIM(path) // "spec_alpha_err_" // TRIM(TabSample_name(i))
      SpecErrorFile_b(i) = TRIM(path) // "spec_beta_err_"  // TRIM(TabSample_name(i)) 
 
   END DO

   KommFile     = TRIM(path) // KommFile
   PreeceFile   = TRIM(path) // PreeceFile
   SternFile    = TRIM(path) // SternFile
   XRFHETE2File = TRIM(path) // XRFHETE2File
   EpGBMFile    = TRIM(path) // EpGBMFile
   eBAT6File    = TRIM(path) // eBAT6File
   
   KommErrorFile     = TRIM(path) // KommErrorFile
   PreeceErrorFile   = TRIM(path) // PreeceErrorFile
   SternErrorFile    = TRIM(path) // SternErrorFile
   XRFHETE2ErrorFile = TRIM(path) // XRFHETE2ErrorFile
   EpGBMErrorFile    = TRIM(path) // EpGBMErrorFile
   eBAT6ErrorFile    = TRIM(path) // eBAT6ErrorFile
   
   eBAT6_EpLFile = TRIM(path) // eBAT6_EpLFile
   
 END SUBROUTINE Generate_Paths
 


 SUBROUTINE Prepare_Histogram_Prop(TabHistlim_inf, TabHistlim_sup, n_call)
   ! Generates the x-axis for luminosity, redshift, Ep, Peak flux histograms
   ! and (re)sets the values of the histograms to 0
   
   REAL(8), DIMENSION(1:N_prop), INTENT(in) :: TabHistlim_inf, TabHistlim_sup
   INTEGER, INTENT(inout)                   :: n_call
   ! Assumptions
   REAL(8) :: LogLmin,LogLmax
   REAL(8) :: zmax
   REAL(8) :: LogEpmin,LogEpmax
   REAL(8) :: LogPmin,LogPmax
   REAL(8) :: alpha_min,alpha_max
   REAL(8) :: beta_min,beta_max
   
   ! Local
   INTEGER :: i
   
   IF (n_call == 0 ) THEN

      IF(Lsave(0)) THEN
         ! Log(L) [erg/s]
         LogLmin = TabHistlim_inf(Prop_LogL)
         LogLmax = TabHistlim_sup(Prop_LogL)
         
         DO i=0, N_L
            TabLogL(i) = LogLmin + (LogLmax-LogLmin) * REAL(i,8)/REAL(N_L,8)
         END DO
      END IF

      IF(zsave(0)) THEN
         ! Redshift z
         zmax = TabHistlim_sup(Prop_z)
         
         DO i=0, N_z
            Tabz(i) = zmax * REAL(i,8)/REAL(N_z,8) 
         END DO
      END IF
      
      IF(Epsave(0)) THEN
         ! Ep [keV]
         LogEpmin = LOG10(TabHistlim_inf(Prop_Ep))
         LogEpmax = LOG10(TabHistlim_sup(Prop_Ep)) 
         
         DO i=0, N_Ep
            TabLogEp(i) = LogEpmin + (LogEpmax-LogEpmin) * REAL(i,8)/REAL(N_Ep,8) 
         END DO
      END IF

      IF(Psave(0)) THEN
         ! Peak Flux [ph/cm2/s between Emin and Emax]
         LogPmin = TabHistlim_inf(Prop_LogP)
         LogPmax = TabHistlim_sup(Prop_LogP)
         
         DO i=0, N_P
            TabLogP(i) = LogPmin + (LogPmax-LogPmin) * REAL(i,8)/REAL(N_P,8)
         END DO
      END IF


      IF(Specsave(0)) THEN
         ! alpha and beta
         alpha_min = TabHistlim_inf(Prop_alpha)
         alpha_max = TabHistlim_sup(Prop_alpha)
         
         DO i=0, N_spec_a
            TabSpec_a(i) = alpha_min + (alpha_max-alpha_min) * REAL(i,8)/REAL(N_spec_a,8)
         END DO
         
         beta_min = TabHistlim_inf(Prop_beta)
         beta_max = TabHistlim_sup(Prop_beta)
         
         DO i=0, N_spec_b
            TabSpec_b(i) = beta_min + (beta_max-beta_min) * REAL(i,8)/REAL(N_spec_b,8)
         END DO
      END IF
      
!!$      ! redshift for eBAT6
!!$      DO i=0, N_eBAT6 
!!$         TabeBAT6_z(i) =  10.d0 * REAL(i,8)/REAL(N_eBAT6,8) 
!!$      END DO

      ! Ep-L plane for eBAT6
      DO i=0, N_eBAT6_EpL ! equally spaced in logscale
         TabeBAT6_EpL(Indice_Ep, i) =  4.d0 * REAL(i,8)/REAL(N_eBAT6_EpL,8) 
         TabeBAT6_EpL(Indice_L, i)  =  48.d0 + (55.d0-48.d0) * REAL(i,8)/REAL(N_eBAT6_EpL,8)
      END DO
      TabeBAT6_EpL = 10**TabeBAT6_EpL ! remove log
      
      n_call = n_call + 1
   END IF
   
   ! --- Reset Histograms --- !

   IF(Lsave(0)) THEN
      TabHistLogL     = 0.d0
      TabHistLogL_master = 0.d0
   END IF

   IF(zsave(0)) THEN
      TabHistz        = 0.d0
      TabHistz_master = 0.d0   
   END IF

   IF(Epsave(0)) THEN
      TabHistLogEp    = 0.d0
      TabHistLogEpobs = 0.d0
      TabHistLogEp_master = 0.d0
      TabHistLogEpobs_master = 0.d0
   END IF

   IF(Psave(0)) THEN
      TabHistLogP     = 0.d0 
      TabHistLogP_master = 0.d0 
   END IF

   IF(Specsave(0)) THEN
      TabHistalpha     = 0.d0 
      TabHistalpha_master = 0.d0
      TabHistbeta      = 0.d0 
      TabHistbeta_master  = 0.d0 
   END IF

   ! Ep-L plane for eBAT6
   TabHisteBAT6_EpL = 0.d0
   TabHisteBAT6_EpL_master = 0.d0
      
 END SUBROUTINE Prepare_Histogram_Prop
 
 SUBROUTINE Reset_Histogram_Chi2()

   IF (Constraint_included(Constraint_Kommers)) THEN
      k_Kommers = 0.d0
      TabHistKomm_P23  = 0.d0
      TabHistKomm_P23_master  = 0.d0
   END IF

   IF (Constraint_included(Constraint_Preece)) THEN
      k_Preece = 0.d0
      TabHistPreece_Ep = 0.d0
      TabHistPreece_Ep_master = 0.d0
   END IF

   IF (Constraint_included(Constraint_Stern)) THEN
      k_Stern = 0.d0
      TabHistStern_P23 = 0.d0 
      TabHistStern_P23_master = 0.d0
   END IF
   
   IF (Constraint_included(Constraint_EpGBM)) THEN
      k_EpGBM = 0.d0
      TabHistEpGBM_Epobs = 0.d0
      TabHistEpGBM_Epobs_master = 0.d0
   END IF
   
   IF (Constraint_included(Constraint_eBAT6)) THEN
      k_eBAT6 = 0.d0
      TabHisteBAT6_z        = 0.d0
      TabHisteBAT6_z_master = 0.d0  
   END IF
   
    
 END SUBROUTINE Reset_Histogram_Chi2

 SUBROUTINE Reset_Histogram_Samples()
   
   IF (Sample_included(Sample_Kommers)) THEN
      k_Kommers = 0.d0
      TabHistKomm_P23  = 0.d0
      TabHistKomm_P23_master  = 0.d0
   END IF

   IF (Sample_included(Sample_Preece)) THEN
      k_Preece = 0.d0
      TabHistPreece_Ep = 0.d0
      TabHistPreece_Ep_master = 0.d0
   END IF

   IF (Sample_included(Sample_Stern)) THEN
      k_Stern = 0.d0
      TabHistStern_P23 = 0.d0 
      TabHistStern_P23_master = 0.d0
   END IF
   
   IF (Sample_included(Sample_EpGBM)) THEN
      k_EpGBM = 0.d0
      TabHistEpGBM_Epobs = 0.d0
      TabHistEpGBM_Epobs_master = 0.d0
   END IF
   
    IF (Sample_included(Sample_eBAT6)) THEN
      ! redshift for eBAT6
      TabHisteBAT6_z = 0.d0
      TabHisteBAT6_z_master = 0.d0
   END IF
   
 END SUBROUTINE Reset_Histogram_Samples

 SUBROUTINE Prepare_Luminosity()
   ! Creates distribution function for luminosity
   INTEGER :: i,j
   REAL(8) :: x1,y1,x2,y2,integ
   REAL(8) :: x_c, xBx, Gamma, gln, Gamma_min
   REAL(8) :: logLbreak, logLmin

   TabFctDistrLum(0) = 0.d0
   
   IF (Model_Lum == Model_LumSch) THEN
      logLmin = LOG10(TabParam_Lum(Param_Lum_Lmin))
      logLbreak = LOG10(TabParam_Lum(Param_Lum_Lbreak))
      TabLogL_CDF(0) = logLmin
      integ = 0.d0
      IF(run_mode == one_run) OPEN(UNIT=909, FILE=TRIM(path)//"probLum.dat")
      WRITE(909,*) 10**logLmin, Calc_Lum(logLmin, Model_Lum, TabParam_Lum) , TabFctDistrLum(0)

      x1 = 10.d0**(TabLogL_CDF(0))
      y1 = Calc_Lum(x1, Model_Lum, TabParam_Lum)
      DO i=1, M_L
         TabLogL_CDF(i) = logLmin + REAL(i,8)/REAL(M_L,8) * (logLbreak-logLmin + 1.5d0)  ! to avoid extremely low probability beyond logLbreak + 1.5
         x2 = 10.d0**(TabLogL_CDF(i))
         y2 = Calc_Lum(x2, Model_Lum, TabParam_Lum)

         integ = 0.5d0 * (x2 - x1) * (y1 + y2)
         x1 = x2
         y1 = y2
         !IF (verbose==2) PRINT*, "integ_Lum =", integ
         
         TabFctDistrLum(i) = TabFctDistrLum(i-1) + integ
         WRITE(909,*) x2, y2, TabFctDistrLum(i)
      END DO
      CLOSE(909)

   ELSE IF (Model_Lum > Model_LumSch) THEN ! CHECK THIS
      IF(rank == master_proc) PRINT*, ' WARNING THIS LUM MODEL HASNT BEEN CHECKED'
      DO i=1, N_L
         
         x1 = 10.d0**(TabLogL(i-1))
         y1 = Calc_Lum(x1, Model_Lum, TabParam_Lum)
         integ = 0.d0
         
         DO j=1, M_L
            x2 = 10.d0**(TabLogL(i-1)+(TabLogL(i)-TabLogL(i-1))*REAL(j,8)/REAL(M_L,8))
            y2 = Calc_Lum(x2, Model_Lum, TabParam_Lum)
            integ = integ + 0.5d0*(x2-x1)*(y2+y1)
            x1 = x2
            y1 = y2
            IF (verbose==2) PRINT*, "integ_Lum =", integ
         END DO
         
         TabFctDistrLum(i) = TabFctDistrLum(i-1) + integ
      END DO
      
   END IF

   IF (verbose==2) PRINT*, 'Final value of F_Lum :', TabFctDistrLum(M_L)
   TabFctDistrLum(0:M_L) = TabFctDistrLum(0:M_L) / TabFctDistrLum(M_L)

 END SUBROUTINE Prepare_Luminosity
 
 SUBROUTINE Prepare_Redshift()
   ! Creates distribution function for redshift
   INTEGER :: i
   REAL(8) :: x1,x2,y1,y2,integ,zmax
   
   
   SELECT CASE(Model_z)
   CASE(Model_zFix)
      zmax = 0.d0
   CASE(Model_zUniform)
      zmax = TabParam_z(Param_z_zmax)
   CASE(Model_zSH)
      zmax = TabParam_z(Param_z_zmax)
   CASE(Model_zDaigne)
      zmax = z_maximum
   CASE(Model_z_evol)
      zmax = TabParam_z(Param_z_zmax)
   CASE(Model_zLi)
      zmax = z_maximum
   CASE(Model_zPesc)
      zmax = z_maximum
   CASE(Model_zBPL)
      zmax = TabParam_z(Param_z_zmax)
   CASE(Model_zBExp)
      zmax = TabParam_z(Param_z_zmax)
   CASE DEFAULT
      STOP "Prepare_Redshift : z model not defined"
   END SELECT
   
   ! indice zmax
   
   IF ( zmax > Tabprecisez(INT(SIZE(Tabprecisez)-1)) ) STOP "Aborting zmax not allowed"
   IF (zmax /= 0.d0 ) THEN
      IF (Model_z /= Model_zPesc) THEN
         izmax = INT(REAL(SIZE(Tabprecisez)-1)*(zmax-Tabprecisez(0))/(Tabprecisez(INT(SIZE(Tabprecisez)-1))-Tabprecisez(0)))
         IF ( ABS(Tabprecisez(izmax)-zmax) > 1.d-4 ) THEN
            WRITE(*,*) "WARNING Moving zmax from ", zmax, " to ",Tabprecisez(izmax)
         END IF
         zmax = Tabprecisez(izmax)
         TabFctDistrz = 0.d0
         IF(run_mode == one_run) OPEN(UNIT=909, FILE=TRIM(path)//"probz.dat")
         IF(run_mode == one_run) OPEN(UNIT=908, FILE=TRIM(path)//"FctDistrz.dat")
         
         DO i=1, izmax
            x1 = Tabprecisez(i-1)
            y1 = Rz(x1, Model_z, Tabparam_z) * TabdVdz(i-1) / (1.d0 + x1)
            
            x2 = Tabprecisez(i)
            y2 = Rz(x2, Model_z, Tabparam_z) * TabdVdz( i ) / (1.d0 + x1)
            
            integ = 0.5d0 * (x2 - x1) * (y1 + y2)
            
            TabFctDistrz(i) = TabFctDistrz(i-1) + integ
            IF(run_mode == one_run)  WRITE(909, '(3ES12.5)') Tabprecisez(i), integ/(Tabprecisez(i)-Tabprecisez(i-1)), TabdVdz(i)/(1.d0 + x1)
            IF(run_mode == one_run)  WRITE(908, '(2ES12.5)') Tabprecisez(i), TabFctDistrz(i)
         END DO
         
         IF(run_mode == one_run) CLOSE(909)
         IF(run_mode == one_run) CLOSE(908)
         pseudo_collapse_rate = TabFctDistrz(izmax)
         !WRITE(*,'(A,ES12.5,A)') ' Pseudo collapse rate :', TabFctDistrz(izmax), ' collapse/yr'
         IF (verbose == 1) WRITE(*,*) 'Final value of F_z (Before normalization) :', TabFctDistrz(izmax)
         TabFctDistrz(0:izmax) = TabFctDistrz(0:izmax) / TabFctDistrz(izmax)
      ELSE
!!$         TabFctDistrz = 0.d0
!!$          OPEN(UNIT=55, FILE='../../catalogs/BAT6_cat/z_cumul_distr_Pesc16.txt')
!!$          READ(55,*) skip ! skip line
!!$          READ(55,*) skip ! skip line
!!$ 
!!$          DO i=0, N_
!!$             READ(55, *) TabGBM_alpha(i), TabFctDistrGBM_alpha(i)
!!$          END DO
!!$          CLOSE(55)
      END IF
      
        
      ELSE ! if zmax = 0
      
   END IF

 END SUBROUTINE Prepare_Redshift
 
 SUBROUTINE Prepare_Spec()
   ! Creates distribution function for alpha and beta
   INTEGER :: i
   CHARACTER :: skip

   TabFctDistrGBM_alpha = 0.d0
   OPEN(UNIT=55, FILE='../../catalogs/GBM_cat/alpha_GBM.txt')
   READ(55,*) skip ! skip line
   READ(55,*) skip ! skip line
   READ(55,*) skip ! skip line
   DO i=0, N_GBM_alpha
      READ(55, *) TabGBM_alpha(i), TabFctDistrGBM_alpha(i)
   END DO    
   CLOSE(55)
   
   TabFctDistrGBM_beta = 0.d0
   OPEN(UNIT=55, FILE='../../catalogs/GBM_cat/beta_GBM.txt')
   READ(55,*) skip ! skip line
   READ(55,*) skip ! skip line
   READ(55,*) skip ! skip line
   DO i=0, N_GBM_beta
      READ(55, *) TabGBM_beta(i), TabFctDistrGBM_beta(i)
   END DO    
   CLOSE(55)
   
 END SUBROUTINE Prepare_Spec

 SUBROUTINE Prepare_Constraints()
   INTEGER   :: i,j
   CHARACTER :: skip
   REAL(8)   :: Global_GRB_rate
   
   ! Kommers (2000) Table 2 : P23 [ph/cm2/s 50-300keV]
   OPEN(FILE='../observational_constraints/Kommers.dat', UNIT=600)
   READ(600, *) skip ! skips the first line of the file 
   j = 1
   DO i=1, 2*N_Komm
      IF (MOD(i,2) /= 0) THEN
         READ(600, '(4ES12.5)') TabKomm_P23(j-1), TabHistKomm_P23obs(j), TabKomm_DRobs(j), TabKomm_DRobserr(j) 
         j = j+1
      ELSE IF (i == 2*N_Komm) THEN
         READ(600,'(ES12.5)') TabKomm_P23(N_Komm)
      ELSE
         READ(600, *) skip 
      END IF
   END DO
   
   CLOSE(600)
   
   TabKomm_LogP23 = LOG10(TabKomm_P23) ! Logscale
   
   !PRINT*,'Global GRB rate from Kommers : ', SUM(TabKomm_DRobs)*4.d0*Pi, 'GRB/year'
   
   DO i=1, N_Komm
      TabKomm_DRDPobs(i)    = TabKomm_DRobs(i)   /(TabKomm_P23(i)-TabKomm_P23(i-1))
      TabKomm_DRDPobserr(i) = TabKomm_DRobserr(i)/(TabKomm_P23(i)-TabKomm_P23(i-1))
   END DO
   
   
   ! Preece (????) : Ep [keV]
   OPEN(FILE='../observational_constraints/preece.eb2.dat', UNIT=601)
   ! READ(601, *) skip ! skips the first line of the file 
   
   DO i=1, N_Preece
      READ(601, *) skip     
      READ(601, *) TabPreece_Ep(i-1), TabHistPreece_Epobs(i)
   END DO
   READ(601, *) TabPreece_Ep(N_Preece)
   !   TabPreece_Ep(N_Preece) = 6309.57d0
   ! DO i=1, N_Preece
   !    WRITE(*,*) " i, TabPreece_Ep(i-1), TabPreece_Ep(i), DN(i)", i, TabPreece_Ep(i-1), TabPreece_Ep(i), TabHistPreece_Epobs(i)
   ! END DO
   ! READ(*,*)
   CLOSE(601)
   
   TabPreece_LogEp = LOG10(TabPreece_Ep)
   TabHistPreece_Epobserr = SQRT(TabHistPreece_Epobs)
   TabHistPreece_Epobs    = TabHistPreece_Epobs      
   
   ! LogNLogP Stern
   !OPEN(FILE='../observational_constraints/lognlogp.stern.dat', UNIT=602)
   OPEN(FILE='../observational_constraints/Stern_lognlogp_rebinned.txt', UNIT=602)
   
   DO i=1, 10   ! skips the first 10 lines of the file 
      READ(602, *) skip
   END DO
   
   DO i=1, N_Stern
      READ(602,*) TabStern_P23(i-1), TabHistStern_P23obs(i), TabHistStern_P23obserr(i)
   END DO
   
   CLOSE(602)
   
   !This was used in old Stern constraint from Daigne+06
   !TabStern_P23(N_Stern) = TabStern_P23(N_Stern-1) + 0.1d0
   TabStern_P23(N_Stern) = 50.d0
   !TabStern_P23 = 10**(TabStern_P23)/0.75d0 ! divide by 0.75 to convert from counts to flux
   
   TabHistStern_P23obs = 10.d0**(TabHistStern_P23obs)
   ! Error propagation from log to linear scale (i.e. the average between the plus and minus errors)
   TabHistStern_P23obserr = TabHistStern_P23obs * ( (10.d0**(TabHistStern_P23obserr)-1.d0) )!+ (1.d0-10.d0**(-TabHistStern_P23obserr)) )
   
   Global_GRB_rate = 0.d0
   DO i=1, N_Stern
      Global_GRB_rate = Global_GRB_rate + LOG10(TabStern_P23(i)/ TabStern_P23(i-1)) * TabHistStern_P23obs(i)
   END DO
   !PRINT*, 'Global GRB rate from Stern+01 : ', Global_GRB_rate, ' GRB/year in 4 pi with pflx in [50-300 keV] above ',TabStern_P23(0),' ph/cm2/s'
   !This was used in old Stern constraint from Daigne+06
   !TabStern_P23 = 10.d0**(TabStern_P23)
   !TabStern_P23 = TabStern_P23 / 0.75d0
   
   
   ! GBM catalog (Gruber 2014) : Ep [keV]
   OPEN(FILE='../observational_constraints/Ep_GBM.txt', UNIT=603)
   READ(603, *) skip ! skips the first lines of the file
   READ(603, *) skip ! skips the first lines of the file
   READ(603, *) skip ! skips the first lines of the file 
   j = 1
   DO i=1, 2*N_EpGBM
      IF (MOD(i,2) /= 0) THEN
         READ(603, *) TabEpGBM_Epobs(j-1), TabHistEpGBM_Epobsobs(j), TabHistEpGBM_Epobsobserr(j)
         j = j+1
      ELSE IF (i == 2*N_EpGBM) THEN
         READ(603,*) TabEpGBM_Epobs(N_EpGBM)
      ELSE
         READ(603, *) skip 
      END IF
   END DO
   TabEpGBM_LogEpobs = LOG10(TabEpGBM_Epobs)
   CLOSE(603)
   
   ! eBAT6 redshift histogram (from Pescalli+16)
   OPEN(FILE='../observational_constraints/eBAT6_constraint.txt', UNIT=605)
   READ(605, *) skip ! skips the first lines of the file
   READ(605, *) skip ! skips the first lines of the file
   READ(605, *) skip ! skips the first lines of the file 
   j = 1
   DO i=1, 2*N_eBAT6
      IF (MOD(i,2) /= 0) THEN
         READ(605, *) TabeBAT6_z(j-1), TabHisteBAT6_zobs(j), TabHisteBAT6_zobserr(j)
         j = j+1
      ELSE IF (i == 2*N_eBAT6) THEN
         READ(605,*) TabeBAT6_z(N_eBAT6)
      ELSE
         READ(605, *) skip 
      END IF
   END DO
   CLOSE(605)


 END SUBROUTINE Prepare_Constraints
  
 SUBROUTINE Fill_Histogram_Prop()
   INTEGER              :: imin, imax, j
      
   ! --- Luminosity [erg/s] --- !
   iBinL = INT( ( LOG10(L)-TabHistlim_inf(Prop_LogL) ) / ( TabHistlim_sup(Prop_LogL)-TabHistlim_inf(Prop_LogL) ) * REAL(N_L,8) ) + 1
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBinL =", iBinL
   
   ! --- Redshift --- !
   iBinz = INT( z / TabHistlim_sup(Prop_z) * REAL(N_z,8) ) + 1
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBinz =", iBinz
   
   ! --- Peak Energy (source frame) [keV] --- !
   IF( (Ep < TabHistlim_inf(Prop_Ep)) .OR. (Ep >= TabHistlim_sup(Prop_Ep) )) THEN
      iBinEp = -1
   ELSE
      iBinEp = INT( ( LOG10(Ep)-LOG10(TabHistlim_inf(Prop_Ep)) ) / ( LOG10(TabHistlim_sup(Prop_Ep))-LOG10(TabHistlim_inf(Prop_Ep)) ) * REAL(N_Ep,8) ) + 1
   END IF
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBinEp =", iBinEp 
   
   ! --- Peak Energy (observer frame) [keV] --- !
   IF( (Epobs < TabHistlim_inf(Prop_Ep)) .OR. (Epobs >= TabHistlim_sup(Prop_Ep) )) THEN
      iBinEpobs = -1
   ELSE
      iBinEpobs = INT( ( LOG10(Epobs)-LOG10(TabHistlim_inf(Prop_Ep)) ) / ( LOG10(TabHistlim_sup(Prop_Ep))-LOG10(TabHistlim_inf(Prop_Ep)) ) * REAL(N_Ep,8) ) + 1
   END IF
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBinEpobs =", iBinEpobs 
   
   ! --- Peak Flux [ph/cm2/s] --- !
   DO i_Sample = 1, N_Samples
      IF(Sample_Included(i_Sample)) THEN
         IF( (LOG10(Peakflux(i_Sample)) < TabHistlim_inf(Prop_LogP)) .OR. (LOG10(Peakflux(i_Sample)) >= TabHistlim_sup(Prop_LogP) )) THEN
            iBinPeakflux(i_Sample) = -1
         ELSE
            iBinPeakflux(i_Sample) = INT( (LOG10(Peakflux(i_Sample))-TabHistlim_inf(Prop_LogP))/( TabHistlim_sup(Prop_LogP)-TabHistlim_inf(Prop_LogP) )*REAL(N_P,8) ) + 1
         END IF
         IF( (verbose == 1) .AND. (rank == master_proc) )WRITE(*,'(A,A,A,I3)') " iBinPeakflux of ",TRIM(TabSample_name(i_Sample))," = ", iBinPeakflux(i_Sample) 
      END IF
   END DO

   ! --- Spec alpha --- !
   IF( (alpha < TabHistlim_inf(Prop_alpha)) .OR. (alpha >= TabHistlim_sup(Prop_alpha) )) THEN
      iBina = -1
   ELSE
      iBina = INT( ( alpha-TabHistlim_inf(Prop_alpha) ) / ( TabHistlim_sup(Prop_alpha)-TabHistlim_inf(Prop_alpha) ) * REAL(N_spec_a,8) ) + 1
   END IF
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBina =", iBina
   
   ! --- Spec beta --- !
   IF( (beta < TabHistlim_inf(Prop_beta)) .OR. (beta >= TabHistlim_sup(Prop_beta) )) THEN
      iBinb = -1
   ELSE
      iBinb = INT( ( beta-TabHistlim_inf(Prop_beta) ) / ( TabHistlim_sup(Prop_beta)-TabHistlim_inf(Prop_beta) ) * REAL(N_spec_b,8) ) + 1
   END IF
   IF (verbose == 2) WRITE(*,'(A,I3)') " iBinb =", iBinb
   

   
   ! ----------- Fill the histograms ------------- !
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample)) THEN        
         IF((iBinL    > 0) &
   & .AND. (iBinL <= N_L)) TabHistLogL(i_Sample,         iBinL) = TabHistLogL(i_Sample,         iBinL) + Prob_det(i_Sample)
         IF(iBinz     > 0) TabHistz(i_Sample,            iBinz) = TabHistz(i_Sample,            iBinz) + Prob_det(i_Sample)
         IF(iBinEp    > 0) TabHistLogEp(i_Sample,       iBinEp) = TabHistLogEp(i_Sample,       iBinEp) + Prob_det(i_Sample)
         IF(iBinEpobs > 0) TabHistLogEpobs(i_Sample, iBinEpobs) = TabHistLogEpobs(i_Sample, iBinEpobs) + Prob_det(i_Sample)
         IF(iBina     > 0) TabHistalpha(i_Sample,        iBina) = TabHistalpha(i_Sample,        iBina) + Prob_det(i_Sample)
         IF(iBinb     > 0) TabHistbeta(i_Sample,         iBinb) = TabHistbeta(i_Sample,         iBinb) + Prob_det(i_Sample)
         IF (i_Sample >= 1) THEN
            IF(iBinPeakflux(i_Sample) > 0) TabHistLogP(i_Sample, iBinPeakflux(i_Sample)) = TabHistLogP(i_Sample, iBinPeakflux(i_Sample)) + Prob_det(i_Sample)
         END IF
      END IF
   END DO
 
 END SUBROUTINE Fill_Histogram_Prop

  
 SUBROUTINE Fill_Histogram_Chi2()
   INTEGER              :: imin, imax, j
   
   ! Comparison with Kommers et al. 2000
   IF(Constraint_Included(Constraint_Kommers)) THEN
      IF ( ( LOG10(Peakflux(Sample_Kommers)) < TabKomm_LogP23(0) ) .OR. ( LOG10(Peakflux(Sample_Kommers)) >= TabKomm_LogP23(N_Komm) ) ) THEN
         iBinKomm = -1
      ELSE
         imin = 1
         imax = N_Komm
         
         DO
            j = (imin+imax)/2
            IF ( LOG10(Peakflux(Sample_Kommers)) >= TabKomm_LogP23(j) ) THEN
               imin = j+1
            ELSE IF ( LOG10(Peakflux(Sample_Kommers)) < TabKomm_LogP23(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinKomm = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinKomm =", iBinKomm
   END IF
   
   ! Comparison with Preece
   IF(Constraint_Included(Constraint_Preece)) THEN
      IF ( ( Epobs < TabPreece_Ep(0) ) .OR. ( Epobs >= TabPreece_Ep(N_Preece) ) ) THEN
         iBinPreece = -1
      ELSE
         imin = 1
         imax = N_Preece
         
         DO
            j = (imin+imax)/2
            IF ( Epobs >= TabPreece_Ep(j) ) THEN
               imin = j+1
            ELSE IF ( Epobs < TabPreece_Ep(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinPreece = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinPreece =", iBinPreece
   END IF
   
   ! Comparison with Stern
   IF(Constraint_Included(Constraint_Stern)) THEN
      IF ( ( Peakflux(Sample_Stern) < TabStern_P23(0) ) .OR. ( Peakflux(Sample_Stern) >= TabStern_P23(N_Stern) ) ) THEN
         iBinStern = -1
      ELSE
         imin = 1
         imax = N_Stern
         
         DO
            j = (imin+imax)/2
            IF ( Peakflux(Sample_Stern) >= TabStern_P23(j) ) THEN
               imin = j+1
            ELSE IF ( Peakflux(Sample_Stern) < TabStern_P23(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinStern = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinStern =", iBinStern
   END IF
   
   ! Comparison with XRFHETE2
   IF(Constraint_Included(Constraint_HETE2)) THEN
      NGRB_HETE2 = NGRB_HETE2 + Prob_det(Sample_HETE2)
      softness =  F_12(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(Instrument_WXM), 30.d0) /&
           & F_12(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(Instrument_FREGATE), TabEmax(Instrument_FREGATE))
      IF (ISNAN(softness)) NaNtest_prop = .TRUE.
      !TabEmax(Instrument_WXM) remember to replace it
      ! PRINT*, "softness =", softness
      IF( softness > 1.d0) THEN
         NGRB_XRFHETE2 = NGRB_XRFHETE2 + Prob_det(Sample_HETE2)
      END IF
   END IF
   
   ! Comparison with GBM Ep distribution from Gruber et al. 2014 catalog
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      IF(  Epobs < TabEpGBM_Epobs(0)  .OR.  Epobs >= TabEpGBM_Epobs(N_EpGBM) ) THEN
         iBinEpGBM = -1
      ELSE
         imin = 1
         imax = N_EpGBM
         
         DO
            j = (imin+imax)/2
            IF ( Epobs >= TabEpGBM_Epobs(j) ) THEN
               imin = j+1
            ELSE IF ( Epobs < TabEpGBM_Epobs(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinEpGBM = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinEpGBM =", iBinEpGBM
   END IF
   
   ! Comparison with eBAT6 redshift distribution from Pescalli et al. 2016
   IF(Constraint_Included(Constraint_eBAT6)) THEN
      IF(  z < TabeBAT6_z(0)  .OR.  z >= TabeBAT6_z(N_eBAT6) ) THEN
         iBineBAT6 = -1
      ELSE
         iBineBAT6 = INT( z / TabeBAT6_z(N_eBAT6) * REAL(N_eBAT6,8) ) + 1
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBineBAT6 =", iBineBAT6
   END IF

   ! --- Simulated observations for comparison with real observations --- !
   ! --------------------------- Fill histograms ------------------------ !

   IF(Constraint_Included(Constraint_Kommers)) THEN
      IF (iBinKomm   > 0) TabHistKomm_P23(iBinKomm)    = TabHistKomm_P23(iBinKomm)    + Prob_det(Sample_Kommers)
   END IF
   IF(Constraint_Included(Constraint_Stern)) THEN
      IF (iBinStern  > 0) TabHistStern_P23(iBinStern)  = TabHistStern_P23(iBinStern)  + Prob_det(Sample_Stern)
   END IF
   IF(Constraint_Included(Constraint_Preece)) THEN
      IF (iBinPreece > 0) TabHistPreece_Ep(iBinPreece) = TabHistPreece_Ep(iBinPreece) + Prob_det(Sample_Preece)
   END IF
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      IF (iBinEpGBM > 0) TabHistEpGBM_Epobs(iBinEpGBM) = TabHistEpGBM_Epobs(iBinEpGBM) + Prob_det(Sample_EpGBM)
   END IF
   IF(Constraint_Included(Constraint_eBAT6)) THEN
      IF (iBineBAT6 > 0) TabHisteBAT6_z(iBineBAT6)     = TabHisteBAT6_z(iBineBAT6)     + Prob_det(Sample_eBAT6)
   END IF
 END SUBROUTINE Fill_Histogram_Chi2

 SUBROUTINE Fill_Samples()
   INTEGER              :: imin, imax, j
   
   ! Comparison with Kommers et al. 2000
   IF(Sample_Included(Sample_Kommers) .AND. .NOT.(Constraint_Included(Constraint_Kommers))) THEN
      IF ( ( LOG10(Peakflux(Sample_Kommers)) < TabKomm_LogP23(0) ) .OR. ( LOG10(Peakflux(Sample_Kommers)) >= TabKomm_LogP23(N_Komm) ) ) THEN
         iBinKomm = -1
      ELSE
         imin = 1
         imax = N_Komm
         
         DO
            j = (imin+imax)/2
            IF ( LOG10(Peakflux(Sample_Kommers)) >= TabKomm_LogP23(j) ) THEN
               imin = j+1
            ELSE IF ( LOG10(Peakflux(Sample_Kommers)) < TabKomm_LogP23(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinKomm = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinKomm =", iBinKomm
   END IF
   
   ! Comparison with Preece
   IF(Sample_Included(Sample_Preece) .AND. .NOT.(Constraint_Included(Constraint_Preece))) THEN
      IF ( ( Epobs < TabPreece_Ep(0) ) .OR. ( Epobs >= TabPreece_Ep(N_Preece) ) ) THEN
         iBinPreece = -1
      ELSE
         imin = 1
         imax = N_Preece
         
         DO
            j = (imin+imax)/2
            IF ( Epobs >= TabPreece_Ep(j) ) THEN
               imin = j+1
            ELSE IF ( Epobs < TabPreece_Ep(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinPreece = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinPreece =", iBinPreece
   END IF
   
   ! Comparison with Stern
   IF(Sample_Included(Sample_Stern) .AND. .NOT.(Constraint_Included(Constraint_Stern))) THEN
      IF ( ( Peakflux(Sample_Stern) < TabStern_P23(0) ) .OR. ( Peakflux(Sample_Stern) >= TabStern_P23(N_Stern) ) ) THEN
         iBinStern = -1
      ELSE
         imin = 1
         imax = N_Stern
         
         DO
            j = (imin+imax)/2
            IF ( Peakflux(Sample_Stern) >= TabStern_P23(j) ) THEN
               imin = j+1
            ELSE IF ( Peakflux(Sample_Stern) < TabStern_P23(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinStern = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinStern =", iBinStern
   END IF
   
   ! Comparison with XRFHETE2
   IF(Sample_Included(Sample_HETE2) .AND. .NOT.(Constraint_Included(Constraint_HETE2))) THEN
      NGRB_HETE2 = NGRB_HETE2 + Prob_det(Sample_HETE2)
      softness =  F_12(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(Instrument_WXM), 30.d0) /&
           & F_12(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(Instrument_FREGATE), TabEmax(Instrument_FREGATE))
      IF (ISNAN(softness)) NaNtest_prop = .TRUE.
      !TabEmax(Instrument_WXM) remember to replace it
      ! PRINT*, "softness =", softness
      IF( softness > 1.d0) THEN
         NGRB_XRFHETE2 = NGRB_XRFHETE2 + Prob_det(Sample_HETE2)
      END IF
   END IF
   
   ! Comparison with GBM Ep distribution from Gruber et al. 2014 catalog
   IF(Sample_Included(Sample_EpGBM) .AND. .NOT.(Constraint_Included(Constraint_EpGBM))) THEN
      IF(  Epobs < TabEpGBM_Epobs(0)  .OR.  Epobs >= TabEpGBM_Epobs(N_EpGBM) ) THEN
         iBinEpGBM = -1
      ELSE
         imin = 1
         imax = N_EpGBM
         
         DO
            j = (imin+imax)/2
            IF ( Epobs >= TabEpGBM_Epobs(j) ) THEN
               imin = j+1
            ELSE IF ( Epobs < TabEpGBM_Epobs(j-1) ) THEN
               imax = j
            ELSE
               EXIT
            END IF
         END DO
         
         iBinEpGBM = j
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBinEpGBM =", iBinEpGBM
   END IF
   
   ! Comparison with z distribution from Pescalli et al. 2016
   IF(Sample_Included(Sample_eBAT6) .AND. .NOT.(Constraint_Included(Constraint_eBAT6))) THEN
      IF(  z < TabeBAT6_z(0)  .OR.  z >= TabeBAT6_z(N_eBAT6) ) THEN
         iBineBAT6 = -1
      ELSE
         iBineBAT6 = INT( z / TabeBAT6_z(N_eBAT6) * REAL(N_eBAT6,8) ) + 1
      END IF
      IF (verbose == 2) WRITE(*,'(A,I3)') " iBineBAT6 =", iBineBAT6
   END IF
    

   ! --- Simulated observations for comparison with real observations --- !
   ! --------------------------- Fill histograms ------------------------ !

   IF(Sample_Included(Sample_Kommers) .AND. .NOT.(Constraint_Included(Constraint_Kommers))) THEN
      IF (iBinKomm   > 0) TabHistKomm_P23(iBinKomm)    = TabHistKomm_P23(iBinKomm)    + Prob_det(Sample_Kommers)
   END IF
   IF(Sample_Included(Sample_Stern) .AND. .NOT.(Constraint_Included(Constraint_Stern))) THEN
      IF (iBinStern  > 0) TabHistStern_P23(iBinStern)  = TabHistStern_P23(iBinStern)  + Prob_det(Sample_Stern)
   END IF
   IF(Sample_Included(Sample_Preece) .AND. .NOT.(Constraint_Included(Constraint_Preece))) THEN
      IF (iBinPreece > 0) TabHistPreece_Ep(iBinPreece) = TabHistPreece_Ep(iBinPreece) + Prob_det(Sample_Preece)
   END IF
   IF(Sample_Included(Sample_EpGBM) .AND. .NOT.(Constraint_Included(Constraint_EpGBM))) THEN
      IF (iBinEpGBM > 0) TabHistEpGBM_Epobs(iBinEpGBM) = TabHistEpGBM_Epobs(iBinEpGBM) + Prob_det(Sample_EpGBM)
   END IF
   IF(Sample_Included(Sample_eBAT6) .AND. .NOT.(Constraint_Included(Constraint_eBAT6))) THEN
      IF (iBineBAT6 > 0) TabHisteBAT6_z(iBineBAT6)     = TabHisteBAT6_z(iBineBAT6) + Prob_det(Sample_eBAT6)
   END IF
   
   CALL Fill_eBAT6()
   
 END SUBROUTINE Fill_Samples


 
 SUBROUTINE Fill_eBAT6()
   INTEGER              :: imin, imax, j
   ! Comparison with Ep-L plane for the eBAT6 sample from Pescalli et al. 2016
   IF(Sample_Included(Sample_eBAT6)) THEN
      IF(hist_flag == 2) THEN
         ! Ep
         IF(  Ep < TabeBAT6_EpL(Indice_Ep,0)  .OR.  Ep >= TabeBAT6_EpL(Indice_Ep,N_eBAT6_EpL) ) THEN
            iBineBAT6_Ep = -1
         ELSE
            iBineBAT6_Ep = INT( ( LOG10(Ep)-LOG10(TabeBAT6_EpL(Indice_Ep,0)) ) / ( LOG10(TabeBAT6_EpL(Indice_Ep,N_eBAT6_EpL))-LOG10(TabeBAT6_EpL(Indice_Ep,0)) ) * REAL(N_eBAT6_EpL,8) ) + 1
         END IF
         IF (verbose == 2) WRITE(*,'(A,I3)') " iBineBAT6_Ep =", iBineBAT6_Ep
         
         ! L
         IF(  L < TabeBAT6_EpL(Indice_L,0)  .OR.  L >= TabeBAT6_EpL(Indice_L,N_eBAT6_EpL) ) THEN
            iBineBAT6_L = -1
         ELSE
            iBineBAT6_L = INT( ( LOG10(L)- LOG10(TabeBAT6_EpL(Indice_L,0)) ) / ( LOG10(TabeBAT6_EpL(Indice_L,N_eBAT6_EpL))-LOG10(TabeBAT6_EpL(Indice_L,0)) ) * REAL(N_eBAT6_EpL,8) ) + 1
         END IF
         IF (verbose == 2) WRITE(*,'(A,I3)') " iBineBAT6_L =", iBineBAT6_L

         IF((iBineBAT6_L > 0) .AND. (iBineBAT6_Ep > 0)) TabHisteBAT6_EpL(iBineBAT6_Ep, iBineBAT6_L) = TabHisteBAT6_EpL(iBineBAT6_Ep, iBineBAT6_L) + Prob_det(Sample_eBAT6)
       
      END IF
   END IF
   

 END SUBROUTINE Fill_eBAT6

 
 SUBROUTINE Normalize_Model()
   ! Calculate normalization coefficients
   REAL(8)  :: x1,x2, x1_forlnL, x2_forlnL
   INTEGER  :: i

   !dof = 0 
   Chi2 = 0.d0
   
   ! Kommers et al. 2000
   IF(Sample_Included(Sample_Kommers)) THEN
      x1 = 0.d0
      x2 = 0.d0
      
      DO i=1, N_Komm
         ! Simulated rate and rate error
         TabHistKomm_DRDP(i)    =      TabHistKomm_P23_master(i)  / ( REAL(Nb_GRB,8) * ( TabKomm_P23(i)-TabKomm_P23(i-1) ) )
         TabHistKomm_DRDPerr(i) = SQRT(TabHistKomm_P23_master(i)) / ( REAL(Nb_GRB,8) * ( TabKomm_P23(i)-TabKomm_P23(i-1) ) )
         ! Note : error is SQRT(N) because of Poisson statistics
         
         x1 = x1 + TabHistKomm_DRDP(i) * TabKomm_DRDPobs(i) / TabKomm_DRDPobserr(i)**2
         x2 = x2 + TabHistKomm_DRDP(i)**2 / TabKomm_DRDPobserr(i)**2
      END DO
     
      IF(ISNAN(x1)) NaNtest_hist = .TRUE.
      IF(ISNAN(x2)) NaNtest_hist = .TRUE.
      k_Kommers = x1 / x2
      IF (ISNAN(k_Kommers)) NaNtest_hist = .TRUE.
      TabHistKomm_DRDP    = TabHistKomm_DRDP    * k_Kommers      ! Normalize
      TabHistKomm_DRDPerr = TabHistKomm_DRDPerr * k_Kommers
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Kommers data) :      k_Kommers =", k_Kommers/(4.d0*Pi), " yr-1                   ]"
   END IF
   
   ! Preece
   IF(Sample_Included(Sample_Preece)) THEN
      x1 = 0.d0
      x2 = 0.d0
      
      DO i=1, N_Preece 
         x1 = x1 + TabHistPreece_Ep_master(i) * TabHistPreece_Epobs(i) / (TabHistPreece_Epobserr(i)**2)
         x2 = x2 + TabHistPreece_Ep_master(i)**2 / (TabHistPreece_Epobserr(i)**2)
      END DO
      IF(x2 .NE. 0.d0) THEN
         k_Preece = x1 / x2
      ELSE
         k_Preece = 0.d0
         Chi2(Constraint_Preece) = 1.d20
      END IF
      IF (ISNAN(k_Preece)) NaNtest_hist = .TRUE.
      TabHistPreece_Ep_master = TabHistPreece_Ep_master * k_Preece      ! Normalize
      TabHistPreece_Eperr     = TabHistPreece_Eperr     * k_Preece
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Preece data)  :       k_Preece =", k_Preece, "                        ]" 
   END IF
   
   ! Stern
   IF(Sample_Included(Sample_Stern)) THEN
      x1 = 0.d0 
      x2 = 0.d0
      x1_forlnL = 0.d0 
      x2_forlnL = 0.d0
      ! Save a record of Stern hist from model as dN to use in likelihood
      TabHistStern_P23_forlnL = TabHistStern_P23_master
      
      DO i=1, N_Stern
         TabHistStern_P23err(i)     = SQRT(TabHistStern_P23_master(i)) / ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) )
         TabHistStern_P23_master(i) =      TabHistStern_P23_master(i)  / ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) )

         IF(TabHistStern_P23obserr(i) > 0.d0) THEN 
            x1 = x1 + (TabHistStern_P23obs(i) * TabHistStern_P23_master(i)) / (TabHistStern_P23obserr(i)**2)
            x2 = x2 + (TabHistStern_P23_master(i)  / TabHistStern_P23obserr(i) )**2
         END IF

         ! Transform Stern from dN/(dlogP * useful_time) to dN for use in likelihood
         TabHistStern_P23obs_forlnL(i) = TabHistStern_P23obs(i) * ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) ) * Delta_t_Stern * Omega_div_4Pi
         TabHistStern_P23obserr_forlnL(i) = TabHistStern_P23obserr(i) * ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) ) * Delta_t_Stern * Omega_div_4Pi
         
        
         IF(TabHistStern_P23obserr_forlnL(i) > 0.d0) THEN 
            x1_forlnL = x1_forlnL + (TabHistStern_P23obs_forlnL(i) * TabHistStern_P23_forlnL(i)) / (TabHistStern_P23obserr_forlnL(i)**2)
            x2_forlnL = x2_forlnL + (TabHistStern_P23_forlnL(i)  / TabHistStern_P23obserr_forlnL(i) )**2
         END IF

      END DO
      k_Stern = x1 / x2
      k_Stern_forlnL = x1_forlnL / x2_forlnL
      !k_Stern = 10.**(x1 / x2)
      !IF (run_mode == one_run) THEN
      IF (ISNAN(k_Stern)) NaNtest_hist = .TRUE.
      IF (ISNAN(k_Stern_forlnL)) NaNtest_hist = .TRUE.
      !END IF
      IF (ISNAN(k_Stern)) k_Stern = 1.d20
      TabHistStern_P23_master = TabHistStern_P23_master * k_Stern
      TabHistStern_P23err     = TabHistStern_P23err * k_Stern
      TabHistStern_P23_forlnL = TabHistStern_P23_forlnL * k_Stern_forlnL
      
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :        k_Stern =", k_Stern ,' yr-1                   ]'
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :   Sim duration =", 1.d0/k_Stern ,' yr                     ]'
     ! IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :          eta_0 =", k_Stern*Nb_GRB/pseudo_collapse_rate, "                        ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :   R_GRB/R_coll =", k_stern*Nb_GRB/collapse_rate_from_SFR , " yr                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :   n_GRB normalization =", k_stern*Nb_GRB / pseudo_collapse_rate, " GRB/yr          ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :  Pseudo collapse rate =", pseudo_collapse_rate, " collapse/yr     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :  Real GRB rate =", k_stern*Nb_GRB, " GRB/yr                 ]"
   END IF
   
   ! XRF HETE2
   IF(Sample_Included(Sample_HETE2)) THEN
      IF(NGRB_HETE2 > 0.d0) THEN
         Frac_XRFHETE2 = NGRB_XRFHETE2 / NGRB_HETE2
      ELSE
         Frac_XRFHETE2 = 0.d0
         Chi2(Constraint_HETE2) = 1.d20
      END IF
      !IF (run_mode == one_run) THEN
      IF (ISNAN(Frac_XRFHETE2)) NaNtest_hist = .TRUE.
      !END IF
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From HETE2                  :  Frac_XRFHETE2 =", Frac_XRFHETE2,'                        ]'   
   END IF

   ! GBM Ep catalog (from Gruber et al. 2014) 
   IF(Sample_Included(Sample_EpGBM)) THEN
      x1 = 0.d0
      x2 = 0.d0
           
      DO i=1, N_EpGBM
         TabHistEpGBM_Epobserr(i) = SQRT(TabHistEpGBM_Epobs_master(i))
         ! Note : error is SQRT(N) because of Poisson statistics
         x1 = x1 + TabHistEpGBM_Epobs_master(i) * TabHistEpGBM_Epobsobs(i) / TabHistEpGBM_Epobsobserr(i)**2
         x2 = x2 + TabHistEpGBM_Epobs_master(i)**2 / TabHistEpGBM_Epobsobserr(i)**2
      END DO
      !PRINT*, " x1, x2 = ",x1,x2
      IF(x2 > 0.d0) THEN
         k_EpGBM = x1 / x2
      ELSE
         !PRINT*, "GBM sample is empty..."
          k_EpGBM = 0.d0
      END IF
      IF (ISNAN(k_EpGBM)) NaNtest_hist = .TRUE.
      TabHistEpGBM_Epobs_master = TabHistEpGBM_Epobs_master * k_EpGBM      ! Normalize
      TabHistEpGBM_Epobserr     = TabHistEpGBM_Epobserr     * k_EpGBM
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From GBM (Gruber data)      :  normalization =", k_EpGBM, "                        ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From GBM (Gruber data)      : GBM efficiency =", k_EpGBM/(k_Stern*Delta_t_EpGBM), "                        ]"
   END IF

   
   ! eBAT6 redshift distribution (from Pescalli et al. 2016) 
   IF(Sample_Included(Sample_eBAT6)) THEN
      x1 = 0.d0
      x2 = 0.d0
           
      DO i=1, N_eBAT6
         ! Note : error is SQRT(N) because of Poisson statistics
         x1 = x1 + TabHisteBAT6_z_master(i) * TabHisteBAT6_zobs(i)
         x2 = x2 + TabHisteBAT6_z_master(i)**2
      END DO
      !PRINT*, " x1, x2 = ",x1,x2
      IF(x2 > 0.d0) THEN
         k_eBAT6 = x1 / x2
      ELSE
         k_eBAT6 = 0.d0
      END IF
      IF (ISNAN(k_eBAT6)) NaNtest_hist = .TRUE.
      TabHisteBAT6_zerr     = SQRT(TabHisteBAT6_z_master)
      TabHisteBAT6_z_master = TabHisteBAT6_z_master * k_eBAT6      ! Normalize
      TabHisteBAT6_zerr     = TabHisteBAT6_zerr     * k_eBAT6
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From eBAT6 (Pescalli data)  :    normalization =", k_eBAT6, "                      ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From eBAT6 (Pescalli data)  : eBAT6 efficiency =", k_eBAT6/(k_Stern*Delta_t_eBAT6), "                      ]"
   END IF
   
   IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
   !IF (run_mode == one_run) THEN
   IF (ISNAN(k_Stern)) NaNtest_hist = .TRUE.
   !END IF
   
   IF (NaNtest_hist) THEN
      WRITE(*,*) "[       k_Kommers      k_Stern      k_Preece      Frac_XRFHETE2      k_EpGBM       k_eBAT6          ]"
      WRITE(*,*) "[       ",   k_Kommers,      k_Stern,      k_Preece,       Frac_XRFHETE2, k_EpGBM, k_eBAT6, "         ]"
      WRITE(*,*) "[                                           TabHistKomm                                             ]"
      WRITE(*,*)  TabHistKomm_DRDP
      WRITE(*,*) "[                                           TabHistStern                                            ]"
      WRITE(*,*)  TabHistStern_P23
      WRITE(*,*) "[                                          TabHistPreece                                            ]"
      WRITE(*,*)  TabHistPreece_Ep
      WRITE(*,*) "[                                          TabHistEpGBM                                             ]"
      WRITE(*,*)  TabHistEpGBM_Epobs
      WRITE(*,*) "[                                          TabHisteBAT6                                             ]"
      WRITE(*,*)  TabHisteBAT6_z
   END IF
   
 END SUBROUTINE Normalize_Model

 SUBROUTINE Calculate_Chi2()
   INTEGER :: i
   
   IF(Constraint_Included(Constraint_Kommers)) THEN
      DO i=1, N_Komm
         Chi2(Constraint_Kommers) = Chi2(Constraint_Kommers) + ( (TabHistKomm_DRDP(i) - TabKomm_DRDPobs(i)) / TabKomm_DRDPobserr(i) )**2
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                   Chi2_Kommers = ", Chi2(Constraint_Kommers), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_Preece)) THEN
      DO i=1, N_Preece
         IF(TabHistPreece_Epobserr(i) > 0.d0) THEN         
            Chi2(Constraint_Preece) = Chi2(Constraint_Preece) + ( (TabHistPreece_Ep_master(i) - TabHistPreece_Epobs(i)) / TabHistPreece_Epobserr(i) )**2
           ! WRITE(*,*) "i, DN-DNobs, err, Chi2 :", i, (TabHistPreece_Ep(i) - TabHistPreece_Epobs(i)), TabHistPreece_Epobserr(i), Chi2(Constraint_Preece)
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                    Chi2_Preece = ", Chi2(Constraint_Preece), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_Stern)) THEN
      DO i=1, N_Stern
         IF(TabHistStern_P23obserr(i) > 0.d0) THEN 
            Chi2(Constraint_Stern) = Chi2(Constraint_Stern) + ( (TabHistStern_P23_master(i)-TabHistStern_P23obs(i)) / TabHistStern_P23obserr(i) )**2
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                     Chi2_Stern = ", Chi2(Constraint_Stern), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_HETE2)) THEN
      Chi2(Constraint_HETE2) = ( (Frac_XRFHETE2 - Frac_XRFHETE2obs) / sigma_XRFHETE2obs )**2
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                  Chi2_XRFHETE2 = ", Chi2(Constraint_HETE2), "                                     ]"  
   END IF
   
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      DO i=1, N_EpGBM
         IF(TabHistEpGBM_Epobsobserr(i) > 0.d0) THEN         
            Chi2(Constraint_EpGBM) = Chi2(Constraint_EpGBM) + ( (TabHistEpGBM_Epobs_master(i) - TabHistEpGBM_Epobsobs(i)) / TabHistEpGBM_Epobsobserr(i) )**2
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                     Chi2_EpGBM = ", Chi2(Constraint_EpGBM), "                                     ]"
   END IF

   IF(Constraint_Included(Constraint_eBAT6)) THEN
      DO i=1, N_eBAT6
         IF(TabHisteBAT6_zobserr(i) > 0.d0) THEN         
            Chi2(Constraint_eBAT6) = Chi2(Constraint_eBAT6) + ( (TabHisteBAT6_z_master(i) - TabHisteBAT6_zobs(i)) / TabHisteBAT6_zobserr(i) )**2
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[     (ignoring empty bins)           Chi2_eBAT6 = ", Chi2(Constraint_eBAT6), "                                     ]"
   END IF
   
   Chi2(0) = SUM(Chi2)
   
 END SUBROUTINE Calculate_Chi2

 SUBROUTINE Calculate_unnormalized_Likelihood()
   INTEGER :: i
 
   lnL = 0.d0
   
   IF(Constraint_Included(Constraint_Kommers)) THEN
      DO i=1, N_Komm
         lnL(Constraint_Kommers) = lnL(Constraint_Kommers) + (TabKomm_DRDPobs(i)+epsilon) * LOG(TabHistKomm_DRDP(i)+epsilon)&
              & - (TabHistKomm_DRDP(i)+epsilon)
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                    lnL_Kommers = ", lnL(Constraint_Kommers), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_Preece)) THEN
      DO i=1, N_Preece
         lnL(Constraint_Preece) = lnL(Constraint_Preece) + (TabHistPreece_Epobs(i)+epsilon) * LOG(TabHistPreece_Ep_master(i)+epsilon) - (TabHistPreece_Ep_master(i)+epsilon)
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                     lnL_Preece = ", lnL(Constraint_Preece), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_Stern)) THEN
      lnL_max_test = 0.d0
      lnL_empty_test = 0.d0
      DO i=1, N_Stern
         lnL(Constraint_Stern) = lnL(Constraint_Stern) + (TabHistStern_P23obs_forlnL(i)+epsilon) * LOG((TabHistStern_P23_forlnL(i)+epsilon)) &
              & - (TabHistStern_P23_forlnL(i)+epsilon)
         lnL_max_test(Constraint_Stern) = lnL_max_test(Constraint_Stern) + (TabHistStern_P23obs_forlnL(i)+epsilon) * LOG((TabHistStern_P23obs_forlnL(i)+epsilon)) &
              & - (TabHistStern_P23obs_forlnL(i)+epsilon)
         lnL_empty_test(Constraint_Stern) = lnL_empty_test(Constraint_Stern) + (TabHistStern_P23obs_forlnL(i)+epsilon) * LOG((epsilon)) &
              & - (epsilon)
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                      lnL_Stern = ", lnL(Constraint_Stern), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                  max lnL_Stern = ", lnL_max_test(Constraint_Stern), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                empty lnL_Stern = ", lnL_empty_test(Constraint_Stern), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_HETE2)) THEN  ! Not modified from Chi2, this has no meaning !!
      lnL(Constraint_HETE2) = ( (Frac_XRFHETE2 - Frac_XRFHETE2obs) / sigma_XRFHETE2obs )**2
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[          no meaning !!!!!         lnL_XRFHETE2 = ", lnL(Constraint_HETE2), "                                     ]"  
   END IF
   
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      DO i=1, N_EpGBM
         lnL(Constraint_EpGBM) = lnL(Constraint_EpGBM) + (TabHistEpGBM_Epobsobs(i)+epsilon) * LOG((TabHistEpGBM_Epobs_master(i)+epsilon))&
              & - (TabHistEpGBM_Epobs_master(i)+epsilon)
         lnL_max_test(Constraint_EpGBM) = lnL_max_test(Constraint_EpGBM) + (TabHistEpGBM_Epobsobs(i)+epsilon) * LOG((TabHistEpGBM_Epobsobs(i)+epsilon))&
              & - (TabHistEpGBM_Epobsobs(i)+epsilon)
         lnL_empty_test(Constraint_EpGBM) = lnL_empty_test(Constraint_EpGBM) + (TabHistEpGBM_Epobsobs(i)+epsilon) * LOG(epsilon)&
              & - (epsilon) 
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                      lnL_EpGBM = ", lnL(Constraint_EpGBM), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                  max lnL_EpGBM = ", lnL_max_test(Constraint_EpGBM), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                empty lnL_EpGBM = ", lnL_empty_test(Constraint_EpGBM), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_eBAT6)) THEN
      lnL_weight(Constraint_eBAT6) = 10.d0
      DO i=1, N_eBAT6
         lnL(Constraint_eBAT6) = lnL(Constraint_eBAT6) + (TabHisteBAT6_zobs(i)+epsilon) * LOG((TabHisteBAT6_z_master(i)+epsilon))&
              & - (TabHisteBAT6_z_master(i)+epsilon)
         lnL_max_test(Constraint_eBAT6) = lnL_max_test(Constraint_eBAT6) + (TabHisteBAT6_zobs(i)+epsilon) * LOG((TabHisteBAT6_zobs(i)+epsilon))&
              & - (TabHisteBAT6_zobs(i)+epsilon)
         lnL_empty_test(Constraint_eBAT6) = lnL_empty_test(Constraint_eBAT6) + (TabHisteBAT6_zobs(i)+epsilon) * LOG(epsilon)&
              & - (epsilon) 
      END DO

      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                      lnL_eBAT6 = ", lnL(Constraint_eBAT6) * lnL_weight(Constraint_eBAT6), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                  max lnL_eBAT6 = ", lnL_max_test(Constraint_eBAT6) * lnL_weight(Constraint_eBAT6), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                                empty lnL_eBAT6 = ", lnL_empty_test(Constraint_eBAT6) * lnL_weight(Constraint_eBAT6), "                                     ]"
   END IF
   
   lnL(0) = SUM(lnL*lnL_weight)
   lnL_max_test(0) = SUM(lnL_max_test*lnL_weight)
   lnL_empty_test(0) = SUM(lnL_empty_test*lnL_weight)
   
 END SUBROUTINE Calculate_unnormalized_Likelihood


 SUBROUTINE Save_Histograms()
   INTEGER :: i, i_Sample
   REAL(8) :: xtemp ! TEMPPPPPP
   ! --- Luminosity --- !

   
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (Lsave(i_Sample)) THEN
            OPEN( UNIT=100, FILE=TRIM(LFile(i_Sample))//".dat" )
            
            DO i=1, N_L
               ! Turn histogram into Lmin * p(L)
               IF (Model_Lum == Model_LumFix) THEN
                  TabHistLogL_master(i_Sample, i) = TabParam_Lum(Param_Lum_L0)   * TabHistLogL_master(i_Sample,i) / ( REAL(Nb_GRB,8) * ( 10**(TabLogL(i)) - 10**(TabLogL(i-1)) ) )
               ELSE 
                  TabHistLogL_master(i_Sample, i) = TabParam_Lum(Param_Lum_Lmin) * TabHistLogL_master(i_Sample,i) / ( REAL(Nb_GRB,8) * ( 10**(TabLogL(i)) - 10**(TabLogL(i-1)) ) )
               END IF
               !    1              2                                     
               ! [Log(L)]     [Lmin*p(L)]  
               WRITE(100, '(2ES12.5)') TabLogL(i-1), TabHistLogL_master(i_Sample, i)
               WRITE(100, '(2ES12.5)') TabLogL(i),   TabHistLogL_master(i_Sample, i)
            END DO
            
            CLOSE(100)
         END IF
      END IF
   END DO
   
   
   ! ----- Redshift ----- !
   
   
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (zsave(i_Sample)) THEN
            OPEN( UNIT=110, FILE=TRIM(zFile(i_Sample))//".dat" )
            OPEN( UNIT=111, FILE=TRIM(zFile_cumul(i_Sample))//".dat" )
            
            xtemp = 0.d0
            WRITE(111, '(2ES12.5)') Tabz(0), xtemp
            DO i=1, N_z
               xtemp = xtemp + TabHistz_master(i_Sample, i)
               WRITE(111, '(2ES12.5)') Tabz(i), xtemp
               ! Normalize
               TabHistz_master(i_Sample, i) = TabHistz_master(i_Sample, i) / ( REAL(Nb_GRB,8) * (Tabz(i) - Tabz(i-1)) )
               !  1           2                         
               ! [z]       [z pdf]     
               WRITE(110, '(2ES12.5)') Tabz(i-1), TabHistz_master(i_Sample, i)
               WRITE(110, '(2ES12.5)') Tabz(i),   TabHistz_master(i_Sample, i)
            END DO
            
            CLOSE(111)
            CLOSE(110)
         END IF
      END IF
   END DO
   
   
   ! ----- Peak Energy ----- !
   
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (Epsave(i_Sample)) THEN
            OPEN( UNIT=120, FILE=TRIM(EpFile(i_Sample))//".dat" )
            
            DO i=1, N_Ep
               ! Normalize
               TabHistLogEp_master(i_Sample,    i) = TabHistLogEp_master(i_Sample,    i) / ( REAL(Nb_GRB,8) * (TabLogEp(i) - TabLogEp(i-1)) )         
               TabHistLogEpobs_master(i_Sample, i) = TabHistLogEpobs_master(i_Sample, i) / ( REAL(Nb_GRB,8) * (TabLogEp(i) - TabLogEp(i-1)) )         
               !  1              2                    3          
               ! [Ep]     [Ep source pdf]        [Ep obs pdf]
               WRITE(120, '(3ES12.5)') TabLogEp(i-1), TabHistLogEp_master(i_Sample, i), TabHistLogEpobs_master(i_Sample, i)
               WRITE(120, '(3ES12.5)') TabLogEp(i),   TabHistLogEp_master(i_Sample, i), TabHistLogEpobs_master(i_Sample, i)
            END DO
            
            CLOSE(120)
         END IF
      END IF
   END DO
   
   ! ----- Peak Flux ----- !
   
   DO i_Sample = 1, N_Samples
      IF (Psave(i_Sample)) THEN
         IF(Sample_Included(i_Sample)) THEN
            OPEN( UNIT=120, FILE=TRIM(PFile(i_Sample))//".dat" )
            
            DO i=1, N_P
               ! Normalize
               TabHistLogP_master(i_Sample, i) = TabHistLogP_master(i_Sample, i) / ( REAL(Nb_GRB,8) * (10**TabLogP(i) - 10**TabLogP(i-1)) )        
               !  1           2                   
               ! [P]       [P pdf]        
               WRITE(120, '(2ES12.5)') TabLogP(i-1), TabHistLogP_master(i_Sample, i)
               WRITE(120, '(2ES12.5)') TabLogP(i),   TabHistLogP_master(i_Sample, i)
            END DO
            
            CLOSE(120)
         END IF
      END IF
   END DO
   
   ! --- Spec alpha --- !
   
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (Specsave(i_Sample)) THEN
            OPEN(UNIT=130, FILE=TRIM(SpecFile_a(i_Sample))//".dat" )
            WRITE(130, '(2ES12.5)') TabSpec_a(0), 0.d0
            DO i=1, N_spec_a
               TabHistalpha_master(i_Sample, i) =  TabHistalpha_master(i_Sample,i) / ( REAL(Nb_GRB,8) * ( TabSpec_a(i) - TabSpec_a(i-1) ) )
               
               !    1              2                                     
               ! [alpha]      [alpha pdf]  
               WRITE(130, '(2ES12.5)') TabSpec_a(i-1), TabHistalpha_master(i_Sample, i)
               WRITE(130, '(2ES12.5)') TabSpec_a(i),   TabHistalpha_master(i_Sample, i)
            END DO
            
            CLOSE(130)
         END IF
      END IF
   END DO
   
   ! --- Spec beta --- !
   
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (Specsave(i_Sample)) THEN
            OPEN(UNIT=140, FILE=TRIM(SpecFile_b(i_Sample))//".dat" )
            WRITE(140, '(2ES12.5)') TabSpec_b(0), 0.d0
            DO i=1, N_spec_b
               TabHistbeta_master(i_Sample, i) =  TabHistbeta_master(i_Sample,i) / ( REAL(Nb_GRB,8) * ( TabSpec_b(i) - TabSpec_b(i-1) ) )
               
               !    1              2                                     
               ! [beta]      [beta pdf]  
               WRITE(140, '(2ES12.5)') TabSpec_b(i-1), TabHistbeta_master(i_Sample, i)
               WRITE(140, '(2ES12.5)') TabSpec_b(i),   TabHistbeta_master(i_Sample, i)
            END DO
            
            CLOSE(140)
         END IF
      END IF
   END DO

 END SUBROUTINE Save_Histograms

 SUBROUTINE Save_eBAT6()
   CHARACTER(len=10) :: eBAT6_EpL_format = '(  ES12.5)'
   INTEGER :: i

   
   WRITE(eBAT6_EpL_format(2:3), '(I2)') N_eBAT6_EpL + 2
     
   IF(Sample_Included(Sample_eBAT6)) THEN
    
         ! redshift
!!$         DO i=2, N_eBAT6 ! start from 2 because ignore first bin
!!$            TabHisteBAT6_z_master(i) = TabHisteBAT6_z_master(i) + TabHisteBAT6_z_master(i-1)
!!$         END DO
!!$         OPEN( UNIT=309, FILE=TRIM(eBAT6File)//".dat" )
!!$
!!$         DO i=1, N_eBAT6
!!$            !  1           2                         
!!$            ! [z]       [number]   
!!$            WRITE(309, '(2ES12.5)') TabeBAT6_z(i-1), TabHisteBAT6_z_master(i)
!!$            WRITE(309, '(2ES12.5)') TabeBAT6_z(i),   TabHisteBAT6_z_master(i)
!!$         END DO
!!$         CLOSE(309)
         
         
         ! Ep-L plane
      OPEN( UNIT=310, FILE=TRIM(eBAT6_EpLFile)//".dat" )
      
      DO i=1, N_eBAT6_EpL
         !      1                   2                         3                
         ! [med_bin Ep]        [med_bin L]       [1D histogram(all_Ep, one_L)]     
         !WRITE(310, eBAT6_EpL_format)  TabeBAT6_EpL(Indice_Ep, i),TabeBAT6_EpL(Indice_L, i), TabHisteBAT6_EpL_master(:,i)
         WRITE(310, eBAT6_EpL_format)  0.5d0*(TabeBAT6_EpL(Indice_Ep, i) + TabeBAT6_EpL(Indice_Ep, i-1) ), 0.5d0*(TabeBAT6_EpL(Indice_L, i) + TabeBAT6_EpL(Indice_L, i-1)), TabHisteBAT6_EpL_master(:,i)
      END DO
      CLOSE(310)
      
   END IF

 END SUBROUTINE Save_eBAT6

 SUBROUTINE Save_Constraints()
   INTEGER :: i
     
   ! ---- P23 Kommers et al. 2000 ---- !
   
   IF (Constraint_save(Constraint_Kommers)) THEN
      OPEN(UNIT=200,  FILE=TRIM(KommFile)//".dat")
      OPEN(UNIT=2000, FILE=TRIM(KommErrorFile)//".dat")
           
      DO i=1, N_Komm
         !     1                   2                       3    
         ! [Log(P23)]   [DRDP hist BATSE23 model]    [DRDP Kommers hist]   
         WRITE(200, '(3ES12.5)') TabKomm_LogP23(i-1), TabHistKomm_DRDP(i), TabKomm_DRDPobs(i)
         WRITE(200, '(3ES12.5)') TabKomm_LogP23(i),   TabHistKomm_DRDP(i), TabKomm_DRDPobs(i)
         
         !     1                   2                        3                         4                           5     
         ! [Log(P23)]   [DRDP hist Kommers model]   [Kommers model err]   [DRDP hist Kommers obs]   [DRDP hist Kommers obs err]
         WRITE(2000, '(5ES12.5)') (TabKomm_LogP23(i-1)+TabKomm_LogP23(i))/2.d0,&
              &                   TabHistKomm_DRDP(i), TabHistKomm_DRDPerr(i), TabKomm_DRDPobs(i), TabKomm_DRDPobserr(i)
      END DO
      
      CLOSE(200)
      CLOSE(2000)
   END IF
   
   IF (Constraint_save(Constraint_Preece)) THEN
      OPEN(UNIT=201,  FILE=TRIM(PreeceFile)//".dat")
      OPEN(UNIT=2001, FILE=TRIM(PreeceErrorFile)//".dat")
      
      DO i=1, N_Preece
         
         !     1                   2                       3    
         ! [Log(Ep)]        [Ep Hist model]         [Ep Preece hist]   
         WRITE(201, '(3ES12.5)') TabPreece_LogEp(i-1), TabHistPreece_Ep_master(i), TabHistPreece_Epobs(i)
         WRITE(201, '(3ES12.5)') TabPreece_LogEp(i),   TabHistPreece_Ep_master(i), TabHistPreece_Epobs(i)
         
         !     1                   2                        3                         4                           5     
         !  [Log(Ep)]       [Ep hist model]        [Ep hist model err]       [Ep Preece hist obs]      [Ep Preece hist obs err]
         WRITE(2001, '(5ES12.5)') (TabPreece_LogEp(i-1)+TabPreece_LogEp(i))/2.d0,&
              &                   TabHistPreece_Ep_master(i), TabHistPreece_Eperr(i), TabHistPreece_Epobs(i), TabHistPreece_Epobserr(i)
      END DO
      
      CLOSE(201)
      CLOSE(2001)
   END IF
   
   IF (Constraint_save(Constraint_Stern)) THEN
      OPEN(UNIT=202,  FILE=TRIM(SternFile)//".dat")
      OPEN(UNIT=2002, FILE=TRIM(SternErrorFile)//".dat")
      
      DO i=1, N_Stern
         
         !     1                   2                       3    
         !   [P23]         [P23 Hist model]         [P23 Stern hist]   
         WRITE(202, '(3ES12.5)') TabStern_P23(i-1), TabHistStern_P23_master(i), TabHistStern_P23obs(i)
         WRITE(202, '(3ES12.5)') TabStern_P23(i),   TabHistStern_P23_master(i), TabHistStern_P23obs(i)
         
         !     1                   2                        3                         4                           5     
         !   [P23]         [P23 hist model]        [P23 hist model err]       [P23 Stern hist obs]      [P23 Stern hist obs err]
         WRITE(2002, '(5ES12.4)') SQRT(TabStern_P23(i-1)*TabStern_P23(i)),&
              &                    TabHistStern_P23_master(i), TabHistStern_P23err(i), TabHistStern_P23obs(i), TabHistStern_P23obserr(i)
      END DO
      
      CLOSE(202)
      CLOSE(2002)
   END IF
   
   IF (Constraint_save(Constraint_EpGBM)) THEN
      OPEN(UNIT=203,  FILE=TRIM(EpGBMFile)//".dat")
      OPEN(UNIT=2003, FILE=TRIM(EpGBMErrorFile)//".dat")
      
      DO i=1, N_EpGBM
         
         !     1                   2                       3    
         ! [Log(Ep)]       [EpGBM model hist]       [EpGBM obs hist]   
         WRITE(203, '(3ES12.5)') TabEpGBM_Epobs(i-1), TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobsobs(i)
         WRITE(203, '(3ES12.5)') TabEpGBM_Epobs(i),   TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobsobs(i)
         
         !     1                   2                        3                         4                           5     
         !  [Log(Ep)]       [Ep hist model]        [Ep hist model err]       [Ep EpGBM hist obs]      [Ep EpGBM hist obs err]
         WRITE(2003, '(5ES12.5)') (TabEpGBM_Epobs(i-1)+TabEpGBM_Epobs(i))/2.d0,&
              &                   TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobserr(i), TabHistEpGBM_Epobsobs(i), TabHistEpGBM_Epobsobserr(i)
      END DO
      
      CLOSE(203)
      CLOSE(2003)
   END IF
   
   IF (Constraint_save(Constraint_eBAT6)) THEN
      OPEN(UNIT=204,  FILE=TRIM(eBAT6File)//".dat")
      OPEN(UNIT=2004, FILE=TRIM(eBAT6ErrorFile)//".dat")
      
      DO i=1, N_eBAT6
         
         !     1                   2                       3    
         !    [z]          [eBAT6 model hist]       [eBAT6 obs hist]   
         WRITE(204, '(3ES12.5)') TabeBAT6_z(i-1), TabHisteBAT6_z_master(i), TabHisteBAT6_zobs(i)
         WRITE(204, '(3ES12.5)') TabeBAT6_z(i),   TabHisteBAT6_z_master(i), TabHisteBAT6_zobs(i)
         
         !     1                   2                        3                         4                           5     
         !    [z]            [z hist model]         [z hist model err]       [z eBAT6 hist obs]      [z eBAT6 hist obs err]
         WRITE(2004, '(5ES12.5)') (TabeBAT6_z(i-1)+TabeBAT6_z(i))/2.d0,&
              &                   TabHisteBAT6_z_master(i), TabHisteBAT6_zerr(i), TabHisteBAT6_zobs(i), TabHisteBAT6_zobserr(i)
      END DO
      
      CLOSE(204)
      CLOSE(2004)
   END IF

   
 END SUBROUTINE Save_Constraints

 SUBROUTINE Post_process_Constraints()
   INTEGER :: i
   CHARACTER(len=10) :: Kommers_format = '(  ES12.5)'
   CHARACTER(len=10) :: Stern_format   = '(  ES12.5)'
   CHARACTER(len=10) :: Preece_format  = '(  ES12.5)'
   CHARACTER(len=10) :: EpGBM_format   = '(  ES12.5)'
   CHARACTER(len=11) :: eBAT6_format   = '(   ES12.5)'
   
   IF (Constraint_save(Constraint_Stern)) THEN
      WRITE(Stern_format(2:3), '(I2)') N_Stern + N_Constraints + 1
      OPEN(UNIT=201,  FILE=TRIM(SternFile)//"_post_proc.dat", POSITION="append")
      WRITE(201, Stern_format) Chi2, TabHistStern_P23_master   
      CLOSE(201)
   END IF
   IF (Constraint_save(Constraint_Preece)) THEN
      WRITE(Preece_format(2:3), '(I2)') N_Preece + N_Constraints + 1
      OPEN(UNIT=202,  FILE=TRIM(PreeceFile)//"_post_proc.dat", POSITION="append")
      WRITE(202, Preece_format) Chi2, TabHistPreece_Ep_master   
      CLOSE(202)
   END IF
   IF (Constraint_save(Constraint_EpGBM)) THEN
      WRITE(EpGBM_format(2:3), '(I2)') N_EpGBM + N_Constraints + 1
      OPEN(UNIT=203,  FILE=TRIM(EpGBMFile)//"_post_proc.dat", POSITION="append")
      WRITE(203, EpGBM_format) Chi2, TabHistEpGBM_Epobs_master   
      CLOSE(203)
   END IF
   IF (Constraint_save(Constraint_eBAT6)) THEN
      WRITE(eBAT6_format(2:4), '(I3)') N_eBAT6 + N_Constraints + 1
      ! make the distribution cumulative
      DO i=2, N_eBAT6 ! start from 2 because ignore first bin
         TabHisteBAT6_z_master(i) = TabHisteBAT6_z_master(i) + TabHisteBAT6_z_master(i-1)
      END DO
      OPEN(UNIT=204,  FILE=TRIM(eBAT6File)//"_post_proc.dat", POSITION="append")
      WRITE(204, eBAT6_format) Chi2, TabHisteBAT6_z_master   
      CLOSE(204)
   END IF
 END SUBROUTINE Post_process_Constraints

 SUBROUTINE Reset_eBAT6_output_files()
   OPEN(UNIT=267, FILE=TRIM(path)//'eBAT6_KS_data_'//str_rank//'.dat')
   WRITE(267,'(A,I2)')"# This is the eBAT6 sample redshift distribution of processor ", rank
   ClOSE(267)
 END SUBROUTINE Reset_eBAT6_output_files
 
 SUBROUTINE WRITE_INFO()
   OPEN(UNIT=83, FILE=TRIM(path)//'info.txt')
   WRITE(83, '(A)') "# This file gathers the information about the run"
   WRITE(83, '(A)') "Directory : "//TRIM(path)
   WRITE(83, '(A,xI2)') "Number of cores :", nb_procs
   IF(RNG == Kiss_rng) THEN
      WRITE(83, '(A)') "RNG : KISS"
   ELSE IF (RNG == MT19937) THEN
      WRITE(83, '(A)') "RNG : MT19937"
   ELSE
      WRITE(83 ,'(A)') "RNG : ERROR"
   END IF
   IF(run_mode == one_run) THEN
      WRITE(83, '(A)') "Run mode : one run"
   ELSE IF(run_mode == param_search) THEN
      WRITE(83, '(A)') "Run mode : parameter search"
   ELSE IF(run_mode == MCMC) THEN
      WRITE(83, '(A)') "Run mode : MCMC"
   ELSE
      WRITE(83,'(A)') "Run mode : ERROR"
   END IF
   WRITE(83, '(A,xES12.5)') 'Number of GRBs :', REAL(Nb_GRB,8)
   WRITE(83, '(A)') 'Parameters explored :'
   IF(lum_explore == 1) WRITE(83, '(A)') " - Luminosity : "
   IF(z_explore == 1) WRITE(83, '(A)') " - Redshift : "
   IF(Ep_explore == 1) WRITE(83, '(A)') " - Ep : "
   IF(spec_explore == 1) WRITE(83, '(A)') " - Spectrum : "
   
   WRITE(83, '(A)') 'Constraints used :'
   
   DO i_Constraint=1, N_Constraints
      IF(Constraint_Included(i_Constraint)) WRITE(83,'(A)') " - "//TRIM(TabConstraint_name(i_Constraint))//" : "
   END DO
   
   WRITE(83, '(A,xF6.3)') 'dof : ', REAL(dof,8)
   WRITE(83, '(A,xF6.3)') '1_sigma : ', delta_chi2_1
   WRITE(83, '(A,xF6.3)') '2_sigma : ', delta_chi2_2
   WRITE(83, '(A,xF6.3)') '3_sigma : ', delta_chi2_3
   CLOSE(83)
   
 END SUBROUTINE WRITE_INFO

 INTEGER(4) FUNCTION Calc_starting_i()
   IF (reprise .EQV. .TRUE.) THEN
      IF (Nb_lines >= N_MCMC_iter) THEN
         starting_i = MOD(Nb_lines , N_MCMC_iter)
      ELSE
         starting_i = Nb_lines
      END IF
   ELSE
      starting_i = 1
   END IF
   Calc_starting_i = starting_i
   
 END FUNCTION Calc_starting_i
 
 
 SUBROUTINE Reprise_mode()
   IF (reprise .EQV. .TRUE.) THEN
      IF(run_mode == MCMC) THEN
         ! Lecture du fichier
         OPEN(UNIT=43, FILE=TRIM(path)//'reprise_MCMC.dat', FORM='unformatted') 
         Nb_lines = 0
         DO
            READ(43, err=996,end=997) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &           
        & Chi2, lnL, k_Kommers, k_Stern, k_Preece, k_EpGBM, k_eBAT6, dof, accepted_rec, tau, Step_Lum, Step_z, Step_Ep,&
        & TabHistLogL_master, TabHistz_master, TabHistLogEp_master, TabHistLogEpobs_master, TabHistLogP_master, &
        & TabHistKomm_DRDP, TabHistPreece_Ep_master, TabHistStern_P23_master, TabHistEpGBM_Epobs_master, TabHisteBAT6_z_master
            Nb_lines = Nb_lines + 1
            ! WRITE(*,*) "[                   reprise #", Nb_lines," : reduced Chi2 = ", Chi2(0)/REAL(dof,8),"              ]"
         END DO
997      IF(rank == master_proc) WRITE(*,*) "[                   Number of lines read in reprise_MCMC.dat : ",Nb_lines,"                              ]"
996      IF(rank == master_proc) WRITE(*,*) "[           [Error] Number of lines read in reprise_MCMC.dat : ",Nb_lines,"                              ]"
         CLOSE(43)
         
      ELSE IF(run_mode == param_search) THEN
         OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted') 
         
         ! Lecture du fichier
         Nb_lines = 0
         DO
            READ(43, err=998,end=999) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
                 & Chi2, k_Kommers, k_Stern, k_Preece, k_EpGBM, Frac_XRFHETE2,dof, chi2_min, Nb_good_models, Nb_models
            
            Nb_lines = Nb_lines + 1
            ! WRITE(*,*) "[                   reprise #", Nb_lines," : reduced Chi2 = ", Chi2(0)/REAL(dof,8),"              ]"
         END DO
999      IF(rank == master_proc) WRITE(*,*) "[                   Number of lines read in reprise.dat : ",Nb_lines,"                              ]"
998      IF(rank == master_proc) WRITE(*,*) "[           [Error] Number of lines read in reprise.dat : ",Nb_lines,"                              ]"
         CLOSE(43)           
      END IF
   END IF
 END SUBROUTINE Reprise_mode

 SUBROUTINE Define_Nb_lines
   ! This part is to find the number of lines in the file
   IF(rank == master_proc) THEN
      OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted')
      Nb_lines = 0
      DO
         READ(43,err=998,end=999) TabSave_ijk_Kiss, TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
              & Chi2, k_Kommers, k_Stern, k_Preece, k_EpGBM, Frac_XRFHETE2, dof, Chi2_min, Nb_good_models, Nb_models
         Nb_lines = Nb_lines + 1
      END DO
998   IF(rank == master_proc)  WRITE(*,*) "[           [Error] Number of lines read in reprise.dat : ",Nb_lines,"                              ]"
999   IF(rank == master_proc)  WRITE(*,*) "[                   Number of lines read in reprise.dat : ",Nb_lines,"                              ]"
      
      CLOSE(43)
      Nb_good_models_to_pp = Nb_good_models
   END IF
   CALL MPI_BCAST(Nb_lines, 1, MPI_INTEGER, master_proc, MPI_COMM_WORLD, code)
   CALL MPI_BCAST(Nb_good_models_to_pp, 1, MPI_INTEGER, master_proc, MPI_COMM_WORLD, code)
   
 END SUBROUTINE Define_Nb_lines

 REAL(8) FUNCTION Calc_kappa_tau(tau_0, epsilon_tau)
   REAL(8), intent(in)  :: tau_0
   REAL(8), intent(in)  :: epsilon_tau

   Calc_kappa_tau = ( epsilon / (tau_0 - 1.d0) )**(0.05d0)
   ! Calculated so that after 20 updates, tau = 1 + epsilon
   
 END FUNCTION Calc_kappa_tau


SUBROUTINE InitECLAIRs(nsigmasECLAIRs)
  REAL(8), intent(in) :: nsigmasECLAIRs
  REAL(8) :: x,y, E1obsECLAIRs, E2obsECLAIRs
  INTEGER :: i
  
  WRITE(*,'(A)')
  WRITE(*,'(A)')   "================================================"
  WRITE(*,'(A)')   "==                                            =="
  WRITE(*,'(A)')   "==          Init ECLAIRs instrument           =="
  WRITE(*,'(A)')   "==                                            =="
  WRITE(*,'(A)')   "================================================"
  WRITE(*,'(A)')
  E1obsECLAIRs = TabEmin(Instrument_ECLAIRs)
  E2obsECLAIRs = TabEmax(Instrument_ECLAIRs)
  WRITE(*,'(2(A,1ES12.5),A)') "ECLAIRs: use energy channel   : ",E1obsECLAIRs," to ",E2obsECLAIRs," keV"
  WRITE(*,'(1(A,1ES12.5),A)') "ECLAIRs: use detection level at ",nsigmasECLAIRs," sigmas"
  WRITE(*,'(A)')
  ExtEclairs=".    .    .    "
  IF (nsigmasECLAIRs<10.D0) THEN
    WRITE(ExtEclairs(2:5),  '(1I4.4)') INT(nsigmasECLAIRS*1000.D0) ! 5.5 -> 5500
  ELSE
    WRITE(ExtEclairs(2:5),  '(1I4.4)') INT(nsigmasECLAIRS*100.D0) ! 10. -> 1000
  END IF
  WRITE(ExtEclairs(7:10), '(1I4.4)') INT(E1obsECLAIRs*10.D0)     ! 4 -> 0040 ; 250 -> 2500
  WRITE(ExtEclairs(12:15),'(1I4.4)') INT(E2obsECLAIRs*10.D0)     ! 4 -> 0040 ; 250 -> 2500

  WRITE(*,'(A,A)')         "ECLAIRs: extension    = ",ExtECLAIRs
  WRITE(*,'(A)')

  ! Effective area

  OPEN(UNIT=100,FILE=TRIM(PathSVOM) // "rf_eff.txt")
  NeffECLAIRS = 0
  DO
    READ(100,*,err=999,end=999) x,y
    NeffECLAIRs = NeffECLAIRs+1
    TABeffECLAIRsE(NeffECLAIRs) = x
    TABeffECLAIRsA(NeffECLAIRs) = y
  END DO
999  CLOSE(100)
  WRITE(*,'(A,I4.4,A)') "ECLAIRs: effective area --> ",NeffECLAIRs," values"
  WRITE(*,'(2(A,1ES12.5),A)') "                       from ",TABeffECLAIRsE(1)," keV -- ",TABeffECLAIRsA(1)," cm2"
  WRITE(*,'(2(A,1ES12.5),A)') "                         to ",TABeffECLAIRsE(NeffECLAIRs)," keV -- ",TABeffECLAIRsA(NeffECLAIRs)," cm2"

  ! Background

  OPEN(UNIT=100,FILE=TRIM(PathSVOM) // "rf_bkg.txt")
  NbkgECLAIRS = 0
  DO
    READ(100,*,err=888,end=888) x,y
    NbkgECLAIRS = NbkgECLAIRS+1
    TABbkgECLAIRsE(NbkgECLAIRs) = x
    TABbkgECLAIRsB(NbkgECLAIRs) = y
  END DO
888  CLOSE(100)
  WRITE(*,'(A,I4.4,A)') "ECLAIRs: background     --> ",NbkgECLAIRs," values"
  WRITE(*,'(2(A,1ES12.5),A)') "                       from ",TABbkgECLAIRsE(1)," keV -- ",TABbkgECLAIRsB(1)," cts/s/keV"
  WRITE(*,'(2(A,1ES12.5),A)') "                         to ",TABbkgECLAIRsE(NbkgECLAIRs)," keV -- ",TABbkgECLAIRsB(NbkgECLAIRs)," cts/s/keV"

  IF (NbkgECLAIRs/=NeffECLAIRs) STOP "Error: different numbers of values for efficiency and background"
  DO i=1, NeffECLAIRs
    IF (TABbkgECLAIRsE(i)/=TABeffECLAIRsE(i)) STOP "Error: different values of energy for efficieny and background"
  END DO

  DO i=1, NeffECLAIRs
    IF (TABeffECLAIRsE(i)<=E1obsECLAIRs) istartECLAIRs=i
    IF (TABeffECLAIRsE(i)<=E2obsECLAIRs) iendECLAIRs=i
  END DO

  WRITE(*,'(2(A,1I4.4))')  "ECLAIRS: energy grid    --> ",istartECLAIRs," to ",iendECLAIRs

  bkgECLAIRsB1 = 0.D0
  DO i=istartECLAIRs+1, iendECLAIRs !2, NbkgECLAIRs
    bkgECLAIRsB1 = bkgECLAIRsB1 + 0.5D0*(TABbkgECLAIRsE(i)-TABbkgECLAIRsE(i-1))*(TABbkgECLAIRsB(i-1)+TABbkgECLAIRsB(i))
  END DO
  WRITE(*,'(A,1ES12.5,A)') "ECLAIRs: background     --> ",bkgECLAIRsB1," cts/s"

  ! offaxis correction

  OPEN(UNIT=100,FILE=TRIM(PathSVOM) // "rf_offaxis.txt")
  DO i=-NoffECLAIRs,NoffECLAIRs
    READ(100,*,err=777) TABoffECLAIRs(i,-NoffECLAIRs:NoffECLAIRs)
!    WRITE(*,*) i,"OK"
  END DO
  CLOSE(100)
  WRITE(*,'(2(A,1ES12.5))') "ECLAIRs: offaxis correction from ",MINVAL(TABoffECLAIRs)," to ",MAXVAL(TABoffECLAIRs)

  ! solid angle associated to each "pixel"

  OPEN(UNIT=100,FILE=TRIM(PathSVOM) // "rf_omega.txt")
  DO i=-NoffECLAIRs,NoffECLAIRs
    READ(100,*,err=777) TABomegaECLAIRs(i,-NoffECLAIRs:NoffECLAIRs)
!    WRITE(*,*) i,"OK"
  END DO
  CLOSE(100)
  WRITE(*,'(2(A,1ES12.5),A)') "ECLAIRs: omega(pixel) = ",MINVAL(TABomegaECLAIRs)," to ",MAXVAL(TABomegaECLAIRs)," sr"
  
  omegaECLAIRs=SUM(TabomegaECLAIRs)
  WRITE(*,'(A,1ES12.5,A)') "ECLAIRs: omega(total) = ",omegaECLAIRs," sr"


  RETURN
777 STOP "Error when reading rf_offaxis.txt"



END SUBROUTINE InitECLAIRs



END PROGRAM Core

