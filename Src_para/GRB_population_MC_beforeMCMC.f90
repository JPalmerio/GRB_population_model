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
  INTEGER                  :: run_mode = 0                                    ! Run mode : 0 (default) is one run mode, 1 is param search
  INTEGER, PARAMETER       :: one_run = 0                                     ! Value for one run mode
  INTEGER, PARAMETER       :: param_search = 1                                ! Value for param search mode
  INTEGER                  :: hist_flag = 0                                    ! Flag to tell MonteCarlo routine to generate histograms or not
  CHARACTER(*), PARAMETER  :: InitFile = '../Input_para/GRB_pop.init'         ! Name of input file
  CHARACTER(Len=255)       :: path                                            ! Path for output
  INTEGER,      PARAMETER  :: verbose    = 0                                  ! Verbosity of the code (0, 1 or 2)
  LOGICAL,      PARAMETER  :: Save_all_GRB = .FALSE.                          ! Writes properties for every GRB in GRB_prop
  LOGICAL                  :: SHOW_READ  = .FALSE.                            ! Verbosity of the ReadInitFile subroutine
  CHARACTER(*), PARAMETER  :: Format_RIF = '(A,A13,A22,A14,ES12.5,A)'         ! Format for ReadInitFile verbosity
  LOGICAL                  :: NaNtest_prop = .FALSE.                          ! To test if any properties are NaN  !! only works during one_run mode
  LOGICAL                  :: NaNtest_hist = .FALSE.                          ! To test if any histograms are NaN  !! only works during one_run mode
  INTEGER                  :: Nb_lines = 0                                    ! Number of lines in reprise file
  INTEGER                  :: Nb_GRB                                          ! Number of GRBs 
  INTEGER                  :: n_call = 0                                      ! If n_call == 0, will generate histograms, else will just reset them to 0
 
  ! ---------- RNG variables ---------- !
  INTEGER, PARAMETER                     :: N_proc_max = 8                       ! Maximum number of processors (8 by default)
  INTEGER, DIMENSION(1:3,0:N_proc_max-1) :: TabInit_ijk_Kiss                     ! Table for initatlizing KISS generator
  INTEGER, DIMENSION(1:3,0:N_proc_max-1) :: TabSave_ijk_Kiss                     ! Table for saving KISS generator
  INTEGER, DIMENSION(0:N_proc_max-1)     :: TabInit_seed_MT                      ! Table for initializing MT19937 generator
  INTEGER                                :: iKiss, jKiss, kKiss                  ! v.a. uniforme : KISS generator
  INTEGER                                :: iKiss_save, jKiss_save, kKiss_save   ! to save v.a.
  INTEGER                                :: MT_seed                              ! MT19937 seed
  INTEGER                                :: Kiss_rng = 0, MT19937 = 1            ! Different cases for RNG
  INTEGER                                :: RNG = 999                            ! Random number generator (Kiss_rng or MT)

  ! ----------- MPI variables ---------- !
  INTEGER :: nb_procs, rank, code, name_length                                ! Variables used to store number of procs, rank of each proc etc..
  INTEGER, PARAMETER :: master_proc = 0                                       ! Rank of the master processor that does the Chi2 calculations


  ! --------- GRB samples and constraints --------- !
  ! ----------------------------------------------- !

  ! --- Samples --- !
  INTEGER, PARAMETER                        :: N_Samples          = 9         ! Number of samples
  INTEGER                                   :: i_Sample                       ! indice for sample loops
  INTEGER, PARAMETER                        :: Sample_Intrinsic   = 0         ! indice for intrinsic sample    (=all GRBs)
  INTEGER, PARAMETER                        :: Sample_Kommers     = 1         ! indice for Kommers sample      (50-300 keV & trigger efficiency from Kommers et al. 2000)
  INTEGER, PARAMETER                        :: Sample_Preece      = 2         ! indice for Preece sample       (50-300 keV & Peak flux > 5     ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_Stern       = 3         ! indice for Stern sample        (50-300 keV)
  INTEGER, PARAMETER                        :: Sample_SWIFTweak   = 4         ! indice for SWIFT sample        (15-150 keV & Peak flux > 0.001 ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_SWIFT       = 5         ! indice for SWIFT sample        (15-150 keV & Peak flux > 0.2   ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_SWIFTbright = 6         ! indice for SWIFT sample        (15-150 keV & Peak flux > 1     ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_HETE2       = 7         ! indice for HETE2 sample        (2-10 keV & 30-400 keV & Peak flux > 1 ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_BAT6ext     = 8         ! indice for BAT6ext sample      (15-150 keV & Peak flux >= 2.6  ph/cm2/s)
  INTEGER, PARAMETER                        :: Sample_GBM         = 9         ! indice for GBM sample          (50-300 keV & Peak flux >= 0.9  ph/cm2/s)
  ! INTEGER, PARAMETER                        :: Sample_HETE2FRE    = 6         ! indice for FREGATE sample      (30-400 keV & Peak flux > 1     ph/cm2/s)
  ! INTEGER, PARAMETER                        :: Sample_HETE2WXM    = 7         ! indice for WXM sample          (2-10   keV & Peak flux > 1     ph/cm2/s)
  ! INTEGER, PARAMETER                        :: Sample_SWIFTBand   = 8         ! formule de band a faire plus tard
  LOGICAL(8), DIMENSION(0:N_Samples)        :: Sample_Included = .FALSE.      ! Table to include sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: TabSample_name                 ! Table for output filename of each sample 
  REAL(8), DIMENSION(0:N_Samples)           :: Threshold                      ! Table for threshold for each sample
  REAL(8), DIMENSION(0:N_Samples)           :: Prob_det                       ! Table for detection probability of instrument for each sample
   
  ! --- Instruments --- !
  INTEGER, PARAMETER                            :: N_Instruments       = 4         ! Number of instruments
  INTEGER                                       :: i_Instrument                    ! indice for instrment loops
  CHARACTER(Len=80), DIMENSION(0:N_Instruments) :: TabInstrument_name              ! Table for output filename of each instrument 
  INTEGER, PARAMETER                            :: Instrument_BATSE    = 1         ! indice for BATSE   instrument (50-300 keV)
  INTEGER, PARAMETER                            :: Instrument_SWIFT    = 2         ! indice for SWIFT   instrument (15-150 keV)
  INTEGER, PARAMETER                            :: Instrument_FREGATE  = 3         ! indice for FREGATE instrument (30-400 keV)
  INTEGER, PARAMETER                            :: Instrument_WXM      = 4         ! indice for WXM     instrument ( 2-10  keV)
  LOGICAL(8), DIMENSION(0:N_Instruments)        :: Instrument_Included = .FALSE.   ! Table to include sample
  REAL(8), DIMENSION(1:N_Instruments)           :: TabEmin                         ! Table for lower energy limit of instrument 
  REAL(8), DIMENSION(1:N_Instruments)           :: TabEmax                         ! Table for higher energy limit of instrument 
  REAL(8), DIMENSION(1:N_Instruments)           :: Peakflux_Instrument             ! Table for Peak fluxes for each instrument 

 

  ! --- Constraints --- !
  INTEGER, PARAMETER                     :: N_Constraints       = 5           ! Number of Constraints
  INTEGER                                :: i_Constraint                      ! indice for Constraint loops
  INTEGER, PARAMETER                     :: Constraint_Kommers  = 1           ! indice for Kommers Constraint (LogN LogP  Kommers et al. 2000, Table 2) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_Preece   = 2           ! indice for Preece  Constraint (LogN LogEp Preece  et al. ????, Table ?) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_Stern    = 3           ! indice for Stern   Constraint (LogN LogP  Stern   et al. 2001, Table ?) Requires BATSE
  INTEGER, PARAMETER                     :: Constraint_XRFHETE2 = 4           ! indice for X-ray Flash fraction constraint (Frac_XRFHETE2, sigma_XRFHETE2) Requires HETE2FRE and HETE2WXM
  INTEGER, PARAMETER                     :: Constraint_EpGBM    = 5           ! indice for Ep constraint (from GBM catalog, Gruber at al. 2014) Requires BATSE
  LOGICAL(8), DIMENSION(0:N_Constraints) :: Constraint_Included = .FALSE.     ! Table to include constraint
  LOGICAL(8), DIMENSION(0:N_Constraints) :: Constraint_save     = .FALSE.     ! Table to save constraint
  CHARACTER(Len=80), DIMENSION(0:N_Constraints) :: TabConstraint_name         ! Table for output filename of each constraint

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
  REAL(8), PARAMETER           :: z_maximum      = 20.d0                      ! maximum redshift



  ! ---------------------------------------------------------- !
  ! --- Source properties / GRB properties in source frame --- !
  ! ---------------------------------------------------------- !
  
  ! --------- Luminosity [erg/s] --------- !
  ! -------------------------------------- !

  ! Assumptions

  INTEGER                   :: Model_Lum        = 0                           ! Luminosity Model indice
  INTEGER, PARAMETER        :: Model_LumFix     = 0                           ! Fixed luminosity. Parameters : L0
  INTEGER, PARAMETER        :: Model_LumPL      = 1                           ! Power-law.        Parameters : Lmin, Lmax, slope
  INTEGER, PARAMETER        :: Model_LumBPL     = 2                           ! Broken power-law. Parameters : Lmin, Lmax, Lbreak, slopeL, slopeH

  INTEGER, PARAMETER               :: NParam_Lum       = 10                   ! Maximum number of parameters for luminosity model
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum     = -999.d0              ! Table of Parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum_min = -999.d0              ! Table of minima for parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabParam_Lum_max = -999.d0              ! Table of maxima for parameters for the Luminosity Models
  REAL(8), DIMENSION(1:NParam_Lum) :: TabBestParam_Lum = -999.d0              ! Table for the best parameters of the Luminosity Models
  INTEGER                          :: lum_explore = 0                         ! Flag for exploration of parameter space of Luminosity model
  
  INTEGER, PARAMETER        :: Param_Lum_L0     = 1                           ! index for L0
  INTEGER, PARAMETER        :: Param_Lum_Lmin   = 2                           ! index for Lmin
  INTEGER, PARAMETER        :: Param_Lum_Lmax   = 3                           ! index for Lmax
  INTEGER, PARAMETER        :: Param_Lum_Lbreak = 4                           ! index for Lbreak
  INTEGER, PARAMETER        :: Param_Lum_slope  = 5                           ! index for slope
  INTEGER, PARAMETER        :: Param_Lum_slopeL = 6                           ! index for slope(low lum)
  INTEGER, PARAMETER        :: Param_Lum_slopeH = 7                           ! index for slope(high lum)

  ! Histogram

  INTEGER, PARAMETER                    :: N_L = 50                           ! Number of bins in TabLogL
  REAL(8), DIMENSION(0:N_L)             :: TabLogL                            ! Table equally spaced in Log(L) 
  REAL(8), DIMENSION(0:N_Samples,1:N_L) :: TabHistLogL = 0.d0                 ! Table for Histogram of L
  REAL(8), DIMENSION(0:N_Samples,1:N_L) :: TabHistLogL_master = 0.d0          ! Table for master Histogram of L
  INTEGER                               :: iBinL                              ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Lsave = .FALSE.                ! Saves (TRUE) Luminosity histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: LFile                          ! Outputfile name for Luminosity for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: LErrorFile                     ! Outputfile name for Luminosity Error for each sample


  ! Random draw

  REAL(8), DIMENSION(0:N_L) :: TabFctDistrLum                                 ! Table for Distribution Function of L
  INTEGER                   :: M_L = 1000                                     ! Precision of Distribution Function integration (reference=1000)(only used if not analytical)


  ! ----------- Redshift ----------- !
  ! -------------------------------- !

  ! Assumptions

  INTEGER                   :: Model_z         = 0                            ! Redshift Model indice
  INTEGER, PARAMETER        :: Model_zFix      = 0                            ! Fixed redshift.       Parameters : z0
  INTEGER, PARAMETER        :: Model_zUniform  = 1                            ! Uniform rate.         Parameters : zmax, R0 
  INTEGER, PARAMETER        :: Model_zSH       = 2                            ! Springel & Hernquist. Parameters : zmax, zm, a, b, R0 
  INTEGER, PARAMETER        :: Model_zDaigne   = 3                            ! Daigne et al. 2006.   Parameters : a, b, c, d, k 

  INTEGER, PARAMETER             :: NParam_z       = 10                       ! Maximum number of parameters for redshift model  
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z     = -999.d0                  ! Table of Parameters for the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z_min = -999.d0                  ! Table of minima for parameters for the Redshift Models
  REAL(8), DIMENSION(1:NParam_z) :: TabParam_z_max = -999.d0                  ! Table of maxima for parameters for the Redshift Models
  INTEGER                        :: redshift_explore = 0                      ! Flag for exploration of parameter space of redshift
  
  INTEGER, PARAMETER        :: Param_z_z0      = 1                            ! index for z0
  INTEGER, PARAMETER        :: Param_z_zmax    = 2                            ! index for zmax
  INTEGER, PARAMETER        :: Param_z_R0      = 3                            ! index for R0
  INTEGER, PARAMETER        :: Param_z_zm      = 4                            ! index for zm (S&H)
  INTEGER, PARAMETER        :: Param_z_a       = 5                            ! index for a  (S&H, Daigne)
  INTEGER, PARAMETER        :: Param_z_b       = 6                            ! index for b  (S&H, Daigne)
  INTEGER, PARAMETER        :: Param_z_c       = 7                            ! index for c  (Daigne)
  INTEGER, PARAMETER        :: Param_z_d       = 8                            ! index for d  (Daigne)

  ! Histogram 

  INTEGER, PARAMETER                    :: N_z = 100                          ! Number of bins in Tabz
  REAL(8), DIMENSION(0:N_z)             :: Tabz                               ! Table equally spaced in z
  REAL(8), DIMENSION(0:N_Samples,1:N_z) :: TabHistz = 0.d0                    ! Table for Histogram of z
  REAL(8), DIMENSION(0:N_Samples,1:N_z) :: TabHistz_master = 0.d0             ! Table for master Histogram of z
  INTEGER                               :: iBinz                              ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: zsave = .FALSE.                ! Saves (TRUE) redshift histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: zFile                          ! Outputfile name for Redshift for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: zFile_cumul                    ! Outputfile name for cumulative Redshift distribution for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: zErrorFile                     ! Outputfile name for Redshift Error for each sample


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
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: SpecFile_a                     ! Outputfile name for alpha for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: SpecFile_b                     ! Outputfile name for beta for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: SpecErrorFile_a                ! Outputfile name for error on alpha for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: SpecErrorFile_b                ! Outputfile name for error on beta for each sample
  



  ! --------- Peak Energy [keV] --------- !
  ! ------------------------------------- !

  ! Assumptions

  INTEGER                    :: Model_Ep          = 0                         ! Peak Energy Model indice
  INTEGER, PARAMETER         :: Model_EpFix       = 0                         ! Fixed Peak Energy (source frame).       Parameters : Ep0
  INTEGER, PARAMETER         :: Model_EpLogNormal = 1                         ! Log-normal distribution (source frame). Parameters : Ep0, sigmaLog
  INTEGER, PARAMETER         :: Model_EpY         = 2                         ! Yonetoku relation TO BE IMPLEMENTED
  INTEGER, PARAMETER         :: Model_EpAmati     = 3                         ! Amati-like relation (source frame)      Parameters : Ep0, sigmaLog, alpha_amati, (Ep = Ep0 * (L/L0)**alpha_amati

  INTEGER, PARAMETER              :: NParam_Ep    = 6                         ! Maximum number of parameters for model of Ep
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep  = -999.d0                   ! Table of Parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep_min  = -999.d0               ! Table of minima for parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabParam_Ep_max  = -999.d0               ! Table of maxima for parameters for the Peak Energy Models
  REAL(8), DIMENSION(1:NParam_Ep) :: TabBestParam_Ep  = -999.d0               ! Table for the best parameters of the Luminosity Models
  INTEGER                         :: Ep_explore = 0                           ! Flag for exploration of parameter space of Peak energy
  
  INTEGER, PARAMETER         :: Param_Ep_Ep0         = 1                      ! index for Ep0
  INTEGER, PARAMETER         :: Param_Ep_sigmaLog    = 2                      ! index for sigmaLog
  INTEGER, PARAMETER         :: Param_Ep_L0          = 3                      ! index for L0 (for Amati)
  INTEGER, PARAMETER         :: Param_Ep_alpha_amati = 4                      ! index for alpha_amati 


  ! Histogram 

  INTEGER, PARAMETER                     :: N_Ep = 50                         ! Number of bins in TabLogEp
  REAL(8), DIMENSION(0:N_Ep)             :: TabLogEp                          ! Table equally spaced in Log(Ep)
  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEp = 0.d0               ! Table for Histogram of Log(Ep)
  REAL(8), DIMENSION(0:N_Samples,1:N_Ep) :: TabHistLogEp_master = 0.d0        ! Table for master Histogram of Log(Ep)
  INTEGER                                :: iBinEp                            ! Indice for the value drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Epsave = .FALSE.               ! Saves (TRUE) Peak Energy histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: EpFile                         ! Outputfile name for Peak Energy for each sample
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: EpErrorFile                    ! Outputfile name for Peak Energy Error for each sample


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

  INTEGER, PARAMETER                    :: N_P = 50                           ! Number of bins in TabP
  REAL(8), DIMENSION(0:N_P)             :: TabLogP                            ! Table equally spaced in Log(P) 
  REAL(8), DIMENSION(1:N_Samples,1:N_P) :: TabHistLogP = 0.d0                 ! Table for Histogram of peak flux P for each sample
  REAL(8), DIMENSION(1:N_Samples,1:N_P) :: TabHistLogP_master = 0.d0          ! Table for Histogram of peak flux P for each sample
  INTEGER, DIMENSION(1:N_Samples)       :: iBinPeakflux                       ! Indices for values drawn

  LOGICAL,           DIMENSION(0:N_Samples) :: Psave = .FALSE.                ! Saves (TRUE) Peak Flux histograms for each sample (Intrinsic, BATSE23...)
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: PFile                          ! Outputfile name for Peak Flux for each sample  
  CHARACTER(Len=80), DIMENSION(0:N_Samples) :: PErrorFile                     ! Outputfile name for Peak Flux Error for each sample


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
  REAL(8)                      :: R_GRB=0.d0                                      ! GRB rate (normalization coefficient) 
  INTEGER                      :: iBinKomm                                        ! Index for histogram

  CHARACTER(Len=80)            :: KommFile       = "Kommers_constraint.dat"       ! Outputfile name for Kommers
  CHARACTER(Len=80)            :: KommErrorFile  = "Kommers_constrainterr.dat"    ! Outputfile name for Kommers Error

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

  CHARACTER(Len=80)              :: PreeceFile       = "Preece_constraint.dat"                    ! Outputfile name for generated Ep histogram
  CHARACTER(Len=80)              :: PreeceErrorFile  = "Preece_constrainterr.dat"                 ! Outputfile name for Error on generated Ep histogram 

  ! Observational data

  REAL(8), DIMENSION(0:N_Preece) :: TabPreece_Ep, TabPreece_LogEp                            ! Table from Preece
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Epobs    = 0.d0                            ! Table for observed BATSE23 Histogram of Ep (data from Preece)
  REAL(8), DIMENSION(1:N_Preece) :: TabHistPreece_Epobserr = 0.d0                            ! Table for observed BATSE23 Histogram error of Ep (data from Preece)

  ! --- BATSE23 : LogN LogP distribution Stern --- !
  ! ---------------------------------------------- !
  
  ! Histogram

  INTEGER, PARAMETER            :: N_Stern                 = 30                              ! Number of bins in TabStern_P23
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23        = 0.d0                            ! Table for Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23_master = 0.d0                            ! Table for Histogram of Peak flux
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23err     = 0.d0                            ! Table for Histogram of Peak flux error
  REAL(8)                       :: k_Stern = 0.d0                                            ! Normalization coefficient
  INTEGER                       :: iBinStern                                                 ! Index for histogram

  CHARACTER(Len=80)              :: SternFile       = "Stern_constraint.dat"                      ! Outputfile name for generated P23 
  CHARACTER(Len=80)              :: SternErrorFile  = "Stern_constrainterr.dat"                   ! Outputfile name for Error on generated P23 

 ! Observational data

  REAL(8), DIMENSION(0:N_Stern) :: TabStern_P23, TabStern_LogP23                             ! Table from Stern
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obs    = 0.d0                             ! Table for LogN LogP from BATSE23 (data from Stern)
  REAL(8), DIMENSION(1:N_Stern) :: TabHistStern_P23obserr = 0.d0                             ! Table for LogN LogP error from BATSE23 (data from Stern)


  ! --- HETE2 : X-ray Flash fraction --- !
  ! ------------------------------------ !
  
  ! Fraction
  
  REAL(8) :: Frac_XRFHETE2
  REAL(8) :: sigma_XRFHETE2
  REAL(8) :: NGRB_XRFHETE2=0.d0, NGRB_HETE2=0.d0
  CHARACTER(Len=80) :: XRFHETE2File      = "XRFHETE2_constraint.dat"                              ! Outputfile name for generated XRF fraction 
  CHARACTER(Len=80) :: XRFHETE2ErrorFile = "XRFHETE2_constrainterr.dat"                           ! Outputfile name for error on generated XRF fraction
  
 ! Observational data

  REAL(8), PARAMETER :: Frac_XRFHETE2obs  = 0.35d0                                           ! X-ray flash fraction detected by HETE2
  REAL(8), PARAMETER :: sigma_XRFHETE2obs = 0.15d0                                           ! Error on X-ray flash fraction detected by HETE2


  ! --- GBM : Ep distribution --- !
  ! ----------------------------- !
  
  INTEGER, PARAMETER            :: N_EpGBM                   = 9                     ! Number of bins in TabEpGBM
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobs        = 0.d0                   ! Table for Histogram of Peak energy
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobs_master = 0.d0                   ! Table for Master Histogram of Peak energy
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobserr     = 0.d0                   ! Table for Histogram of error on peak energy
  REAL(8)                       :: norm_EpGBM                = 0.d0                   ! Normalization coefficient
  INTEGER                       :: iBinEpGBM                                          ! Index for histogram

  CHARACTER(Len=80)             :: EpGBMFile       = "EpGBM_constraint.dat"           ! Outputfile name for EpGBM
  CHARACTER(Len=80)             :: EpGBMErrorFile  = "EpGBM_constrainterr.dat"        ! Outputfile name for EpGBM Error

  ! Observational data

  REAL(8), DIMENSION(0:N_EpGBM) :: TabEpGBM_Epobs, TabEpGBM_LogEpobs                  ! Extracted histogram from GBM catalog
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobsobs    = 0.d0                    ! Table for observed GBM Histogram of Epobs (data from Gruber et al 2014 catalog)
  REAL(8), DIMENSION(1:N_EpGBM) :: TabHistEpGBM_Epobsobserr = 0.d0                    ! Table for observed error of GBM Histogram of Epobs (data from Gruber et al 2014 catalog)



  ! -------- Chi2 -------- !
  ! ---------------------- !

  REAL(8), DIMENSION(0:N_Constraints) :: Chi2                                                ! Table for Chi2 adjusment of the model for each constraint
  INTEGER :: dof                                                                             ! Numbers of degrees of freedom
  REAL(8) :: delta_chi2_1                                                                    ! Chi2 interval in which "good models" are included at 1 sigma (68.3%)
  REAL(8) :: delta_chi2_2                                                                    ! Chi2 interval in which "good models" are included at 2 sigma (95.5%)
  REAL(8) :: delta_chi2_3                                                                    ! Chi2 interval in which "good models" are included at 3 sigma (99.7%)
  REAL(8) :: Chi2_min = 1.d21                                                                ! Lowest Chi2 found
  INTEGER :: good_model=0, Nb_good_models=0, Nb_models=0                                     ! Number of models and good models (defined above)

  REAL(8) :: t
  INTEGER :: ii !temp
  REAL(8) :: Collapse_rate = 0.d0
    


  ! ----- Variables for the GRB draws ----- !

  ! Intrinsic
  REAL(8) :: L,z,D_L,Ep,alpha,beta              ! variables used for each GRB draw
  REAL(8) :: ktild                              ! Normalization of spectrum
  REAL(8) :: t_star                             ! Luminosity BPL draw
  
  ! Assumptions
  REAL(8) :: z0,zmax,a,b,zm,R0
  REAL(8) :: Ep0, sigmaLog
  
  ! Observed
  REAL(8) :: softness
  REAL(8) :: Epobs
  
  
   

  ! ------ MCMC stuff ----- !

  INTEGER :: Niter = 1000000, MCMC_run_nb
  REAL(8) :: old_slope, old_Lmin, old_Lmax
  REAL(8) :: slope_step, Lmin_step, slope_start, Lmin_start, Lmax_start, Lmax_step
  
  ! --------------------------------- !
  ! ------- End of declaration ------ !
  ! --------------------------------- !


  ! ----- Initialization ----- !

  CALL MPI_INIT(code)                                     ! MPI environment starts here
  !CALL PTIM_start(user_label)                             ! times stuff
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)          ! gives the rank of the processor
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)      ! gives the number of processors being used
  
  IF(rank == master_proc) SHOW_READ = .TRUE.         ! Only master proc has verbosity (otherwise unreadable)
  CALL InitPhysics(.FALSE.)                          ! Initializes Omega_M, Omega_Lambda etc... and filters
  CALL InitCosmology(SHOW_READ)                      ! Creates precise tables for z, dVdz, D_L and Vz
  CALL Generate_Names()                              ! Creates TabSample_name
  CALL ReadInitFile(InitFile)                        ! Reads the input file and sets the models according to the user's choice
  IF(run_mode == param_search) CALL ReadParamSearchFile()   ! Reads the input file for parameter search mode
  IF(run_mode == one_run) CALL Calculate_dof()
  CALL Generate_Paths()                              ! Creates the paths for the save files
  CALL Init_RNG()                                    ! Initializes random number generator  
  CALL InitBoundaries()                              ! Creates various boundaries for the samples 
  !IF(rank == master_proc)  WRITE(*,*) " Initializion OK"
  ! CALL TestKiss()

  ! --- Observational constraints --- !

  CALL Prepare_Constraints()                        ! Generate the tables with observational data
  !IF(rank == master_proc)  WRITE(*,*) " Prepare_Constraints OK"
  delta_chi2_1 = Calc_DeltaChi2(0.6827d0, REAL(dof,8)) 
  delta_chi2_2 = Calc_DeltaChi2(0.9545d0, REAL(dof,8)) 
  delta_chi2_3 = Calc_DeltaChi2(0.9973d0, REAL(dof,8)) 

  IF(rank == master_proc) THEN
     CALL WRITE_INFO()
  END IF
  
  ! Reprise mode 
  CALL Reprise_mode()
  !IF(rank == master_proc)  WRITE(*,*) " Reprise_mode OK"

  ! ----- MonteCarlo routine ----- !

  IF(run_mode == param_search) THEN
    
     DO
        good_model = 0

        ! Save random generator indexes for reproductibility
        CALL SaveKISS()
        !IF(rank == master_proc)  WRITE(*,*) " Save_KISS OK"

        ! Random draw of parameters
        IF(rank == master_proc) THEN
           CALL Draw_model_parameters()
        END IF
        !IF(rank == master_proc)  WRITE(*,*) " Draw_model_parameters OK"
        CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
        CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)

        CALL MonteCarlo(hist_flag)
        !IF(rank == master_proc)  WRITE(*,*) " MonteCarlo  OK"

      
        IF(Chi2(0) < Chi2_min) THEN
           ! Rerun model with correct seeds and save histograms
           iKiss = TabSave_ijk_KISS(1,rank)
           jKiss = TabSave_ijk_KISS(2,rank)
           kKiss = TabSave_ijk_KISS(3,rank)
           IF(rank == master_proc) THEN
              CALL Draw_model_parameters()
           END IF
           CALL MPI_BCAST(TabParam_Lum,  NParam_Lum,  MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_z,    NParam_z,    MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_Ep,   NParam_Ep,   MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           CALL MPI_BCAST(TabParam_spec, NParam_spec, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
           hist_flag = 1
           CALL MonteCarlo(hist_flag)
           hist_flag = 0
           IF(rank == master_proc) THEN
              WRITE(*,*) "New best model :", Chi2(0)," (compared to :",Chi2_min,"). Histograms saved."
              CALL Update_best_parameters()
           END IF
           Chi2_min = Chi2(0)
        END IF
        
        IF(Chi2(0) <= dof + delta_chi2_3) THEN ! 3 sigma interval for now 
           good_model = 1
           Nb_good_models = Nb_good_models + 1      
        END IF
       
        Nb_models = Nb_models + 1
       
        ! Write in chi2_models
        
        IF(rank == master_proc) THEN    
           OPEN(UNIT=56, FILE=TRIM(path)//'chi2_models.dat', POSITION='append')  
        END IF
        ! This needs to be automatized
        IF(rank == master_proc) WRITE(56,*) good_model, TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_Lmax),TabParam_Lum(Param_Lum_slope),&
             & TabParam_Ep(Param_Ep_Ep0), TabParam_Ep(Param_Ep_sigmaLog), Chi2(0), Chi2(Constraint_Stern), Chi2(Constraint_EpGBM), k_Stern, norm_EpGBM
        CLOSE(56)
        
        ! Write in reprise
        OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted', POSITION='append')
        IF(rank == master_proc) WRITE(43) TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
             & Chi2, R_GRB, k_Stern, k_Preece, norm_EpGBM, Frac_XRFHETE2, dof, chi2_min, Nb_good_models, Nb_models
        CLOSE(43) ! Close reprise file

        IF(Nb_good_models == 500) THEN 
           IF(rank == master_proc)  WRITE(*,*) "Reached 500 good models"
           EXIT
        END IF
        IF(Nb_models >= 50000) THEN
           IF(rank == master_proc)   WRITE(*,*) "Reached 50000 models"
           EXIT
        END IF
        
        IF(rank == master_proc)THEN  ! This needs to be automatized
           IF(MOD(Nb_models, 5) == 0) WRITE(*,'(I6,I4,A,F5.1,A,5ES12.5)') Nb_models, Nb_good_models, ' Best Chi2 : ',Chi2_min,' for Lmin, Lmax, slope, Ep0, sigma :',&
             & TabBestParam_Lum(Param_Lum_Lmin), TabBestParam_Lum(Param_Lum_Lmax), TabBestParam_Lum(Param_Lum_slope), TabBestParam_Ep(Param_Ep_Ep0), TabBestParam_Ep(Param_Ep_sigmaLog)
        END IF
     END DO

    
    

  ! Markov Chain Monte Carlo
!!$  ELSE IF(Mode_MCMC) THEN
!!$     Lmin_min = 49.d0
!!$     Lmin_max = 52.d0
!!$     Lmax_min = 52.d0
!!$     Lmax_max = 55.d0
!!$     slope_min = 0.5d0
!!$     slope_max = 3.d0
!!$     slope_step = 0.125d0
!!$     Lmin_step  = 0.2d0
!!$     Lmax_step  = 0.2d0
!!$    
!!$     ! Random starting point :
!!$     Lmin_start = 1.60*10**50.d0
!!$     Lmax_start = 1.98*10**53.d0
!!$     slope_start = 1.70d0
!!$     MCMC_run_nb = 0
!!$
!!$
!!$     CALL One_MCMC_run(Lmin_start, slope_start, Lmin_step, slope_step, Niter, MCMC_run_nb)
!!$     t = uniform()
!!$     Lmin_start = 10.d0**( Lmin_min + t * (Lmin_max - Lmin_min) )
!!$     t = uniform()
!!$     Lmax_start = 10.d0**( Lmax_min + t * (Lmax_max - Lmax_min) )
!!$     t = uniform()
!!$     slope_start = slope_min + t * (slope_max - slope_min)
!!$     WRITE(*,*) " Lmin_start, slope_start = ", Lmin_start, slope_start
!!$ 
 
 
  ELSE

  IF (rank == master_proc) WRITE(*,*) " You chose one run mode"   
     
  DO ii = 1, 1

     CALL SaveKISS()
     
     CALL MonteCarlo(hist_flag) 
 
     ! Write in reprise
     OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted', POSITION='append')
     IF(rank == master_proc) WRITE(43) TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
             & Chi2, R_GRB, k_Stern, k_Preece, norm_EpGBM, Frac_XRFHETE2, dof, chi2_min, Nb_good_models, Nb_models
      CLOSE(43) ! Close reprise file
  END DO
  
END IF
 
  
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
    CHARACTER(*), PARAMETER                         :: GRB_PropFile = "GRB_Properties.dat"    ! Filename for saving the properties
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
    CALL Prepare_Luminosity()                  ! Prepare distribution function for Luminosity draws  (essentially unused for now : 23/09/2016, because Lum func is analytical)
    ! CAREFUL THIS HAS BEEN CHANGED 03/05 (because not used if Lum Funct is analytical) ----> 23/09/2016 Note : not a big deal for now, no needs for caps...
    CALL Prepare_Spec()                        ! Prepare distribution function for alpha and beta draws (currently only from GBM catalog)

  !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Select model  OK"
    
    ! ---------------------------------- !
    ! --- (1.bis) Prepare histograms --- !
    ! ---------------------------------- !
    
    IF(hist_flag == 1) THEN
       CALL Prepare_Histogram_Prop(TabHistlim_inf, TabHistlim_sup, n_call)        ! Create histogram limits and resets them
    END IF
    CALL Reset_Histogram_Chi2()                                                 ! Resets histograms used in Chi2 calculation
  !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Prepare histgromas  OK"
    ! -------------------------------------- !
    ! ---------- 2. Generate GRBs ---------- !
    ! -------------------------------------- !
    
    IF(run_mode == one_run) THEN
       IF(Save_all_GRB .EQV. .TRUE.) THEN
          OPEN(UNIT=66, FILE=GRB_PropFile)                ! Save all GRB properties, be careful the file is heavy (several GB)
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
       !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo z  OK : ", z
       CALL Draw_Luminosity(L)
       !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Lum  OK : ", L
       CALL Draw_alpha_beta(alpha, beta)
       !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo alpha, beta  OK : ", alpha, beta
       CALL Draw_Ep(Ep)
       !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Ep  OK  : ", Ep
       ! ------------ Observed properties ----------- !      
       Epobs = Calc_Epobs(Ep,z)  
       CALL Calculate_Peakflux_Instrument()
     !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Caluclate_Peakflux_Instrument  OK"
       CALL Calculate_Detection_Probability()
     !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Calculate_Detection_Probability  OK"
       ! --------------- Histograms --------------- !
       ! ------------------------------------------ !

       IF(hist_flag == 1) THEN
          CALL Fill_Histogram_Prop()                       ! Find the bins and add the probability of detection to the histograms of L, z, Ep etc...
       END IF

       CALL Fill_Histogram_Chi2()                          ! Find the bins and add the probability of detection to the histograms to use in Chi2 calculation

     !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Histograms  OK"

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
          
          IF(Save_all_GRB .EQV. .TRUE.)  WRITE(66,'(6ES12.5)') REAL(i,8), L, z, Epobs, alpha, beta
          
          
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
                  
          IF(Save_all_GRB .EQV. .TRUE.) CLOSE(66)     
                
          ! Track progress
          IF(rank == master_proc) THEN
             IF( MOD(100.d0 * REAL(i,8)*REAL(nb_procs,8)/REAL(Nb_GRB,8), 5.d0) == 0 ) THEN
                WRITE(*,'(I3,1x)', advance='no') INT(100.d0 * REAL(i,8)*REAL(nb_procs,8)/REAL(Nb_GRB))
             END IF
          END IF
       END IF
       
    END DO

    IF(run_mode == one_run) THEN
       IF(rank == master_proc) WRITE(*, '(A)') "      ]"
    END IF
   ! --------------------------------- !
   ! ----- End of : Generate GRB ----- !
   ! --------------------------------- !
   
   ! ----------------------------------------------------------------------------------- END OF MAIN LOOP ----------------------------------------------------------------------------------- !

   ! -------- Combine histograms of each procs -------- !
   IF(Constraint_Included(Constraint_Kommers)) CALL MPI_REDUCE(TabHistKomm_P23,    TabHistKomm_P23_master,    N_Komm,   MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
   IF(Constraint_Included(Constraint_Stern))   CALL MPI_REDUCE(TabHistStern_P23,   TabHistStern_P23_master,   N_Stern,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
   IF(Constraint_Included(Constraint_Preece))  CALL MPI_REDUCE(TabHistPreece_Ep,   TabHistPreece_Ep_master,   N_Preece, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)
   IF(Constraint_Included(Constraint_EpGBM))   CALL MPI_REDUCE(TabHistEpGBM_Epobs, TabHistEpGBM_Epobs_master, N_EpGBM,  MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD, code)

   IF(hist_flag == 1) THEN
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
            CALL MPI_REDUCE(    GRB_Prop(i_Sample,:),     GRB_Prop_master(i_Sample,:), N_saved_Prop+1, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD,code)
            CALL MPI_REDUCE(average_Prop(i_Sample,:), average_Prop_master(i_Sample,:), N_saved_Prop+1, MPI_REAL8, MPI_SUM, master_proc, MPI_COMM_WORLD,code)
         END IF
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
      END IF

      ! ---------- 3. Normalizations ---------- !
      ! --------------------------------------- !

      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                                                                                   ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                    Normalization coefficients                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      CALL Normalize_Model()
      !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Normalize_Model  OK"
      
      ! ----- 4. Chi squared calculation ----- !
      ! -------------------------------------- !
      
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                                                                                   ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[                                          Chi squared                                              ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      CALL Calculate_Chi2()
      !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo Calculate_Chi2  OK"
      
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                           Chi2 = ", Chi2(0), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                           Delta Chi2 (3 sigma) = ", dof + delta_chi2_3, "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A,I4.0,A)')   "[                             Degrees of freedom =   ", dof, "                                           ]"
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                   Reduced Chi2 = ", Chi2(0)/REAL(dof,8), "                                     ]"
      IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
      
      ! ----- (5.) Save histograms ----- !
      ! -------------------------------- !
      
      IF(hist_flag == 1) CALL Save_Histograms()
      
      
      ! ----- (6.) Save constraints ----- !
      ! --------------------------------- !
      
      IF(hist_flag == 1)  CALL Save_Constraints()
      
      IF(run_mode == one_run)  WRITE(*,*) " "
      IF(run_mode == one_run)  WRITE(*,'(A)') " ==================================== Monte Carlo done ============================================== "
      IF(run_mode == one_run)  WRITE(*,*) " "
      
      
   END IF ! Ends master rank computation
   
  ! Broadcast chi2
   CALL MPI_BCAST(Chi2, N_Constraints, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
   

 END SUBROUTINE MonteCarlo
 

 
 SUBROUTINE ReadInitFile(Name)
   INTEGER :: s_L,s_z,s_Spec,s_Ep,s_P,s_Constraint,i,incl,temp
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
   IF(SHOW_READ) WRITE(*,'(A,A37,A)')           "[            Output path : ",ADJUSTL(path)," ]"

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
      hist_flag = 1
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : one run                   ]"
   ELSE IF (run_mode == param_search) THEN
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : parameter search          ]"
   ELSE
      IF(SHOW_READ) WRITE(*,'(A)') "[                    run mode chosen : INVALID                   ]"
   END IF


   ! Number of GRBs 
   
   CALL ReadInteger(50, Nb_GRB) 
   Nb_GRB = 6 * 10**(Nb_GRB)
   IF(SHOW_READ) WRITE(*,'(A,1ES12.5,A)') "[            Number of GRBs simulated : ", REAL(Nb_GRB,8), "             ]"
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
      
   CASE(Model_LumBPL) ! Broken Power Law distribution
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmin))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmin =",TabParam_Lum(Param_Lum_Lmin),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lmax))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmax =",TabParam_Lum(Param_Lum_Lmax),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_Lbreak))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lbreak =",TabParam_Lum(Param_Lum_Lbreak),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slopeL))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeL =",TabParam_Lum(Param_Lum_slopeL),"  ]"
      CALL ReadReal(50, TabParam_Lum(Param_Lum_slopeH))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeH =",TabParam_Lum(Param_Lum_slopeH),"  ]"
      
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
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Log Normal ->", "mu =", TabParam_Ep(Param_Ep_Ep0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_sigmaLog))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Log Normal ->", " sigma =", TabParam_Ep(Param_Ep_sigmaLog),"  ]"
   CASE(Model_EpY) ! Yonetoku 
      STOP "Model_Ep : Yonetoku not coded yet"     
   CASE(Model_EpAmati)
      CALL ReadReal(50, TabParam_Ep(Param_Ep_Ep0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", "Ep0 =", TabParam_Ep(Param_Ep_Ep0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_sigmaLog))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " sigma =", TabParam_Ep(Param_Ep_sigmaLog),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_L0))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " L0 =", TabParam_Ep(Param_Ep_L0),"  ]"
      CALL ReadReal(50, TabParam_Ep(Param_Ep_alpha_amati))
      IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Ep :", "Amati-like ->", " alpha_amati =", TabParam_Ep(Param_Ep_alpha_amati),"  ]"
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
            Instrument_Included(Instrument_SWIFT) = .TRUE.
         CASE(Sample_SWIFT)
            Instrument_Included(Instrument_SWIFT) = .TRUE.
         CASE(Sample_SWIFTbright)
            Instrument_Included(Instrument_SWIFT) = .TRUE.
         CASE(Sample_HETE2)
            Instrument_Included(Instrument_FREGATE) = .TRUE.
            Instrument_Included(Instrument_WXM) = .TRUE.
         CASE(Sample_BAT6ext)
            Instrument_Included(Instrument_SWIFT) = .TRUE.
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
         CASE(Constraint_XRFHETE2)
            Sample_Included(Sample_HETE2) = .TRUE.
            Instrument_Included(Instrument_FREGATE) = .TRUE.
            Instrument_Included(Instrument_WXM) = .TRUE.
         CASE(Constraint_EpGBM)
            Sample_Included(Sample_GBM) = .TRUE.
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
   OPEN(UNIT=51, FILE='../Input_para/lum_param_search.init')
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

          
      CASE(Model_LumBPL) ! Broken Power Law distribution
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmin_min =",TabParam_Lum_min(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmin))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmin_max =",TabParam_Lum_max(Param_Lum_Lmin),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmax_min =",TabParam_Lum_min(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lmax))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lmax_max =",TabParam_Lum_max(Param_Lum_Lmax),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lbreak_min =",TabParam_Lum_min(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_Lbreak))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "Lbreak_max =",TabParam_Lum_max(Param_Lum_Lbreak),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slopeL))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeL_min =",TabParam_Lum_min(Param_Lum_slopeL),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slopeL))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeL_max =",TabParam_Lum_max(Param_Lum_slopeL),"  ]"
         CALL ReadReal(51, TabParam_Lum_min(Param_Lum_slopeH))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeH_min =",TabParam_Lum_min(Param_Lum_slopeH),"  ]"
         CALL ReadReal(51, TabParam_Lum_max(Param_Lum_slopeH))
         IF(SHOW_READ) WRITE(*,Format_RIF) "[ ","Model_Lum :", "Broken Power Law ->", "slopeH_max =",TabParam_Lum_max(Param_Lum_slopeH),"  ]"

   END SELECT
   END IF
         
   CLOSE(51)

   ! Redshift
   OPEN(UNIT=52, FILE='../Input_para/redshift_param_search.init')
   CALL ReadInteger(52, redshift_explore)

   IF(redshift_explore == 1) THEN
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
      CASE DEFAULT
         IF(SHOW_READ) WRITE(*,'(A)') "[           No parameter search for this redshift model          ]"
         IF(SHOW_READ) WRITE(*,'(A)') "[ -------------------------------------------------------------- ]"   
      END SELECT
   END IF
   CLOSE(52)
   
   ! Peak energy
   OPEN(UNIT=53, FILE='../Input_para/Ep_param_search.init')
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
      CASE(Model_LumBPL)
         IF(TabParam_Lum_min(Param_Lum_Lmin)   /= TabParam_Lum_max(Param_Lum_Lmin))   dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lmax)   /= TabParam_Lum_max(Param_Lum_Lmax))   dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_Lbreak) /= TabParam_Lum_max(Param_Lum_Lbreak)) dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slopeL) /= TabParam_Lum_max(Param_Lum_slopeL)) dof = dof - 1
         IF(TabParam_Lum_min(Param_Lum_slopeH) /= TabParam_Lum_max(Param_Lum_slopeH)) dof = dof - 1
      END SELECT
   END IF
   
   IF(redshift_explore == 1) THEN
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
         IF(i_Constraint == Constraint_XRFHETE2) dof = dof + 1
         IF(i_Constraint == Constraint_EpGBM)    dof = dof + N_EpGBM  - 1   
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
   IF(Model_Lum == Model_LumBPL) t_star = Calc_t_star(TabParam_Lum)
   
   
   ! --- Spectral Model --- !
   
   IF(Model_Spec == Model_SpecBPLFix .OR. Model_Spec == Model_SpecBandFix) THEN
      alpha = TabParam_Spec(Param_spec_alpha)
      beta  = TabParam_Spec(Param_spec_beta)
      ktild = Calc_ktild(alpha, beta, Model_Spec)
   END IF
   IF(verbose == 1) PRINT*, "ktild =", ktild
   
   
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
         CASE(Model_LumBPL)
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
 END SUBROUTINE Draw_model_parameters
 
 SUBROUTINE Update_best_parameters()

   
   OPEN(UNIT=876, FILE=TRIM(path)//'best_chi2.txt')
   WRITE(876,'(A)') "# This file contains the information of the best model"
   WRITE(876,'(A)') "# "
   WRITE(876,'(A,ES12.5)') "# Best Chi2 : ", Chi2_min
   WRITE(876,'(A,I3)') "# degrees of freedom : ", dof
   WRITE(876,'(A,F6.3)') "# 3 sigma interval : ", delta_chi2_3
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
   WRITE(876,'(A)') "# Luminosity function : (2) Broken Power Law (5 args : Lmin, Lmax, Lbreak, slopeL, slopeH)"
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
   CASE(Model_LumBPL)
      TabBestParam_Lum(Param_Lum_Lmin) = TabParam_Lum(Param_Lum_Lmin)
      TabBestParam_Lum(Param_Lum_Lmax) = TabParam_Lum(Param_Lum_Lmax)
      TabBestParam_Lum(Param_Lum_Lbreak) = TabParam_Lum(Param_Lum_Lbreak)
      TabBestParam_Lum(Param_Lum_slopeL) = TabParam_Lum(Param_Lum_slopeL)
      TabBestParam_Lum(Param_Lum_slopeH) = TabParam_Lum(Param_Lum_slopeH)
      WRITE(876,'(I1)') Model_LumBPL
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmin)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lmax)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_Lbreak)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slopeL)
      WRITE(876,'(ES12.5)') TabBestParam_Lum(Param_Lum_slopeH)
   END SELECT
   
   WRITE(876,'(A)') "# "
   WRITE(876,'(A)') "# Peak energy : (0) Fix (1 arg : Ep0)"
   WRITE(876,'(A)') "# Peak energy : (1) LogNormal (2 args : Ep0[keV], sigmaLog)"
   WRITE(876,'(A)') "# Peak energy : (3) LogNormal (4 args : Ep0[keV], sigmaLog, L0[erg/s], alpha_amati)"
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

   CLOSE(876)

 END SUBROUTINE Update_best_parameters

 SUBROUTINE Draw_Redshift(z)
   ! --- Redshift ---!
   REAL(8), INTENT(out) :: z
   INTEGER              :: imin, imax, j

   IF(Model_z /= Model_zFix) THEN
      DO
         t = uniform()
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
   
   IF (verbose==1 .OR. verbose==2) PRINT*, "z =", z
   IF (verbose==1 .OR. verbose==2) PRINT*, "D_L = ", D_L, "[Mpc]"
   IF (ISNAN(z)) NaNtest_prop = .TRUE.
   
 END SUBROUTINE Draw_Redshift
 
 SUBROUTINE Draw_Luminosity(L)
   ! -- Luminosity [erg/s] -- !
   REAL(8), INTENT(out) :: L
  
   SELECT CASE(Model_Lum)
   CASE(Model_LumPL) 
      t = uniform()
      IF(TabParam_Lum(Param_Lum_slope) == 1.d0) THEN
         L = TabParam_Lum(Param_Lum_Lmin) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**t
      ELSE
         L = TabParam_Lum(Param_Lum_Lmin) * &
              & ( 1.d0 - t*( 1.d0 - (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**(1.d0-TabParam_Lum(Param_Lum_slope)) )&
              & )**(1.d0/ (1.d0-TabParam_Lum(Param_Lum_slope)) )
      END IF
      
   CASE(Model_LumBPL)
      t = uniform()
      IF(TabParam_Lum(Param_Lum_slopeL) == 1) THEN
         IF(TabParam_Lum(Param_Lum_slopeH) == 1) THEN
            L = TabParam_Lum(Param_Lum_Lmin) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lmin))**t
         ELSE
            IF(t < t_star) THEN 
               L = TabParam_Lum(Param_Lum_Lmin) * (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))**(t/t_star)
            ELSE
               L = TabParam_Lum(Param_Lum_Lbreak) * ( 1.d0 - (t-t_star) * (1.d0 -&
                    & ( TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak) )**(1.d0-TabParam_Lum(Param_Lum_slopeH)) )  / (1.d0-t_star) )&
                    & **( 1.d0/(1.d0 - TabParam_Lum(Param_Lum_slopeH)) )
            END IF
         END IF
      ELSE
         IF(TabParam_Lum(Param_Lum_slopeH) == 1) THEN
            t = uniform()
            IF(t < t_star) THEN
               L = TabParam_Lum(Param_Lum_Lmin) * (1.d0 - t/t_star * (1.d0 - (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))&
                    & **(1.d0 -TabParam_Lum(Param_Lum_slopeL))) )**(1.d0/(1.d0-TabParam_Lum(Param_Lum_slopeL)))
            ELSE 
               L = TabParam_Lum(Param_Lum_Lbreak) * (TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak))&
                    & **((t-t_star)/(1.d0-t_star))
            END IF
         ELSE
            t = uniform()
            IF(t < t_star) THEN
               L = TabParam_Lum(Param_Lum_Lmin) * (1.d0 - t/t_star * (1.d0 - (TabParam_Lum(Param_Lum_Lbreak)/TabParam_Lum(Param_Lum_Lmin))&
                    & **(1.d0 -TabParam_Lum(Param_Lum_slopeL))) )**(1.d0/(1.d0-TabParam_Lum(Param_Lum_slopeL)))
            ELSE
               L = TabParam_Lum(Param_Lum_Lbreak) * ( 1.d0 - (t-t_star) * (1.d0 -&
                    & ( TabParam_Lum(Param_Lum_Lmax)/TabParam_Lum(Param_Lum_Lbreak) )**(1.d0-TabParam_Lum(Param_Lum_slopeH)) )  / (1.d0-t_star) )&
                    & **( 1.d0/(1.d0 - TabParam_Lum(Param_Lum_slopeH)) )
            END IF
         END IF
      END IF
      
   CASE DEFAULT
      STOP "Monte Carlo (2) : ERROR Lum, not implemented yet"
      
   END SELECT
   IF (verbose==1 .OR. verbose==2) PRINT*, "L =", L, "[erg/s]"
   IF (ISNAN(L)) NaNtest_prop = .TRUE.
   
 END SUBROUTINE Draw_Luminosity

 SUBROUTINE Draw_alpha_beta(alpha, beta)
   ! --- Spectral Model --- ! 
   REAL(8), INTENT(out) :: alpha, beta
   INTEGER              :: imin, imax, j
   
   IF(Model_Spec == Model_SpecBandD) THEN
      
      DO 
         alpha = gaussian(1.d0, 0.5d0)
         IF (alpha < 2.0d0) EXIT
      END DO
      
      DO 
         t = uniform()
         IF( t <= 0.15625d0) THEN
            beta =  2.d0 + SQRT(t/2.5d0)
         ELSE
            beta = (9.d0 - SQRT(13.5d0*(1.d0-t)) )/2.5d0
         END IF
         IF (beta > 2.0d0) EXIT
      END DO
      
   ELSE IF (Model_Spec == Model_SpecBandGBM) THEN
      DO
         t = uniform()
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
      !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, alpha   OK : ", alpha
      DO
         t = uniform()
         IF((t > 0.d0) .AND. (t < 1.0d0)) EXIT
      END DO
      !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, t for beta   OK : ", t
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
        !IF(rank == master_proc)  WRITE(*,*) " in MonteCarlo, in Draw_alpha_beta, beta   OK : ", beta
   END IF
   ktild = Calc_ktild(alpha, beta, Model_Spec) 
   IF (ISNAN(ktild)) NaNtest_prop = .TRUE.
   IF (ISNAN(alpha)) NaNtest_prop = .TRUE.
   IF (ISNAN(beta)) NaNtest_prop = .TRUE.
   IF(ISNAN(ktild)) WRITE(*,'(A, 3ES12.5)')  "param(ktild, alpha, beta-alpha)  :  ", ktild, alpha, (beta-alpha) 
   IF(verbose == 1) THEN
      WRITE(*,'(A,F5.2)') 'alpha = ', alpha
      WRITE(*,'(A,F5.2)') 'beta = ', beta        
   END IF
   
   
 END SUBROUTINE Draw_alpha_beta
 
 SUBROUTINE Draw_Ep(Ep)
   ! --- Peak Energy [keV] --- !
   REAL(8), INTENT(out) :: Ep
   
   SELECT CASE(Model_Ep)
   CASE(Model_EpLogNormal)
      DO
         Ep = 10**( gaussian(LOG10(TabParam_Ep(Param_Ep_Ep0)), TabParam_Ep(Param_Ep_sigmaLog)) )
         IF(Ep > 0.d0) EXIT
      END DO
   CASE(Model_EpAmati)
      DO 
         t  = gaussian(0.d0, TabParam_Ep(Param_Ep_sigmaLog))
         Ep = TabParam_Ep(Param_Ep_Ep0) * (L/TabParam_Ep(Param_Ep_L0))**TabParam_Ep(Param_Ep_alpha_amati) * &
              & 10.d0**( SQRT(1.d0 + TabParam_Ep(Param_Ep_alpha_amati)**2) * t )
        
         IF(Ep > 0.d0) EXIT
      END DO
   END SELECT
   IF (ISNAN(Ep)) NaNtest_prop = .TRUE.
   
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
   TabEmin(Instrument_SWIFT) = 15.d0
   TabEmax(Instrument_SWIFT) = 150.d0

   ! FREGATE
   TabEmin(Instrument_FREGATE) = 30.d0
   TabEmax(Instrument_FREGATE) = 400.d0

   ! WXM
   TabEmin(Instrument_WXM) = 2.d0
   TabEmax(Instrument_WXM) = 10.d0

   ! --- Detection Threshold [ph/cm2/s] --- !
   ! -------------------------------------- !

   Threshold(Sample_Intrinsic)   = 0.d0
   Threshold(Sample_Kommers)     = 0.d0   ! 50-300 keV
   Threshold(Sample_Preece)      = 5.d0   ! 50-300 keV
   Threshold(Sample_Stern)       = 0.d0   ! 50-300 keV
   Threshold(Sample_SWIFTweak)   = 0.01d0 ! 15-150 keV
   Threshold(Sample_SWIFT)       = 0.2d0  ! 15-150 keV
   Threshold(Sample_SWIFTbright) = 1.d0   ! 15-150 keV
   Threshold(Sample_HETE2)       = 1.d0   ! 2-10 keV or 30-400 keV
   Threshold(Sample_BAT6ext)     = 2.6d0  ! 15-150 keV
   Threshold(Sample_GBM)         = 0.9d0  ! 50-300 keV
 
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
 
 SUBROUTINE TestKiss()
  INTEGER :: i
  REAL(8) :: x
  
  DO i=1, 1000000000
     x = uniform()
     IF (MODULO(i,100000000)==0) WRITE(*,'(4(A,I12))') "i=",i,"  Kiss : i=",iKiss," j=",jKiss," k=",kKiss
  END DO
END SUBROUTINE TestKiss

SUBROUTINE SaveKISS()
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
     WRITE(325, *) '# Last saved KISS seeds i,j,k for ',nb_procs,' cores' 
     DO i=0, nb_procs-1
        WRITE(325, *) TabSave_ijk_KISS(1,i),TabSave_ijk_KISS(2,i),TabSave_ijk_KISS(3,i)
     END DO
     CLOSE(325)
  END IF
END SUBROUTINE SaveKISS

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
      call RANDOM_NUMBER(uniform)
   END IF
 END FUNCTION uniform



 FUNCTION Kiss()
   INTEGER :: Kiss
   iKiss=69069*iKiss+32606797
   jKiss=Melange(Melange(jKiss,17),-15)
   kKiss=Melange(IAND(Melange(kKiss,18),2147483647),-13)
   Kiss=iKiss+jKiss+kKiss
 END FUNCTION Kiss


 
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
   gaussian = mu + sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   !sigma * ( COS(theta) * SQRT(-LOG(t)) ) 
   ! for comparison with Daigne 2006 add sqrt(2)
 END FUNCTION gaussian
 


 REAL(8) FUNCTION Rz(z, Model_z, Tabparam_z)
   ! GRB rate as a function of redshift [GRB/yr/Mpc^3]
   ! b   : slope at low z
   ! b-a : slope at high z
   
   REAL(8), INTENT(in)                  :: z
   INTEGER, INTENT(in)                  :: Model_z
   REAL(8), DIMENSION(1:10), INTENT(in) :: Tabparam_z
   REAL(8)                              :: a, b, c, d, zm
   
   SELECT CASE(Model_z)
   CASE(Model_zFix)
      Rz = 1.d0
   CASE(Model_zUniform)
      Rz = 1.d0
   CASE(Model_zSH)
      zm = TabParam_z(Param_z_zm)
      a  = TabParam_z(Param_z_a)
      b  = TabParam_z(Param_z_b)
      Rz = a * EXP(b*(z-zm)) / ( a-b + b*EXP(a*(z-zm)) )
   CASE(Model_zDaigne)
      a = TabParam_z(Param_z_a)
      b = TabParam_z(Param_z_b)
      c = TabParam_z(Param_z_c)
      d = TabParam_z(Param_z_d)
      Rz = 0.0122d0 * a * EXP(b*z) / ( EXP(c*z) + d )
   END SELECT
   
 END FUNCTION Rz
 
 REAL(8) FUNCTION Calc_Epobs(Ep,z)
   REAL(8), INTENT(in) :: Ep, z
   ! --- Peak Energy (obs) [keV] --- !
   Calc_Epobs = Ep / (1.d0 + z)
   IF (ISNAN(Epobs)) NaNtest_prop = .TRUE.
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

   CASE(Model_LumBPL)
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

   CASE DEFAULT
      STOP "Calc_Lum : error, model not defined."
      
   END SELECT
  
 END FUNCTION Calc_Lum
 
 REAL(8) FUNCTION Calc_ktild(alpha, beta, Model_Spec)
   ! Calculates the normalization factor for the spectral shape
   ! Uses the Gamma functions for optimized speed. See stats.f90 for details
   REAL(8), INTENT(in) :: alpha, beta 
   INTEGER, INTENT(in) :: Model_Spec
   REAL(8)             :: epsilon,x1,x2,x_c,xmax, xmin, dlogx, xBx1, xBx2, xBx, Gamma, gln, xBx_test
   INTEGER             :: j, M

   IF(Model_Spec == Model_SpecBPLFix) THEN
      Calc_ktild = (2.d0-alpha) * (beta-2.d0) / (beta-alpha)
   ELSE
      x_c     = (beta - alpha) / (2.d0 - alpha)
      CALL GammaSeries(Gamma, 2.d0-alpha, beta-alpha, gln)
      xBx = Gamma * EXP(gln) /( (2.d0 - alpha)**(2.d0-alpha) )
      xBx = xBx +  x_c**(2.d0-alpha) * EXP(alpha-beta) / (beta - 2.d0)
      Calc_ktild = 1.d0/xBx
   END IF

 END FUNCTION Calc_ktild
 
 SUBROUTINE Calculate_Peakflux_Instrument()
   ! --- Number of Photons between Emin and Emax [ph/cm^2/s] --- !

   DO i_Instrument=1, N_Instruments
      IF(Instrument_Included(i_Instrument)) THEN
         Peakflux_Instrument(i_Instrument) = Nb_ph(L, z, Ep, D_L, alpha, beta, ktild, TabEmin(i_Instrument), TabEmax(i_Instrument))
         IF (verbose==1 .OR. verbose==2) THEN 
            WRITE(*,'(A,A20,ES12.5,A,F3.0,A,F4.0,A)')&
                 & " Instrument : " , TRIM(TabInstrument_name(i_Instrument)) // " Peakflux =",Peakflux_Instrument(i_Instrument),&
                 & " [ph/cm2/s between ",TabEmin(i_Instrument)," and ", TabEmax(i_Instrument), " keV]"
         END IF
         IF (ISNAN(Peakflux_Instrument(i_Instrument))) NaNtest_prop = .TRUE.
      END IF
   END DO
   
 END SUBROUTINE Calculate_Peakflux_Instrument

 SUBROUTINE Calculate_Detection_Probability()
   
   ! ---------------- Detection Probability [0,1] ---------------- !
   
   DO i_Sample=0, N_Samples
      IF(i_Sample >= 1) THEN
         IF(Sample_Included(i_Sample)) THEN      
            SELECT CASE(i_Sample)
            CASE(Sample_Kommers)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               ! BATSE23 (Kommers et al. 2000, eq (4) ) :
               Prob_det(Sample_Kommers) = 0.5d0 * ( 1.d0 + ERF(-4.801d0 + 29.868d0 * Peakflux_Instrument(Instrument_BATSE) ))
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
               Prob_det(Sample_Stern) = 1.d0
               
            CASE(Sample_SWIFTweak)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_SWIFT)
               IF(Peakflux_Instrument(Instrument_SWIFT) >= Threshold(Sample_Swiftweak)) THEN
                  Prob_det(Sample_Swiftweak) = 1.d0
               ELSE
                  Prob_det(Sample_Swiftweak) = 0.d0
               END IF
               
            CASE(Sample_SWIFT)
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_SWIFT)
               IF(Peakflux_Instrument(Instrument_SWIFT) >= Threshold(Sample_Swift)) THEN
                  Prob_det(Sample_Swift) = 1.d0
               ELSE
                  Prob_det(Sample_Swift) = 0.d0
               END IF
               
            CASE(Sample_SWIFTbright) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_SWIFT)
               IF(Peakflux_Instrument(Instrument_SWIFT) >= Threshold(Sample_Swiftbright)) THEN
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
            CASE(Sample_BAT6ext) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_SWIFT)
               IF(Peakflux_Instrument(Instrument_SWIFT) >= Threshold(Sample_BAT6ext)) THEN
                  Prob_det(Sample_BAT6ext) = 1.d0
               ELSE
                  Prob_det(Sample_BAT6ext) = 0.d0
               END IF
            CASE(Sample_GBM) 
               Peakflux(i_Sample) = Peakflux_Instrument(Instrument_BATSE)
               IF(Peakflux_Instrument(Instrument_BATSE) >= Threshold(Sample_GBM)) THEN
                  Prob_det(Sample_GBM) = 1.d0
               ELSE
                  Prob_det(Sample_GBM) = 0.d0
               END IF  
            CASE DEFAULT
               STOP "Error, undefined sample."
               
            END SELECT
            
            IF (ISNAN(Peakflux(i_Sample))) NaNtest_prop = .TRUE.
            
            IF (verbose==1 .OR. verbose==2) PRINT*, "Prob_det of ",TRIM(TabSample_name(i_Sample)) ," = ", Prob_det(i_Sample)       
         END IF
      ELSE
         Prob_det(i_Sample) = 1.d0
      END IF
      !         IF(verbose == 1 ) PRINT*, "Prob_det of "//TRIM(TabSample_name(i_Sample))//" :", Prob_det(i_Sample)
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
   IF(verbose == 1) WRITE(*,*) " Number of steps for Nb_ph integration : ", M2
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
   Nb_ph = (1.d0+z) * L /(4.d0*Pi*(Ep*keV)*(D_L*Mpc)**2) * B
 END FUNCTION Nb_ph


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
   IF(verbose==1) PRINT*, "B =", B
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
   TabSample_name(Sample_BAT6ext)     = "BAT6ext_sample"
   TabSample_name(Sample_GBM)         = "GBM_sample"
   
   TabInstrument_name(Instrument_BATSE)   = "BATSE"
   TabInstrument_name(Instrument_SWIFT)   = "SWIFT"
   TabInstrument_name(Instrument_FREGATE) = "FREGATE"
   TabInstrument_name(Instrument_WXM)     = "WXM"

   TabConstraint_name(Constraint_Kommers)  = "Kommers"
   TabConstraint_name(Constraint_Preece)   = "Preece"
   TabConstraint_name(Constraint_Stern)    = "Stern"
   TabConstraint_name(Constraint_XRFHETE2) = "XRFHETE2"
   TabConstraint_name(Constraint_EpGBM)    = "EpGBM"
   
 END SUBROUTINE Generate_Names

 SUBROUTINE Generate_Paths()
   INTEGER :: i
   DO i=0,N_Samples
      LFile(i)      = TRIM(path) // "luminosity_"  // TRIM(TabSample_name(i)) // ".dat"
      zFile(i)      = TRIM(path) // "redshift_"    // TRIM(TabSample_name(i)) // ".dat"
      EpFile(i)     = TRIM(path) // "peakenergy_"  // TRIM(TabSample_name(i)) // ".dat"
      PFile(i)      = TRIM(path) // "peakflux_"    // TRIM(TabSample_name(i)) // ".dat"
      SpecFile_a(i) = TRIM(path) // "spec_alpha_"  // TRIM(TabSample_name(i)) // ".dat"
      SpecFile_b(i) = TRIM(path) // "spec_beta_"   // TRIM(TabSample_name(i)) // ".dat"

      zFile_cumul(i)      = TRIM(path) // "cumul_redshift_"    // TRIM(TabSample_name(i)) // ".dat"
      
      LErrorFile(i)      = TRIM(path) // "luminosity_err_" // TRIM(TabSample_name(i)) // ".dat"
      zErrorFile(i)      = TRIM(path) // "redshift_err_"   // TRIM(TabSample_name(i)) // ".dat"
      EpErrorFile(i)     = TRIM(path) // "peakenergy_err_" // TRIM(TabSample_name(i)) // ".dat"
      PErrorFile(i)      = TRIM(path) // "peakflux_err_"   // TRIM(TabSample_name(i)) // ".dat"
      SpecErrorFile_a(i) = TRIM(path) // "spec_alpha_err_" // TRIM(TabSample_name(i)) // ".dat"
      SpecErrorFile_b(i) = TRIM(path) // "spec_beta_err_"  // TRIM(TabSample_name(i)) // ".dat"  
 
   END DO

   KommFile     = TRIM(path) // KommFile
   PreeceFile   = TRIM(path) // PreeceFile
   SternFile    = TRIM(path) // SternFile
   XRFHETE2File = TRIM(path) // XRFHETE2File
   EpGBMFile    = TRIM(path) // EpGBMFile

   KommErrorFile     = TRIM(path) // KommErrorFile
   PreeceErrorFile   = TRIM(path) // PreeceErrorFile
   SternErrorFile    = TRIM(path) // SternErrorFile
   XRFHETE2ErrorFile = TRIM(path) // XRFHETE2ErrorFile
   EpGBMErrorFile    = TRIM(path) // EpGBMErrorFile
   
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
   
 END SUBROUTINE Prepare_Histogram_Prop
 
 SUBROUTINE Reset_Histogram_Chi2()

   IF (Constraint_included(Constraint_Kommers)) THEN
      R_GRB = 0.d0
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
      norm_EpGBM = 0.d0
      TabHistEpGBM_Epobs = 0.d0
      TabHistEpGBM_Epobs_master = 0.d0
   END IF
   
 END SUBROUTINE Reset_Histogram_Chi2
 

 SUBROUTINE Prepare_Luminosity()
   ! Creates distribution function for luminosity
   INTEGER :: i,j
   REAL(8) :: x1,y1,x2,y2,integ
   
   TabFctDistrLum(0) = 0.d0
   
   IF (Model_Lum > 2) THEN ! CHANGE THIS
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
      
      IF (verbose==1 .OR. verbose==2) PRINT*, 'Final value of F_Lum :', TabFctDistrLum(N_L)
      IF (verbose==1 .OR. verbose==2) READ(*,*)
   END IF
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
   CASE DEFAULT
      STOP "Prepare_Redshift : z model not defined"
   END SELECT
   
   ! indice zmax
   
   IF ( zmax > Tabprecisez(INT(SIZE(Tabprecisez)-1)) ) STOP "Aborting zmax not allowed"
   IF (zmax /= 0.d0 ) THEN
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
      Collapse_rate = TabFctDistrz(izmax)
      !WRITE(*,'(A,ES12.5,A)') ' Collapse rate :', TabFctDistrz(izmax), ' collapse/yr'
      IF (verbose == 1 .OR. verbose == 2) WRITE(*,*) 'Final value of F_z (Before normalization) :', TabFctDistrz(izmax)   
      TabFctDistrz(0:izmax) = TabFctDistrz(0:izmax) / TabFctDistrz(izmax)
   ELSE
      TabFctDistrz = 0.d0
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
    OPEN(FILE='../observational_constraints/lognlogp.stern.dat', UNIT=602)

    DO i=1, 11   ! skips the first 11 lines of the file 
       READ(602, *) skip
    END DO
    
    DO i=1, N_Stern
       READ(602,*) TabStern_P23(i-1), TabHistStern_P23obs(i), TabHistStern_P23obserr(i)
    END DO

    CLOSE(602)


    TabStern_P23(N_Stern) = TabStern_P23(N_Stern-1) + 0.2d0
    
    Global_GRB_rate = 0.d0
    DO i=1, N_Stern
       Global_GRB_rate = Global_GRB_rate + (TabStern_P23(i) - TabStern_P23(i-1)) * 10d0**TabHistStern_P23obs(i)
    END DO
    !PRINT*, 'Global GRB rate from Stern   : ', Global_GRB_rate, ' GRB/year'
    TabStern_P23 = 10.d0**(TabStern_P23)
    TabStern_P23 = TabStern_P23 / 0.75d0
    TabHistStern_P23obs = 10.d0**(TabHistStern_P23obs)
    
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

 END SUBROUTINE Prepare_Constraints
 
 SUBROUTINE Fill_Histogram_Prop()
   INTEGER              :: imin, imax, j

   ! --- Luminosity [erg/s] --- !
   iBinL = INT( ( LOG10(L)-TabHistlim_inf(Prop_LogL) ) / ( TabHistlim_sup(Prop_LogL)-TabHistlim_inf(Prop_LogL) ) * REAL(N_L,8) ) + 1
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinL =", iBinL
   
   ! --- Redshift --- !
   iBinz = INT( z / TabHistlim_sup(Prop_z) * REAL(N_z,8) ) + 1
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinz =", iBinz
   
   ! --- Peak Energy (source frame) [keV] --- !
   IF( (Ep < TabHistlim_inf(Prop_Ep)) .OR. (Ep >= TabHistlim_sup(Prop_Ep) )) THEN
      iBinEp = -1
   ELSE
      iBinEp = INT( ( LOG10(Ep)-LOG10(TabHistlim_inf(Prop_Ep)) ) / ( LOG10(TabHistlim_sup(Prop_Ep))-LOG10(TabHistlim_inf(Prop_Ep)) ) * REAL(N_Ep,8) ) + 1
   END IF
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinEp =", iBinEp 
   
   ! --- Peak Energy (observer frame) [keV] --- !
   IF( (Epobs < TabHistlim_inf(Prop_Ep)) .OR. (Epobs >= TabHistlim_sup(Prop_Ep) )) THEN
      iBinEpobs = -1
   ELSE
      iBinEpobs = INT( ( LOG10(Epobs)-LOG10(TabHistlim_inf(Prop_Ep)) ) / ( LOG10(TabHistlim_sup(Prop_Ep))-LOG10(TabHistlim_inf(Prop_Ep)) ) * REAL(N_Ep,8) ) + 1
   END IF
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinEpobs =", iBinEpobs 
   
   ! --- Peak Flux [ph/cm2/s] --- !
   DO i_Sample = 1, N_Samples
      IF(Sample_Included(i_Sample)) THEN
         IF( (LOG10(Peakflux(i_Sample)) < TabHistlim_inf(Prop_LogP)) .OR. (LOG10(Peakflux(i_Sample)) >= TabHistlim_sup(Prop_LogP) )) THEN
            iBinPeakflux(i_Sample) = -1
         ELSE
            iBinPeakflux(i_Sample) = INT( (LOG10(Peakflux(i_Sample))-TabHistlim_inf(Prop_LogP))/( TabHistlim_sup(Prop_LogP)-TabHistlim_inf(Prop_LogP) )*REAL(N_P,8) ) + 1
         END IF
         IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,A,A,I3)') " iBinPeakflux of ",TRIM(TabSample_name(i_Sample))," = ", iBinPeakflux(i_Sample) 
      END IF
   END DO

   ! --- Spec alpha --- !
   IF( (alpha < TabHistlim_inf(Prop_alpha)) .OR. (alpha >= TabHistlim_sup(Prop_alpha) )) THEN
      iBina = -1
   ELSE
      iBina = INT( ( alpha-TabHistlim_inf(Prop_alpha) ) / ( TabHistlim_sup(Prop_alpha)-TabHistlim_inf(Prop_alpha) ) * REAL(N_spec_a,8) ) + 1
   END IF
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBina =", iBina
   
   ! --- Spec beta --- !
   IF( (beta < TabHistlim_inf(Prop_beta)) .OR. (beta >= TabHistlim_sup(Prop_beta) )) THEN
      iBinb = -1
   ELSE
      iBinb = INT( ( beta-TabHistlim_inf(Prop_beta) ) / ( TabHistlim_sup(Prop_beta)-TabHistlim_inf(Prop_beta) ) * REAL(N_spec_b,8) ) + 1
   END IF
   IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinb =", iBinb
   

   
   ! ----------- Fill the histograms ------------- !
   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample)) THEN        
         IF(iBinL     > 0) TabHistLogL(i_Sample,         iBinL) = TabHistLogL(i_Sample,         iBinL) + Prob_det(i_Sample)
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
      IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinKomm =", iBinKomm
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
      IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinPreece =", iBinPreece
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
      IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinStern =", iBinStern
   END IF
   
   ! Comparison with XRFHETE2
   IF(Constraint_Included(Constraint_XRFHETE2)) THEN
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
      IF (verbose == 1 .OR. verbose == 2) WRITE(*,'(A,I3)') " iBinEpGBM =", iBinEpGBM
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
      IF (iBinEpGBM > 0) TabHistEpGBM_Epobs(iBinEpGBM) = TabHistEpGBM_Epobs(iBinEpGBM) + Prob_det(Sample_GBM)
     ! PRINT*,"iBinEpGBM, Prob_det :", iBinEpGBM, Prob_det(Sample_GBM)
   END IF
 END SUBROUTINE Fill_Histogram_Chi2

 SUBROUTINE Normalize_Model()
   ! Calculate normalization coefficients
   REAL(8)  :: x1,x2
   INTEGER  :: i

   !dof = 0 
   Chi2 = 0.d0
   
   ! Kommers et al. 2000
   IF(Constraint_Included(Constraint_Kommers)) THEN
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
      
      R_GRB = x1 / x2
      IF (ISNAN(R_GRB)) NaNtest_hist = .TRUE.
      TabHistKomm_DRDP    = TabHistKomm_DRDP    * R_GRB      ! Normalize
      TabHistKomm_DRDPerr = TabHistKomm_DRDPerr * R_GRB
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Kommers data) :          R_GRB =", R_GRB, "                        ]"
   END IF
   
   ! Preece
   IF(Constraint_Included(Constraint_Preece)) THEN
      x1 = 0.d0
      x2 = 0.d0
      
      DO i=1, N_Preece
         !TabHistPreece_Eperr(i) = SQRT(TabHistPreece_Ep(i)) / (TabPreece_Ep(i)-TabPreece_Ep(i-1))        
         !TabHistPreece_Ep(i)    =      TabHistPreece_Ep(i)  / (TabPreece_Ep(i)-TabPreece_Ep(i-1)) 
         ! Note : error is SQRT(N) for the moment
         
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
   IF(Constraint_Included(Constraint_Stern)) THEN
      x1 = 0.d0 
      x2 = 0.d0
      
      DO i=1, N_Stern
         TabHistStern_P23err(i)     = SQRT(TabHistStern_P23_master(i)) / ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) )
         TabHistStern_P23_master(i) =      TabHistStern_P23_master(i)  / ( LOG10(TabStern_P23(i)/TabStern_P23(i-1)) )
         
         IF(TabHistStern_P23_master(i) > 0.d0) THEN 
            x1 = x1 + LOG10(TabHistStern_P23obs(i) / TabHistStern_P23_master(i)) / (TabHistStern_P23obserr(i)**2)
            x2 = x2 + 1.d0 / TabHistStern_P23obserr(i)**2
         ELSE
            Chi2(Constraint_Stern) = 1.d20
         END IF
      END DO
      k_Stern = 10.d0**(x1 / x2)
      !IF (run_mode == one_run) THEN
      IF (ISNAN(k_Stern)) NaNtest_hist = .TRUE.
      !END IF
      IF (ISNAN(k_Stern)) k_Stern = 1.d20
      TabHistStern_P23_master = TabHistStern_P23_master * k_Stern
      TabHistStern_P23err     = TabHistStern_P23err * k_Stern
      
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :        k_Stern =", k_Stern ,'                        ]'
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From BATSE23 (Stern data)   :         GRB/SN =", k_Stern*Nb_GRB/Collapse_rate, "                        ]"
      
   END IF
   
   ! XRF HETE2
   IF(Constraint_Included(Constraint_XRFHETE2)) THEN
      IF(NGRB_HETE2 > 0.d0) THEN
         Frac_XRFHETE2 = NGRB_XRFHETE2 / NGRB_HETE2
      ELSE
         Frac_XRFHETE2 = 0.d0
         Chi2(Constraint_XRFHETE2) = 1.d20
      END IF
      !IF (run_mode == one_run) THEN
      IF (ISNAN(Frac_XRFHETE2)) NaNtest_hist = .TRUE.
      !END IF
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From HETE2                  :  Frac_XRFHETE2 =", Frac_XRFHETE2,'                        ]'   
   END IF

   ! GBM Ep catalog (from Gruber et al. 2014) 
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      x1 = 0.d0
      x2 = 0.d0
      !PRINT *, "R_GRB = ", R_GRB
      
      DO i=1, N_EpGBM
         TabHistEpGBM_Epobserr(i) = SQRT(TabHistEpGBM_Epobs_master(i))
         ! Note : error is SQRT(N) because of Poisson statistics
         x1 = x1 + TabHistEpGBM_Epobs_master(i) * TabHistEpGBM_Epobsobs(i) / TabHistEpGBM_Epobsobserr(i)**2
         x2 = x2 + TabHistEpGBM_Epobs_master(i)**2 / TabHistEpGBM_Epobsobserr(i)**2
      END DO
      !PRINT*, " x1, x2 = ",x1,x2
      IF(x2 > 0.d0) THEN
         norm_EpGBM = x1 / x2
      ELSE
         !PRINT*, "GBM sample is empty..."
          norm_EpGBM = 0.d0
      END IF
      IF (ISNAN(norm_EpGBM)) NaNtest_hist = .TRUE.
      TabHistEpGBM_Epobs_master = TabHistEpGBM_Epobs_master * norm_EpGBM      ! Normalize
      TabHistEpGBM_Epobserr     = TabHistEpGBM_Epobserr     * norm_EpGBM
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.5,A)') "[                 From GBM (Gruber data)      :    normalization =", norm_EpGBM, "                      ]"
   END IF

   
   IF(run_mode == one_run)  WRITE(*,'(A)') "[ ------------------------------------------------------------------------------------------------- ]"
      
   !IF (run_mode == one_run) THEN
   IF (ISNAN(k_Stern)) NaNtest_hist = .TRUE.
   !END IF
   
   IF (NaNtest_hist) THEN
      WRITE(*,*) "[       R_GRB       k_Stern           k_Preece         Frac_XRFHETE2,     norm_EpGBM                ]"
      WRITE(*,*) "[       ",   R_GRB,      k_Stern,      k_Preece,       Frac_XRFHETE2, norm_EpGBM, "         ]"
      WRITE(*,*) "[                                           TabHistKomm                                             ]"
      WRITE(*,*)  TabHistKomm_DRDP
      WRITE(*,*) "[                                           TabHistStern                                            ]"
      WRITE(*,*)  TabHistStern_P23
      WRITE(*,*) "[                                          TabHistPreece                                            ]"
      WRITE(*,*)  TabHistPreece_Ep
      WRITE(*,*) "[                                          TabHistEpGBM                                             ]"
      WRITE(*,*)  TabHistEpGBM_Epobs
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
         IF(TabHistStern_P23_master(i) > 0.d0) THEN 
            Chi2(Constraint_Stern) = Chi2(Constraint_Stern) + ( LOG10(TabHistStern_P23_master(i)/TabHistStern_P23obs(i)) / TabHistStern_P23obserr(i) )**2
            
         ELSE
            IF(run_mode == one_run)  WRITE(*,'(A)') "[                 !!! WARNING : empty bin in Stern histogram, setting Chi2 to 10^20 !!!             ]"
            Chi2(Constraint_Stern) = 1.d20
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                     Chi2_Stern = ", Chi2(Constraint_Stern), "                                     ]"
   END IF
   
   IF(Constraint_Included(Constraint_XRFHETE2)) THEN
      Chi2(Constraint_XRFHETE2) = ( (Frac_XRFHETE2 - Frac_XRFHETE2obs) / sigma_XRFHETE2obs )**2
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                  Chi2_XRFHETE2 = ", Chi2(Constraint_XRFHETE2), "                                     ]"  
   END IF
   
   IF(Constraint_Included(Constraint_EpGBM)) THEN
      DO i=1, N_EpGBM
         IF(TabHistEpGBM_Epobsobserr(i) > 0.d0) THEN         
            Chi2(Constraint_EpGBM) = Chi2(Constraint_EpGBM) + ( (TabHistEpGBM_Epobs_master(i) - TabHistEpGBM_Epobsobs(i)) / TabHistEpGBM_Epobsobserr(i) )**2
         END IF
      END DO
      IF(run_mode == one_run)  WRITE(*,'(A,ES12.2,A)') "[                                     Chi2_EpGBM = ", Chi2(Constraint_EpGBM), "                                     ]"
   END IF
   Chi2(0) = SUM(Chi2)
   
 END SUBROUTINE Calculate_Chi2



 SUBROUTINE Save_Histograms()
   INTEGER :: i, i_Sample
   REAL(8) :: xtemp ! TEMPPPPPP
   ! --- Luminosity --- !

   DO i_Sample = 0, N_Samples
      IF(Sample_Included(i_Sample))THEN
         IF (Lsave(i_Sample)) THEN
            OPEN( UNIT=100, FILE=TRIM(LFile(i_Sample)) )

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
            OPEN( UNIT=110, FILE=TRIM(zFile(i_Sample)) )
            OPEN( UNIT=111, FILE=TRIM(zFile_cumul(i_Sample)) )

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
            OPEN( UNIT=120, FILE=TRIM(EpFile(i_Sample)) )
      
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
            OPEN( UNIT=120, FILE=TRIM(PFile(i_Sample)) )
            
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
            OPEN(UNIT=130, FILE=TRIM(SpecFile_a(i_Sample)) )
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
            OPEN(UNIT=140, FILE=TRIM(SpecFile_b(i_Sample)) )
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


 SUBROUTINE Save_Constraints()
   INTEGER :: i

   ! ---- P23 Kommers et al. 2000 ---- !

   IF (Constraint_save(Constraint_Kommers)) THEN
      OPEN(UNIT=200,  FILE=KommFile)
      OPEN(UNIT=2000, FILE=KommErrorFile)
      
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
      OPEN(UNIT=201,  FILE=PreeceFile)
      OPEN(UNIT=2001, FILE=PreeceErrorFile)

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
      OPEN(UNIT=202,  FILE=SternFile)
      OPEN(UNIT=2002, FILE=SternErrorFile)

      DO i=1, N_Stern
        
         !     1                   2                       3    
         !   [P23]         [P23 Hist model]         [P23 Stern hist]   
         WRITE(202, '(3ES12.5)') TabStern_P23(i-1), TabHistStern_P23_master(i), TabHistStern_P23obs(i)
         WRITE(202, '(3ES12.5)') TabStern_P23(i),   TabHistStern_P23_master(i), TabHistStern_P23obs(i)

         !     1                   2                        3                         4                           5     
         !   [P23]         [P23 hist model]        [P23 hist model err]       [P23 Stern hist obs]      [P23 Stern hist obs err]
         WRITE(2002, '(5ES12.5)') SQRT(TabStern_P23(i-1)*TabStern_P23(i)),&
              &                    TabHistStern_P23_master(i), TabHistStern_P23err(i), TabHistStern_P23obs(i), TabHistStern_P23obserr(i)
      END DO
      
      CLOSE(202)
      CLOSE(2002)
   END IF
   
   IF (Constraint_save(Constraint_EpGBM)) THEN
      OPEN(UNIT=203,  FILE=EpGBMFile)
      OPEN(UNIT=2003, FILE=EpGBMErrorFile)

      DO i=1, N_EpGBM
        
         !     1                   2                       3    
         ! [Log(Ep)]       [EpGBM model hist]       [EpGBM obs hist]   
         WRITE(203, '(3ES12.5)') TabEpGBM_LogEpobs(i-1), TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobsobs(i)
         WRITE(203, '(3ES12.5)') TabEpGBM_LogEpobs(i),   TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobsobs(i)

         !     1                   2                        3                         4                           5     
         !  [Log(Ep)]       [Ep hist model]        [Ep hist model err]       [Ep EpGBM hist obs]      [Ep EpGBM hist obs err]
         WRITE(2003, '(5ES12.5)') (TabEpGBM_LogEpobs(i-1)+TabEpGBM_LogEpobs(i))/2.d0,&
              &                   TabHistEpGBM_Epobs_master(i), TabHistEpGBM_Epobserr(i), TabHistEpGBM_Epobsobs(i), TabHistEpGBM_Epobsobserr(i)
      END DO
      
      CLOSE(203)
      CLOSE(2003)
   END IF

 END SUBROUTINE Save_Constraints

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
     ELSE
        WRITE(83,'(A)') "Run mode : ERROR"
     END IF
     WRITE(83, '(A,xES12.5)') 'Number of GRBs :', REAL(Nb_GRB,8)
     WRITE(83, '(A)') 'Parameters explored :'
     IF(lum_explore == 1) WRITE(83, '(A)') " - Luminosity"
     IF(redshift_explore == 1) WRITE(83, '(A)') " - Redshift"
     IF(Ep_explore == 1) WRITE(83, '(A)') " - Ep"
     IF(spec_explore == 1) WRITE(83, '(A)') " - Spectrum"
     
     WRITE(83, '(A)') 'Constraints used :'
     
     DO i_Constraint=1, N_Constraints
        IF(Constraint_Included(i_Constraint)) WRITE(83,'(A)') " - "//TRIM(TabConstraint_name(i_Constraint))
     END DO

     WRITE(83, '(A,xF6.3)') 'dof', REAL(dof,8)
     WRITE(83, '(A,xF6.3)') '1_sigma', delta_chi2_1
     WRITE(83, '(A,xF6.3)') '2_sigma', delta_chi2_2
     WRITE(83, '(A,xF6.3)') '3_sigma', delta_chi2_3
     CLOSE(83)

   END SUBROUTINE WRITE_INFO
   
 SUBROUTINE Reprise_mode()
   IF (reprise .EQV. .TRUE.) THEN
      OPEN(UNIT=43, FILE=TRIM(path)//'reprise.dat', FORM='unformatted') 
  
      ! Lecture du fichier
      Nb_lines = 0
      DO
        READ(43, err=999,end=999) TabParam_Lum, TabParam_z, TabParam_Spec, TabParam_Ep, &
            & Chi2, R_GRB, k_Stern, k_Preece, norm_EpGBM, Frac_XRFHETE2
         
         Nb_lines = Nb_lines + 1
         WRITE(*,*) "[                   reprise #", Nb_lines," : reduced Chi2 = ", Chi2(0)/REAL(dof,8),"              ]"
      END DO
999   WRITE(*,*) "[                   Number of lines read in reprise.dat : ",Nb_lines,"                              ]"
      
      CLOSE(43)           
   END IF
 END SUBROUTINE Reprise_mode
 
 
!!$
!!$ SUBROUTINE One_MCMC_run(Lmin_start, slope_start, Lmin_step, slope_step, Niter, MCMC_run_nb)
!!$   ! This will change to be more general, just testing it now
!!$   CHARACTER(Len=30)      :: MCMC_path='MCMC_path_run_'
!!$   REAL(8), INTENT(in)    :: Lmin_start, slope_start, Lmin_step, slope_step
!!$   INTEGER, INTENT(in)    :: Niter
!!$   INTEGER, INTENT(inout) :: MCMC_run_nb
!!$
!!$
!!$   IF(rank == master_proc) THEN
!!$      IF(MCMC_run_nb < 10) WRITE(MCMC_path(15:20),'(A,I1,A)') 'v',MCMC_run_nb,'.dat'
!!$      IF(MCMC_run_nb >= 10) WRITE(MCMC_path(15:21),'(A,I2,A)') 'v',MCMC_run_nb,'.dat'
!!$      WRITE(*,*) TRIM(path)//TRIM(MCMC_path)
!!$      OPEN(FILE=TRIM(path)//TRIM(MCMC_path), UNIT=57)
!!$   END IF
!!$
!!$   TabParam_Lum(Param_Lum_Lmin) = Lmin_start
!!$   TabParam_Lum(Param_Lum_Lmax) = Lmax_start
!!$   TabParam_Lum(Param_Lum_slope) = slope_start
!!$
!!$  
!!$   CALL MonteCarlo(hist_flag)
!!$
!!$   IF(rank == master_proc) Chi2_min = Chi2(0)
!!$   IF(rank == master_proc) WRITE(*,*) Chi2(0)
!!$   IF(rank == master_proc) WRITE(57,*) TabParam_Lum(Param_Lum_Lmin), TabParam_Lum(Param_Lum_Lmax), TabParam_Lum(Param_Lum_slope), Chi2_min
!!$   
!!$   DO ii = 1, Niter 
!!$      IF(rank == master_proc) THEN
!!$         WRITE(*,*) "Currently at : ", 100*REAL(ii,8)/REAL(Niter,8) ," %"
!!$         k_Stern    = 0.d0
!!$         k_Preece   = 0.d0
!!$         ! Save current spot
!!$         old_Lmin  = TabParam_Lum(Param_Lum_Lmin) 
!!$         old_Lmax  = TabParam_Lum(Param_Lum_Lmax)
!!$         old_slope = TabParam_Lum(Param_Lum_slope)
!!$
!!$         DO
!!$            ! make small step in parameter space
!!$            t = uniform()
!!$            TabParam_Lum(Param_Lum_Lmin)  = 10.**(LOG10(TabParam_Lum(Param_Lum_Lmin))  + Lmin_step  * (t-0.5d0))
!!$            IF (LOG10(TabParam_Lum(Param_Lum_Lmin)) > Lmin_min .AND. LOG10(TabParam_Lum(Param_Lum_Lmin)) < Lmin_max) EXIT
!!$         END DO
!!$        
!!$         DO
!!$            ! make small step in parameter space
!!$            t = uniform()
!!$            TabParam_Lum(Param_Lum_Lmax)  = 10.**(LOG10(TabParam_Lum(Param_Lum_Lmax))  + Lmax_step  * (t-0.5d0))
!!$            IF (LOG10(TabParam_Lum(Param_Lum_Lmax)) > Lmax_min .AND. LOG10(TabParam_Lum(Param_Lum_Lmax)) < Lmax_max) EXIT
!!$         END DO
!!$         
!!$         DO
!!$            t = uniform()
!!$            TabParam_Lum(Param_Lum_slope) = TabParam_Lum(Param_Lum_slope) + slope_step * (t-0.5d0)
!!$            IF (TabParam_Lum(Param_Lum_slope) > slope_min .AND. TabParam_Lum(Param_Lum_slope) < slope_max) EXIT
!!$         END DO
!!$         
!!$         WRITE(*,*) "Lmin, Lmax, slope :", TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_Lmax), TabParam_Lum(Param_Lum_slope) 
!!$      END IF
!!$
!!$      ! Broadcast Luminosity parameters for all procs
!!$      CALL MPI_BCAST(TabParam_Lum, NParam_Lum, MPI_REAL8, master_proc, MPI_COMM_WORLD, code)
!!$       !WRITE(*,*) " Rank", rank, " Param_lum :", TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_Lmax),TabParam_Lum(Param_Lum_slope)
!!$
!!$      CALL MonteCarlo(hist_flag)
!!$      
!!$      IF(rank == master_proc) THEN
!!$         IF(Chi2(0) <= Chi2_min) THEN
!!$            WRITE(*,*) "New best model :", Chi2(0)," (compared to :",Chi2_min,")"
!!$            Chi2_min = Chi2(0)
!!$            WRITE(57,*) TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_Lmax), TabParam_Lum(Param_Lum_slope), Chi2_min
!!$         ELSE 
!!$            t = uniform()
!!$            ! not sure I understand this line, need to go back and think about it
!!$            IF( t >= Chi2_min/Chi2(0)) THEN 
!!$               Chi2_min = Chi2(0)
!!$               WRITE(57,*) TabParam_Lum(Param_Lum_Lmin),TabParam_Lum(Param_Lum_Lmax), TabParam_Lum(Param_Lum_slope), Chi2_min
!!$            ELSE
!!$               TabParam_Lum(Param_Lum_Lmin)  =  old_Lmin
!!$               TabParam_Lum(Param_Lum_Lmax)  =  old_Lmax
!!$               TabParam_Lum(Param_Lum_slope) =  old_slope 
!!$            END IF
!!$         END IF
!!$      END IF
!!$      
!!$      
!!$   END DO
!!$
!!$  IF(rank == master_proc)THEN
!!$     WRITE(*,*) 'Finished MCMC run number : ', MCMC_run_nb
!!$     CLOSE(57)
!!$     MCMC_run_nb = MCMC_run_nb + 1
!!$  END IF
!!$
!!$ END SUBROUTINE One_MCMC_run
!!$
!!$


 SUBROUTINE Post_processing()

   

 END SUBROUTINE Post_processing

END PROGRAM Core

