MODULE Physics
  IMPLICIT None
  PRIVATE

  PUBLIC                     :: InitPhysics ! Module initialization
  PUBLIC                     :: CalcRandom  ! Random number 0-1

  INTEGER, PARAMETER, Public :: U = KIND(1.0D0)

  REAL(U), PARAMETER, Public :: Pi           = 3.141592653589793238462643383279502884197_U ! sans unite
  REAL(U), PARAMETER, Public :: cLight       = REAL(2.99792458D10,U)    ! cm/s
  REAL(U), PARAMETER, Public :: hPlanck      = REAL(6.6260755D-27,U)    ! erg.s
  REAL(U), PARAMETER, Public :: hbar         = REAL(1.05457266D-27,U)   ! erg.s (h/2pi)
  REAL(U), PARAMETER, Public :: eElectron    = REAL(4.8032068D-10,U)    ! esu
  REAL(U), PARAMETER, Public :: mElectron    = REAL(9.1093897D-28,U)    ! g (i.e. 510.99906 keV/c2)
  REAL(U), PARAMETER, Public :: mProton      = REAL(1.6726231D-24,U)    ! g (i.e. 938.27231 MeV/c2)
  REAL(U), PARAMETER, Public :: uamu         = REAL(1.6605402D-24,U)    ! g (unite de mass atomique)
  REAL(U), PARAMETER, Public :: SigmaThomson = REAL(6.6524616D-25,U)    ! cm2
  REAL(U), PARAMETER, Public :: GNewton      = REAL(6.67259D-8,U)       ! cgs
  REAL(U), PARAMETER, Public :: NAvogadro    = REAL(6.0221367D23,U)     ! mol
  REAL(U), PARAMETER, Public :: kBoltzmann   = REAL(1.380658D-16)       ! erg/K
  REAL(U), PARAMETER, Public :: SigmaStefan  = REAL(5.67051D-5,U)       ! erg/s/cm2/K4 
  REAL(U), PARAMETER, Public :: aStefan      = REAL(7.56591D-15,U)      ! erg/cm3/K4 (recalcule = 4Sigma/c)

  REAL(U), PARAMETER, Public :: au           = REAL(1.4959787066D13,U)  ! cm
  REAL(U), PARAMETER, Public :: pc           = REAL(3.0856775807D18,U)  ! cm
  REAL(U), PARAMETER, Public :: kpc          = REAL(3.0856775807D21,U)  ! cm
  REAL(U), PARAMETER, Public :: Mpc          = REAL(3.0856775807D24,U)  ! cm
  REAL(U), PARAMETER, Public :: Gpc          = REAL(3.0856775807D27,U)  ! cm
  REAL(U), PARAMETER, Public :: Rsun         = REAL(6.96D10,U)          ! cm
  REAL(U), PARAMETER, Public :: Msun         = REAL(1.98892D33,U)       ! g
  REAL(U), PARAMETER, Public :: Lsun         = REAL(3.846D33,U)         ! erg/s
  REAL(U), PARAMETER, Public :: day          = 86400._U                 ! s       (24 h)
  REAL(U), PARAMETER, Public :: yr           = REAL(3.156D7,U)          ! s
  REAL(U), PARAMETER, Public :: ly           = REAL(0.9461D18,U)        ! cm

  REAL(U), PARAMETER, Public :: eV           = REAL(1.60217733D-12,U)   ! erg
  REAL(U), PARAMETER, Public :: keV          = REAL(1.60217733D-9,U)    ! erg
  REAL(U), PARAMETER, Public :: MeV          = REAL(1.60217733D-6,U)    ! erg
  REAL(U), PARAMETER, Public :: GeV          = REAL(1.60217733D-3,U)    ! erg
  REAL(U), PARAMETER, Public :: TeV          = REAL(1.60217733D0,U)     ! erg
  
  REAL(U), PARAMETER, Public :: erg          = REAL(6.24150636D5,U)     ! MeV


  REAL(U), PARAMETER, Public :: MeVperc2     = REAL(1.78266270D-27,U)   ! g

  REAL(U), PARAMETER, Public :: Jy           = REAL(1.D-23,U)           ! erg/cm2/s/Hz
  REAL(U), PARAMETER, Public :: mJy          = REAL(1.D-26,U)           ! erg/cm2/s/Hz
  REAL(U), PARAMETER, Public :: microJy      = REAL(1.D-29,U)           ! erg/cm2/s/Hz

  REAL(U), PARAMETER, Public :: micron       = REAL(1.D-4,U)            ! cm

! "Standard" cosmology :
! ----------------------

  REAL(U), Public            :: OmegaM,OmegaL,OmegaK,h

! Standard filters :
! ------------------

  REAL(U), Public            :: LambdaU,nuU,DeltaLambdaU,DeltaNuU,FlambdaRefU,FnuRefU
  REAL(U), Public            :: LambdaB,nuB,DeltaLambdaB,DeltaNuB,FlambdaRefB,FnuRefB
  REAL(U), Public            :: LambdaV,nuV,DeltaLambdaV,DeltaNuV,FlambdaRefV,FnuRefV
  REAL(U), Public            :: LambdaR,nuR,DeltaLambdaR,DeltaNuR,FlambdaRefR,FnuRefR
  REAL(U), Public            :: LambdaI,nuI,DeltaLambdaI,DeltaNuI,FlambdaRefI,FnuRefI
  REAL(U), Public            :: LambdaJ,nuJ,DeltaLambdaJ,DeltaNuJ,FlambdaRefJ,FnuRefJ
  REAL(U), Public            :: LambdaH,nuH,DeltaLambdaH,DeltaNuH,FlambdaRefH,FnuRefH
  REAL(U), Public            :: LambdaK,nuK,DeltaLambdaK,DeltaNuK,FlambdaRefK,FnuRefK

CONTAINS

SUBROUTINE InitPhysics(Show)
  LOGICAL, INTENT(in) :: Show
  CALL InitStandardCosmology(Show)
  CALL InitStandardFilters(Show)
END SUBROUTINE InitPhysics

SUBROUTINE InitStandardCosmology(Show)
  LOGICAL, INTENT(in) :: Show
  ! Old values WMAP
  OmegaM = 0.27_U
  OmegaL = 0.73_U
  OmegaK = 0.00_U
  h      = 0.71_U
  ! New values PLANCK
  OmegaM = 0.286_U
  OmegaL = 0.714_U
  OmegaK = 0.000_U
  h      = 0.696_U

  IF (Show) THEN
    WRITE(*,'(4(A,1ES9.2))') "Physics:InitStandardCosmology:   Omega(matter) = ",OmegaM,&
                                                        &" Omega(Lambda) = ",OmegaL,&
                                                        &" Omega(total) = ",(1._U-OmegaK),&
                                                        &" h = ",h
  END IF
END SUBROUTINE InitStandardCosmology

SUBROUTINE InitStandardFilters(Show)
  LOGICAL, INTENT(in) :: Show
  LambdaU      = 0.365_U*micron ! cm
  LambdaB      = 0.44_U *micron ! cm
  LambdaV      = 0.55_U *micron ! cm
  LambdaR      = 0.70_U *micron ! cm
  LambdaI      = 0.90_U *micron ! cm
  LambdaJ      = 1.25_U *micron ! cm
  LambdaH      = 1.65_U *micron ! cm
  LambdaK      = 2.2_U  *micron ! cm

  DeltaLambdaU = 0.068_U*micron ! cm
  DeltaLambdaB = 0.098_U*micron ! cm
  DeltaLambdaV = 0.089_U*micron ! cm
  DeltaLambdaR = 0.22_U *micron ! cm
  DeltaLambdaI = 0.24_U *micron ! cm
  DeltaLambdaJ = 0.3_U  *micron ! cm
  DeltaLambdaH = 0.4_U  *micron ! cm
  DeltaLambdaK = 0.6_U  *micron ! cm

  nuU          = cLight/LambdaU ! Hz
  nuB          = cLight/LambdaB ! Hz
  nuV          = cLight/LambdaV ! Hz
  nuR          = cLight/LambdaR ! Hz
  nuI          = cLight/LambdaI ! Hz
  nuJ          = cLight/LambdaJ ! Hz
  nuH          = cLight/LambdaH ! Hz
  nuK          = cLight/LambdaK ! Hz

  DeltanuU     = cLight*DeltaLambdaU/(LambdaU-DeltaLambdaU/2._U)/(LambdaU+DeltaLambdaU/2._U) ! Hz
  DeltanuB     = cLight*DeltaLambdaB/(LambdaB-DeltaLambdaB/2._U)/(LambdaB+DeltaLambdaB/2._U) ! Hz 
  DeltanuV     = cLight*DeltaLambdaV/(LambdaV-DeltaLambdaV/2._U)/(LambdaV+DeltaLambdaV/2._U) ! Hz 
  DeltanuR     = cLight*DeltaLambdaR/(LambdaR-DeltaLambdaR/2._U)/(LambdaR+DeltaLambdaR/2._U) ! Hz 
  DeltanuI     = cLight*DeltaLambdaI/(LambdaI-DeltaLambdaI/2._U)/(LambdaI+DeltaLambdaI/2._U) ! Hz 
  DeltanuJ     = cLight*DeltaLambdaJ/(LambdaJ-DeltaLambdaJ/2._U)/(LambdaJ+DeltaLambdaJ/2._U) ! Hz 
  DeltanuH     = cLight*DeltaLambdaH/(LambdaH-DeltaLambdaH/2._U)/(LambdaH+DeltaLambdaH/2._U) ! Hz 
  DeltanuK     = cLight*DeltaLambdaK/(LambdaK-DeltaLambdaK/2._U)/(LambdaK+DeltaLambdaK/2._U) ! Hz 

  FlambdaRefU  = 0.427_U           ! erg/cm2/s/cm
  FlambdaRefB  = 0.661_U           ! erg/cm2/s/cm
  FlambdaRefV  = 0.364_U           ! erg/cm2/s/cm
  FlambdaRefR  = 0.174_U           ! erg/cm2/s/cm
  FlambdaRefI  = 0.0832_U          ! erg/cm2/s/cm
  FlambdaRefJ  = 0.0318_U          ! erg/cm2/s/cm
  FlambdaRefH  = 0.0118_U          ! erg/cm2/s/cm
  FlambdaRefK  = 0.00417_U         ! erg/cm2/s/cm

  FnuRefU      = REAL(1.90D-20 ,U) ! erg/cm2/s/Hz
  FnuRefB      = REAL(4.27D-20 ,U) ! erg/cm2/s/Hz
  FnuRefV      = REAL(3.67D-20 ,U) ! erg/cm2/s/Hz
  FnuRefR      = REAL(2.84D-20 ,U) ! erg/cm2/s/Hz
  FnuRefI      = REAL(2.25D-20 ,U) ! erg/cm2/s/Hz
  FnuRefJ      = REAL(1.65D-20 ,U) ! erg/cm2/s/Hz
  FnuRefH      = REAL(1.07D-20 ,U) ! erg/cm2/s/Hz
  FnuRefK      = REAL(0.673D-20,U) ! erg/cm2/s/Hz
  
  IF (Show) THEN
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         U band nu = ",nuU," Hz  Delta nu = ",DeltanuU," Hz --> E = ",(hPlanck*nuU/eV)," eV Delta E = ",(DeltanuU*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         B band nu = ",nuB," Hz  Delta nu = ",DeltanuB," Hz --> E = ",(hPlanck*nuB/eV)," eV Delta E = ",(DeltanuB*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         V band nu = ",nuV," Hz  Delta nu = ",DeltanuV," Hz --> E = ",(hPlanck*nuV/eV)," eV Delta E = ",(DeltanuV*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         R band nu = ",nuR," Hz  Delta nu = ",DeltanuR," Hz --> E = ",(hPlanck*nuR/eV)," eV Delta E = ",(DeltanuR*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         I band nu = ",nuI," Hz  Delta nu = ",DeltanuI," Hz --> E = ",(hPlanck*nuI/eV)," eV Delta E = ",(DeltanuI*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         J band nu = ",nuJ," Hz  Delta nu = ",DeltanuJ," Hz --> E = ",(hPlanck*nuJ/eV)," eV Delta E = ",(DeltanuJ*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         H band nu = ",nuH," Hz  Delta nu = ",DeltanuH," Hz --> E = ",(hPlanck*nuH/eV)," eV Delta E = ",(DeltanuH*hPlanck/eV)," eV"
    WRITE(*,'(4(A,1ES9.2),A)') "Physics:InitStandardfilters:         K band nu = ",nuK," Hz  Delta nu = ",DeltanuK," Hz --> E = ",(hPlanck*nuK/eV)," eV Delta E = ",(DeltanuK*hPlanck/eV)," eV"
  END IF
END SUBROUTINE InitStandardFilters

FUNCTION CalcRandom()
  REAL(U) :: CalcRandom
  REAL(U), DIMENSION(1:1) :: Tab

  CALL Random_Number(Tab)
  CalcRandom = Tab(1)
END FUNCTION

END MODULE Physics
