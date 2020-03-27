MODULE Cosmology
  USE Physics
  IMPLICIT None
  PRIVATE
  
  PUBLIC                           :: E                                 ! E(z) = SQRT(OmegaM*(1+z)**3 + OmegaL) for OmegaK = 0
  PUBLIC                           :: InitCosmology                     ! Creates precise z, D_L and dVdz tables
  
  INTEGER, PARAMETER               :: zlim = 100                        ! upper z used to create tables
  INTEGER, PARAMETER               :: NN = 100000                       ! Size of Table for Luminosity Distance
  REAL(8), DIMENSION(0:NN), Public :: TabD_L                            ! Table for Luminosity Distance [Mpc]
  REAL(8), DIMENSION(0:NN), Public :: Tabprecisez                       ! Table for z (very precise)
  REAL(8), DIMENSION(0:NN), Public :: TabdVdz                           ! Table for dVdz [Mpc^3]
  REAL(8), DIMENSION(0:NN), Public :: TabVz                             ! Table for V [Mpc^3]
  INTEGER                          :: i                                 ! Loop indice
  CHARACTER(*), PARAMETER          :: Format_IC1 ='(A,A25,A20,A)'       ! Format for InitCosmology output on screen
  CHARACTER(*), PARAMETER          :: Format_IC2 ='(A,A45,A20,A)'       ! Format for InitCosmology output on screen
  CHARACTER(*), PARAMETER          :: precisezfile = "../data/cosmology/Tabprecisez.dat"  ! Name of precise z file
  CHARACTER(*), PARAMETER          :: D_Lfile      = "../data/cosmology/TabD_L.dat"       ! Name of D_L file
  CHARACTER(*), PARAMETER          :: dVdzfile     = "../data/cosmology/TabdVdz.dat"      ! Name of dVdz file
  CHARACTER(*), PARAMETER          :: Vzfile       = "../data/cosmology/TabVz.dat"        ! Name of Vz file      
  LOGICAL                          :: D_Lexist                          ! Does the D_L file exist ? yes = .TRUE.
  LOGICAL                          :: precisezexist                     ! Does the precise z file exist ?
  LOGICAL                          :: dVdzexist                         ! Does the dVdz file exist ?
  LOGICAL                          :: Vzexist                           ! Does the Vz file exist ?

CONTAINS

SUBROUTINE InitCosmology(SHOW_READ)
  LOGICAL, intent(in) :: SHOW_READ
  INQUIRE(FILE = precisezfile, EXIST = precisezexist)
  INQUIRE(FILE = D_Lfile,      EXIST = D_Lexist)
  INQUIRE(FILE = dVdzfile,     EXIST = dVdzexist)
  INQUIRE(FILE = Vzfile,       EXIST = Vzexist)

  IF(SHOW_READ) WRITE(*,'(A)') " "
  IF(SHOW_READ) WRITE(*,'(A)') " ====================== Initializing cosmology ====================== "
  IF(SHOW_READ) WRITE(*,'(A)') " "

  OmegaM = 0.27d0
  OmegaL = 0.73d0
  h      = 0.71d0

  IF(SHOW_READ) WRITE(*,'(A)')            "[                            Parameters                             ]"
  IF(SHOW_READ) WRITE(*,'(A)')            "[ ----------------------------------------------------------------- ]"
  IF(SHOW_READ) WRITE(*,'(A,A10,F5.3,A)') "[                         ", "OmegaM = ", OmegaM,   "                           ]"
  IF(SHOW_READ) WRITE(*,'(A,A10,F5.3,A)') "[                         ", "OmegaL = ", OmegaL,   "                           ]"
  IF(SHOW_READ) WRITE(*,'(A,A10,F5.1,A)') "[                         ", "H0     = ", 100.d0*h, " (km/s/Mpc)                ]"
  IF(SHOW_READ) WRITE(*,'(A)')            "[ ----------------------------------------------------------------- ]"
  IF(SHOW_READ) WRITE(*,'(A)') " "
  IF(SHOW_READ) WRITE(*,'(A)')            "[                              Files                                ]"
  IF(SHOW_READ) WRITE(*,'(A)')            "[ ----------------------------------------------------------------- ]"

  ! Create precise z grid
  IF(precisezexist) THEN 
     IF(SHOW_READ) WRITE(*,Format_IC1) "[         ", "precise z file exists : ", precisezfile, "             ]"
     OPEN(UNIT=1000, FILE = precisezfile, FORM="unformatted")
     READ(1000) Tabprecisez
     CLOSE(UNIT=1000)
  ELSE
     IF(SHOW_READ) WRITE(*,Format_IC2) "[ ", "precise z file does not exist, creating : ", precisezfile, " ]"
     Tabprecisez(0) = 0.d0
     DO i=1, NN
        Tabprecisez(i) = REAL(zlim,8) * REAL(i,8)/REAL(NN,8)
     END DO
     OPEN(UNIT=1000, FILE = precisezfile, FORM="unformatted")
     WRITE(1000) Tabprecisez
     CLOSE(UNIT=1000)
  END IF

  ! Create grid for Luminosity Distance
  IF(D_Lexist) THEN
     IF(SHOW_READ) WRITE(*,Format_IC1) "[         ", "D_L file exists : ", D_Lfile, "             ]"
     OPEN(UNIT=1001, FILE = D_Lfile, FORM="unformatted")
     READ(1001) TabD_L
     CLOSE(UNIT=1001)
  ELSE
     IF(SHOW_READ) WRITE(*,Format_IC2) "[ ", "D_L file does not exist, creating : ", D_Lfile, " ]"
     TabD_L(0) = 0.d0
     DO i=1, NN
        TabD_L(i) = (1.d0 + Tabprecisez(i))/(1.d0+Tabprecisez(i-1)) * TabD_L(i-1) +&
                  & (1.d0 + Tabprecisez(i)) * cLight/(10**7*h) * &
                  & 0.5d0 * (REAL(zlim,8)/REAL(NN,8)) * (1.d0/E(Tabprecisez(i)) + 1.d0/E(Tabprecisez(i-1))) 
     END DO
     OPEN(UNIT=1001, FILE = D_Lfile, FORM="unformatted")
     WRITE(1001) TabD_L
     CLOSE(UNIT=1001)
  END IF

  ! Create grid for dVdz
  IF(dVdzexist) THEN
     IF(SHOW_READ) WRITE(*,Format_IC1) "[         ", "dVdz file exists : ", dVdzfile, "             ]"
     OPEN(UNIT=1002, FILE = dVdzfile, FORM="unformatted")
     READ(1002) TabdVdz
     CLOSE(UNIT=1002)
  ELSE
     IF(SHOW_READ) WRITE(*,Format_IC2) "[ ", "dVdz file does not exist, creating : ", dVdzfile, " ]"
     TabdVdz(0) = 0.d0
     DO i=1, NN
        TabdVdz(i) = 4.d0 * Pi * cLight/(10**7*h) * &
                   & TabD_L(i)**2/(1.d0 + Tabprecisez(i))**2 * 1.d0/E(Tabprecisez(i))
     END DO
     OPEN(UNIT=1002, FILE = dVdzfile, FORM="unformatted")
     WRITE(1002) TabdVdz
     CLOSE(UNIT=1002)
  END IF

 ! Create grid for Vz
  IF(Vzexist) THEN
     IF(SHOW_READ) WRITE(*,Format_IC1) "[         ", "Vz file exists : ", Vzfile, "             ]"
     OPEN(UNIT=1003, FILE = Vzfile, FORM="unformatted")
     READ(1003) TabVz
     CLOSE(UNIT=1003)
  ELSE
     IF(SHOW_READ) WRITE(*,Format_IC2) "[ ", "Vz file does not exist, creating : ", Vzfile, " ]"
     TabVz(0) = 0.d0
     DO i=1, NN
        TabVz(i) = TabVz(i-1) + 0.5d0 * (REAL(zlim,8)/REAL(NN,8)) *&
                 & (TabdVdz(i)/(1.d0+Tabprecisez(i)) + TabdVdz(i-1)/(1.d0+Tabprecisez(i-1)))
     END DO
     OPEN(UNIT=1003, FILE = Vzfile, FORM="unformatted")
     WRITE(1003) TabVz
     CLOSE(UNIT=1003)
  END IF
  IF(SHOW_READ) WRITE(*,'(A)') "[ ----------------------------------------------------------------- ]"
  IF(SHOW_READ) WRITE(*,'(A)') " "
  IF(SHOW_READ) WRITE(*,'(A)') " ================== Done initializing cosmology ====================="
  IF(SHOW_READ) WRITE(*,'(A)') " "

END SUBROUTINE InitCosmology



  
! Used by LuminDist
REAL(8) FUNCTION E(z) ! True for OmegaK = 0
  REAL(8), INTENT(in) :: z
  E = SQRT(OmegaM*(1+z)**3 + OmegaL)
END FUNCTION E





END MODULE Cosmology
