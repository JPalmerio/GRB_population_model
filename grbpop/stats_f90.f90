MODULE stats_f90
  IMPLICIT None
  PRIVATE
  PUBLIC :: ERF
  PUBLIC :: ERFC
  PUBLIC :: GammaP
  PUBLIC :: GammaQ
  PUBLIC :: GammaSeries
  PUBLIC :: GammaContinuedFraction
  PUBLIC :: GammaLn

CONTAINS

   FUNCTION ERF(x)
     REAL(8), INTENT(in) :: x
     REAL(8)             :: ERF
     IF (x<0.d0) THEN
        ERF = -GammaP(0.5d0,(x**2.d0))
     ELSE
        ERF = GammaP(0.5d0,(x**2.d0))
     END IF
   END FUNCTION ERF
   
   FUNCTION ERFC(x)
     REAL(8), INTENT(in) :: x
     REAL(8) :: ERFC
     IF (x<0.d0) THEN
        ERFC = 1.d0+GammaP(0.5d0,(x**2.d0))
     ELSE
        ERFC = GammaQ(0.5d0,(x**2.d0))
     END IF
   END FUNCTION ERFC
   
   FUNCTION GammaP(a,x)
     REAL(8), INTENT(in) :: a,x
     REAL(8)             :: GammaP
     REAL(8)             :: Gamma,gln
     
     IF ((x<0.d0).or.(a<=0.d0)) THEN
        WRITE(*,*) "x=",x,"a=",a
        GammaP = 0.5d0
        WRITE(*,*) "Bad arguments in GammaP !!!"
        RETURN
        !    STOP "Bad arguments in GammaP"
     END IF
     IF (x<(a+1.d0)) THEN
        CALL GammaSeries(Gamma,a,x,gln)
        GammaP = Gamma
     ELSE
        CALL GammaContinuedFraction(Gamma,a,x,gln)
        GammaP = 1.d0-Gamma
     END IF
   END FUNCTION GammaP
   
   FUNCTION GammaQ(a,x)
     REAL(8), INTENT(in) :: a,x
     REAL(8)             :: GammaQ
     REAL(8)             :: Gamma,gln
     
     IF ((x<0.d0).or.(a<=0.d0)) THEN
        WRITE(*,*) "x=",x,"a=",a
        GammaQ = 0.5d0
        WRITE(*,*) "Bad arguments in GammaQ"
        RETURN
        !STOP "Bad arguments in GammaQ"
     END IF
     IF (x<(a+1.d0)) THEN
        CALL GammaSeries(Gamma,a,x,gln)
        GammaQ = 1.d0-Gamma
     ELSE
        CALL GammaContinuedFraction(Gamma,a,x,gln)
        GammaQ = Gamma
     END IF
   END FUNCTION GammaQ
   
   SUBROUTINE GammaSeries(Gamma,a,x,gln)
     REAL(8), INTENT(out) :: Gamma
     REAL(8), INTENT(in)  :: a,x
     REAL(8), INTENT(out) :: gln
     
     INTEGER, PARAMETER   :: iterMax = 100
     REAL(8), PARAMETER   :: eps = REAL(3.D-7,8)
     
     INTEGER              :: n
     REAL(8)              :: ap,del,sum
     
     gln = GammaLn(a)
     IF (x <= 0.d0) THEN
        IF (x < 0.d0) STOP "x < 0 in GammaSeries !!!"
        Gamma = 0.d0
        RETURN
     END IF
     ap  = a
     sum = 1.d0/a
     del = sum
     DO n=1, itermax
        ap  = ap  + 1.d0
        del = del * x/ap
        sum = sum + del
        IF (ABS(del) < ABS(sum)*Eps) EXIT
     END DO
     IF (ABS(del) >= ABS(sum)*Eps) STOP "a too large, itermax too small in GammaSeries"
     
     Gamma = sum*EXP(-x+a*LOG(x)-gln)
   END SUBROUTINE GammaSeries
   
   SUBROUTINE GammaContinuedFraction(Gamma,a,x,gln)
     REAL(8), INTENT(in)  :: a,x
     REAL(8), INTENT(out) :: Gamma,gln
     
     INTEGER, PARAMETER   :: iterMax = 100
     REAL(8), PARAMETER   :: eps     = REAL(3.D-7,8)
     REAL(8), PARAMETER   :: FPMin   = REAL(1.D-30,8)
     
     INTEGER :: i
     REAL(8) :: an,b,c,d,del,h
     
     gln = GammaLn(a)
     
     b   = x  + 1.d0-a
     c   = 1.d0/FPMIN
     d   = 1.d0/b
     h   = d
     
     DO i=1, itermax
        an = -i*(i-a)
        b  = b + 2.d0
        d  = an*d+b

        IF (ABS(d) < FPMin) d = FPMin
        c  = b+an/c
        
        IF (ABS(c) < FPMin) c = FPMIN
        d  = 1.d0/d
        del=d*c
        h  =h*del
        
        IF (ABS(del-1.d0) < Eps) EXIT
     END DO
     IF (ABS(del-1.d0) >= Eps) STOP "a too large, ITMAX too small in GammaContinuedFraction"
     
     Gamma = EXP(-x+a*LOG(x)-gln)*h
   END SUBROUTINE GammaContinuedFraction
   
   FUNCTION GammaLn(xx)
     REAL(8), INTENT(in) :: xx
     REAL(8) :: GammaLn
     
     INTEGER :: j
     REAL(8) :: x,y,ser,stp,tmp
     REAL(8), DIMENSION(1:6) :: cof
     
     cof(1) = REAL(76.18009172947146d0,8)
     cof(2) = REAL(-86.50532032941677d0,8)
     cof(3) = REAL(24.01409824083091d0,8)
     cof(4) = REAL(-1.231739572450155d0,8)
     cof(5) = REAL(.1208650973866179d-2,8)
     cof(6) = REAL(-.5395239384953d-5,8)
     
     stp    = REAL(2.5066282746310005d0,8)
     
     x   = xx
     y   = x
     tmp = x+5.5d0
     tmp = (x+0.5d0)*LOG(tmp)-tmp
     ser = REAL(1.000000000190015d0,8)
     DO j=1, 6
        y   = y + 1.d0
        ser = ser+cof(j)/y
     END DO
     GammaLn = tmp+LOG(stp*ser/x)
   END FUNCTION GammaLn

END MODULE stats_f90
