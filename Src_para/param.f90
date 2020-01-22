MODULE ParameterFile
  USE Physics
  IMPLICIT None
  PRIVATE
  PUBLIC :: ReadLine    ! Read a line in a parameter file
  PUBLIC :: ReadReal    ! Read a real number in a parameter file
  PUBLIC :: ReadInteger ! Read an integer in a parameter file

  LOGICAL :: ShowSkip = .FALSE.

CONTAINS

SUBROUTINE ReadReal(n,x)
  INTEGER, INTENT(in)  :: n    ! file unit
  REAL(U), INTENT(out) :: x    ! real to be read
  CHARACTER(Len=255)   :: Line
  CALL ReadLine(n,Line)
  READ(Line,'(1ES12.5)') x
END SUBROUTINE ReadReal

SUBROUTINE ReadInteger(n,i)
  INTEGER, INTENT(in)  :: n    ! file unit
  INTEGER, INTENT(out) :: i    ! integer to be read
  CHARACTER(Len=255)   :: Line
  CALL ReadLine(n,Line)
  READ(Line,'(1I8)') i
END SUBROUTINE ReadInteger

! ---------------------------------------------------------------------------------------------
! Readline : skip all lines starting with '#' in a parameter file and read the next 'real' line
! ---------------------------------------------------------------------------------------------

SUBROUTINE ReadLine(n,line)
  INTEGER, INTENT(in)             :: n    ! file unit
  CHARACTER(Len=255), INTENT(out) :: line ! line read

  DO
     READ(n,'(A)',err=1,end=1) Line
     IF (Line(1:1)/='#') THEN 
        EXIT
     ELSE
        IF (ShowSkip) WRITE(*,'(A,A)') "Skip line ",TRIM(Line)
     END IF
  END DO
  RETURN
  
1 WRITE(*,'(A,I4)') "ERROR (Readline) ; unit = ",n

  STOP
END SUBROUTINE ReadLine

END MODULE ParameterFile
