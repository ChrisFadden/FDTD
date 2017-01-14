!> Prints Ez output to a .csv file
MODULE PrintField
  USE Grid, ONLY : SizeX, SizeY, Ez
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0);
  REAL(dp) :: Min;
  REAL(dp) :: Max;
  
  CONTAINS

    !> Prints output
    !! \param Ez field
    
    SUBROUTINE Print_Initialize()
      IMPLICIT NONE
      Min = 0.0
      Max = 0.0
    END SUBROUTINE

    SUBROUTINE PrintEz(file)
      IMPLICIT NONE
      INTEGER :: i, j
      CHARACTER (LEN=*), PARAMETER :: dir="./build/output/";
      CHARACTER(LEN=*), PARAMETER :: ext=".csv";
      CHARACTER(LEN=3), INTENT(IN) :: file;
      
      CHARACTER(LEN=22) :: filename;
      filename = dir//file//ext;

      OPEN(UNIT=1, FILE=filename, FORM="FORMATTED", &
                   STATUS="REPLACE", ACTION="WRITE");
      
      DO i = 1, SizeX
        DO j = 1, SizeY

          IF(Ez(i,j) > Max) THEN
            Max = Ez(i,j)
          ELSE IF(Ez(i,j) < Min) THEN
            Min = Ez(i,j)
          END IF
          WRITE(UNIT=1, FMT="(F12.6)",ADVANCE="NO") Ez(i,j);
          WRITE(UNIT=1, FMT="(A2)",ADVANCE="NO") ", ";
        END DO
          WRITE(UNIT=1,FMT=*) 
      END DO
             
      CLOSE(UNIT=1);

    END SUBROUTINE
    
    SUBROUTINE Print_Limits()
      IMPLICIT NONE

      OPEN(UNIT=1, FILE="./build/output/limits.txt", FORM="FORMATTED", &
                   STATUS="REPLACE", ACTION="WRITE");
      
      WRITE(UNIT=1, FMT=*) Min;
      WRITE(UNIT=1, FMT=*) Max;
  
    END SUBROUTINE

END MODULE
