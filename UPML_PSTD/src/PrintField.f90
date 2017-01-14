!> Prints Ez output to a .csv file
MODULE PrintField
  USE Grid, ONLY: Hz
  USE Parameters, ONLY: SizeX, SizeY 
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0);

  CONTAINS

    !> Prints output
    !! \param Ez field

    SUBROUTINE PrintHz(file)
      IMPLICIT NONE
      INTEGER :: i, j
      CHARACTER (LEN=*), PARAMETER :: dir="./build/output/";
      CHARACTER(LEN=*), PARAMETER :: ext=".csv";
      CHARACTER(LEN=3), INTENT(IN) :: file;
      
      CHARACTER(LEN=22) :: filename;
      filename = dir//file//ext;
      
      OPEN(UNIT=1, FILE=filename, FORM="FORMATTED", &
                   STATUS="REPLACE", ACTION="WRITE");
      
      DO i = 1, SizeX-1
        DO j = 1, SizeY-1
          WRITE(UNIT=1, FMT="(ES14.7)",ADVANCE="NO") Hz(i,j);
          WRITE(UNIT=1, FMT="(A2)",ADVANCE="NO") ", ";
        END DO
          WRITE(UNIT=1,FMT=*) 
      END DO
             
      CLOSE(UNIT=1);

    END SUBROUTINE

END MODULE
