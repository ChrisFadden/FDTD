PROGRAM Main    
  USE Parameters 
  USE Grid
  USE Update
  USE UPML
  USE PrintField
  USE MPI_utils
  USE OMP_LIB 
  IMPLICIT NONE 
  INCLUDE "fftw3.f"

  INTEGER :: t, mystart, myend;
  CALL MPI_UTILS_INITIALIZE();

  CALL GRID_Initialize();
  CALL UPML_Initialize();
  
  DO t = 1,MaxTime
    CALL UpdateH();
    CALL UpdateE(t);
  END DO
  
  CALL PrintHz("800");
 
  CALL MPI_UTILS_FINALIZE();
END PROGRAM
