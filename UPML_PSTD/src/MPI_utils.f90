MODULE MPI_UTILS
  USE MPI 
  IMPLICIT NONE
  PRIVATE 
  
  !*******************
  ! Private Variables
  !*******************
  INTEGER :: rank;
  INTEGER :: numProc;
  
  !*****************
  ! Public Functions
  !*****************
  PUBLIC :: MPI_UTILS_INITIALIZE
  PUBLIC :: MPI_UTILS_FINALIZE
  PUBLIC :: GetIndex

  CONTAINS
    SUBROUTINE MPI_UTILS_INITIALIZE()
      IMPLICIT NONE 
      INTEGER :: ierr;

      CALL MPI_Init(ierr);
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr);
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProc,ierr);

    END SUBROUTINE

    SUBROUTINE MPI_UTILS_FINALIZE()
      IMPLICIT NONE
      INTEGER :: ierr;
      
      CALL MPI_FINALIZE(ierr);

    END SUBROUTINE
    
    SUBROUTINE GetIndex(start,finish,mystart,myend)
      IMPLICIT NONE

      INTEGER :: mystart, myend, start, finish;
      INTEGER :: size;

      size = (finish - start);

      mystart = (size/numProc) * rank + start;
      
      IF(MOD(size,numProc) > numProc) THEN
        mystart = mystart + rank;
        myend = mystart + (size/numProc);
      ELSE
        mystart = mystart + MOD(size,numProc);
        myend = mystart + (size / numProc) - 1;
      END IF
      
      IF(rank .EQ. 0) THEN
        mystart = start;
      END IF
      
      IF(rank .EQ. (numProc-1)) THEN
        myend = finish;
      END IF
      
      IF(numProc >= size+1) THEN
        mystart = start + rank;
        myend = mystart;
        PRINT*, "NumProc is close to size, beware of undefined behavior"
      END IF

     !PRINT*, "PROC ", rank," START ", mystart," END ", myend;

    END SUBROUTINE

END MODULE MPI_UTILS





