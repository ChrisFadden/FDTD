MODULE FFTW_utils
  !USE Grid
  IMPLICIT NONE 
!  PRIVATE
  
  INCLUDE "fftw3.f"
  
  INTEGER, PRIVATE, PARAMETER :: dp = KIND(1.d0);
  
  CONTAINS
    SUBROUTINE Take_FFT(inArray, outArray)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:) :: inArray;
      COMPLEX(dp), DIMENSION(:) :: outArray;
      INTEGER :: Plan;
      
      !Create Plan
      CALL DFFTW_PLAN_DFT_R2C_1D(Plan,SIZE(inArray),inArray,outArray,FFTW_ESTIMATE);
      !Execute FFT
      CALL DFFTW_EXECUTE_DFT_R2C(Plan,inArray,outArray);
      !Destroy Plan
      CALL DFFTW_DESTROY_PLAN(Plan);
 
    END SUBROUTINE
    
    SUBROUTINE TAKE_IFFT(inArray, outArray)
      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(:) :: inArray;
      REAL(dp), DIMENSION(:) :: outArray;
      INTEGER :: Plan;
      
      !Create Plan
      CALL DFFTW_PLAN_DFT_C2R_1D(Plan,SIZE(outArray),inArray,outArray,FFTW_ESTIMATE);
      !Execute FFT
      CALL DFFTW_EXECUTE_DFT_C2R(Plan,inArray,outArray);
      !Destroy Plan
      CALL DFFTW_DESTROY_PLAN(Plan);
    END SUBROUTINE 
    
    SUBROUTINE FFT_PARTIAL_X(inArray, lenX, outputArray)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:) :: inArray;
      COMPLEX(dp), DIMENSION(SIZE(inArray)/2 + 1) :: outArray;
      REAL(dp), DIMENSION(SIZE(inArray)/2 + 1) :: kx;
      REAL(dp), DIMENSION(:) :: outputArray;
      INTEGER :: i,j; 
      INTEGER :: lenX;
      REAL(dp) :: pi; 
      COMPLEX(dp), PARAMETER :: jj = (0.d0,1.d0);
      
      pi = ACOS(-1.d0);

      DO i = 1,SIZE(inArray)/2 
          kx(i) = 2.d0*pi * REAL(i-1) / lenX; 
      END DO
      
      kx(SIZE(inArray)/2 + 1) = 0;
     
      CALL TAKE_FFT(inArray,outArray); 
      outArray = kx * jj * outArray;       
      CALL TAKE_IFFT(outArray,outputArray);
      
      outputArray = outputArray / SIZE(inArray); 
    END SUBROUTINE
    
    !This is basically equivalent to FFT_PARTIAL_X
    SUBROUTINE FFT_PARTIAL_Y(inArray, lenY, outputArray)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:) :: inArray;
      COMPLEX(dp), DIMENSION(SIZE(inArray)/2 + 1) :: outArray;
      REAL(dp), DIMENSION(SIZE(inArray)/2 + 1) :: ky;
      REAL(dp), DIMENSION(:) :: outputArray;
      INTEGER :: i,j;
      INTEGER :: lenY;
      REAL(dp) :: pi;
      COMPLEX(dp), PARAMETER :: jj = (0.d0,1.d0);

      pi = ACOS(-1.d0);
      
      DO j = 1,SIZE(inArray)/2
        ky(j) = 2.d0*pi * REAL(j-1) / lenY;
      END DO
      
      ky(SIZE(inArray)/2 + 1) = 0;
      CALL TAKE_FFT(inArray,outArray);
      outArray = ky * jj * outArray;
      CALL TAKE_IFFT(outArray,outputArray);

      outputArray = outputArray / SIZE(inArray);
    END SUBROUTINE

END MODULE FFTW_utils
