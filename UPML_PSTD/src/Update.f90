!> PSTD Update Equations
MODULE Update
  USE PARAMETERS
  USE Grid
  USE FFTW_utils
  USE omp_lib
  IMPLICIT NONE
  PRIVATE
  !------------------
  ! PRIVATE variables
  !------------------
  INTEGER, PARAMETER :: dp = KIND(1.d0);
  !-------------------
  ! Public Subroutines
  !-------------------
  PUBLIC :: UpdateH
  PUBLIC :: UpdateE

  CONTAINS
    SUBROUTINE UpdateH()
      IMPLICIT NONE
      
      INTEGER :: i,j; !Loop Variables
      BzStore = Bz;
      
      !$OMP PARALLEL DO
      DO j = 1,SizeY-1
       CALL FFT_PARTIAL_X(Ey(:,j),SizeX,Ey_dx(:,j)); 
      END DO
      !$OMP END PARALLEL DO 

      Ey_dx = Ey_dx / dx; 
      
      !$OMP PARALLEL DO
      DO i = 1,SizeX-1
       CALL FFT_PARTIAL_Y(Ex(i,:),SizeY,Ex_dy(i,:)); 
      END DO
      !$OMP END PARALLEL DO
      Ex_dy = Ex_dy / dx;
      
      !$OMP PARALLEL DO
      DO j = 1,SizeY-1
        DO i = 1,SizeX-1
          Bz(i,j) = D1hz(i,j) * Bz(i,j) - D2hz(i,j) * &
                              (Ey_dx(i,j) - Ex_dy(i,j));          
                   
          Hz(i,j) = D3hz(i,j) * Hz(i,j) + D4hz(i,j) * &
                    (D5hz(i,j) * Bz(i,j) - D6hz(i,j)*BzStore(i,j));
        END DO !end j loop
      END DO !end i loop
      !$OMP END PARALLEL DO

    END SUBROUTINE !H Updates
    
    SUBROUTINE UpdateE(t)
      IMPLICIT NONE
      INTEGER :: i,j,t;
       
      DxStore = DxField;
      
      !$OMP PARALLEL DO
      DO i = 1,SizeX-1
        CALL FFT_PARTIAL_Y(Hz(i,:),SizeY-1,Hz_dy(i,:));
      END DO
      !$OMP END PARALLEL DO
      Hz_dy = Hz_dy / dx;
      
      !$OMP PARALLEL DO
      DO j = 1,SizeY-1
        CALL FFT_PARTIAL_X(Hz(:,j),SizeX-1,Hz_dx(:,j));
      END DO
      !$OMP END PARALLEL DO
      Hz_dx = Hz_dx / dx;
      
      !$OMP PARALLEL DO
      DO j = 1,SizeY-1
        DO i = 1, SizeX-1        
          DxField(i,j) = C1ex(i,j) * DxField(i,j) + &
                        C2ex(i,j) * Hz_dy(i,j); 

          Ex(i,j) = C3ex(i,j)*Ex(i,j)+C4ex(i,j) * &
                           (C5ex(i,j)*DxField(i,j) - &
                            C6ex(i,j)*DxStore(i,j));
        END DO !j loop
      END DO !i loop
      !$OMP END PARALLEL DO 

      DyStore = DyField;
      
      !$OMP PARALLEL DO
      DO j = 1,SizeY-1
        DO i = 1,SizeX-1      
          DyField(i,j) = C1ey(i,j) * DyField(i,j) - &
                         C2ey(i,j) * Hz_dx(i,j); 
                                                                          
          Ey(i,j) = C3ey(i,j)*Ey(i,j)+C4ey(i,j)* &
                       (C5ey(i,j)*DyField(i,j) - &
                        C6ey(i,j)*DyStore(i,j));
        END DO !j loop
      END DO !i loop
      !$OMP END PARALLEL DO

      !*************
      !Apply Source
      !*************
      
      !Sine Excitation
      Ey(is,js) = Ey(is,js) + 20000*sin(2*pi*freq * dt*t); 
      Ey(is,js+1) = Ey(is,js+1) + 20000*sin(2*pi*freq * dt*t);
      !Gaussian Excitation
!      Ey(is,js) = Ey(is,js) - 20000.d0 * ((t*dt - t0) / tw) * &
!                              DEXP(-1.d0 * ((t*dt - t0) / tw)**2);
!      Ey(is,js+1) = Ey(is,js+1) - 20000.d0 * ((t*dt - t0) / tw) * &
!                              DEXP(-1.d0 * ((t*dt - t0) / tw)**2);
 
    END SUBROUTINE !E Updates

END MODULE Update

































