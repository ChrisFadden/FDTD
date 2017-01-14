

!> Main FDTD Update Equations
MODULE FDTD
  USE FDTD_Constants, ONLY : pi
  USE Grid 
  USE TFSF, ONLY : EzInc, HxInc, HyInc

  IMPLICIT NONE
  
  !>  Electric and Magnetic Fields
  !!
  !! \details  Because of the staggering of the Yee grid,
  !! the fields are offset from each other.  Therefore,
  !! Hx(m,n+1/2), Hy(m+1/2,n), and Ez(m,n) share the 
  !! same index, [m, n].

  !==================
  ! PRIVATE variables
  !==================
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.d0);

  !=====================
  ! PUBLIC variables
  !=====================
 
  CONTAINS

    !>  H-field FDTD Update Loops 
    !! \param Hx field
    !! \param Hy field

    SUBROUTINE UpdateH()
      IMPLICIT NONE
      
      INTEGER :: i, j;

      !==========
      ! Hx Update
      !==========
      DO j = 1, SizeY-1
        DO i = 1, SizeX-1
          
          !Chxh = 1.d0 
          !Chxe = dt / (mu0 * dx);
          Hx(i,j) = Chxh(i,j) * Hx(i,j) &
            + Chxe(i,j) * den_hy(j) * ( Ez(i,j) - Ez(i,j+1) ); 
        END DO 
      END DO 
        
      !==========
      ! Hy Update
      !==========
      DO j = 1, SizeY-1
        DO i = 1, SizeX-1
          Hy(i,j) = Chyh(i,j) * Hy(i,j) &
            + Chye(i,j) * den_hx(i) * ( Ez(i+1, j) - Ez(i,j) );
        END DO 
      END DO 
      
    END SUBROUTINE 
    
    !> Ez FDTD Update Loop
    !! \param Ez
    SUBROUTINE UpdateE()
      IMPLICIT NONE
      INTEGER :: i, j;
      
      !Ceze = 1.d0
      !Cezh = dt / (eps0 * dx)
      DO j = 2, SizeY - 1
        DO i = 2, SizeX - 1
          Ez(i,j) = Ceze(i,j) * Ez(i,j) &
            + Cezh(i,j)* den_ex(i)*( (Hy(i,j) - Hy(i-1,j)) &
            + den_ey(j)*(Hx(i,j-1) - Hx(i,j)) );
        END DO
      END DO
      
    END SUBROUTINE 
    
    !>  Hx, Hy TF/SF Updates 
    !! \param Hx
    !! \param Hy
    !! \param EzInc
    SUBROUTINE TFSF_Hupdate()
      IMPLICIT NONE
      INTEGER :: i, j;

      j = TFSF_y0; 
      DO i = TFSF_x0, TFSF_x1
        Hx(i,j-1) = Hx(i,j-1) + Chxe(i,j)*EzInc(i,j);
      END DO
      
      j = TFSF_y1;
      DO i = TFSF_x0, TFSF_x1
        Hx(i,j) = Hx(i,j) - Chxe(i,j)*EzInc(i,j);
      END DO

      i = TFSF_x0;
      DO j = TFSF_y0, TFSF_y1
        Hy(i-1,j) = Hy(i-1, j) - Chye(i,j)*EzInc(i,j);
      END DO 
      
      i = TFSF_x1;
      DO j = TFSF_y0, TFSF_y1
        Hy(i,j) = Hy(i,j) + Chye(i,j)*EzInc(i,j);
      END DO

    END SUBROUTINE
    
    !>  Ez TF/SF Updates 
    !! \param Ez
    !! \param HxInc
    !! \param HyInc
    SUBROUTINE TFSF_Eupdate()
    IMPLICIT NONE
    INTEGER :: i, j;
    
    j = TFSF_y0;
    DO i = TFSF_x0, TFSF_x1
      Ez(i,j) = Ez(i,j) + Cezh(i,j)*HxInc(i,j-1);
    END DO

    j = TFSF_y1;
    DO i = TFSF_x0, TFSF_x1
      Ez(i,j) = Ez(i,j) - Cezh(i,j)*HxInc(i,j);
    END DO

    i = TFSF_x0;
    DO j = TFSF_y0, TFSF_y1
      Ez(i,j) = Ez(i,j) - Cezh(i,j)*HyInc(i-1,j);
    END DO
    
    i = TFSF_x1;
    DO j = TFSF_y0, TFSF_y1
      Ez(i,j) = Ez(i,j) + Cezh(i,j)*HyInc(i,j);
    END DO 

    END SUBROUTINE
    
    !> Source function for FDTD
    !! \param freq
    !! \param time
    FUNCTION Source(n)
      REAL(dp) :: Source
      REAL(dp) :: t0;
      REAL(dp) :: tw;
      INTEGER, INTENT(IN) :: n 
      
      IF(SRC_GAUSS .EQ. 0) THEN
        Source = 20.d0 * DSIN(2.d0 * pi * freq * dt * n);   
      ELSE
        tw = 0.5 / freq;
        t0 = 4.0 * tw;
        Source = -2.d0 * (( n * dt - t0) / tw) * &
                 DEXP(-((n * dt - t0) / tw) * ((n * dt - t0) / tw));
      END IF
    END FUNCTION 

END MODULE FDTD









