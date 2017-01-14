
!> Convolutional Perfectly Matched Layer Implementation   
MODULE CPML
  USE FDTD_Constants, ONLY: eps0, mu0
  USE Grid 
   
  IMPLICIT NONE
  PRIVATE !make variables private by default 
  
  !==================
  ! PRIVATE Variables
  !==================
  INTEGER, PARAMETER :: m = 3, ma = 1;
  INTEGER, PARAMETER :: dp = KIND(1.d0);

  REAL(dp) :: sig_x_max;
  REAL(dp) :: sig_y_max;
  
  REAL(dp), PARAMETER :: alpha_x_max = 0.24;
  REAL(dp), PARAMETER :: alpha_y_max = alpha_x_max;

  REAL(dp), PARAMETER :: kappa_x_max = 15.0;
  REAL(dp), PARAMETER :: kappa_y_max = kappa_x_max;

  !**************
  ! Arrays
  !**************
   
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_x, ce_x, alphae_x, sige_x;
  REAL(dp), DIMENSION(:), ALLOCATABLE :: kappae_x; 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: bh_x, ch_x, alphah_x, sigh_x; 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: kappah_x;
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_y, ce_y, alphae_y, sige_y;
  REAL(dp), DIMENSION(:), ALLOCATABLE :: kappae_y;
  REAL(dp), DIMENSION(:), ALLOCATABLE :: bh_y, ch_y, alphah_y, sigh_y;
  REAL(dp), DIMENSION(:), ALLOCATABLE :: kappah_y;

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Ezx1;  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Ezx2;  
  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Ezy1;
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Ezy2;
   
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Hyx1;
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Hyx2;
  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Hxy1;
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psi_Hxy2;
  
  !===================
  ! Public Subroutines
  !===================
  PUBLIC :: CPML_Initialize
  PUBLIC :: CPML_Finalize
  PUBLIC :: CPML_UpdateHx
  PUBLIC :: CPML_UpdateHy
  PUBLIC :: CPML_UpdateEz

  CONTAINS
    SUBROUTINE CPML_INITIALIZE()
      IMPLICIT NONE 
      INTEGER :: i, j, ii, jj;
      
      !*******************
      !Allocate arrays
      !*******************
      ALLOCATE(be_x(PML_Size),ce_x(PML_Size),alphae_x(PML_Size));
      ALLOCATE(sige_x(PML_Size),kappae_x(PML_Size));
      ALLOCATE(bh_x(PML_Size-1),ch_x(PML_Size-1),alphah_x(PML_Size-1));
      ALLOCATE(sigh_x(PML_Size-1),kappah_x(PML_Size-1));

      ALLOCATE(be_y(PML_Size),ce_y(PML_Size),alphae_y(PML_Size));
      ALLOCATE(sige_y(PML_Size),kappae_y(PML_Size));
      ALLOCATE(bh_y(PML_Size-1),ch_y(PML_Size-1),alphah_y(PML_Size-1));
      ALLOCATE(sigh_y(PML_Size-1),kappah_y(PML_Size-1));

      ALLOCATE(psi_Ezx1(PML_Size,SizeY),psi_Ezx2(PML_Size,SizeY));
      ALLOCATE(psi_Ezy1(SizeX,PML_Size),psi_Ezy2(SizeX,PML_Size));

      ALLOCATE(psi_Hyx1(PML_Size-1,SizeY),psi_Hyx2(PML_Size-1,SizeY));
      ALLOCATE(psi_Hxy1(SizeX,PML_Size-1),psi_Hxy2(SizeX,PML_Size-1));
      
      kappae_x(:) = 1.d0;
      kappah_x(:) = 1.d0;

      kappae_y(:) = 1.d0;
      kappah_y(:) = 1.d0;

      psi_Ezx1(:,:) = 0.d0;
      psi_Ezx2(:,:) = 0.d0;
      psi_Ezy1(:,:) = 0.d0;
      psi_Ezy2(:,:) = 0.d0;
      
      psi_Hyx1(:,:) = 0.d0;
      psi_Hyx2(:,:) = 0.d0;
      psi_Hxy1(:,:) = 0.d0;
      psi_Hxy2(:,:) = 0.d0;
        
      sig_x_max = 0.75 * (0.8*(m+1)/(dx * (mu0/eps0)**0.5));
      sig_y_max = 0.75 * (0.8*(m+1)/(dy * (mu0/eps0)**0.5));
           
      DO i = 1, PML_Size 
        sige_x(i) = sig_x_max * ( (PML_Size - i) / (PML_Size - 1.0))**m;
        alphae_x(i) = alpha_x_max * ( (i-1.0) / (PML_Size-1.0) )**ma;
        kappae_x(i) = 1.0 + (kappa_x_max - 1.0) * &
                        ((PML_Size - i) / (PML_Size - 1.0))**m;
        be_x(i) = DEXP(-(sige_x(i) / kappae_x(i) + &
                    alphae_x(i)) * dt / eps0);
        IF((sige_x(i) == 0.0d0) .AND. &
          (alphae_x(i) == 0.0d0) .AND. (i == PML_Size)) THEN
            ce_x(i) = 0.0d0
        ELSE
          ce_x(i) = sige_x(i) * (be_x(i) - 1.0) / & 
            (sige_x(i) + kappae_x(i) * alphae_x(i)) / kappae_x(i)
        END IF
      END DO
            
      DO i = 1, PML_Size - 1
        sigh_x(i) = sig_x_max * ( (PML_Size - i - 0.5) / (PML_Size - 1.0))**m;
        alphah_x(i) = alpha_x_max * ((i - 0.5) / (PML_Size - 1.0))**ma;
        kappah_x(i) = 1.0 + (kappa_x_max - 1.0) * &
                      ((PML_Size - i - 0.5) / (PML_Size - 1.0))**m;
        bh_x(i) = DEXP(-(sigh_x(i) / kappah_x(i) + &
                          alphah_x(i)) * dt / eps0);
        ch_x(i) = sigh_x(i) * (bh_x(i) - 1.0) / &
                   (sigh_x(i) + kappah_x(i) * alphah_x(i)) / kappah_x(i);

      END DO 

      DO j = 1, PML_Size
        sige_y(j) = sig_y_max * ( (PML_Size - j) / (PML_Size-1.0))**m;
        alphae_y(j) = alpha_y_max * ( (j-1) / (PML_Size-1.0))**ma;
        kappae_y(j) = 1.0 + (kappa_y_max - 1.0) * &
                        ((PML_Size - j) / (PML_Size-1.0))**m;
        be_y(j) = DEXP(-(sige_y(j) / kappae_y(j) + &
                    alphae_y(j)) * dt / eps0);

        IF((sige_y(j) == 0.0) .AND. &
          (alphae_y(j) == 0.0) .AND. (j == PML_Size)) THEN
          ce_y(j) = 0.0d0
        ELSE
          ce_y(j) = sige_y(j) * (be_y(j) - 1.0) / &
            (sige_y(j) + kappae_y(j) * alphae_y(j)) / kappae_y(j)
        END IF
             
      END DO

      DO j = 1, PML_Size-1
        sigh_y(j) = sig_y_max * ( (PML_Size - j - 0.5) / (PML_Size-1.0))**m;
        alphah_y(j) = alpha_y_max * ( (j - 0.5) / (PML_Size-1.0))**ma;
        kappah_y(j) = 1.0 + (kappa_y_max - 1.0) * &
                      ((PML_Size - j - 0.5) / (PML_Size - 1.0))**m;
        bh_y(j) = DEXP(-(sigh_y(j) / kappah_y(j) + &
                    alphah_y(j)) * dt /eps0);
        ch_y(j) = sigh_y(j) * (bh_y(j) - 1.0) / &
                    (sigh_y(j) + kappah_y(j)*alphah_y(j)) / kappah_y(j);                
      END DO
      
      !need hx, hy, ex, ey
      DO i = 1, PML_Size-1
          den_hx(i) = den_hx(i) / kappah_x(i);
      END DO
       
      ii = PML_Size-1;
      DO i = (SizeX + 1) - PML_Size, SizeX-1 
        den_hx(i) = den_hx(i) / kappah_x(ii);
        ii = ii - 1;
      END DO
      
      DO j = 1, PML_Size-1
        den_hy(j) = den_hy(j) / kappah_y(j);
      END DO

      jj = PML_Size-1;
      DO j = (SizeY+1) - PML_Size, SizeY-1
        den_hy(j) = den_hy(j) / kappah_y(jj);
        jj = jj - 1;
      END DO
      
      DO j = 1, PML_Size
        den_ey(j) = den_ey(j) / kappae_y(j);
      END DO
         
      DO i = 1, PML_Size  
        den_ex(i) = den_ex(i) / kappae_x(i);
      END DO
      
      jj = PML_Size;
      DO j = (SizeY+1) - PML_Size, SizeY-1
        den_ey(j) = den_ey(j) / kappae_y(jj);
        jj = jj-1;
      END DO

      ii = PML_Size;
      DO i = (SizeX+1) - PML_Size, SizeX-1
        den_ex(i) = den_ex(i) / kappae_x(ii);
        ii = ii-1;
      END DO
      
    END SUBROUTINE
    
    SUBROUTINE CPML_UpdateHx()
      IMPLICIT NONE
      INTEGER :: i, j, jj;
      
      DO i = 1, SizeX-1
        
        !***********************
        !bottom Hx, j-direction
        !***********************
        DO j = 1, PML_Size-1
          psi_Hxy1(i,j) = bh_y(j) * psi_Hxy1(i,j) + &
                          ch_y(j) * (Ez(i,j) - Ez(i,j+1)) / dy;
          Hx(i,j) = Hx(i,j) + (dt/mu0) * psi_Hxy1(i,j); 
        END DO !j loop
        
        !*******************
        !top Hx, j-direction
        !*******************
        jj = PML_Size-1;
        DO j = (SizeY + 1) - PML_Size, SizeY-1
          psi_Hxy2(i,jj) = bh_y(jj) * psi_Hxy2(i,jj) + &
                          ch_y(jj) * (Ez(i,j) - Ez(i,j+1)) / dy;
          Hx(i,j) = Hx(i,j) + (dt/mu0) * psi_Hxy2(i,jj);  
          jj = jj-1;
        END DO !j loop
      
      END DO !i loop

    END SUBROUTINE
    
    SUBROUTINE CPML_UpdateHy()
      IMPLICIT NONE
      INTEGER :: i, j, ii;
      
      DO j = 1, SizeY-1
        
        !**********************
        !bottom Hy, i-direction
        !**********************
        DO i = 1, PML_Size-1
          psi_Hyx1(i,j) = bh_x(i) * psi_Hyx1(i,j) + &
                           ch_x(i) * (Ez(i+1,j) - Ez(i,j))/dx;
          Hy(i,j) = Hy(i,j) + (dt/mu0) * psi_Hyx1(i,j);
        END DO !i loop
        
        !*******************
        !top Hy, i-direction
        !*******************
        ii = PML_Size-1;
        DO i = (SizeX+1) - PML_Size, SizeX-1
          psi_Hyx2(ii,j) = bh_x(ii) * psi_Hyx2(ii,j) + &
                            ch_x(ii) * (Ez(i+1,j) - Ez(i,j))/dx;
          Hy(i,j) = Hy(i,j) + (dt/mu0) * psi_Hyx2(ii,j);
          ii = ii-1;
        END DO !i loop 
      END DO !j loop

    END SUBROUTINE

    SUBROUTINE CPML_UpdateEz()
      IMPLICIT NONE
      INTEGER :: i, j, ii, jj;
      
      DO j = 2, SizeY-1
        !**********************
        !bottom Ez, i-direction
        !**********************
        DO i = 2, PML_Size
          psi_Ezx1(i,j) = be_x(i) * psi_Ezx1(i,j) + &
                           ce_x(i) * (Hy(i,j) - Hy(i-1,j))/dx;
          Ez(i,j) = Ez(i,j) + (dt / eps0) * psi_Ezx1(i,j); 
        END DO !i loop
        
        !*******************
        !top Ez, i-direction
        !*******************
        ii = PML_Size;
        DO i = (SizeX+1) - PML_Size, SizeX-1
          psi_Ezx2(ii,j) = be_x(ii) * psi_Ezx2(ii,j) + &
                            ce_x(ii) * (Hy(i,j) - Hy(i-1,j))/dx;
          Ez(i,j) = Ez(i,j) + (dt / eps0) * psi_Ezx2(ii,j);
          ii = ii-1;
        END DO !i loop
      END DO !j loop
      
      DO i = 2, SizeX-1
        !**********************
        !bottom Ez, j-direction
        !**********************
        DO j = 2, PML_Size
          psi_Ezy1(i,j) = be_y(j) * psi_Ezy1(i,j) + &
                           ce_y(j) * (Hx(i,j) - Hx(i,j-1))/dy;
          Ez(i,j) = Ez(i,j) + (dt / eps0)*psi_Ezy1(i,j);
        END DO
        
        !*******************
        !top Ez, j-direction
        !*******************
        jj = PML_Size; 
        DO j = (SizeY+1) - PML_Size, SizeY-1
          psi_Ezy2(i,jj) = be_y(jj) * psi_Ezy2(i,jj) + &
                            ce_y(jj) * (Hx(i,j) - Hx(i,j-1))/dy;
          Ez(i,j) = Ez(i,j) + (dt / eps0) * psi_Ezy2(i,jj);
          jj = jj-1
        END DO
      END DO !i loop

    END SUBROUTINE
    
    SUBROUTINE CPML_Finalize()
      DEALLOCATE(be_x,ce_x,alphae_x);
      DEALLOCATE(sige_x,kappae_x);
      DEALLOCATE(bh_x,ch_x,alphah_x);
      DEALLOCATE(sigh_x,kappah_x);

      DEALLOCATE(be_y,ce_y,alphae_y);
      DEALLOCATE(sige_y,kappae_y);
      DEALLOCATE(bh_y,ch_y,alphah_y);
      DEALLOCATE(sigh_y,kappah_y);

      DEALLOCATE(psi_Ezx1,psi_Ezx2);
      DEALLOCATE(psi_Ezy1,psi_Ezy2);

      DEALLOCATE(psi_Hyx1,psi_Hyx2);
      DEALLOCATE(psi_Hxy1,psi_Hxy2);
      
    END SUBROUTINE

END MODULE CPML














