!> PSTD Constants
MODULE UPML
  USE PARAMETERS
  USE Grid
  IMPLICIT NONE
  PRIVATE
  !------------------
  ! PRIVATE variables
  !------------------
  INTEGER, PARAMETER :: dp = KIND(1.d0);

  REAL(dp), PARAMETER :: rmax = EXP(-16.d0);
  REAL(dp), PARAMETER :: orderbc = 4.0;
  REAL(dp) :: delbc, kmax;
  REAL(dp) :: sigmam, sigfactor, kfactor;
  
  REAL(dp) :: x1, x2, y1, y2;
  REAL(dp) :: sigma, ki;
  REAL(dp) :: facm, facp;

  !-------------------
  ! Public Subroutines
  !-------------------
  PUBLIC :: UPML_Initialize

  CONTAINS
    SUBROUTINE UPML_Initialize()
      IMPLICIT NONE
      INTEGER :: i, j; !Loop Variables
      !--------------
      ! X Boundary
      !--------------
      delbc = PML_Size * dx;
      sigmam = -LOG(rmax)*(orderbc+1.d0) / (2.d0 * eta * delbc);
      sigfactor = sigmam / (dx * (delbc**orderbc) * (orderbc+1.d0));
      kmax = 11;
      kfactor = (kmax - 1.d0) / dx / (orderbc + 1.d0) / (delbc**orderbc);
      
      DO i = 1,PML_Size
        x1 = (PML_Size - i + 1) * dx;
        x2 = (PML_Size - i) * dx;
        sigma = sigfactor * (x1**(orderbc+1) - x2**(orderbc+1));
        ki = 1 + kfactor*(x1**(orderbc+1) - x2**(orderbc+1));
        facm = (2*epsr*eps0*ki - sigma*dt);
        facp = (2*epsr*eps0*ki + sigma*dt);

        C5ex(i,:) = facp;
        C5ex(SizeX-i,:) = facp;
        C6ex(i,:) = facm;
        C6ex(SizeX-i,:) = facm;

        D1hz(i,:) = facm / facp;
        D1hz(SizeX-i,:) = facm/facp;
        D2hz(i,:) = 2.d0 * epsr * eps0 * dt / facp;
        D2hz(SizeX-i,:) = 2.d0 * epsr * eps0 *dt / facp;
        
        !Grid Cell Boundary
        x1 = (PML_Size - i + 1.5) * dx;
        x2 = (PML_Size - i + 0.5) * dx;
        sigma = sigfactor*(x1**(orderbc+1) - x2**(orderbc+1));
        ki = 1.d0 + kfactor * (x1**(orderbc+1) - x2**(orderbc+1));
        facm = (2.d0 * epsr * eps0 * ki - sigma*dt);
        facp = (2.d0 * epsr * eps0 * ki + sigma*dt);
        
        C3ey(i,:) = facm / facp;
        C3ey(SizeX-i+1,:) = facm / facp;
        C4ey(i,:) = 1.d0 / facp / epsr / eps0;
        C4ey(SizeX-i+1,:) = 1.0 / facp / epsr /eps0; 
      END DO

      !PEC Boudary Walls
      C3ey(1,:) = -1.d0;
      C3ey(SizeX,:) = -1.d0;
      C4ey(1,:) = 0.d0;
      C4ey(SizeX,:) = 0.d0;
      
      !--------------------
      !y-varying parameters
      !--------------------
      delbc = PML_Size * dx;
      sigmam = -LOG(rmax) * epsr * eps0 * cc * (orderbc+1.d0) / (2.d0 * delbc);
      sigfactor = sigmam / (dx * (delbc**orderbc) * (orderbc + 1.d0));
      kmax = 11;
      kfactor = (kmax - 1.d0) / dx / (orderbc + 1.d0) / (delbc**orderbc);
      
      DO j = 1, PML_Size
        y1 = (PML_Size - j + 1) * dx;
        y2 = (PML_Size - j) * dx;
        sigma = sigfactor * (y1**(orderbc+1) - y2**(orderbc+1));
        ki = 1.d0 + kfactor * (y1**(orderbc+1) - y2**(orderbc+1));
        facm = (2*epsr*eps0*ki - sigma*dt);
        facp = (2*epsr*eps0*ki + sigma*dt);

        C5ey(:,j) = facp;
        C5ey(:,SizeY-j) = facp;
        C6ey(:,j) = facm;
        C6ey(:,SizeY-j) = facm;

        D3hz(:,j) = facm / facp;
        D3hz(:,SizeY-j) = facm / facp;
        D4hz(:,j) = 1.d0 / facp / mur / mu0;
        D4hz(:,SizeY-j) = 1.d0 / facp / mur / mu0;
        
        !Grid Cell Boundary
        y1 = (PML_Size - j + 1.5) * dx;
        y2 = (PML_Size - j + 0.5) * dx;
        sigma = sigfactor * (y1**(orderbc+1) - y2**(orderbc+1));
        ki = 1.d0 + kfactor*(y1**(orderbc+1) - y2**(orderbc+1));
        facm = (2*epsr*eps0*ki - sigma*dt);
        facp = (2*epsr*eps0*ki + sigma*dt);

        C1ex(:,j) = facm / facp;
        C1ex(:,SizeY-j + 1) = facm / facp;
        C2ex(:,j) = 2.d0 * epsr * eps0 * dt/ facp;
        C2ex(:,SizeY-j + 1) = 2.d0 * epsr * eps0 * dt / facp;

      END DO
      
      C1ex(:,1) = -1;
      C1ex(:,SizeY) = -1;
      C2ex(:,1) = 0;
      C2ex(:,SizeY) = 0;

    END SUBROUTINE !Initialize

END MODULE UPML

































