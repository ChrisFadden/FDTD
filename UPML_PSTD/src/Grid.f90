!> PSTD Constants
MODULE Grid
  USE PARAMETERS
  IMPLICIT NONE

  !------------------
  ! PRIVATE variables
  !------------------
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.d0);

  !-----------------
  ! PUBLIC variables
  !-----------------
 
  !-----------------
  ! Field Arrays
  !-----------------
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: Ex;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: Ey;
  
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: Ex_dy; 
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: Ey_dx; 

  REAL(dp), DIMENSION(SizeX-1,SizeY) :: DxField;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: DyField;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: DxStore;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: DyStore; 
  
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: Hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: Hz_dx;
  REAL(dp), DIMENSION(SizeY-1,SizeY-1) :: Hz_dy;

  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: Bz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: BzStore;

  !-------------------
  ! Coefficient Arrays
  !-------------------
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C1ex;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C2ex;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C3ex;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C4ex;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C5ex;
  REAL(dp), DIMENSION(SizeX-1,SizeY) :: C6ex;
  
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C1ey;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C2ey;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C3ey;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C4ey;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C5ey;
  REAL(dp), DIMENSION(SizeX,SizeY-1) :: C6ey;
  
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D1hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D2hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D3hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D4hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D5hz;
  REAL(dp), DIMENSION(SizeX-1,SizeY-1) :: D6hz;
  
  !-------------------
  ! Public Subroutines
  !-------------------
  PUBLIC :: Grid_Initialize

  CONTAINS
    SUBROUTINE Grid_Initialize()
      IMPLICIT NONE
      REAL(dp), PARAMETER :: C1 = 1.d0;
      REAL(dp), PARAMETER :: C2 = dt;
      REAL(dp), PARAMETER :: C3 = 1.0;
      REAL(dp), PARAMETER :: C4 = 1.d0 / 2.d0 / epsr / epsr / eps0 / eps0;
      REAL(dp), PARAMETER :: C5 = 2.0 * epsr * eps0;
      REAL(dp), PARAMETER :: C6 = 2.0 * epsr * eps0;

      REAL(dp), PARAMETER :: D1 = 1.d0;
      REAL(dp), PARAMETER :: D2 = dt;
      REAL(dp), PARAMETER :: D3 = 1.d0;
      REAL(dp), PARAMETER :: D4 = 1.d0 / 2.d0 / epsr / eps0 / mur / mu0;
      REAL(dp), PARAMETER :: D5 = 2.0 * epsr * eps0;
      REAL(dp), PARAMETER :: D6 = 2.0 * epsr * eps0;
      
      C1ex(:,PML_Size+1:SizeY-PML_Size) = C1;
      C2ex(:,PML_Size+1:SizeY-PML_Size) = C2;
      C3ex(:,:) = C3;
      C4ex(:,:) = C4;
      C5ex(PML_Size+1:SizeX-1-PML_Size,:) = C5;
      C6ex(PML_Size+1:SizeX-1-PML_Size,:) = C6;
      
      C1ey(:,:) = C1;
      C2ey(:,:) = C2;
      C3ey(PML_Size+1:SizeX-PML_Size,:) = C3;
      C4ey(PML_Size+1:SizeX-PML_Size,:) = C4;
      C5ey(:,PML_Size+1:SizeY-1-PML_Size) = C5;
      C6ey(:,PML_Size+1:SizeY-1-PML_Size) = C6;

      D1hz(PML_Size+1:SizeX-1-PML_Size,:) = D1;
      D2hz(PML_Size+1:SizeX-1-PML_Size,:) = D2;
      D3hz(:,PML_Size+1:SizeY-1-PML_Size) = D3;
      D4hz(:,PML_Size+1:SizeY-1-PML_Size) = D4;
      D5hz(:,:) = D5;
      D6hz(:,:) = D6;

      Ex(:,:) = 0;
      Ey(:,:) = 0;

      DxField(:,:) = 0;
      DyField(:,:) = 0;
      DxStore(:,:) = 0;
      DyStore(:,:) = 0; 
  
      Hz(:,:) = 0;
      Bz(:,:) = 0;
      BzStore(:,:) = 0;

    END SUBROUTINE !Initialize

END MODULE Grid

































