!> PSTD Constants
MODULE PARAMETERS
  IMPLICIT NONE

  !------------------
  ! PRIVATE variables
  !------------------
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.d0);

  !-----------------
  ! PUBLIC variables
  !-----------------
  
  REAL(dp), PARAMETER :: cc = 299792458.d0;
  REAL(dp), PARAMETER :: mu0 = 4.d0 * acos(-1.d0) * 10.d0**(-7.d0);
  REAL(dp), PARAMETER :: eps0 = 1.d0 / (cc*cc*mu0);
  REAL(dp), PARAMETER :: eta0 = SQRT(mu0 / eps0); 
  
  REAL(dp), PARAMETER :: mur = 1.d0;
  REAL(dp), PARAMETER :: epsr = 1.d0;
  REAL(dp), PARAMETER :: eta = eta0 * SQRT(mur / epsr);
  
  REAL(dp), PARAMETER :: pi = ACOS(-1.d0);

      !---------------
      ! Grid Constants
      !---------------
  INTEGER, PARAMETER :: SizeX_NoPML = 100;
  INTEGER, PARAMETER :: SizeY_NoPML = 100;

  INTEGER, PARAMETER :: PML_Size = 10;
  
  INTEGER, PARAMETER :: SizeX = SizeX_NoPML + 2*PML_Size;
  INTEGER, PARAMETER :: SizeY = SizeY_NoPML + 2*PML_Size;
      
  INTEGER, PARAMETER :: is = SizeX / 2;
  INTEGER, PARAMETER :: js = SizeY / 2;

      !--------------------
      ! Grid Parameters
      !--------------------
  REAL(dp), PARAMETER :: dx = 0.002;
  REAL(dp), PARAMETER :: dt = dx / (SQRT(2.d0) * pi * cc);
  INTEGER, PARAMETER ::  MaxTime = 300;

      !---------------------
      ! Source Excitation
      !---------------------
  REAL(dp), PARAMETER :: freq = 10.d0 * 10**9;
  REAL(dp), PARAMETER :: tw = 0.5 / freq;
  REAL(dp), PARAMETER :: t0 = 4.d0 * tw;
 
 END MODULE PARAMETERS

































