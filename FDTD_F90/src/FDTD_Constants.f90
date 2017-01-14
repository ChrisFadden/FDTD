! MODULE: FDTD Constants Module

!>  \brief Contains various constants for an FDTD simulation, including sizes of various domains, and electric/magnetic constants.\n  
!>
! Description
!> dp = FORTRAN constant for double precision floating point on any OS \n
!> pi = Standard constant (3.14159...) \n
!> cc = The speed of light in m/s \n
!> mu0 = magnetic permeability in a vacuum \n
!> eps0 = electric permittivity in a vacuum \n
!> freq = The highest frequency for the simulation.  Used for Fourier Transform  when calculating power. \n
!> Nlambda = Number of points per wavelength \n
!> dx = The spatial step in the x direction \n
!> dy = The spatial step in the y direction \n
!> dt = The time step increment\n
!> TotalTime = The total number of time steps in the simulation \n
!> TFSF_Size = The size of the Total Field region, in cells, which must encompass the scatterer \n
!> PML_Size = The size, in cells, of the Perfectly Matched Layer, to prevent reflections from boundaries. \n
!> SizeX = The size in cells, in the x direction, of the total simulation \n
!> SizeY = The size in cells, in the y direction, of the total simualtion \n
!> TFSF_x0 = The first cell of the Total Field region in the x direction \n
!> TFSF_x1 = The last cell of the Total Field region in the x direction \n
!> TFSF_y0 = The first cell of the Total Field region in the y direction \n
!> TFSF_y1 = The last cell fo the Total Field region in the y direction \n
!> phi = The direction of the plane wave incident on the total field box \n
!> cosphi = The cosine of phi\n
!> sinphi = The sine of phi\n

MODULE FDTD_Constants
  IMPLICIT NONE
  
  !==================
  ! PRIVATE variables
  !==================
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.d0);
  
  !=================
  ! PUBLIC variables
  !=================
  
  !*******************
  !Constants
  !*******************
  REAL(dp), PARAMETER  :: pi   = 4.0d0 * DATAN(1.0d0);
  
  REAL(dp), PARAMETER  :: cc   = 299792458.d0;
  REAL(dp), PARAMETER  :: mu0  = 4.0d0 * pi * 10.0d0**(-7.0d0);
  REAL(dp), PARAMETER  :: eps0 = 1.0d0 / (cc * cc * mu0);

END MODULE FDTD_Constants






















