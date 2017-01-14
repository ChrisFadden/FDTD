/** @file Grid.h
 *  @brief Class definition for global electromagnetics grid
 *
 *  @author Chris Fadden
 *
 */

#ifndef _GRID_H
#define _GRID_H

#include <cmath>
#include <vector>
#include <string>

class Grid {
 public:
  /** @brief Constructor of the Grid class
   *
   *  The grid is initialized using values taken
   *  from the file created by the Java GUI.  Please
   *  take note of the main python script used to run
   *  this, and keep in mind the filename is hardcoded
   *  at this time.
   */
  Grid();

  /** @brief returns the maximum time of simulation
   *
   *   This ensures protection of the Grid class variables,
   *   and acts as a get function for the main program
   */
  int getMaxTime();

  /** @brief Returns value of the harmonic source
   *
   *  This evaluates a harmonic source, a sine wave,
   *  at the specified time, and returns that value.
   *  This has a single frequency content, and is mostly
   *  used for testing purposes.
   *
   *  @param t the time at which to evaluate the function
   *  @return the value at the specified time
   */

  double HarmonicSource(int);

  /** @brief Returns value of the gaussian source
   *
   *    This evaluates a gaussian source, an exponential raised
   *    to the $x^2$.  This has the benefit in a time domain
   *    simulation of having multiple frequency components,
   *    and therefore can give a more broadband response for
   *    the simulation
   *
   *  @param t the time at which to evaluate the function
   *  @return the value at the specified time
   */
  double GaussianSource(int);

  double Source(int);

  int SizeX;
  int SizeY;
  int PML_size;
  double rError;  // reflection error

  int t = 0;
  int MaxTime;

  int isrc; /**< x-coordinate of the source on the grid */
  int jsrc; /**< y-coordinate of the source on the grid */

  const int cc = 299792458;
  const double mu0 = 16 * atan(1) * 1.0e-7;
  const double eps0 = 1.0 / (cc * cc * mu0);

  double dx;
  double dt;
  double CFL;

  enum SourceType { HARMONIC, GAUSSIAN };
  SourceType src;
  double freq;

  double tw;
  double t0;

  std::vector<double> epsr;
  std::vector<double> mur;
};
#endif
