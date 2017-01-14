/** @file err.cpp
 *  @brief print out error messages
 *
 *  This file checks the user defined
 *  parameters, and ensures they make
 *  sense for a simulation.
 *
 *  Anything that is impossible to simulate,
 *  such as a source located at a
 *  negative location, will print out
 *  an error.
 *
 *  Recommended good practice, such as
 *  having a wavelength in between the
 *  scatter and source, will print out
 *  a warning.
 */
#include "error.h"
#include <iostream>

#define RED "\033[1m\033[31m"
#define YELLOW "\033[1m\033[33m"
#define WHITE "\033[1m"
#define RESET "\033[0m"

/** @brief red error message printer
 *
 *  This simply prints out the passed
 *  string in bold with a red WARNING
 *  in front of it, to alert the user
 *  that there was an error, and the
 *  simulation was not completed
 */
void error_msg(std::string s) {
  std::cout << RED << "ERROR: " << RESET;
  std::cout << WHITE << s << RESET << std::endl;
}

/** @brief yellow warning message printer
 *
 *  This prints out the passed string in
 *  bold with a yellow WARNING in front.
 *  This alerts the user that the parameters
 *  they have chosen do not follow good FDTD
 *  practice, and could have erroneous results
 */
void warning_msg(std::string s) {
  std::cout << YELLOW << "WARNING: " << RESET;
  std::cout << WHITE << s << RESET << std::endl;
}

/** @brief identifies bad input parameters
 *
 *  This is a series of if statements checking
 *  user input parameters.  This only checks
 *  the most obvious bad parameters, and should
 *  not be absolutely relied on for good
 *  simulation results
 *
 *  @return true if there was an error, false otherwise
 */
bool check_input_parameters(grid *g) {
  bool error = false;

  if (dx < 0) {
    error_msg("Negative frequencies not allowed");
    error = true;
  }
  if (!sizeX) {
    error_msg("size X never set");
    error = true;
  }
  if (!sizeY) {
    error_msg("size Y never set");
    error = true;
  }
  if (!(i_src || tfsf_x0)) {
    error_msg("Source not present on grid (x index never set)");
    error = true;
  }
  if (!(j_src || tfsf_y0)) {
    error_msg("Source not present on grid (y index never set)");
    error = true;
  }
  if (!(k_src || tfsf_z0)) {
    error_msg("Source not present on grid (z index never set)");
    error = true;
  }

  if (i_src > sizeX) {
    error_msg("Source not present on grid (x index too big)");
    error = true;
  }
  if (j_src > sizeY) {
    error_msg("Source not present on grid (y index too big)");
    error = true;
  }

  if (k_src > sizeZ) {
    error_msg("Source not present on grid (z index too big)");
    error = true;
  }
  if (tfsf_x1 < tfsf_x0) {
    error_msg("TFSF x1 less than x0");
    error = true;
  }

  if (tfsf_y1 < tfsf_y0) {
    error_msg("TFSF y1 less than y0");
    error = true;
  }

  if (tfsf_z1 < tfsf_z0) {
    error_msg("TFSF z1 less than z0");
    error = true;
  }

  if (tfsf_x1 > sizeX) {
    error_msg("TFSF larger than grid (x-direction)");
    error = true;
  }

  if (tfsf_y1 > sizeY) {
    error_msg("TFSF larger than grid (y-direction)");
    error = true;
  }

  if (tfsf_z1 > sizeZ) {
    error_msg("TFSF larger than grid (z-direction)");
    error = true;
  }

  if ((tfsf_x0 - pml_size) && tfsf_x0 < NLambda) {
    std::string s = "Recommend tfsf_x0 greater than NLambda\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if ((tfsf_y0 - pml_size) && tfsf_y0 < NLambda) {
    std::string s = "Recommend tfsf_y0 greater than NLambda\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if ((tfsf_z0 - pml_size) && tfsf_z0 < NLambda) {
    std::string s = "Recommend tfsf_z0 greater than NLambda\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if ((sizeX - tfsf_x1) < NLambda) {
    std::string s = "Recommend tfsf_x1 spaced more than NLambda from sizeX\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if ((sizeY - tfsf_y1) < NLambda) {
    std::string s = "Recommend tfsf_y1 spaced more than NLambda from sizeY\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if ((sizeZ - tfsf_z1) < NLambda) {
    std::string s = "Recommend tfsf_z1 spaced more than NLambda from sizeZ\n";
    s += "         This is for a useful scatter field response";
    warning_msg(s);
  }

  if (pml_size < 10) {
    warning_msg("Recommend PML size > 10 for good performance");
  }

  return error;
}
