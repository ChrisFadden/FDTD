/** @file Field.h
 *  @brief Class definition of the wave equation field
 *
 *  @author Chris Fadden
 *
 */

#ifndef _FIELD_H
#define _FIELD_H

#include "Grid.h"

class Field {
 public:
  /** @brief Constructor of the Field class
   *
   *  The field is initialized with values
   *  from the  inheritance of the Grid.  Its
   *  values therefore depend on the default
   *  constructor of the Grid.
   */
  Field(Grid*);
  /** @brief Updates the field value
   *
   *    This function is in the time loop,
   *  and updates the field value using finite differences.
   *
   *  @param n the current number of iterations
   *  @return The field value is updated
   */
  void Update(int);
  void ReplaceFields();

  /** @brief Outputs the field values
   *
   *   This function prints the values of the field at all points
   *   to a csv file, so the results can be seen graphically.
   */
  void Print(std::string);
  void PrintLimits();

  // These replace macros, and ease the indexing of 2-D arrays
  double Dim2(int, int);
  double& UOld(int, int);
  double& UNew(int, int);
  double& U(int, int);
  double& Ca(int, int);

  int SizeX();
  int SizeY();
  int isrc();
  int jsrc();
  Grid::SourceType src();

  double dt();
  double cc();
  double dx();
  double Epsr(int);
  double Mur(int);

  std::vector<double> uOld;
  std::vector<double> uNew;
  std::vector<double> u;

  std::vector<double> ca; /**< array of coefficients for the update equation*/

  Grid* g;
  double MaxValue;
  double MinValue;
  int rank;
  int size;
  int start;
  int end;
};
#endif
