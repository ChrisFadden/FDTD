/** @file ABC.h
 *  @brief ABC implementation for wave equation
 *
 *  @author Chris Fadden
 */

#ifndef _ABC_H
#define _ABC_H

#include "Grid.h"
#include "Field.h"

class ABC {
 public:
  ABC(Grid*, Field*);

  Field* f;
  Grid* g;
  
  double ABC_C1;
  double ABC_C2;
  double ABC_C3;

  void ApplyBC();

};

#endif
