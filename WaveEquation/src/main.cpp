#include <iostream>
#include "Grid.h"
#include "Field.h"
#include "ABC.h"

int main(int argc, char** argv) {
  Grid g;
  Field f(&g);
  ABC mur(&g, &f);

  // Movie Variables
  int frameSkip = g.getMaxTime() / 100;
  int p = 0;

  for (int n = 0; n < g.getMaxTime(); n++) {
    
    f.Update(n);
    mur.ApplyBC();
    f.ReplaceFields();

    if(n % frameSkip == 0)
    {
      f.Print(std::to_string(p));
      p++;
    }
    
  }

   f.PrintLimits();
   
   return 0;
}
