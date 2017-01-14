#ifndef SOURCE_H
#define SOURCE_H

#include "constants.h"
#include "field.h"
#include "grid.h"
#include <cstddef>

// ABSTRACT BASE CLASS
class SourceFunc {
public:
  virtual ~SourceFunc() {}
  virtual real calcSource(const std::size_t t, const real offset) = 0;
};

// DERIVED CLASSES
class HarmonicSource : public SourceFunc {
public:
  // Constructor
  HarmonicSource(real Freq, grid *Grid);

  // Function
  real calcSource(const std::size_t t, const real offset = 0);
  // Variables
private:
  real freq;
  grid *g;
};

class GaussianSource : public SourceFunc {
public:
  // Constructor
  GaussianSource(real Freq, grid *Grid);

  // Function
  real calcSource(const std::size_t t, const real offset = 0);
  // Variables
private:
  real tw, t0;
  grid *g;
};

// BASE CLASS
class Source {
public:
  // Constructor
  Source() : f(NULL), g(NULL), source(NULL) {}
  Source(field *F, grid *G, SourceFunc *S) : f(F), g(G), source(S) {}
  virtual ~Source() { delete source; }
  // Functions
  virtual void applySource(const std::size_t t);
  virtual void initialize();
  // Variables
  field *f;
  grid *g;
  SourceFunc *source;
};
// DERIVED CLASSES
class TFSF_Source : public Source {
public:
  TFSF_Source(field *F, grid *G, SourceFunc *S) : f(F), g(G), source(S) {}
  ~TFSF_Source() { delete source; }
  void applySource(const std::size_t t);
  int calcDistance(int i, int j, int k);
  void initialize();
  void calcIncH();
  void calcIncE();
  // variables
  field *f;
  grid *g;
  SourceFunc *source;
  real kinc[3];
  real rcomp[3];
  std::size_t dmax;
  array Hinc;
  array Einc;
};

#endif
