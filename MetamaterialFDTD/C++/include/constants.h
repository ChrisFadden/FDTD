#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <vector>

// Numeric Constants
#define cc (299792458)
#define cc2 (89875517873681764)
#define mu0 (1.2566370614359172953057e-6)
#define eps0 (8.854187817620389850536e-12)
#define eta0 (376.73031346177)

#define pi (3.14159265358979323846)
#define ns (1e-9)
#define us (1e-6)
#define ms (1e-3)
#define GHz (1e9)
#define MHz (1e6)
#define kHz (1e3)

// Simulation Constants
#define CFL 0.57735026919
#define NLambda 5

// PML constants
#define pml_size 10
#define pml_m 3
#define pml_ma 1
#define alpha_max 0.24
#define kappa_max 15.0

// Custom Numeric Types
typedef float real;
typedef std::complex<real> cplx;
typedef std::vector<real> array;
typedef std::size_t loop;
#endif
