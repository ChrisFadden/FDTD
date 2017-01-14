/** @file calcCurl.cpp
 *  @brief Implements the Curl of the update equations
 */

#include "update.h"
#include <cstddef>
#include <iostream>

void calcCurlEx(field *f, grid *g) {

  for (loop i = 0; i != nx - 1; ++i) {
    for (loop j = 1; j != ny - 1; ++j) {
      for (loop k = 0; k != nz - 1; ++k) {
        curlA(i, j, k) = (Hz(i, j, k) - Hz(i, j - 1, k)) * denEy(j);
        curlB(i, j, k) = (Hy(i, j, k + 1) - Hy(i, j, k)) * denEz(k);
      }
    }
  }

  return;
}

void calcCurlEy(field *f, grid *g) {

  for (loop i = 1; i != nx - 1; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 0; k != nz - 1; ++k) {
        curlA(i, j, k) = (Hz(i - 1, j, k) - Hz(i, j, k)) * denEx(i);
        curlB(i, j, k) = (Hx(i, j, k) - Hx(i, j, k + 1)) * denEz(k);
      }
    }
  }

  return;
}

void calcCurlEz(field *f, grid *g) {
  for (loop i = 1; i != nx - 1; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 0; k != nz - 1; ++k) {
        curlA(i, j, k) = (Hy(i, j, k) - Hy(i - 1, j, k)) * denEx(i);
        curlB(i, j, k) = (Hx(i, j, k) - Hx(i, j - 1, k)) * denEy(j);
      }
    }
  }

  return;
}

void calcCurlHx(field *f, grid *g) {
  for (loop i = 0; i != nx; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 1; k != nz - 1; ++k) {
        curlA(i, j, k) = (Ez(i, j, k) - Ez(i, j + 1, k)) * denHy(j);
        curlB(i, j, k) = (Ey(i, j, k - 1) - Ey(i, j, k)) * denHz(k);
      }
    }
  }
  return;
}
void calcCurlHy(field *f, grid *g) {
  for (loop i = 0; i != nx - 1; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 1; k != nz - 1; ++k) {
        curlA(i, j, k) = (Ez(i + 1, j, k) - Ez(i, j, k)) * denHx(i);
        curlB(i, j, k) = (Ex(i, j, k) - Ex(i, j, k - 1)) * denHz(k);
      }
    }
  }
  return;
}

void calcCurlHz(field *f, grid *g) {
  for (loop i = 0; i != nx - 1; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 0; k != nz - 1; ++k) {
        curlA(i, j, k) = (Ey(i, j, k) - Ey(i + 1, j, k)) * denHx(i);
        curlB(i, j, k) = (Ex(i, j, k) - Ex(i, j + 1, k)) * denHy(j);
      }
    }
  }

  return;
}
