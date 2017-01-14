#include "constants.h"
#include "source.h"
#include <iostream>
#include <math.h>

#define TFSF_offset 5

#define HxInc(i, j, k)                                                         \
  Hinc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (sin(TFSFpsi) * sin(phi) + cos(TFSFpsi) * cos(theta) * cos(phi))

#define HyInc(i, j, k)                                                         \
  Hinc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (-sin(TFSFpsi) * cos(phi) + cos(TFSFpsi) * cos(theta) * sin(phi))

#define HzInc(i, j, k)                                                         \
  Hinc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (-cos(TFSFpsi) * sin(theta))

#define ExInc(i, j, k)                                                         \
  Einc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (cos(TFSFpsi) * sin(phi) - sin(TFSFpsi) * cos(theta) * cos(phi))

#define EyInc(i, j, k)                                                         \
  Einc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (-cos(TFSFpsi) * cos(phi) - sin(TFSFpsi) * cos(theta) * sin(phi))

#define EzInc(i, j, k)                                                         \
  Einc.at(calcDistance((i), (j), (k)) + TFSF_offset) *                         \
      (sin(TFSFpsi) * sin(theta))

#define E_dx curlA(i, 0, 0)
#define H_dy curlB(i, 0, 0)

void TFSF_Source::initialize() {
  kinc[0] = sin(theta) * cos(phi);
  kinc[1] = sin(theta) * cos(theta);
  kinc[2] = cos(theta);

  std::size_t x, y, z;
  x = tfsf_x1 - tfsf_x0;
  y = tfsf_y1 - tfsf_y0;
  z = tfsf_z1 - tfsf_z0;

  // k * r
  dmax = abs(ceil(x * kinc[0] + y * kinc[1] + z * kinc[2]));
  dmax = dmax + 2 * (TFSF_offset + 1);
  Hinc.resize(dmax);
  Einc.resize(dmax);
}

int TFSF_Source::calcDistance(int i, int j, int k) {
  if (theta <= pi / 2.0) {
    if (phi <= pi / 2.0) {
      rcomp[0] = i - (int)tfsf_x0;
      rcomp[1] = j - (int)tfsf_y0;
      rcomp[2] = k - (int)tfsf_z0;
    } else if (phi <= pi) {
      rcomp[0] = i - (int)tfsf_x1;
      rcomp[1] = j - (int)tfsf_y0;
      rcomp[2] = k - (int)tfsf_z0;
    } else if (phi <= 3.0 * pi / 2.0) {
      rcomp[0] = i - (int)tfsf_x1;
      rcomp[1] = j - (int)tfsf_y1;
      rcomp[2] = k - (int)tfsf_z0;
    } else {
      rcomp[0] = i - (int)tfsf_x0;
      rcomp[1] = j - (int)tfsf_y1;
      rcomp[2] = k - (int)tfsf_z0;
    }
  } else if (theta < pi) {
    if (phi <= pi / 2.0) {
      rcomp[0] = i - (int)tfsf_x0;
      rcomp[1] = j - (int)tfsf_y0;
      rcomp[2] = k - (int)tfsf_z1;
    } else if (phi <= pi) {
      rcomp[0] = i - (int)tfsf_x1;
      rcomp[1] = j - (int)tfsf_y0;
      rcomp[2] = k - (int)tfsf_z1;
    } else if (phi <= 3.0 * pi / 2.0) {
      rcomp[0] = i - (int)tfsf_x1;
      rcomp[1] = j - (int)tfsf_y1;
      rcomp[2] = k - (int)tfsf_z1;
    } else {
      rcomp[0] = i - (int)tfsf_x0;
      rcomp[1] = j - (int)tfsf_y1;
      rcomp[2] = k - (int)tfsf_z1;
    }
  }

  return (kinc[0] * rcomp[0] + kinc[1] * rcomp[1] + kinc[2] * rcomp[2]);
}

void TFSF_Source::calcIncE() {

  for (loop i = 0; i != dmax - 1; ++i) {
    curlA(i, 0, 0) = (Einc[i] - Einc[i + 1]) / dx;
  }
}

void TFSF_Source::calcIncH() {

  for (loop i = 1; i != dmax; ++i) {
    curlB(i, 0, 0) = (Hinc[i - 1] - Hinc[i]) / dx;
  }

  return;
}

void TFSF_Source::applySource(const std::size_t t) {
  // Update Fields
  loop i, j, k;

  // Hx
  j = tfsf_y0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Hx(i, j - 1, k) += (Db / dy) * EzInc(i, j, k);
    }
  }

  j = tfsf_y1;
  for (i = tfsf_x0; i != tfsf_x1 + 1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Hx(i, j, k) -= (Db / dy) * EzInc(i, j, k);
    }
  }

  k = tfsf_z0;
  for (i = tfsf_x0; i != tfsf_x1 + 1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1; ++j) {
      Hx(i, j, k - 1) -= (Db / dz) * EyInc(i, j, k);
    }
  }

  k = tfsf_z1;
  for (i = tfsf_x0; i != tfsf_x1 + 1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1; ++j) {
      Hx(i, j, k) += (Db / dz) * EyInc(i, j, k);
    }
  }

  // Hy
  i = tfsf_x0;
  for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Hy(i - 1, j, k) -= (Db / dx) * EzInc(i, j, k);
    }
  }

  i = tfsf_x1;
  for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Hy(i, j, k) += (Db / dx) * EzInc(i, j, k);
    }
  }

  k = tfsf_z0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y1; j != tfsf_y1 + 1; ++j) {
      Hy(i, j, k - 1) += (Db / dz) * ExInc(i, j, k);
    }
  }

  k = tfsf_z1;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
      Hy(i, j, k) -= (Db / dz) * ExInc(i, j, k);
    }
  }

  // Hz
  i = tfsf_x0;
  for (j = tfsf_y0; j != tfsf_y1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Hz(i - 1, j, k) += (Db / dx) * EyInc(i, j, k);
    }
  }

  i = tfsf_x1;
  for (j = tfsf_y0; j != tfsf_y1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Hz(i, j, k) -= (Db / dx) * EyInc(i, j, k);
    }
  }

  j = tfsf_y0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Hz(i, j - 1, k) -= (Db / dy) * ExInc(i, j, k);
    }
  }

  j = tfsf_y1;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Hz(i, j, k) += (Db / dy) * ExInc(i, j, k);
    }
  }

  // Update Incident Fields
  calcIncE();
  Hinc[dmax - 1] = Hinc[dmax - 2];
  for (i = 0; i != dmax - 1; i++) {
    Hinc[i] = Hinc[i] + Db * E_dx;
  }

  Einc[0] = Einc[1];
  Einc[2] = source->calcSource(t, 0);

  calcIncH();
  for (i = 1; i != dmax; i++) {
    Einc[i] = Einc[i] + Cb * H_dy;
  }

  // Ex
  j = tfsf_y0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Ex(i, j, k) -= (Cb / dy) * HzInc(i, j - 1, k);
    }
  }
  j = tfsf_y1;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Ex(i, j, k) += (Cb / dy) * HzInc(i, j, k);
    }
  }

  k = tfsf_z0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
      Ex(i, j, k) += (Cb / dz) * HyInc(i, j, k - 1);
    }
  }

  k = tfsf_z1;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
      Ex(i, j, k) -= (Cb / dz) * HyInc(i, j, k);
    }
  }

  // Ey
  i = tfsf_x0;
  for (j = tfsf_y0; j != tfsf_y1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Ey(i, j, k) += (Cb / dx) * HzInc(i - 1, j, k);
    }
  }

  i = tfsf_x1;
  for (j = tfsf_y0; j != tfsf_y1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1 + 1; ++k) {
      Ey(i, j, k) -= (Cb / dx) * HzInc(i, j, k);
    }
  }

  k = tfsf_z0;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
      Ey(i, j, k) -= (Cb / dz) * HxInc(i, j, k - 1);
    }
  }

  k = tfsf_z1;
  for (i = tfsf_x0; i != tfsf_x1; ++i) {
    for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
      Ey(i, j, k) += (Cb / dz) * HxInc(i, j, k);
    }
  }

  // Ez
  i = tfsf_x0;
  for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
    for (k = tfsf_z0; k != tfsf_z0; ++k) {
      Ez(i, j, k) -= (Cb / dx) * HyInc(i - 1, j, k);
    }
  }

  i = tfsf_x1;
  for (j = tfsf_y0; j != tfsf_y1 + 1; ++j) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Ez(i, j, k) += (Cb / dx) * HyInc(i, j, k);
    }
  }

  j = tfsf_y0;
  for (i = tfsf_x0; i != tfsf_x1 + 1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Ez(i, j, k) += (Cb / dy) * HxInc(i, j - 1, k);
    }
  }

  j = tfsf_y1;
  for (i = tfsf_x0; i != tfsf_x1 + 1; ++i) {
    for (k = tfsf_z0; k != tfsf_z1; ++k) {
      Ez(i, j, k) -= (Cb / dy) * HxInc(i, j, k);
    }
  }
}
