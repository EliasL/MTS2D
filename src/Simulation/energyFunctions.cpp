#include "energyFunctions.h"
#include "Eigen/src/Core/Matrix.h"
// #include "itensor/index.h"
// #include "itensor/itensor.h"

// Defined in
// A Variational Model for Reconstructive Phase Transformations in Crystals,
// and their Relation to Dislocations and Plasticity
// Sergio Conti & Giovanni Zanzotto

// Generated by python script
namespace ContiPotential {

double energyDensity(double C11, double C22, double C12, double beta, double K,
                     double noise) {
  double x0 = C11 * C22 - pow(C12, 2);
  double x1 = noise * x0;
  double x2 = pow(x0, -1.0 / 2.0);
  double x3 = C11 * x2;
  double x4 = C22 * x2;
  double x5 = pow(x3 - x4, 2);
  double x6 = C12 * x2;
  double x7 = x3 + x4;
  double x8 = -4 * x6 + x7;
  double x9 = pow(x8, 2);
  double x10 = pow(x5 + 0.33333333333333331 * x9, 3);
  double x11 = x5 * x8 - 0.1111111111111111 * pow(x8, 3);
  double x12 = pow(x11, 2);
  double x13 = -x6 + x7;
  double x14 = 0.25 * x5 + 0.083333333333333329 * x9;
  double x15 = x11 * x14;

  double phi =
      K * (x1 - log(x1)) +
      beta * (-0.0064709595959595969 * x10 + 0.000946969696969697 * x12 +
              0.012345679012345677 * pow(x13, 4) * x14 +
              0.035353535353535352 * x15 * (x3 + x4 - x6)) +
      0.005681818181818182 * x10 + 0.037037037037037028 * x11 * pow(x13, 3) +
      0.032196969696969696 * x12 - 0.24242424242424243 * x15 * (x3 + x4 - x6);
  return phi;
}

Matrix2d stress(double C11, double C22, double C12, double beta, double K,
                double noise) {
  double x0 = C11 * C22;
  double x1 = pow(C12, 2);
  double x2 = x0 - x1;
  double x3 = 1.0 / x2;
  double x4 = pow(x2, -3.0 / 2.0);
  double x5 = pow(C22, 2) * x4;
  double x6 = pow(x2, -1.0 / 2.0);
  double x7 = -3 * x6;
  double x8 = x0 * x4;
  double x9 = x7 + (3.0 / 2.0) * x8;
  double x10 = (3.0 / 2.0) * x5 + x9;
  double x11 = C12 * x6;
  double x12 = C11 * x6;
  double x13 = C22 * x6;
  double x14 = x12 + x13;
  double x15 = -4 * x11 + x14;
  double x16 = x12 - x13;
  double x17 = pow(x16, 2);
  double x18 = -0.1111111111111111 * pow(x15, 3) + x15 * x17;
  double x19 = -x11 + x14;
  double x20 = 0.037037037037037028 * x18 * pow(x19, 2);
  double x21 = -2.6666666666666665 * x6;
  double x22 = x21 + 1.3333333333333333 * x8;
  double x23 = pow(x15, 2);
  double x24 = 0.25 * x17 + 0.083333333333333329 * x23;
  double x25 = 0.090909090909090912 * x24;
  double x26 = x18 * x25;
  double x27 = pow(x17 + 0.33333333333333331 * x23, 2);
  double x28 = 2 * x6;
  double x29 = x28 - x8;
  double x30 = x16 * (x29 + x5);
  double x31 = -x28 + x8;
  double x32 = x15 * (4 * C12 * C22 * x4 - x31 - x5);
  double x33 = x27 * (3 * x30 + 1.0 * x32);
  double x34 = 0.25 * x30 + 0.083333333333333329 * x32;
  double x35 = x18 * x34;
  double x36 = -2.6666666666666665 * x11 + 2.6666666666666665 * x12 +
               2.6666666666666665 * x13;
  double x37 = 0.090909090909090912 * x36;
  double x38 = 2 * C12;
  double x39 = x38 * x4;
  double x40 = C22 * x39;
  double x41 = -x40;
  double x42 = -x6 + (1.0 / 2.0) * x8;
  double x43 = x17 * (-x41 - x42 - 1.0 / 2.0 * x5);
  double x44 = x23 * (6 * C12 * C22 * x4 - x10);
  double x45 = x15 * x30;
  double x46 = x43 - 0.1111111111111111 * x44 + x45;
  double x47 = pow(x19, 3);
  double x48 = 0.037037037037037028 * x47;
  double x49 = x25 * x36;
  double x50 = x18 * (2 * x43 - 0.22222222222222221 * x44 + 2 * x45);
  double x51 = -4 * x6;
  double x52 = x51 + 2 * x8;
  double x53 = 0.012345679012345677 * x24 * x47;
  double x54 = 0.012345679012345677 * pow(x19, 4);
  double x55 = -2.333333333333333 * x6;
  double x56 = x55 + 1.1666666666666665 * x8;
  double x57 = 0.015151515151515152 * x24;
  double x58 = x18 * x57;
  double x59 = -2.333333333333333 * x11 + 2.333333333333333 * x12 +
               2.333333333333333 * x13;
  double x60 = 0.015151515151515152 * x59;
  double x61 = x57 * x59;
  double x62 = pow(C11, 2) * x4;
  double x63 = (3.0 / 2.0) * x62 + x9;
  double x64 = x16 * (-x29 - x62);
  double x65 = x15 * (4 * C11 * C12 * x4 - x31 - x62);
  double x66 = x27 * (3 * x64 + 1.0 * x65);
  double x67 = 0.25 * x64 + 0.083333333333333329 * x65;
  double x68 = x18 * x67;
  double x69 = -2 * C11 * C12 * x4;
  double x70 = x17 * (-x42 - 1.0 / 2.0 * x62 - x69);
  double x71 = x23 * (6 * C11 * C12 * x4 - x63);
  double x72 = x15 * x64;
  double x73 = x70 - 0.1111111111111111 * x71 + x72;
  double x74 = x18 * (2 * x70 - 0.22222222222222221 * x71 + 2 * x72);
  double x75 = x1 * x4;
  double x76 = C12 * x4;
  double x77 = C11 * x76;
  double x78 = C22 * x76;
  double x79 = 3 * x77 + 3 * x78;
  double x80 = C11 * x39;
  double x81 = x16 * (x41 + x80);
  double x82 = x15 * (x40 - 8 * x6 - 8 * x75 + x80);
  double x83 = x27 * (3 * x81 + 1.0 * x82);
  double x84 = 0.25 * x81 + 0.083333333333333329 * x82;
  double x85 = x18 * x84;
  double x86 = x51 - 4 * x75;
  double x87 = x17 * (x77 + x78 + x86);
  double x88 = x15 * x81;
  double x89 = x23 * (-12 * x6 - 12 * x75 + x79);
  double x90 = x87 + x88 - 0.1111111111111111 * x89;
  double x91 = x18 * (2 * x87 + 2 * x88 - 0.22222222222222221 * x89);

  double dPhi_dC11 =
      C22 * K * (noise - x3) +
      beta * (-0.0064709595959595969 * x33 + x34 * x54 + x35 * x60 + x46 * x61 +
              0.000946969696969697 * x50 - x53 * (x41 + 2 * x5 + x52) -
              x58 * (-1.1666666666666665 * C12 * C22 * x4 +
                     1.1666666666666665 * x5 + x56)) +
      (1.0 / 2.0) * x20 * (3 * C12 * C22 * x4 - 2 * x10) +
      x26 * (-1.3333333333333333 * C12 * C22 * x4 + x22 +
             1.3333333333333333 * x5) +
      0.005681818181818182 * x33 - x35 * x37 + x46 * x48 - x46 * x49 +
      0.032196969696969696 * x50;
  double dPhi_dC22 =
      C11 * K * (noise - x3) +
      beta * (-x53 * (x52 + 2 * x62 + x69) + x54 * x67 -
              x58 * (-1.1666666666666665 * C11 * C12 * x4 + x56 +
                     1.1666666666666665 * x62) +
              x60 * x68 + x61 * x73 - 0.0064709595959595969 * x66 +
              0.000946969696969697 * x74) +
      (1.0 / 2.0) * x20 * (3 * C11 * C12 * x4 - 2 * x63) +
      x26 * (-1.3333333333333333 * C11 * C12 * x4 + x22 +
             1.3333333333333333 * x62) -
      x37 * x68 + x48 * x73 - x49 * x73 + 0.005681818181818182 * x66 +
      0.032196969696969696 * x74;
  double dPhi_dC12 =
      K * (2 * C12 * x3 - noise * x38) +
      beta * (x53 * (4 * x77 + 4 * x78 + x86) + x54 * x84 +
              x58 * (x55 - 2.333333333333333 * x75 + 2.333333333333333 * x77 +
                     2.333333333333333 * x78) +
              x60 * x85 + x61 * x90 - 0.0064709595959595969 * x83 +
              0.000946969696969697 * x91) +
      x20 * (x7 - 3 * x75 + x79) -
      x26 * (x21 - 2.6666666666666665 * x75 + 2.6666666666666665 * x77 +
             2.6666666666666665 * x78) -
      x37 * x85 + x48 * x90 - x49 * x90 + 0.005681818181818182 * x83 +
      0.032196969696969696 * x91;

  return Matrix2d{{dPhi_dC11, dPhi_dC12 / 2}, {dPhi_dC12 / 2, dPhi_dC22}};
}
} // namespace ContiPotential