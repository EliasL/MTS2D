#include "energyFunctions.h"
#include "Eigen/src/Core/Matrix.h"
#include "itensor/index.h"
#include "itensor/itensor.h"

// Defined in
// A Variational Model for Reconstructive Phase Transformations in Crystals,
// and their Relation to Dislocations and Plasticity
// Sergio Conti & Giovanni Zanzotto
namespace ContiPotential {

double I1(double c11, double c22, double c12) {
  return (1.0 / 3.0) * (c11 + c22 - c12);
}

double I2(double c11, double c22, double c12) {
  return (1.0 / 4.0) * pow(c11 - c22, 2) +
         (1.0 / 12.0) * pow(c11 + c22 - 4 * c12, 2);
}

double I3(double c11, double c22, double c12) {
  return pow(c11 - c22, 2) * (c11 + c22 - 4 * c12) -
         (1.0 / 9.0) * pow(c11 + c22 - 4 * c12, 3);
}

double psi1(double I1, double I2, double I3) {
  return pow(I1, 4) * I2 - (41.0 * pow(I2, 3) / 99.0) +
         (7 * I1 * I2 * I3 / 66.0) + pow(I3, 2) / 1056.0;
}

double psi2(double I1, double I2, double I3) {
  return 4.0 * pow(I2, 3) / 11.0 + pow(I1, 3) * I3 -
         (8.0 * I1 * I2 * I3 / 11.0) + (17.0 * pow(I3, 2) / 528.0);
}

double phi_d(double c11, double c22, double c12, double beta) {
  double sqrtDet = sqrt(c11 * c22 - c12 * c12);
  c11 /= sqrtDet;
  c22 /= sqrtDet;
  c12 /= sqrtDet;

  double _I1 = I1(c11, c22, c12);
  double _I2 = I2(c11, c22, c12);
  double _I3 = I3(c11, c22, c12);

  return beta * psi1(_I1, _I2, _I3) + psi2(_I1, _I2, _I3);
}

double phi_v(double detC, double K, double noise) {
  return K * (detC - log(detC)) * noise;
}

double energyDensity(double c11, double c22, double c12, double beta, double K,
                     double noise) {
  double detC = c11 * c22 - c12 * c12;
  return phi_d(c11, c22, c12, beta) + phi_v(detC, K, noise);
}

Matrix2d stress(double c11, double c22, double c12, double beta, double K,
                double noise) {
  // Generated from python script

  double x0 = c11 * c22;
  double x1 = pow(c12, 2);
  double x2 = x0 - x1;
  double x3 = 1.0 / x2;
  double x4 = K * noise;
  double x5 = pow(x2, -3.0 / 2.0);
  double x6 = pow(c22, 2) * x5;
  double x7 = pow(x2, -1.0 / 2.0);
  double x8 = -3 * x7;
  double x9 = x0 * x5;
  double x10 = x8 + (3.0 / 2.0) * x9;
  double x11 = x10 + (3.0 / 2.0) * x6;
  double x12 = c12 * x7;
  double x13 = c11 * x7;
  double x14 = c22 * x7;
  double x15 = x13 + x14;
  double x16 = -4 * x12 + x15;
  double x17 = x13 - x14;
  double x18 = pow(x17, 2);
  double x19 = -0.1111111111111111 * pow(x16, 3) + x16 * x18;
  double x20 = -x12 + x15;
  double x21 = 0.037037037037037028 * x19 * pow(x20, 2);
  double x22 = -2.6666666666666665 * x7;
  double x23 = x22 + 1.3333333333333333 * x9;
  double x24 = pow(x16, 2);
  double x25 = 0.25 * x18 + 0.083333333333333329 * x24;
  double x26 = 0.090909090909090912 * x25;
  double x27 = x19 * x26;
  double x28 = pow(x18 + 0.33333333333333331 * x24, 2);
  double x29 = 2 * x7;
  double x30 = x29 - x9;
  double x31 = x17 * (x30 + x6);
  double x32 = -x29 + x9;
  double x33 = x16 * (4 * c12 * c22 * x5 - x32 - x6);
  double x34 = x28 * (3 * x31 + x33);
  double x35 = 0.25 * x31 + 0.083333333333333329 * x33;
  double x36 = x19 * x35;
  double x37 = -2.6666666666666665 * x12 + 2.6666666666666665 * x13 +
               2.6666666666666665 * x14;
  double x38 = 0.090909090909090912 * x37;
  double x39 = 2 * c12;
  double x40 = x39 * x5;
  double x41 = c22 * x40;
  double x42 = -x41;
  double x43 = -x7 + (1.0 / 2.0) * x9;
  double x44 = x18 * (-x42 - x43 - 1.0 / 2.0 * x6);
  double x45 = x24 * (6 * c12 * c22 * x5 - x11);
  double x46 = x16 * x31;
  double x47 = x44 - 0.1111111111111111 * x45 + x46;
  double x48 = pow(x20, 3);
  double x49 = 0.037037037037037028 * x48;
  double x50 = x26 * x37;
  double x51 = x19 * (2 * x44 - 0.22222222222222221 * x45 + 2 * x46);
  double x52 = -4 * x7;
  double x53 = x52 + 2 * x9;
  double x54 = 0.012345679012345677 * x25 * x48;
  double x55 = 0.012345679012345677 * pow(x20, 4);
  double x56 = -2.333333333333333 * x7;
  double x57 = x56 + 1.1666666666666665 * x9;
  double x58 = 0.015151515151515152 * x25;
  double x59 = x19 * x58;
  double x60 = -2.333333333333333 * x12 + 2.333333333333333 * x13 +
               2.333333333333333 * x14;
  double x61 = 0.015151515151515152 * x60;
  double x62 = x58 * x60;
  double x63 = pow(c11, 2) * x5;
  double x64 = x10 + (3.0 / 2.0) * x63;
  double x65 = x17 * (-x30 - x63);
  double x66 = x16 * (4 * c11 * c12 * x5 - x32 - x63);
  double x67 = x28 * (3 * x65 + x66);
  double x68 = 0.25 * x65 + 0.083333333333333329 * x66;
  double x69 = x19 * x68;
  double x70 = -2 * c11 * c12 * x5;
  double x71 = x18 * (-x43 - 1.0 / 2.0 * x63 - x70);
  double x72 = x24 * (6 * c11 * c12 * x5 - x64);
  double x73 = x16 * x65;
  double x74 = x71 - 0.1111111111111111 * x72 + x73;
  double x75 = x19 * (2 * x71 - 0.22222222222222221 * x72 + 2 * x73);
  double x76 = x1 * x5;
  double x77 = c12 * x5;
  double x78 = c11 * x77;
  double x79 = c22 * x77;
  double x80 = 3 * x78 + 3 * x79;
  double x81 = c11 * x40;
  double x82 = x17 * (x42 + x81);
  double x83 = x16 * (x41 - 8 * x7 - 8 * x76 + x81);
  double x84 = x28 * (3 * x82 + x83);
  double x85 = 0.25 * x82 + 0.083333333333333329 * x83;
  double x86 = x19 * x85;
  double x87 = x52 - 4 * x76;
  double x88 = x18 * (x78 + x79 + x87);
  double x89 = x16 * x82;
  double x90 = x24 * (-12 * x7 - 12 * x76 + x80);
  double x91 = x88 + x89 - 0.1111111111111111 * x90;
  double x92 = x19 * (2 * x88 + 2 * x89 - 0.22222222222222221 * x90);

  double dPhi_dC11 =
      beta * (-0.0064709595959595969 * x34 + x35 * x55 + x36 * x61 + x47 * x62 +
              0.000946969696969697 * x51 - x54 * (x42 + x53 + 2 * x6) -
              x59 * (-1.1666666666666665 * c12 * c22 * x5 + x57 +
                     1.1666666666666665 * x6)) -
      c22 * x4 * (x3 - 1) + (1.0 / 2.0) * x21 * (3 * c12 * c22 * x5 - 2 * x11) +
      x27 * (-1.3333333333333333 * c12 * c22 * x5 + x23 +
             1.3333333333333333 * x6) +
      0.005681818181818182 * x34 - x36 * x38 + x47 * x49 - x47 * x50 +
      0.032196969696969696 * x51;
  double dPhi_dC22 =
      beta * (-x54 * (x53 + 2 * x63 + x70) + x55 * x68 -
              x59 * (-1.1666666666666665 * c11 * c12 * x5 + x57 +
                     1.1666666666666665 * x63) +
              x61 * x69 + x62 * x74 - 0.0064709595959595969 * x67 +
              0.000946969696969697 * x75) -
      c11 * x4 * (x3 - 1) + (1.0 / 2.0) * x21 * (3 * c11 * c12 * x5 - 2 * x64) +
      x27 * (-1.3333333333333333 * c11 * c12 * x5 + x23 +
             1.3333333333333333 * x63) -
      x38 * x69 + x49 * x74 - x50 * x74 + 0.005681818181818182 * x67 +
      0.032196969696969696 * x75;
  double dPhi_dC12 =
      beta * (x54 * (4 * x78 + 4 * x79 + x87) + x55 * x85 +
              x59 * (x56 - 2.333333333333333 * x76 + 2.333333333333333 * x78 +
                     2.333333333333333 * x79) +
              x61 * x86 + x62 * x91 - 0.0064709595959595969 * x84 +
              0.000946969696969697 * x92) +
      x21 * (-3 * x76 + x8 + x80) -
      x27 * (x22 - 2.6666666666666665 * x76 + 2.6666666666666665 * x78 +
             2.6666666666666665 * x79) -
      x38 * x86 + x4 * (2 * c12 * x3 - x39) + x49 * x91 - x50 * x91 +
      0.005681818181818182 * x84 + 0.032196969696969696 * x92;

  return Matrix2d{{dPhi_dC11, dPhi_dC12 / 2}, {dPhi_dC12 / 2, dPhi_dC22}};
}

itensor::ITensor hessian(double c11, double c22, double c12, double beta,
                         double K, double noise) {
  double x0 = c11 * c22;
  double x1 = pow(c12, 2);
  double x2 = x0 - x1;
  double x3 = 1.0 / x2;
  double x4 = K * noise;
  double x5 = pow(x2, -1.0 / 2.0);
  double x6 = c12 * x5;
  double x7 = c11 * x5;
  double x8 = c22 * x5;
  double x9 = x7 + x8;
  double x10 = -x6 + x9;
  double x11 = pow(x10, 2);
  double x12 = pow(x2, -3.0 / 2.0);
  double x13 = pow(c22, 2);
  double x14 = x12 * x13;
  double x15 = -3 * x5;
  double x16 = x0 * x12;
  double x17 = x15 + (3.0 / 2.0) * x16;
  double x18 = (3.0 / 2.0) * x14 + x17;
  double x19 = (3.0 / 2.0) * c12 * c22 * x12 - x18;
  double x20 = x11 * x19;
  double x21 = -4 * x6 + x9;
  double x22 = x7 - x8;
  double x23 = pow(x22, 2);
  double x24 = -0.1111111111111111 * pow(x21, 3) + x21 * x23;
  double x25 = 0.037037037037037028 * x24;
  double x26 = -2.6666666666666665 * x5;
  double x27 = 1.3333333333333333 * x16 + x26;
  double x28 =
      1.3333333333333333 * c12 * c22 * x12 - 1.3333333333333333 * x14 - x27;
  double x29 = pow(x21, 2);
  double x30 = 0.25 * x23 + 0.083333333333333329 * x29;
  double x31 = 0.090909090909090912 * x30;
  double x32 = x24 * x31;
  double x33 = x23 + 0.33333333333333331 * x29;
  double x34 = pow(x33, 2);
  double x35 = 2 * x5;
  double x36 = -x16 + x35;
  double x37 = x14 + x36;
  double x38 = x22 * x37;
  double x39 = -x35;
  double x40 = x16 + x39;
  double x41 = x14 + x40;
  double x42 = 4 * c12 * c22 * x12 - x41;
  double x43 = x21 * x42;
  double x44 = 3 * x38 + x43;
  double x45 = x34 * x44;
  double x46 = 0.25 * x38 + 0.083333333333333329 * x43;
  double x47 = x24 * x46;
  double x48 = -2.6666666666666665 * x6 + 2.6666666666666665 * x7 +
               2.6666666666666665 * x8;
  double x49 = 0.090909090909090912 * x48;
  double x50 = 2 * c12;
  double x51 = c22 * x12;
  double x52 = x50 * x51;
  double x53 = -x52;
  double x54 = (1.0 / 2.0) * x14;
  double x55 = (1.0 / 2.0) * x16;
  double x56 = -x5 + x55;
  double x57 = -x53 - x54 - x56;
  double x58 = x23 * x57;
  double x59 = 6 * c12 * c22 * x12 - x18;
  double x60 = x29 * x59;
  double x61 = x21 * x38;
  double x62 = x58 - 0.1111111111111111 * x60 + x61;
  double x63 = pow(x10, 3);
  double x64 = 0.037037037037037028 * x63;
  double x65 = x31 * x48;
  double x66 = 2 * x58 - 0.22222222222222221 * x60 + 2 * x61;
  double x67 = x24 * x66;
  double x68 = -4 * x5;
  double x69 = 2 * x16 + x68;
  double x70 = -2 * x14 - x53 - x69;
  double x71 = x63 * x70;
  double x72 = 0.012345679012345677 * x30;
  double x73 = 0.012345679012345677 * pow(x10, 4);
  double x74 = 0.015151515151515152 * x30;
  double x75 = -2.333333333333333 * x5;
  double x76 = 1.1666666666666665 * x16 + x75;
  double x77 =
      1.1666666666666665 * c12 * c22 * x12 - 1.1666666666666665 * x14 - x76;
  double x78 = x24 * x77;
  double x79 =
      -2.333333333333333 * x6 + 2.333333333333333 * x7 + 2.333333333333333 * x8;
  double x80 = 0.015151515151515152 * x79;
  double x81 = x74 * x79;
  double x82 = pow(c11, 2);
  double x83 = x12 * x82;
  double x84 = x17 + (3.0 / 2.0) * x83;
  double x85 = (3.0 / 2.0) * c11 * c12 * x12 - x84;
  double x86 = x11 * x25;
  double x87 =
      1.3333333333333333 * c11 * c12 * x12 - x27 - 1.3333333333333333 * x83;
  double x88 = -x36 - x83;
  double x89 = x22 * x88;
  double x90 = x40 + x83;
  double x91 = 4 * c11 * c12 * x12 - x90;
  double x92 = x21 * x91;
  double x93 = 3 * x89 + x92;
  double x94 = x34 * x93;
  double x95 = 0.25 * x89 + 0.083333333333333329 * x92;
  double x96 = x24 * x95;
  double x97 = (1.0 / 2.0) * x83;
  double x98 = -2 * c11 * c12 * x12;
  double x99 = -x56 - x97 - x98;
  double x100 = x23 * x99;
  double x101 = 6 * c11 * c12 * x12 - x84;
  double x102 = x101 * x29;
  double x103 = x21 * x89;
  double x104 = x100 - 0.1111111111111111 * x102 + x103;
  double x105 = 2 * x100 - 0.22222222222222221 * x102 + 2 * x103;
  double x106 = x105 * x24;
  double x107 = -x69 - 2 * x83 - x98;
  double x108 = x63 * x72;
  double x109 =
      1.1666666666666665 * c11 * c12 * x12 - x76 - 1.1666666666666665 * x83;
  double x110 = x24 * x74;
  double x111 = x1 * x12;
  double x112 = c11 * x12;
  double x113 = c12 * x112;
  double x114 = c12 * x51;
  double x115 = 3 * x113 + 3 * x114;
  double x116 = -3 * x111 + x115 + x15;
  double x117 = x112 * x50;
  double x118 = x117 + x53;
  double x119 = x118 * x22;
  double x120 = x117 + x52;
  double x121 = -8 * x111 + x120 - 8 * x5;
  double x122 = x121 * x21;
  double x123 = 3 * x119 + x122;
  double x124 = x123 * x34;
  double x125 = -2.6666666666666665 * x111 + 2.6666666666666665 * x113 +
                2.6666666666666665 * x114 + x26;
  double x126 = 0.25 * x119 + 0.083333333333333329 * x122;
  double x127 = x126 * x24;
  double x128 = -4 * x111 + x68;
  double x129 = x113 + x114 + x128;
  double x130 = x129 * x23;
  double x131 = x119 * x21;
  double x132 = -12 * x111 + x115 - 12 * x5;
  double x133 = x132 * x29;
  double x134 = x130 + x131 - 0.1111111111111111 * x133;
  double x135 = 2 * x130 + 2 * x131 - 0.22222222222222221 * x133;
  double x136 = x135 * x24;
  double x137 = 4 * x113 + 4 * x114 + x128;
  double x138 = -2.333333333333333 * x111 + 2.333333333333333 * x113 +
                2.333333333333333 * x114 + x75;
  double x139 = pow(x2, -2);
  double x140 = x139 * x4;
  double x141 = pow(x2, -5.0 / 2.0);
  double x142 = x13 * x141;
  double x143 = c12 * x142;
  double x144 = pow(c22, 3) * x141;
  double x145 = 3 * x51;
  double x146 = c11 * x142;
  double x147 = -x145 + (9.0 / 4.0) * x146;
  double x148 = (9.0 / 4.0) * x144 + x147;
  double x149 = 2.6666666666666665 * x51;
  double x150 = 2.0 * x146 - x149;
  double x151 = -x114;
  double x152 = x10 * x25;
  double x153 = x152 * x19;
  double x154 = 0.18181818181818182 * x28;
  double x155 = 2 * x38;
  double x156 = x44 * (x155 + 0.66666666666666663 * x43);
  double x157 = 0.005681818181818182 * x33;
  double x158 = (3.0 / 2.0) * x144;
  double x159 = 2 * x51;
  double x160 = (3.0 / 2.0) * x146;
  double x161 = x159 - x160;
  double x162 = x22 * (-x158 - x161);
  double x163 = x5 - x55;
  double x164 = x37 * (x163 + x54);
  double x165 = -6 * x143;
  double x166 = -x159;
  double x167 = x160 + x166;
  double x168 = x21 * (x158 + x165 + x167);
  double x169 = x42 * x57;
  double x170 = x34 * (3 * x162 + 3 * x164 + x168 + x169);
  double x171 = 0.18181818181818182 * x48;
  double x172 = 0.25 * x162 + 0.25 * x164 + 0.083333333333333329 * x168 +
                0.083333333333333329 * x169;
  double x173 = x172 * x24;
  double x174 = x62 * x66;
  double x175 = 3 * x143;
  double x176 = -x175;
  double x177 = (3.0 / 4.0) * x146 - x51;
  double x178 = x23 * ((3.0 / 4.0) * x144 + x176 + x177);
  double x179 = x29 * (-9 * x143 + x148);
  double x180 = x162 * x21;
  double x181 = x164 * x21;
  double x182 = x43 * x59;
  double x183 = x155 * x57 + x178 - 0.1111111111111111 * x179 + x180 + x181 -
                0.1111111111111111 * x182;
  double x184 = x24 * (2 * x178 - 0.22222222222222221 * x179 + 2 * x180 +
                       2 * x181 - 0.22222222222222221 * x182 + 4 * x38 * x57);
  double x185 = 4 * x51;
  double x186 = 3 * x146 - x185;
  double x187 = x70 * x72;
  double x188 = 2.333333333333333 * x51;
  double x189 = 1.7499999999999998 * x146 - x188;
  double x190 = 0.0064709595959595969 * x33;
  double x191 = 0.030303030303030304 * x62;
  double x192 = x141 * x82;
  double x193 = c12 * x192;
  double x194 = pow(c11, 3) * x141;
  double x195 = 3 * x112;
  double x196 = c22 * x192;
  double x197 = -x195 + (9.0 / 4.0) * x196;
  double x198 = (9.0 / 4.0) * x194 + x197;
  double x199 = 2.6666666666666665 * x112;
  double x200 = 2.0 * x196 - x199;
  double x201 = c11 * c12 * x12 - x90;
  double x202 = 0.18181818181818182 * x87;
  double x203 = 0.074074074074074056 * x11;
  double x204 = 2 * x89;
  double x205 = x204 + 0.66666666666666663 * x92;
  double x206 = x205 * x93;
  double x207 = x104 * x30;
  double x208 = 2 * x112;
  double x209 = (3.0 / 2.0) * x194;
  double x210 = (3.0 / 2.0) * x196;
  double x211 = x22 * (x208 + x209 - x210);
  double x212 = -x163 - x97;
  double x213 = x212 * x88;
  double x214 = -6 * x193;
  double x215 = -x208 + x210;
  double x216 = x21 * (x209 + x214 + x215);
  double x217 = x91 * x99;
  double x218 = x34 * (3 * x211 + 3 * x213 + x216 + x217);
  double x219 = x104 * x95;
  double x220 = 0.25 * x211 + 0.25 * x213 + 0.083333333333333329 * x216 +
                0.083333333333333329 * x217;
  double x221 = x220 * x24;
  double x222 = x104 * x105;
  double x223 = 3 * x193;
  double x224 = -x223;
  double x225 = -x112 + (3.0 / 4.0) * x196;
  double x226 = x23 * ((3.0 / 4.0) * x194 + x224 + x225);
  double x227 = x29 * (-9 * x193 + x198);
  double x228 = x21 * x211;
  double x229 = x21 * x213;
  double x230 = x101 * x92;
  double x231 = x204 * x99 + x226 - 0.1111111111111111 * x227 + x228 + x229 -
                0.1111111111111111 * x230;
  double x232 = x24 * (2 * x226 - 0.22222222222222221 * x227 + 2 * x228 +
                       2 * x229 - 0.22222222222222221 * x230 + 4 * x89 * x99);
  double x233 = 4 * x112;
  double x234 = 3 * x196 - x233;
  double x235 = x11 * x72;
  double x236 = x107 * x235;
  double x237 = 2.333333333333333 * x112;
  double x238 = 1.7499999999999998 * x196 - x237;
  double x239 = 0.024691358024691353 * x63;
  double x240 = 0.030303030303030304 * x109;
  double x241 = 0.030303030303030304 * x79;
  double x242 = c12 * x12;
  double x243 = pow(c12, 3) * x141;
  double x244 = x1 * x141;
  double x245 = 9 * x244;
  double x246 = c11 * x245 + c22 * x245 + x145 + x195;
  double x247 = -2 * x111 + x120 + x39;
  double x248 = 8.0 * x244;
  double x249 = 0.18181818181818182 * x125;
  double x250 = 2 * x119;
  double x251 = 0.66666666666666663 * x122 + x250;
  double x252 = x123 * x251;
  double x253 = x134 * x30;
  double x254 = x113 + x151;
  double x255 = x118 * x254;
  double x256 = 6 * x244;
  double x257 = c22 * x256;
  double x258 = c11 * x256 + x208;
  double x259 = x22 * (x166 - x257 + x258);
  double x260 = x121 * x129;
  double x261 = x159 + x257;
  double x262 = x21 * (-24 * x242 - 24 * x243 + x258 + x261);
  double x263 = x34 * (3 * x255 + 3 * x259 + x260 + x262);
  double x264 = x126 * x134;
  double x265 = 0.25 * x255 + 0.25 * x259 + 0.083333333333333329 * x260 +
                0.083333333333333329 * x262;
  double x266 = x24 * x265;
  double x267 = x134 * x135;
  double x268 = x21 * x255;
  double x269 = 3 * x244;
  double x270 = -12 * x242 - 12 * x243;
  double x271 = x23 * (c11 * x269 + c22 * x269 + x112 + x270 + x51);
  double x272 = x21 * x259;
  double x273 = x29 * (-36 * x242 - 36 * x243 + x246);
  double x274 = x122 * x132;
  double x275 = x129 * x250 + x268 + x271 + x272 - 0.1111111111111111 * x273 -
                0.1111111111111111 * x274;
  double x276 = x24 * (4 * x119 * x129 + 2 * x268 + 2 * x271 + 2 * x272 -
                       0.22222222222222221 * x273 - 0.22222222222222221 * x274);
  double x277 = 12 * x244;
  double x278 = c22 * x277 + x185;
  double x279 = c11 * x277 + x233;
  double x280 = 6.9999999999999991 * x244;
  double x281 = 0.030303030303030304 * x138;
  double x282 = c12 * x0 * x141;
  double x283 = x147 + x197;
  double x284 = 0.090909090909090912 * x28;
  double x285 = 0.090909090909090912 * x47;
  double x286 = x205 * x44;
  double x287 = x28 * x31;
  double x288 = x31 * x62;
  double x289 = x212 * x37;
  double x290 = x22 * (x161 + x215);
  double x291 = x42 * x99;
  double x292 = 4 * x242 - 6 * x282;
  double x293 = x21 * (x167 + x215 + x292);
  double x294 = x34 * (3 * x289 + 3 * x290 + x291 + x293);
  double x295 = x62 * x95;
  double x296 = x104 * x46;
  double x297 = 0.25 * x289 + 0.25 * x290 + 0.083333333333333329 * x291 +
                0.083333333333333329 * x293;
  double x298 = x24 * x297;
  double x299 = x104 * x66;
  double x300 = x12 * x50 - 3 * x282;
  double x301 = x23 * (x177 + x225 + x300);
  double x302 = x21 * x289;
  double x303 = x21 * x290;
  double x304 = x29 * (6 * x242 - 9 * x282 + x283);
  double x305 = x59 * x92;
  double x306 = x301 + x302 + x303 - 0.1111111111111111 * x304 -
                0.1111111111111111 * x305 + x38 * x99 + x57 * x89;
  double x307 =
      x24 * (x155 * x99 + x204 * x57 + 2 * x301 + 2 * x302 + 2 * x303 -
             0.22222222222222221 * x304 - 0.22222222222222221 * x305);
  double x308 = x11 * x187;
  double x309 = 0.012345679012345677 * x71;
  double x310 = 0.012345679012345677 * x63;
  double x311 = x310 * x46;
  double x312 = 0.015151515151515152 * x78;
  double x313 = 0.015151515151515152 * x47;
  double x314 = x74 * x77;
  double x315 = x62 * x74;
  double x316 = x140 * x50;
  double x317 = c22 * x244;
  double x318 = 3 * x242 - 9.0 / 2.0 * x282;
  double x319 = -9.0 / 2.0 * x143 + x318;
  double x320 = 2.6666666666666665 * x242 - 4.0 * x282;
  double x321 = x254 * x37;
  double x322 = x22 * (x175 + x300);
  double x323 = x129 * x42;
  double x324 = x21 * (x176 + x278 + x300);
  double x325 = 3 * x321 + 3 * x322 + x323 + x324;
  double x326 = x126 * x62;
  double x327 = x134 * x46;
  double x328 = 0.25 * x321 + 0.25 * x322 + 0.083333333333333329 * x323 +
                0.083333333333333329 * x324;
  double x329 = x24 * x328;
  double x330 = x242 - 3.0 / 2.0 * x282;
  double x331 = x23 * (-3.0 / 2.0 * x143 + x261 + x330);
  double x332 = x21 * x321;
  double x333 = x21 * x322;
  double x334 = x29 * (18 * x317 + x319 + 6 * x51);
  double x335 = x122 * x59;
  double x336 = x119 * x57 + x129 * x38 + x331 + x332 + x333 -
                0.1111111111111111 * x334 - 0.1111111111111111 * x335;
  double x337 = x129 * x155 + x250 * x57 + 2 * x331 + 2 * x332 + 2 * x333 -
                0.22222222222222221 * x334 - 0.22222222222222221 * x335;
  double x338 = 2.333333333333333 * x242 - 3.4999999999999996 * x282;
  double x339 = c11 * x244;
  double x340 = -9.0 / 2.0 * x193 + x318;
  double x341 = x254 * x88;
  double x342 = x22 * (-x223 - x300);
  double x343 = x129 * x91;
  double x344 = x21 * (x224 + x279 + x300);
  double x345 = 3 * x341 + 3 * x342 + x343 + x344;
  double x346 = x104 * x126;
  double x347 = x134 * x95;
  double x348 = 0.25 * x341 + 0.25 * x342 + 0.083333333333333329 * x343 +
                0.083333333333333329 * x344;
  double x349 = x24 * x348;
  double x350 = x23 * (-3.0 / 2.0 * x193 + x258 + x330);
  double x351 = x21 * x341;
  double x352 = x21 * x342;
  double x353 = x29 * (6 * x112 + 18 * x339 + x340);
  double x354 = x101 * x122;
  double x355 = x119 * x99 + x129 * x89 + x350 + x351 + x352 -
                0.1111111111111111 * x353 - 0.1111111111111111 * x354;
  double x356 = x129 * x204 + x250 * x99 + 2 * x350 + 2 * x351 + 2 * x352 -
                0.22222222222222221 * x353 - 0.22222222222222221 * x354;

  double dPhi_dC11 =
      beta * (-0.0064709595959595969 * x45 + x46 * x73 + x47 * x80 + x62 * x81 +
              0.000946969696969697 * x67 + x71 * x72 + x74 * x78) -
      c22 * x4 * (x3 - 1) + x20 * x25 - x28 * x32 + 0.005681818181818182 * x45 -
      x47 * x49 + x62 * x64 - x62 * x65 + 0.032196969696969696 * x67;
  double dPhi_dC22 = beta * (x104 * x81 + 0.000946969696969697 * x106 +
                             x107 * x108 + x109 * x110 + x73 * x95 + x80 * x96 -
                             0.0064709595959595969 * x94) -
                     c11 * x4 * (x3 - 1) + x104 * x64 - x104 * x65 +
                     0.032196969696969696 * x106 - x32 * x87 - x49 * x96 +
                     x85 * x86 + 0.005681818181818182 * x94;
  double dPhi_dC12 =
      beta *
          (x108 * x137 + x110 * x138 - 0.0064709595959595969 * x124 +
           x126 * x73 + x127 * x80 + x134 * x81 + 0.000946969696969697 * x136) +
      x116 * x86 + 0.005681818181818182 * x124 - x125 * x32 - x127 * x49 +
      x134 * x64 - x134 * x65 + 0.032196969696969696 * x136 +
      x4 * (2 * c12 * x3 - x50);
  double dPhi_dC11_dC11 =
      beta * (x108 * (3 * x144 + x176 + x186) +
              x110 * (-1.7499999999999998 * x143 + 1.7499999999999998 * x144 +
                      x189) -
              x156 * x190 - 0.0064709595959595969 * x170 + x172 * x73 +
              x173 * x80 + 0.000946969696969697 * x174 + x183 * x81 +
              0.000946969696969697 * x184 + x187 * x20 + x191 * x30 * x77 +
              x191 * x46 * x79 + 0.024691358024691353 * x46 * x71 +
              0.030303030303030304 * x46 * x78) +
      x13 * x140 - x153 * (x151 + x41) - x154 * x30 * x62 - x154 * x47 +
      x156 * x157 + 0.005681818181818182 * x170 - x171 * x46 * x62 -
      x173 * x49 + 0.032196969696969696 * x174 + x183 * x64 - x183 * x65 +
      0.032196969696969696 * x184 + 0.074074074074074056 * x20 * x62 -
      x32 * (-2.0 * x143 + 2.0 * x144 + x150) -
      1.0 / 4.0 * x86 * (9 * x143 - 4 * x148);
  double dPhi_dC22_dC22 =
      beta *
          (x107 * x239 * x95 + x108 * (3 * x194 + x224 + x234) +
           x110 *
               (-1.7499999999999998 * x193 + 1.7499999999999998 * x194 + x238) -
           x190 * x206 + x207 * x240 - 0.0064709595959595969 * x218 +
           x219 * x241 + x220 * x73 + x221 * x80 + 0.000946969696969697 * x222 +
           x231 * x81 + 0.000946969696969697 * x232 + x236 * x85 + x240 * x96) +
      x104 * x203 * x85 + x140 * x82 + x152 * x201 * x85 + x157 * x206 -
      x171 * x219 - x202 * x207 - x202 * x96 + 0.005681818181818182 * x218 -
      x221 * x49 + 0.032196969696969696 * x222 + x231 * x64 - x231 * x65 +
      0.032196969696969696 * x232 - x32 * (-2.0 * x193 + 2.0 * x194 + x200) -
      1.0 / 4.0 * x86 * (9 * x193 - 4 * x198);
  double dPhi_dC12_dC12 =
      beta *
          (x108 * (x270 + x278 + x279) +
           x110 * (c11 * x280 + c22 * x280 + x188 + x237 -
                   6.9999999999999991 * x242 - 6.9999999999999991 * x243) +
           x116 * x137 * x235 + x126 * x137 * x239 + x127 * x281 - x190 * x252 +
           x241 * x264 + x253 * x281 - 0.0064709595959595969 * x263 +
           x265 * x73 + x266 * x80 + 0.000946969696969697 * x267 + x275 * x81 +
           0.000946969696969697 * x276) +
      x116 * x134 * x203 + x116 * x152 * x247 - x127 * x249 + x157 * x252 -
      x171 * x264 - x249 * x253 + 0.005681818181818182 * x263 - x266 * x49 +
      0.032196969696969696 * x267 + x275 * x64 - x275 * x65 +
      0.032196969696969696 * x276 -
      x32 * (c11 * x248 + c22 * x248 + x149 + x199 - 8.0 * x242 - 8.0 * x243) +
      2 * x4 * (2 * x1 * x139 + x3 - 1) - x86 * (9 * x242 + 9 * x243 - x246);
  double dPhi_dC11_dC22 =
      beta *
          (x104 * x314 + x107 * x311 + x108 * (x186 + x234 + x300) +
           x109 * x313 + x109 * x315 +
           x110 * (x189 + x238 + 1.1666666666666665 * x242 -
                   1.7499999999999998 * x282) -
           x190 * x286 - 0.0064709595959595969 * x294 + x295 * x80 +
           x296 * x80 + x297 * x73 + x298 * x80 + 0.000946969696969697 * x299 +
           x306 * x81 + 0.000946969696969697 * x307 + x308 * x85 + x309 * x95 +
           x312 * x95) +
      0.037037037037037028 * x104 * x20 - x104 * x287 +
      0.037037037037037028 * x11 * x62 * x85 + x153 * x201 + x157 * x286 -
      x284 * x96 - x285 * x87 - x288 * x87 + 0.005681818181818182 * x294 -
      x295 * x49 - x296 * x49 - x298 * x49 + 0.032196969696969696 * x299 +
      x306 * x64 - x306 * x65 + 0.032196969696969696 * x307 -
      x32 * (x150 + x200 + 1.3333333333333333 * x242 - 2.0 * x282) +
      x4 * (x0 * x139 - x3 + 1) +
      (1.0 / 4.0) * x86 * (6 * x242 - 9 * x282 + 4 * x283);
  double dPhi_dC11_dC12 =
      beta *
          (x108 * (x165 + x261 + x292) +
           x110 * (-3.4999999999999996 * x143 + 3.4999999999999996 * x317 +
                   x338 + 1.1666666666666665 * x51) +
           x116 * x308 + x126 * x309 + x126 * x312 + x134 * x314 +
           0.000946969696969697 * x134 * x66 + x137 * x311 + x138 * x313 +
           x138 * x315 - x190 * x251 * x44 + 0.000946969696969697 * x24 * x337 -
           0.0064709595959595969 * x325 * x34 + x326 * x80 + x327 * x80 +
           x328 * x73 + x329 * x80 + x336 * x81) -
      c22 * x316 + 0.037037037037037028 * x10 * x19 * x24 * x247 +
      0.037037037037037028 * x11 * x116 * x62 +
      0.037037037037037028 * x11 * x134 * x19 +
      0.018518518518518514 * x11 * x24 * (9 * x317 + 2 * x319 + 3 * x51) -
      x125 * x285 - x125 * x288 - x127 * x284 - x134 * x287 +
      0.032196969696969696 * x134 * x66 + 0.032196969696969696 * x24 * x337 +
      0.005681818181818182 * x251 * x33 * x44 -
      x32 * (-4.0 * x143 + 4.0 * x317 + x320 + 1.3333333333333333 * x51) +
      0.005681818181818182 * x325 * x34 - x326 * x49 - x327 * x49 - x329 * x49 +
      0.037037037037037028 * x336 * x63 - x336 * x65;
  double dPhi_dC22_dC12 =
      beta *
          (x104 * x138 * x74 + 0.000946969696969697 * x105 * x134 +
           x107 * x126 * x310 + x108 * (x214 + x258 + x292) +
           0.015151515151515152 * x109 * x127 + x109 * x134 * x74 +
           x110 * (1.1666666666666665 * x112 - 3.4999999999999996 * x193 +
                   x338 + 3.4999999999999996 * x339) +
           x116 * x236 + x137 * x310 * x95 + 0.015151515151515152 * x138 * x96 -
           x190 * x251 * x93 + 0.000946969696969697 * x24 * x356 -
           0.0064709595959595969 * x34 * x345 + x346 * x80 + x347 * x80 +
           x348 * x73 + x349 * x80 + x355 * x81) -
      c11 * x316 + 0.037037037037037028 * x10 * x24 * x247 * x85 +
      0.037037037037037028 * x104 * x11 * x116 - x104 * x125 * x31 +
      0.032196969696969696 * x105 * x134 +
      0.037037037037037028 * x11 * x134 * x85 +
      0.018518518518518514 * x11 * x24 * (3 * x112 + 9 * x339 + 2 * x340) -
      0.090909090909090912 * x125 * x96 - 0.090909090909090912 * x127 * x87 -
      x134 * x31 * x87 + 0.032196969696969696 * x24 * x356 +
      0.005681818181818182 * x251 * x33 * x93 -
      x32 * (1.3333333333333333 * x112 - 4.0 * x193 + x320 + 4.0 * x339) +
      0.005681818181818182 * x34 * x345 - x346 * x49 - x347 * x49 - x349 * x49 +
      0.037037037037037028 * x355 * x63 - x355 * x65;

  // First derivative
  Matrix2d stress = Matrix2d{{dPhi_dC11, dPhi_dC12 / 2}, //
                             {dPhi_dC12 / 2, dPhi_dC22}};

  // Second derivative
  Matrix2d topLeft = Matrix2d{{dPhi_dC11_dC11, dPhi_dC11_dC12 / 2},
                              {dPhi_dC11_dC12 / 2, dPhi_dC11_dC22}};

  Matrix2d topRight = Matrix2d{{dPhi_dC11_dC12 / 2, dPhi_dC12_dC12 / 4},
                               {dPhi_dC12_dC12 / 4, dPhi_dC22_dC12 / 2}};

  Matrix2d bottomRight = Matrix2d{{dPhi_dC11_dC22, dPhi_dC11_dC12 / 2},
                                  {dPhi_dC11_dC12 / 2, dPhi_dC22_dC22}};

  std::vector<std::vector<Matrix2d>> matrices = {{topLeft, topRight},
                                                 {stress, bottomRight}};

  // Create a 2x2x2x2 ITensor
  auto i = itensor::Index(2, "i");
  auto j = itensor::Index(2, "j");
  auto k = itensor::Index(2, "k");
  auto l = itensor::Index(2, "l");

  itensor::ITensor tensor(i, j, k, l);

  for (int i1 = 1; i1 <= 2; ++i1) {
    for (int j1 = 1; j1 <= 2; ++j1) {
      for (int k1 = 1; k1 <= 2; ++k1) {
        for (int l1 = 1; l1 <= 2; ++l1) {
          tensor.set(i(i1), j(j1), k(k1), l(l1),
                     matrices[i1 - 1][j1 - 1](k1 - 1, l1 - 1));
        }
      }
    }
  }
  return tensor;
}

} // namespace ContiPotential
