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
  return K * (detC - log(detC) - noise);
}

double energyDensity(double c11, double c22, double c12, double beta, double K,
                     double noise) {
  double detC = c11 * c22 - c12 * c12;
  return phi_d(c11, c22, c12, beta) + phi_v(detC, K, noise);
}

Matrix2d stress(double c11, double c22, double c12, double beta, double K) {
  // Generated from python script

  double x0 = c11 * c22;
  double x1 = pow(c12, 2);
  double x2 = x0 - x1;
  double x3 = 1.0 / x2;
  double x4 = pow(x2, -3.0 / 2.0);
  double x5 = pow(c22, 2) * x4;
  double x6 = pow(x2, -1.0 / 2.0);
  double x7 = -3 * x6;
  double x8 = x0 * x4;
  double x9 = x7 + (3.0 / 2.0) * x8;
  double x10 = (3.0 / 2.0) * x5 + x9;
  double x11 = c12 * x6;
  double x12 = c11 * x6;
  double x13 = c22 * x6;
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
  double x32 = x15 * (4 * c12 * c22 * x4 - x31 - x5);
  double x33 = x27 * (3 * x30 + x32);
  double x34 = 0.25 * x30 + 0.083333333333333329 * x32;
  double x35 = x18 * x34;
  double x36 = -2.6666666666666665 * x11 + 2.6666666666666665 * x12 +
               2.6666666666666665 * x13;
  double x37 = 0.090909090909090912 * x36;
  double x38 = 2 * c12;
  double x39 = x38 * x4;
  double x40 = c22 * x39;
  double x41 = -x40;
  double x42 = -x6 + (1.0 / 2.0) * x8;
  double x43 = x17 * (-x41 - x42 - 1.0 / 2.0 * x5);
  double x44 = x23 * (6 * c12 * c22 * x4 - x10);
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
  double x62 = pow(c11, 2) * x4;
  double x63 = (3.0 / 2.0) * x62 + x9;
  double x64 = x16 * (-x29 - x62);
  double x65 = x15 * (4 * c11 * c12 * x4 - x31 - x62);
  double x66 = x27 * (3 * x64 + x65);
  double x67 = 0.25 * x64 + 0.083333333333333329 * x65;
  double x68 = x18 * x67;
  double x69 = -2 * c11 * c12 * x4;
  double x70 = x17 * (-x42 - 1.0 / 2.0 * x62 - x69);
  double x71 = x23 * (6 * c11 * c12 * x4 - x63);
  double x72 = x15 * x64;
  double x73 = x70 - 0.1111111111111111 * x71 + x72;
  double x74 = x18 * (2 * x70 - 0.22222222222222221 * x71 + 2 * x72);
  double x75 = x1 * x4;
  double x76 = c12 * x4;
  double x77 = c11 * x76;
  double x78 = c22 * x76;
  double x79 = 3 * x77 + 3 * x78;
  double x80 = c11 * x39;
  double x81 = x16 * (x41 + x80);
  double x82 = x15 * (x40 - 8 * x6 - 8 * x75 + x80);
  double x83 = x27 * (3 * x81 + x82);
  double x84 = 0.25 * x81 + 0.083333333333333329 * x82;
  double x85 = x18 * x84;
  double x86 = x51 - 4 * x75;
  double x87 = x17 * (x77 + x78 + x86);
  double x88 = x15 * x81;
  double x89 = x23 * (-12 * x6 - 12 * x75 + x79);
  double x90 = x87 + x88 - 0.1111111111111111 * x89;
  double x91 = x18 * (2 * x87 + 2 * x88 - 0.22222222222222221 * x89);

  double dPhi_dC11 =
      beta * (-0.0064709595959595969 * x33 + x34 * x54 + x35 * x60 + x46 * x61 +
              0.000946969696969697 * x50 - x53 * (x41 + 2 * x5 + x52) -
              x58 * (-1.1666666666666665 * c12 * c22 * x4 +
                     1.1666666666666665 * x5 + x56)) -
      c22 * K * (x3 - 1) + (1.0 / 2.0) * x20 * (3 * c12 * c22 * x4 - 2 * x10) +
      x26 * (-1.3333333333333333 * c12 * c22 * x4 + x22 +
             1.3333333333333333 * x5) +
      0.005681818181818182 * x33 - x35 * x37 + x46 * x48 - x46 * x49 +
      0.032196969696969696 * x50;
  double dPhi_dC22 =
      beta * (-x53 * (x52 + 2 * x62 + x69) + x54 * x67 -
              x58 * (-1.1666666666666665 * c11 * c12 * x4 + x56 +
                     1.1666666666666665 * x62) +
              x60 * x68 + x61 * x73 - 0.0064709595959595969 * x66 +
              0.000946969696969697 * x74) -
      c11 * K * (x3 - 1) + (1.0 / 2.0) * x20 * (3 * c11 * c12 * x4 - 2 * x63) +
      x26 * (-1.3333333333333333 * c11 * c12 * x4 + x22 +
             1.3333333333333333 * x62) -
      x37 * x68 + x48 * x73 - x49 * x73 + 0.005681818181818182 * x66 +
      0.032196969696969696 * x74;
  double dPhi_dC12 =
      beta * (x53 * (4 * x77 + 4 * x78 + x86) + x54 * x84 +
              x58 * (x55 - 2.333333333333333 * x75 + 2.333333333333333 * x77 +
                     2.333333333333333 * x78) +
              x60 * x85 + x61 * x90 - 0.0064709595959595969 * x83 +
              0.000946969696969697 * x91) +
      K * (2 * c12 * x3 - x38) + x20 * (x7 - 3 * x75 + x79) -
      x26 * (x21 - 2.6666666666666665 * x75 + 2.6666666666666665 * x77 +
             2.6666666666666665 * x78) -
      x37 * x85 + x48 * x90 - x49 * x90 + 0.005681818181818182 * x83 +
      0.032196969696969696 * x91;

  return Matrix2d{{dPhi_dC11, dPhi_dC12 / 2}, {dPhi_dC12 / 2, dPhi_dC22}};
}

ITensor hessian(double c11, double c22, double c12, double beta, double K) {
  double x0 = c11 * c22;
  double x1 = pow(c12, 2);
  double x2 = x0 - x1;
  double x3 = 1.0 / x2;
  double x4 = pow(x2, -1.0 / 2.0);
  double x5 = c12 * x4;
  double x6 = c11 * x4;
  double x7 = c22 * x4;
  double x8 = x6 + x7;
  double x9 = -x5 + x8;
  double x10 = pow(x9, 2);
  double x11 = pow(x2, -3.0 / 2.0);
  double x12 = pow(c22, 2);
  double x13 = x11 * x12;
  double x14 = -3 * x4;
  double x15 = x0 * x11;
  double x16 = x14 + (3.0 / 2.0) * x15;
  double x17 = (3.0 / 2.0) * x13 + x16;
  double x18 = (3.0 / 2.0) * c12 * c22 * x11 - x17;
  double x19 = x10 * x18;
  double x20 = -4 * x5 + x8;
  double x21 = x6 - x7;
  double x22 = pow(x21, 2);
  double x23 = -0.1111111111111111 * pow(x20, 3) + x20 * x22;
  double x24 = 0.037037037037037028 * x23;
  double x25 = -2.6666666666666665 * x4;
  double x26 = 1.3333333333333333 * x15 + x25;
  double x27 =
      1.3333333333333333 * c12 * c22 * x11 - 1.3333333333333333 * x13 - x26;
  double x28 = pow(x20, 2);
  double x29 = 0.25 * x22 + 0.083333333333333329 * x28;
  double x30 = 0.090909090909090912 * x29;
  double x31 = x23 * x30;
  double x32 = x22 + 0.33333333333333331 * x28;
  double x33 = pow(x32, 2);
  double x34 = 2 * x4;
  double x35 = -x15 + x34;
  double x36 = x13 + x35;
  double x37 = x21 * x36;
  double x38 = -x34;
  double x39 = x15 + x38;
  double x40 = x13 + x39;
  double x41 = 4 * c12 * c22 * x11 - x40;
  double x42 = x20 * x41;
  double x43 = 3 * x37 + x42;
  double x44 = x33 * x43;
  double x45 = 0.25 * x37 + 0.083333333333333329 * x42;
  double x46 = x23 * x45;
  double x47 = -2.6666666666666665 * x5 + 2.6666666666666665 * x6 +
               2.6666666666666665 * x7;
  double x48 = 0.090909090909090912 * x47;
  double x49 = 2 * c12;
  double x50 = c22 * x11;
  double x51 = x49 * x50;
  double x52 = -x51;
  double x53 = (1.0 / 2.0) * x13;
  double x54 = (1.0 / 2.0) * x15;
  double x55 = -x4 + x54;
  double x56 = -x52 - x53 - x55;
  double x57 = x22 * x56;
  double x58 = 6 * c12 * c22 * x11 - x17;
  double x59 = x28 * x58;
  double x60 = x20 * x37;
  double x61 = x57 - 0.1111111111111111 * x59 + x60;
  double x62 = pow(x9, 3);
  double x63 = 0.037037037037037028 * x62;
  double x64 = x30 * x47;
  double x65 = 2 * x57 - 0.22222222222222221 * x59 + 2 * x60;
  double x66 = x23 * x65;
  double x67 = -4 * x4;
  double x68 = 2 * x15 + x67;
  double x69 = -2 * x13 - x52 - x68;
  double x70 = x62 * x69;
  double x71 = 0.012345679012345677 * x29;
  double x72 = 0.012345679012345677 * pow(x9, 4);
  double x73 = 0.015151515151515152 * x29;
  double x74 = -2.333333333333333 * x4;
  double x75 = 1.1666666666666665 * x15 + x74;
  double x76 =
      1.1666666666666665 * c12 * c22 * x11 - 1.1666666666666665 * x13 - x75;
  double x77 = x23 * x76;
  double x78 =
      -2.333333333333333 * x5 + 2.333333333333333 * x6 + 2.333333333333333 * x7;
  double x79 = 0.015151515151515152 * x78;
  double x80 = x73 * x78;
  double x81 = pow(c11, 2);
  double x82 = x11 * x81;
  double x83 = x16 + (3.0 / 2.0) * x82;
  double x84 = (3.0 / 2.0) * c11 * c12 * x11 - x83;
  double x85 = x10 * x24;
  double x86 =
      1.3333333333333333 * c11 * c12 * x11 - x26 - 1.3333333333333333 * x82;
  double x87 = -x35 - x82;
  double x88 = x21 * x87;
  double x89 = x39 + x82;
  double x90 = 4 * c11 * c12 * x11 - x89;
  double x91 = x20 * x90;
  double x92 = 3 * x88 + x91;
  double x93 = x33 * x92;
  double x94 = 0.25 * x88 + 0.083333333333333329 * x91;
  double x95 = x23 * x94;
  double x96 = (1.0 / 2.0) * x82;
  double x97 = -2 * c11 * c12 * x11;
  double x98 = -x55 - x96 - x97;
  double x99 = x22 * x98;
  double x100 = 6 * c11 * c12 * x11 - x83;
  double x101 = x100 * x28;
  double x102 = x20 * x88;
  double x103 = -0.1111111111111111 * x101 + x102 + x99;
  double x104 = -0.22222222222222221 * x101 + 2 * x102 + 2 * x99;
  double x105 = x104 * x23;
  double x106 = -x68 - 2 * x82 - x97;
  double x107 = x62 * x71;
  double x108 =
      1.1666666666666665 * c11 * c12 * x11 - x75 - 1.1666666666666665 * x82;
  double x109 = x23 * x73;
  double x110 = x1 * x11;
  double x111 = c11 * x11;
  double x112 = c12 * x111;
  double x113 = c12 * x50;
  double x114 = 3 * x112 + 3 * x113;
  double x115 = -3 * x110 + x114 + x14;
  double x116 = x111 * x49;
  double x117 = x116 + x52;
  double x118 = x117 * x21;
  double x119 = x116 + x51;
  double x120 = -8 * x110 + x119 - 8 * x4;
  double x121 = x120 * x20;
  double x122 = 3 * x118 + x121;
  double x123 = x122 * x33;
  double x124 = -2.6666666666666665 * x110 + 2.6666666666666665 * x112 +
                2.6666666666666665 * x113 + x25;
  double x125 = 0.25 * x118 + 0.083333333333333329 * x121;
  double x126 = x125 * x23;
  double x127 = -4 * x110 + x67;
  double x128 = x112 + x113 + x127;
  double x129 = x128 * x22;
  double x130 = x118 * x20;
  double x131 = -12 * x110 + x114 - 12 * x4;
  double x132 = x131 * x28;
  double x133 = x129 + x130 - 0.1111111111111111 * x132;
  double x134 = 2 * x129 + 2 * x130 - 0.22222222222222221 * x132;
  double x135 = x134 * x23;
  double x136 = 4 * x112 + 4 * x113 + x127;
  double x137 = -2.333333333333333 * x110 + 2.333333333333333 * x112 +
                2.333333333333333 * x113 + x74;
  double x138 = pow(x2, -2);
  double x139 = K * x138;
  double x140 = pow(x2, -5.0 / 2.0);
  double x141 = x12 * x140;
  double x142 = c12 * x141;
  double x143 = pow(c22, 3) * x140;
  double x144 = 3 * x50;
  double x145 = c11 * x141;
  double x146 = -x144 + (9.0 / 4.0) * x145;
  double x147 = (9.0 / 4.0) * x143 + x146;
  double x148 = 2.6666666666666665 * x50;
  double x149 = 2.0 * x145 - x148;
  double x150 = -x113;
  double x151 = x24 * x9;
  double x152 = x151 * x18;
  double x153 = 0.18181818181818182 * x27;
  double x154 = 2 * x37;
  double x155 = x43 * (x154 + 0.66666666666666663 * x42);
  double x156 = 0.005681818181818182 * x32;
  double x157 = (3.0 / 2.0) * x143;
  double x158 = 2 * x50;
  double x159 = (3.0 / 2.0) * x145;
  double x160 = x158 - x159;
  double x161 = x21 * (-x157 - x160);
  double x162 = x4 - x54;
  double x163 = x36 * (x162 + x53);
  double x164 = -6 * x142;
  double x165 = -x158;
  double x166 = x159 + x165;
  double x167 = x20 * (x157 + x164 + x166);
  double x168 = x41 * x56;
  double x169 = x33 * (3 * x161 + 3 * x163 + x167 + x168);
  double x170 = 0.18181818181818182 * x47;
  double x171 = 0.25 * x161 + 0.25 * x163 + 0.083333333333333329 * x167 +
                0.083333333333333329 * x168;
  double x172 = x171 * x23;
  double x173 = x61 * x65;
  double x174 = 3 * x142;
  double x175 = -x174;
  double x176 = (3.0 / 4.0) * x145 - x50;
  double x177 = x22 * ((3.0 / 4.0) * x143 + x175 + x176);
  double x178 = x28 * (-9 * x142 + x147);
  double x179 = x161 * x20;
  double x180 = x163 * x20;
  double x181 = x42 * x58;
  double x182 = x154 * x56 + x177 - 0.1111111111111111 * x178 + x179 + x180 -
                0.1111111111111111 * x181;
  double x183 = x23 * (2 * x177 - 0.22222222222222221 * x178 + 2 * x179 +
                       2 * x180 - 0.22222222222222221 * x181 + 4 * x37 * x56);
  double x184 = 4 * x50;
  double x185 = 3 * x145 - x184;
  double x186 = x69 * x71;
  double x187 = 2.333333333333333 * x50;
  double x188 = 1.7499999999999998 * x145 - x187;
  double x189 = 0.0064709595959595969 * x32;
  double x190 = 0.030303030303030304 * x61;
  double x191 = x140 * x81;
  double x192 = c12 * x191;
  double x193 = pow(c11, 3) * x140;
  double x194 = 3 * x111;
  double x195 = c22 * x191;
  double x196 = -x194 + (9.0 / 4.0) * x195;
  double x197 = (9.0 / 4.0) * x193 + x196;
  double x198 = 2.6666666666666665 * x111;
  double x199 = 2.0 * x195 - x198;
  double x200 = c11 * c12 * x11 - x89;
  double x201 = x151 * x84;
  double x202 = 0.18181818181818182 * x86;
  double x203 = 0.074074074074074056 * x10;
  double x204 = 2 * x88;
  double x205 = x204 + 0.66666666666666663 * x91;
  double x206 = x205 * x92;
  double x207 = x103 * x29;
  double x208 = 2 * x111;
  double x209 = (3.0 / 2.0) * x193;
  double x210 = (3.0 / 2.0) * x195;
  double x211 = x21 * (x208 + x209 - x210);
  double x212 = -x162 - x96;
  double x213 = x212 * x87;
  double x214 = -6 * x192;
  double x215 = -x208 + x210;
  double x216 = x20 * (x209 + x214 + x215);
  double x217 = x90 * x98;
  double x218 = x33 * (3 * x211 + 3 * x213 + x216 + x217);
  double x219 = x103 * x94;
  double x220 = 0.25 * x211 + 0.25 * x213 + 0.083333333333333329 * x216 +
                0.083333333333333329 * x217;
  double x221 = x220 * x23;
  double x222 = x103 * x104;
  double x223 = 3 * x192;
  double x224 = -x223;
  double x225 = -x111 + (3.0 / 4.0) * x195;
  double x226 = x22 * ((3.0 / 4.0) * x193 + x224 + x225);
  double x227 = x28 * (-9 * x192 + x197);
  double x228 = x20 * x211;
  double x229 = x20 * x213;
  double x230 = x100 * x91;
  double x231 = x204 * x98 + x226 - 0.1111111111111111 * x227 + x228 + x229 -
                0.1111111111111111 * x230;
  double x232 = x23 * (2 * x226 - 0.22222222222222221 * x227 + 2 * x228 +
                       2 * x229 - 0.22222222222222221 * x230 + 4 * x88 * x98);
  double x233 = 4 * x111;
  double x234 = 3 * x195 - x233;
  double x235 = x10 * x71;
  double x236 = x106 * x235;
  double x237 = 2.333333333333333 * x111;
  double x238 = 1.7499999999999998 * x195 - x237;
  double x239 = 0.024691358024691353 * x62;
  double x240 = 0.030303030303030304 * x108;
  double x241 = 0.030303030303030304 * x78;
  double x242 = c12 * x11;
  double x243 = pow(c12, 3) * x140;
  double x244 = x1 * x140;
  double x245 = 9 * x244;
  double x246 = c11 * x245 + c22 * x245 + x144 + x194;
  double x247 = -2 * x110 + x119 + x38;
  double x248 = 8.0 * x244;
  double x249 = 0.18181818181818182 * x124;
  double x250 = 2 * x118;
  double x251 = 0.66666666666666663 * x121 + x250;
  double x252 = x122 * x251;
  double x253 = x133 * x29;
  double x254 = x112 + x150;
  double x255 = x117 * x254;
  double x256 = 6 * x244;
  double x257 = c22 * x256;
  double x258 = c11 * x256 + x208;
  double x259 = x21 * (x165 - x257 + x258);
  double x260 = x120 * x128;
  double x261 = x158 + x257;
  double x262 = x20 * (-24 * x242 - 24 * x243 + x258 + x261);
  double x263 = x33 * (3 * x255 + 3 * x259 + x260 + x262);
  double x264 = x125 * x133;
  double x265 = 0.25 * x255 + 0.25 * x259 + 0.083333333333333329 * x260 +
                0.083333333333333329 * x262;
  double x266 = x23 * x265;
  double x267 = x133 * x134;
  double x268 = x20 * x255;
  double x269 = 3 * x244;
  double x270 = -12 * x242 - 12 * x243;
  double x271 = x22 * (c11 * x269 + c22 * x269 + x111 + x270 + x50);
  double x272 = x20 * x259;
  double x273 = x28 * (-36 * x242 - 36 * x243 + x246);
  double x274 = x121 * x131;
  double x275 = x128 * x250 + x268 + x271 + x272 - 0.1111111111111111 * x273 -
                0.1111111111111111 * x274;
  double x276 = x23 * (4 * x118 * x128 + 2 * x268 + 2 * x271 + 2 * x272 -
                       0.22222222222222221 * x273 - 0.22222222222222221 * x274);
  double x277 = 12 * x244;
  double x278 = c22 * x277 + x184;
  double x279 = c11 * x277 + x233;
  double x280 = 6.9999999999999991 * x244;
  double x281 = 0.030303030303030304 * x137;
  double x282 = c12 * x0 * x140;
  double x283 = x146 + x196;
  double x284 = 0.090909090909090912 * x27;
  double x285 = 0.090909090909090912 * x46;
  double x286 = 0.037037037037037028 * x19;
  double x287 = 0.037037037037037028 * x10;
  double x288 = x287 * x61;
  double x289 = x205 * x43;
  double x290 = x27 * x30;
  double x291 = x30 * x61;
  double x292 = x212 * x36;
  double x293 = x21 * (x160 + x215);
  double x294 = x41 * x98;
  double x295 = 4 * x242 - 6 * x282;
  double x296 = x20 * (x166 + x215 + x295);
  double x297 = x33 * (3 * x292 + 3 * x293 + x294 + x296);
  double x298 = x61 * x94;
  double x299 = x103 * x45;
  double x300 = 0.25 * x292 + 0.25 * x293 + 0.083333333333333329 * x294 +
                0.083333333333333329 * x296;
  double x301 = x23 * x300;
  double x302 = x103 * x65;
  double x303 = x11 * x49 - 3 * x282;
  double x304 = x22 * (x176 + x225 + x303);
  double x305 = x20 * x292;
  double x306 = x20 * x293;
  double x307 = x28 * (6 * x242 - 9 * x282 + x283);
  double x308 = x58 * x91;
  double x309 = x304 + x305 + x306 - 0.1111111111111111 * x307 -
                0.1111111111111111 * x308 + x37 * x98 + x56 * x88;
  double x310 =
      x23 * (x154 * x98 + x204 * x56 + 2 * x304 + 2 * x305 + 2 * x306 -
             0.22222222222222221 * x307 - 0.22222222222222221 * x308);
  double x311 = x10 * x186;
  double x312 = 0.012345679012345677 * x70;
  double x313 = 0.012345679012345677 * x62;
  double x314 = x313 * x45;
  double x315 = 0.015151515151515152 * x77;
  double x316 = 0.015151515151515152 * x46;
  double x317 = x73 * x76;
  double x318 = x61 * x73;
  double x319 = x139 * x49;
  double x320 = c22 * x244;
  double x321 = 3 * x242 - 9.0 / 2.0 * x282;
  double x322 = -9.0 / 2.0 * x142 + x321;
  double x323 = 2.6666666666666665 * x242 - 4.0 * x282;
  double x324 = x251 * x43;
  double x325 = x254 * x36;
  double x326 = x21 * (x174 + x303);
  double x327 = x128 * x41;
  double x328 = x20 * (x175 + x278 + x303);
  double x329 = x33 * (3 * x325 + 3 * x326 + x327 + x328);
  double x330 = x125 * x61;
  double x331 = x133 * x45;
  double x332 = 0.25 * x325 + 0.25 * x326 + 0.083333333333333329 * x327 +
                0.083333333333333329 * x328;
  double x333 = x23 * x332;
  double x334 = x133 * x65;
  double x335 = x242 - 3.0 / 2.0 * x282;
  double x336 = x22 * (-3.0 / 2.0 * x142 + x261 + x335);
  double x337 = x20 * x325;
  double x338 = x20 * x326;
  double x339 = x28 * (18 * x320 + x322 + 6 * x50);
  double x340 = x121 * x58;
  double x341 = x118 * x56 + x128 * x37 + x336 + x337 + x338 -
                0.1111111111111111 * x339 - 0.1111111111111111 * x340;
  double x342 =
      x23 * (x128 * x154 + x250 * x56 + 2 * x336 + 2 * x337 + 2 * x338 -
             0.22222222222222221 * x339 - 0.22222222222222221 * x340);
  double x343 = 2.333333333333333 * x242 - 3.4999999999999996 * x282;
  double x344 = c11 * x244;
  double x345 = -9.0 / 2.0 * x192 + x321;
  double x346 = x251 * x92;
  double x347 = x254 * x87;
  double x348 = x21 * (-x223 - x303);
  double x349 = x128 * x90;
  double x350 = x20 * (x224 + x279 + x303);
  double x351 = x33 * (3 * x347 + 3 * x348 + x349 + x350);
  double x352 = x103 * x125;
  double x353 = x133 * x94;
  double x354 = 0.25 * x347 + 0.25 * x348 + 0.083333333333333329 * x349 +
                0.083333333333333329 * x350;
  double x355 = x23 * x354;
  double x356 = x104 * x133;
  double x357 = x22 * (-3.0 / 2.0 * x192 + x258 + x335);
  double x358 = x20 * x347;
  double x359 = x20 * x348;
  double x360 = x28 * (6 * x111 + 18 * x344 + x345);
  double x361 = x100 * x121;
  double x362 = x118 * x98 + x128 * x88 + x357 + x358 + x359 -
                0.1111111111111111 * x360 - 0.1111111111111111 * x361;
  double x363 =
      x23 * (x128 * x204 + x250 * x98 + 2 * x357 + 2 * x358 + 2 * x359 -
             0.22222222222222221 * x360 - 0.22222222222222221 * x361);

  double dPhi_dC11 =
      beta * (-0.0064709595959595969 * x44 + x45 * x72 + x46 * x79 + x61 * x80 +
              0.000946969696969697 * x66 + x70 * x71 + x73 * x77) -
      c22 * K * (x3 - 1) + x19 * x24 - x27 * x31 + 0.005681818181818182 * x44 -
      x46 * x48 + x61 * x63 - x61 * x64 + 0.032196969696969696 * x66;
  double dPhi_dC22 = beta * (x103 * x80 + 0.000946969696969697 * x105 +
                             x106 * x107 + x108 * x109 + x72 * x94 + x79 * x95 -
                             0.0064709595959595969 * x93) -
                     c11 * K * (x3 - 1) + x103 * x63 - x103 * x64 +
                     0.032196969696969696 * x105 - x31 * x86 - x48 * x95 +
                     x84 * x85 + 0.005681818181818182 * x93;
  double dPhi_dC12 =
      beta *
          (x107 * x136 + x109 * x137 - 0.0064709595959595969 * x123 +
           x125 * x72 + x126 * x79 + x133 * x80 + 0.000946969696969697 * x135) +
      K * (2 * c12 * x3 - x49) + x115 * x85 + 0.005681818181818182 * x123 -
      x124 * x31 - x126 * x48 + x133 * x63 - x133 * x64 +
      0.032196969696969696 * x135;

  double dPhi_dC11_dC11 =
      beta * (x107 * (3 * x143 + x175 + x185) +
              x109 * (-1.7499999999999998 * x142 + 1.7499999999999998 * x143 +
                      x188) -
              x155 * x189 - 0.0064709595959595969 * x169 + x171 * x72 +
              x172 * x79 + 0.000946969696969697 * x173 + x182 * x80 +
              0.000946969696969697 * x183 + x186 * x19 + x190 * x29 * x76 +
              x190 * x45 * x78 + 0.024691358024691353 * x45 * x70 +
              0.030303030303030304 * x45 * x77) +
      x12 * x139 - x152 * (x150 + x40) - x153 * x29 * x61 - x153 * x46 +
      x155 * x156 + 0.005681818181818182 * x169 - x170 * x45 * x61 -
      x172 * x48 + 0.032196969696969696 * x173 + x182 * x63 - x182 * x64 +
      0.032196969696969696 * x183 + 0.074074074074074056 * x19 * x61 -
      x31 * (-2.0 * x142 + 2.0 * x143 + x149) -
      1.0 / 4.0 * x85 * (9 * x142 - 4 * x147);
  double dPhi_dC22_dC22 =
      beta *
          (x106 * x239 * x94 + x107 * (3 * x193 + x224 + x234) +
           x109 *
               (-1.7499999999999998 * x192 + 1.7499999999999998 * x193 + x238) -
           x189 * x206 + x207 * x240 - 0.0064709595959595969 * x218 +
           x219 * x241 + x220 * x72 + x221 * x79 + 0.000946969696969697 * x222 +
           x231 * x80 + 0.000946969696969697 * x232 + x236 * x84 + x240 * x95) +
      x103 * x203 * x84 + x139 * x81 + x156 * x206 - x170 * x219 + x200 * x201 -
      x202 * x207 - x202 * x95 + 0.005681818181818182 * x218 - x221 * x48 +
      0.032196969696969696 * x222 + x231 * x63 - x231 * x64 +
      0.032196969696969696 * x232 - x31 * (-2.0 * x192 + 2.0 * x193 + x199) -
      1.0 / 4.0 * x85 * (9 * x192 - 4 * x197);
  double dPhi_dC12_dC12 =
      beta *
          (x107 * (x270 + x278 + x279) +
           x109 * (c11 * x280 + c22 * x280 + x187 + x237 -
                   6.9999999999999991 * x242 - 6.9999999999999991 * x243) +
           x115 * x136 * x235 + x125 * x136 * x239 + x126 * x281 - x189 * x252 +
           x241 * x264 + x253 * x281 - 0.0064709595959595969 * x263 +
           x265 * x72 + x266 * x79 + 0.000946969696969697 * x267 + x275 * x80 +
           0.000946969696969697 * x276) +
      2 * K * (2 * x1 * x138 + x3 - 1) + x115 * x133 * x203 +
      x115 * x151 * x247 - x126 * x249 + x156 * x252 - x170 * x264 -
      x249 * x253 + 0.005681818181818182 * x263 - x266 * x48 +
      0.032196969696969696 * x267 + x275 * x63 - x275 * x64 +
      0.032196969696969696 * x276 -
      x31 * (c11 * x248 + c22 * x248 + x148 + x198 - 8.0 * x242 - 8.0 * x243) -
      x85 * (9 * x242 + 9 * x243 - x246);
  double dPhi_dC11_dC22 =
      beta *
          (x103 * x317 + x106 * x314 + x107 * (x185 + x234 + x303) +
           x108 * x316 + x108 * x318 +
           x109 * (x188 + x238 + 1.1666666666666665 * x242 -
                   1.7499999999999998 * x282) -
           x189 * x289 - 0.0064709595959595969 * x297 + x298 * x79 +
           x299 * x79 + x300 * x72 + x301 * x79 + 0.000946969696969697 * x302 +
           x309 * x80 + 0.000946969696969697 * x310 + x311 * x84 + x312 * x94 +
           x315 * x94) +
      K * (x0 * x138 - x3 + 1) + x103 * x286 - x103 * x290 + x152 * x200 +
      x156 * x289 - x284 * x95 - x285 * x86 + x288 * x84 - x291 * x86 +
      0.005681818181818182 * x297 - x298 * x48 - x299 * x48 - x301 * x48 +
      0.032196969696969696 * x302 + x309 * x63 - x309 * x64 -
      x31 * (x149 + x199 + 1.3333333333333333 * x242 - 2.0 * x282) +
      0.032196969696969696 * x310 +
      (1.0 / 4.0) * x85 * (6 * x242 - 9 * x282 + 4 * x283);
  double dPhi_dC11_dC12 =
      beta * (x107 * (x164 + x261 + x295) +
              x109 * (-3.4999999999999996 * x142 + 3.4999999999999996 * x320 +
                      x343 + 1.1666666666666665 * x50) +
              x115 * x311 + x125 * x312 + x125 * x315 + x133 * x317 +
              x136 * x314 + x137 * x316 + x137 * x318 - x189 * x324 -
              0.0064709595959595969 * x329 + x330 * x79 + x331 * x79 +
              x332 * x72 + x333 * x79 + 0.000946969696969697 * x334 +
              x341 * x80 + 0.000946969696969697 * x342) -
      c22 * x319 + x115 * x288 - x124 * x285 - x124 * x291 - x126 * x284 +
      x133 * x286 - x133 * x290 + x152 * x247 + x156 * x324 -
      x31 * (-4.0 * x142 + 4.0 * x320 + x323 + 1.3333333333333333 * x50) +
      0.005681818181818182 * x329 - x330 * x48 - x331 * x48 - x333 * x48 +
      0.032196969696969696 * x334 + x341 * x63 - x341 * x64 +
      0.032196969696969696 * x342 +
      (1.0 / 2.0) * x85 * (9 * x320 + 2 * x322 + 3 * x50);
  double dPhi_dC22_dC12 =
      beta *
          (x103 * x137 * x73 + x106 * x125 * x313 +
           x107 * (x214 + x258 + x295) + 0.015151515151515152 * x108 * x126 +
           x108 * x133 * x73 +
           x109 * (1.1666666666666665 * x111 - 3.4999999999999996 * x192 +
                   x343 + 3.4999999999999996 * x344) +
           x115 * x236 + x136 * x313 * x94 + 0.015151515151515152 * x137 * x95 -
           x189 * x346 - 0.0064709595959595969 * x351 + x352 * x79 +
           x353 * x79 + x354 * x72 + x355 * x79 + 0.000946969696969697 * x356 +
           x362 * x80 + 0.000946969696969697 * x363) -
      c11 * x319 + x103 * x115 * x287 - x103 * x124 * x30 -
      0.090909090909090912 * x124 * x95 - 0.090909090909090912 * x126 * x86 +
      x133 * x287 * x84 - x133 * x30 * x86 + x156 * x346 + x201 * x247 -
      x31 * (1.3333333333333333 * x111 - 4.0 * x192 + x323 + 4.0 * x344) +
      0.005681818181818182 * x351 - x352 * x48 - x353 * x48 - x355 * x48 +
      0.032196969696969696 * x356 + x362 * x63 - x362 * x64 +
      0.032196969696969696 * x363 +
      (1.0 / 2.0) * x85 * (3 * x111 + 9 * x344 + 2 * x345);

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

  // Create a 2x2x2 ITensor
  auto i = itensor::Index(2, "i");
  auto j = itensor::Index(2, "j");
  auto k = itensor::Index(2, "k");

  ITensor tensor(i, j, k);

  // Assign values from matrices to the tensor
  tensor.set(i(1), j(1), k(1), topLeft(0, 0));
  tensor.set(i(1), j(1), k(2), topLeft(0, 1));
  tensor.set(i(1), j(2), k(1), topLeft(1, 0));
  tensor.set(i(1), j(2), k(2), topLeft(1, 1));

  tensor.set(i(2), j(1), k(1), topRight(0, 0));
  tensor.set(i(2), j(1), k(2), topRight(0, 1));
  tensor.set(i(2), j(2), k(1), topRight(1, 0));
  tensor.set(i(2), j(2), k(2), topRight(1, 1));

  tensor.set(i(1), j(2), k(1), bottomRight(0, 0));
  tensor.set(i(1), j(2), k(2), bottomRight(0, 1));
  tensor.set(i(2), j(1), k(1), bottomRight(1, 0));
  tensor.set(i(2), j(1), k(2), bottomRight(1, 1));
  return tensor;
}

} // namespace ContiPotential
