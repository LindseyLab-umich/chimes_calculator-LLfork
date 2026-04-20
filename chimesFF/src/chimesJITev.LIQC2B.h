
#include <stdio.h>
#include <stdlib.h>
#include <string>

inline double Power(double X, int pow)
{
  switch (pow) {
  case 0: return 1.0 ;
  case 1: return X ;
  case 2: return X*X ;
  case 3: return X*X*X ;
  case 4: return X*X*X*X ;
  case 5: return X*X*X*X*X ;
  case 6: return X*X*X*X*X*X ;
  case 7: return X*X*X*X*X*X*X ;
  case 8: return X*X*X*X*X*X*X*X ;        
  default:
     return X * Power(X,pow-1) ;
  }
}

struct ChimesJITev {

   void poly_2B(double x0, int pairidx, double *e, double *f0) ;
   void CC(double x0, double *e, double *f0) ;
   void poly_3B(double x0, double x1, double x2, int tripidx, double *e,
                double *f0, double *f1, double *f2) ;
   void poly_4B(double x0, double x1, double x2, double x3, double x4, double x5,
                int quadidx, double *e,
                double *f0, double *f1, double *f2, double *f3, double *f4, double *f5) ;
} ; 
static std::string jit_commands = \
"! https://doi.org/10.1021/acs.jctc.7b00867\n"
"!\n"
"! 2-body only liquid carbon model fit (no stress tensor in fit)\n"
"!\n"
"! WARNING: This force field is only intended for liquid carbon from 5000 and 2.43 gcc\n"
"! WARNING: This force field comes with no guarantees\n"
"!\n"
"USECOUL: false\n"
"FITCOUL: false\n"
"USEPOVR: false\n"
"FITPOVR: false\n"
"USE3BCH: false\n"
"USE4BCH: false\n"
"\n"
"PAIRTYP: CHEBYSHEV  12 0 0 -1 1\n"
"\n"
"ATOM TYPES: 1\n"
"\n"
"# TYPEIDX #	# ATM_TYP #	# ATMCHRG #	# ATMMASS #\n"
"0		C		0		12.011\n"
"\n"
"ATOM PAIRS: 1\n"
"\n"
"!# PAIRIDX #	# ATM_TY1 #	# ATM_TY1 #	# S_MINIM #	# S_MAXIM #	# S_DELTA #	# CHBDIST #	# MORSE_LAMBDA #\n"
"	0               C               C               1.0         3.15               0.01            MORSE           1.25\n"
"\n"
"FCUT TYPE: CUBIC\n"
"\n"
"PAIR CHEBYSHEV PENALTY DIST: 0.01\n"
"PAIR CHEBYSHEV PENALTY SCALING: 1E8\n"
"\n"
"ATOM PAIR TRIPLETS: 0\n"
"\n"
"PAIR CHEBYSHEV PARAMS\n"
"\n"
"PAIRTYPE PARAMS: 0 C C\n"
"\n"
"0  285.73883308072\n"
"1  -213.71388752372\n"
"2  358.53331099031\n"
"3  -172.12400486549\n"
"4  44.775023503150\n"
"5  -34.154784921509\n"
"6  30.632345544482\n"
"7  -33.336059893072\n"
"8  11.483163813684\n"
"9  -0.99086720791180\n"
"10 -3.3830138904188\n"
"11 1.2480108628453\n"
"\n"
"PAIRMAPS: 1\n"
"0 CC\n"
"\n"
"ENDFILE\n"
"\n"
;
static std::string jit_file = "../../serial_interface/tests/force_fields/published_params.liqC.2b.cubic.txt";
static std::string jit_date = "2025-12-22" ;
