
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
   void OO(double x0, double *e, double *f0) ;
   void HH(double x0, double *e, double *f0) ;
   void OH(double x0, double *e, double *f0) ;
   void poly_3B(double x0, double x1, double x2, int tripidx, double *e,
                double *f0, double *f1, double *f2) ;
   void poly_4B(double x0, double x1, double x2, double x3, double x4, double x5,
                int quadidx, double *e,
                double *f0, double *f1, double *f2, double *f3, double *f4, double *f5) ;
} ; 
static std::string jit_commands = \
"! Date  2021-06-09\n"
"!\n"
"! Number of variables            =  27\n"
"! Number of equations            =  28514\n"
"! svd algorithm used\n"
"! eps (= args.eps*dmax)          =   6.4810e-02\n"
"! SVD regularization factor      =  1.0000e-05\n"
"! RMS force error                =  1.1392e+01\n"
"! max abs variable               =  1.2402e+02\n"
"! number of fitting vars         =  25\n"
"! Bayesian Information Criterion =  1.3900e+05\n"
"!\n"
"USECOUL: false\n"
"FITCOUL: false\n"
"USE3BCH: false\n"
"USE4BCH: false\n"
"\n"
"PAIRTYP: CHEBYSHEV  8 0 0 -1 1\n"
"\n"
"ATOM TYPES: 2\n"
"\n"
"# TYPEIDX #	# ATM_TYP #	# ATMCHRG #	# ATMMASS #\n"
"0		O		0.0		15.9994\n"
"1		H		0.0		1.0079\n"
"\n"
"ATOM PAIRS: 3\n"
"\n"
"# PAIRIDX #	# ATM_TY1 #	# ATM_TY1 #	# S_MINIM #	# S_MAXIM #	# CHBDIST #	# MORSE_LAMBDA #\n"
"	0               O               O               0.6             6               MORSE           1.25\n"
"	1               H               H               0.6             6               MORSE           1.25\n"
"	2               O               H               0.6             6               MORSE           1.25\n"
"\n"
"FCUT TYPE: CUBIC\n"
"\n"
"ATOM PAIR TRIPLETS: 0\n"
"ATOM PAIR QUADRUPLETS: 0\n"
"\n"
"PAIR CHEBYSHEV PARAMS\n"
"\n"
"PAIRTYPE PARAMS: 0 O O\n"
"\n"
"  0  -1.2402135452659e+02\n"
"  1  -1.2194250391651e+02\n"
"  2   7.9233371953572e+01\n"
"  3   9.3768161098852e+01\n"
"  4  -3.2951551982120e+01\n"
"  5  -1.0932950044247e+02\n"
"  6  -5.0309959010988e+01\n"
"  7  -1.9491633485732e+01\n"
"\n"
"PAIRTYPE PARAMS: 1 H H\n"
"\n"
"  0  -2.1709290361960e+01\n"
"  1  -9.7664540386521e+00\n"
"  2   1.3757279364043e+01\n"
"  3   3.6426444952133e+01\n"
"  4   4.2229224485115e+01\n"
"  5   3.6175513360451e+01\n"
"  6   1.8620903475289e+01\n"
"  7   6.0796989873361e+00\n"
"\n"
"PAIRTYPE PARAMS: 2 O H\n"
"\n"
"  0   1.4732200706827e+01\n"
"  1   4.4018649095207e+01\n"
"  2   1.1438706910015e+01\n"
"  3   8.9040002428182e+00\n"
"  4  -1.9788436460739e+01\n"
"  5  -6.0378477953627e+00\n"
"  6  -8.8614335450873e+00\n"
"  7   1.6158044249025e-01\n"
"\n"
"TRIPLET CHEBYSHEV PARAMS\n"
"\n"
"QUADRUPLET CHEBYSHEV PARAMS\n"
"\n"
"\n"
"PAIRMAPS: 4\n"
"1 HH\n"
"2 HO\n"
"2 OH\n"
"0 OO\n"
"\n"
"ENDFILE\n"
;
static std::string jit_file = "../../serial_interface/tests/force_fields/test_params.h2o_2bcheby.txt";
static std::string jit_date = "2025-12-22" ;
