#include "chimesJITev.h"
void ChimesJITev::poly_2B(double x0, int pairidx, double *e, double *f0)
{
  switch(pairidx) {
   case 0:
   ChimesJITev::XX(x0, e, f0) ;
   break ;
   case 1:
   ChimesJITev::YY(x0, e, f0) ;
   break ;
   case 2:
   ChimesJITev::XY(x0, e, f0) ;
   break ;
  default: printf("Bad pair index in ChimesJIT") ; exit(1) ; 
  }
}

void ChimesJITev::poly_3B(double x0, double x1, double x2, int tripidx, double *e,
                        double *f0, double *f1, double *f2)
{
  switch(tripidx) {
   case 0:
   ChimesJITev::XXX(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 1:
   ChimesJITev::XXY(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 2:
   ChimesJITev::XYY(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 3:
   ChimesJITev::YYY(x0, x1, x2, e, f0, f1, f2) ;
   break ;
  default: printf("Bad triplet index in ChimesJIT") ; exit(1) ; 
  }
}
 
void ChimesJITev::poly_4B(double x0, double x1, double x2, double x3, double x4, double x5,
                        int quadidx, double *e,
                        double *f0, double *f1, double *f2, double *f3, double *f4, double *f5)
{
  switch(quadidx) {
   case 0:
   ChimesJITev::XXXX(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 1:
   ChimesJITev::XXXY(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 2:
   ChimesJITev::XXYY(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 3:
   ChimesJITev::XYYY(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 4:
   ChimesJITev::YYYY(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
  default: printf("Bad quaduplet index in ChimesJIT") ; exit(1) ; 
  }
}
