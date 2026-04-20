#include "chimesJITev.h"
void ChimesJITev::poly_2B(double x0, int pairidx, double *e, double *f0)
{
  switch(pairidx) {
   case 0:
   ChimesJITev::NN(x0, e, f0) ;
   break ;
   case 1:
   ChimesJITev::HH(x0, e, f0) ;
   break ;
   case 2:
   ChimesJITev::NH(x0, e, f0) ;
   break ;
  default: printf("Bad pair index in ChimesJIT") ; exit(1) ; 
  }
}

void ChimesJITev::poly_3B(double x0, double x1, double x2, int tripidx, double *e,
                        double *f0, double *f1, double *f2)
{
  switch(tripidx) {
   case 0:
   ChimesJITev::NNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 1:
   ChimesJITev::HNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 2:
   ChimesJITev::HHN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 3:
   ChimesJITev::HHH(x0, x1, x2, e, f0, f1, f2) ;
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
   ChimesJITev::NNNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 1:
   ChimesJITev::HNNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 2:
   ChimesJITev::HHNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 3:
   ChimesJITev::HHHN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 4:
   ChimesJITev::HHHH(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
  default: printf("Bad quaduplet index in ChimesJIT") ; exit(1) ; 
  }
}
