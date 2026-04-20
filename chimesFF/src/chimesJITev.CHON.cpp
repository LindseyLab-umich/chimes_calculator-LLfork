#include "chimesJITev.h"
void ChimesJITev::poly_2B(double x0, int pairidx, double *e, double *f0)
{
  switch(pairidx) {
   case 0:
   ChimesJITev::CC(x0, e, f0) ;
   break ;
   case 1:
   ChimesJITev::HH(x0, e, f0) ;
   break ;
   case 2:
   ChimesJITev::OO(x0, e, f0) ;
   break ;
   case 3:
   ChimesJITev::NN(x0, e, f0) ;
   break ;
   case 4:
   ChimesJITev::CH(x0, e, f0) ;
   break ;
   case 5:
   ChimesJITev::CO(x0, e, f0) ;
   break ;
   case 6:
   ChimesJITev::CN(x0, e, f0) ;
   break ;
   case 7:
   ChimesJITev::HO(x0, e, f0) ;
   break ;
   case 8:
   ChimesJITev::HN(x0, e, f0) ;
   break ;
   case 9:
   ChimesJITev::ON(x0, e, f0) ;
   break ;
  default: printf("Bad pair index in ChimesJIT") ; exit(1) ; 
  }
}

void ChimesJITev::poly_3B(double x0, double x1, double x2, int tripidx, double *e,
                        double *f0, double *f1, double *f2)
{
  switch(tripidx) {
   case 0:
   ChimesJITev::CCC(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 1:
   ChimesJITev::CCH(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 2:
   ChimesJITev::CCO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 3:
   ChimesJITev::CCN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 4:
   ChimesJITev::CHH(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 5:
   ChimesJITev::CHO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 6:
   ChimesJITev::CHN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 7:
   ChimesJITev::COO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 8:
   ChimesJITev::CNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 9:
   ChimesJITev::CNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 10:
   ChimesJITev::HHH(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 11:
   ChimesJITev::HHO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 12:
   ChimesJITev::HHN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 13:
   ChimesJITev::HOO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 14:
   ChimesJITev::HNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 15:
   ChimesJITev::HNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 16:
   ChimesJITev::OOO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 17:
   ChimesJITev::NOO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 18:
   ChimesJITev::NNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 19:
   ChimesJITev::NNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
  default: printf("Bad triplet index in ChimesJIT") ; exit(1) ; 
  }
}
 
void ChimesJITev::poly_4B(double x0, double x1, double x2, double x3, double x4, double x5,
                        int quadidx, double *e,
                        double *f0, double *f1, double *f2, double *f3, double *f4, double *f5)
{
  switch(quadidx) {
  default: printf("Bad quaduplet index in ChimesJIT") ; exit(1) ; 
  }
}
