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
   ChimesJITev::NN(x0, e, f0) ;
   break ;
   case 3:
   ChimesJITev::OO(x0, e, f0) ;
   break ;
   case 4:
   ChimesJITev::CH(x0, e, f0) ;
   break ;
   case 5:
   ChimesJITev::CN(x0, e, f0) ;
   break ;
   case 6:
   ChimesJITev::CO(x0, e, f0) ;
   break ;
   case 7:
   ChimesJITev::HN(x0, e, f0) ;
   break ;
   case 8:
   ChimesJITev::HO(x0, e, f0) ;
   break ;
   case 9:
   ChimesJITev::NO(x0, e, f0) ;
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
   ChimesJITev::CCN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 3:
   ChimesJITev::CCO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 4:
   ChimesJITev::CHH(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 5:
   ChimesJITev::CHN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 6:
   ChimesJITev::CHO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 7:
   ChimesJITev::CNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 8:
   ChimesJITev::CNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 9:
   ChimesJITev::COO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 10:
   ChimesJITev::HHH(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 11:
   ChimesJITev::HHN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 12:
   ChimesJITev::HHO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 13:
   ChimesJITev::HNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 14:
   ChimesJITev::HNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 15:
   ChimesJITev::HOO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 16:
   ChimesJITev::NNN(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 17:
   ChimesJITev::NNO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 18:
   ChimesJITev::NOO(x0, x1, x2, e, f0, f1, f2) ;
   break ;
   case 19:
   ChimesJITev::OOO(x0, x1, x2, e, f0, f1, f2) ;
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
   ChimesJITev::CCCC(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 1:
   ChimesJITev::CCCH(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 2:
   ChimesJITev::CCCN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 3:
   ChimesJITev::CCCO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 4:
   ChimesJITev::CCHH(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 5:
   ChimesJITev::CCHN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 6:
   ChimesJITev::CCHO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 7:
   ChimesJITev::CCNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 8:
   ChimesJITev::CCNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 9:
   ChimesJITev::CCOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 10:
   ChimesJITev::CHHH(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 11:
   ChimesJITev::CHHN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 12:
   ChimesJITev::CHHO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 13:
   ChimesJITev::CHNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 14:
   ChimesJITev::CHNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 15:
   ChimesJITev::CHOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 16:
   ChimesJITev::CNNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 17:
   ChimesJITev::CNNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 18:
   ChimesJITev::CNOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 19:
   ChimesJITev::COOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 20:
   ChimesJITev::HHHH(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 21:
   ChimesJITev::HHHN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 22:
   ChimesJITev::HHHO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 23:
   ChimesJITev::HHNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 24:
   ChimesJITev::HHNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 25:
   ChimesJITev::HHOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 26:
   ChimesJITev::HNNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 27:
   ChimesJITev::HNNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 28:
   ChimesJITev::HNOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 29:
   ChimesJITev::HOOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 30:
   ChimesJITev::NNNN(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 31:
   ChimesJITev::NNNO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 32:
   ChimesJITev::NNOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 33:
   ChimesJITev::NOOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
   case 34:
   ChimesJITev::OOOO(x0, x1, x2, x3, x4, x5, e, f0, f1, f2, f3, f4, f5) ;
   break ;
  default: printf("Bad quaduplet index in ChimesJIT") ; exit(1) ; 
  }
}
