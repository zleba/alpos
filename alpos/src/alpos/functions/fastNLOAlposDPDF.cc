// Author: Daniel Britzger
// DESY, 02.03.2017

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO interface for alpos                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>
//#include "alpos/functions/fastNLOAlposDiffH12006FitB.h"
#include "alpos/functions/fastNLOAlposDPDF.h"
#include "alpos/ATheory.h"

using namespace std;

//______________________________________________________________________________

extern "C" {
   // tint is the 'maximum t'
//   void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs, double *tint);  
//   void diffpdferr_(double* xpom, double*  zpom, double*  muf, int* ifit, int* ierr,  double *pdfs); 
   void diffpdferr_(double* xpom, double*  zpom, double*  muf, int* ifit, int* ierr,  int* ireg, double *pdfs);
}  


//______________________________________________________________________________



// vector<double> fastNLOAlposDPDF::GetXFX(double x, double Q) const {
//    //
//    //  GetXFX is used to get the parton array from the
//    //  pre-defined pdf-interface.
//    //

//    // Example access (1)
//    // Mind: Call PAR(PDF) somewhere in advance to guarantee correctly updated values
//    return QUICK(PDF, ({x,Q}) );

//    // Example access (2)
//    //return QUICK_VAR(PDF, 2, x,Q );

//    // Example access (3)
//    // vector<double> xpmuf = {x,Q};
//    // return QUICK_VEC(PDF, xpmuf );

//    // Example access (4)
//    // SET(PDF.xp,x,0);
//    // SET(PDF.muf,Q,0);
//    // vector<double> vals = VALUES(PDF);
//    // return vals;

// }




fastNLOAlposDPDF::fastNLOAlposDPDF(std::string filename) : fastNLODiffReader(filename) , ftint(-1.) {

}

//______________________________________________________________________________

bool fastNLOAlposDPDF::InitPDF() {
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // LHAPDF interface:
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinitialized the set PDF-set.

   // nothing todo
   // PAR(Alpha_s);
   // PAR(PDF);
   return true;
}


//______________________________________________________________________________

double fastNLOAlposDPDF::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   
   return QUICK(Alpha_s, ({Q}) )[0];

   // static const int NF=4; // from h12006B_wrapper.h
   // static const double b0 = (11. - 2./3.*NF);  // The beta coefficients of the QCD beta function
   // static const  double b1 = (51. - 19./3.*NF);

   // //      double lmd = 0.399; // according to matthias
   // static const double lmd = 0.3395; // according to matthias
   // double t = log(Q/lmd);
   // double asMz = 1.0/(b0*t);

   // return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI ;
}


//______________________________________________________________________________



std::vector<double> fastNLOAlposDPDF::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //

   return QUICK(DPDF, ({xpom,zpom,muf}) );

   // -- directly h12006Fit from fortran
   {
      std::vector < double > xfx(13);

      // ----
      int ifit = 2;
      int ierr = 0;
      int ireg = 0;
   
      diffpdferr_(&xpom, &zpom, &muf, &ifit, &ierr,  &ireg, &xfx[0]); 
      // ----
      // double tint = ftint;
      // diffpdf_(&xpom,&zpom,&muf,&xfx[0],&tint);
      // ---
      //debug<<"xpom="<<xpom<<"\tzpom="<<zpom<<"\tmuf="<<muf<<"\tgluon = "<<xfx[6]<<endl;
      return xfx;
   }
}

//______________________________________________________________________________

