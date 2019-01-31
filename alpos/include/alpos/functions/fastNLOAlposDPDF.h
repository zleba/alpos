// DB 02/01/2019

#ifndef FASTNLOALPOSDPDF
#define FASTNLOALPOSDPDF


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "fastnlotk/fastNLODiffReader.h"

//______________________________________________________________________________


// extern "C" {
//    // tint is the 'maximum t'
// //   void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs, double *tint);
//    void diffpdferr_(double* xpom, double*  zpom, double*  muf, int* ifit, int* ierr,  double *pdfs);
// }


//______________________________________________________________________________


class fastNLOAlposDPDF : public fastNLODiffReader {

public:

   fastNLOAlposDPDF(std::string filename);// : fastNLODiffReader(filename) , ftint(-1.) { };
   ~fastNLOAlposDPDF(void) {
      ;
   };

   void SetTIntegratedRange(double tmax) { ftint = tmax;}
   double GetTIntegratedRange() const {return  ftint;}

   void SetAlposName(const std::string& name) {fAlposName = name;}
   const std::string& GetAlposName() const { return fAlposName;}

   std::string fAlposName;
protected:

   double ftint;
   // inherited functions
   double EvolveAlphas(double Q) const;
   bool InitPDF() ;
   std::vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;
};



//______________________________________________________________________________




// fastNLOAlposDPDF::fastNLOAlposDPDF(std::string filename) : fastNLODiffReader(filename) , ftint(-1.) {
// }


//______________________________________________________________________________

/*
double fastNLOAlposDPDF::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //

   static const int NF=4; // from h12006B_wrapper.h
   static const double b0 = (11. - 2./3.*NF);  // The beta coefficients of the QCD beta function
   static const  double b1 = (51. - 19./3.*NF);

   //      double lmd = 0.399; // according to matthias
   static const double lmd = 0.3395; // according to matthias
   double t = log(Q/lmd);
   double asMz = 1.0/(b0*t);

   return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI ;
}


//______________________________________________________________________________


bool fastNLOAlposDPDF::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // nothing todo!
   return true;
}


//______________________________________________________________________________



std::vector<double> fastNLOAlposDPDF::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //
   std::vector < double > xfx(13);
   double tint = ftint;
   diffpdf_(&xpom,&zpom,&muf,&xfx[0],&tint);
   //debug<<"xpom="<<xpom<<"\tzpom="<<zpom<<"\tmuf="<<muf<<"\tgluon = "<<xfx[6]<<endl;
   return xfx;
}


//______________________________________________________________________________

*/

#endif
