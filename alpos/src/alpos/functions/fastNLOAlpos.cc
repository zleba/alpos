// Author: Daniel Britzger
// DESY, 15.01.2015

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO interface for alpos                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "alpos/functions/fastNLOAlpos.h"
#include "alpos/ATheory.h"

using namespace std;



//______________________________________________________________________________


fastNLOAlpos::fastNLOAlpos(string name) : fastNLOReader(name) {
}


//______________________________________________________________________________


fastNLOAlpos::fastNLOAlpos(string name, string LHAPDFFile, int PDFMember) : fastNLOReader(name) {
   // Everything set. Do cross sections calculation.
}


//______________________________________________________________________________


double fastNLOAlpos::EvolveAlphas(double Q) const {
   //debug<<"EvolveAlphas with Q="<<Q<<endl;
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //

   // working example
   // SET(Alpha_s.mur,Q,0);
   // double as = PAR(Alpha_s); // since alpha_s returns a single parameter
   // return as;

   // avoid repeated lookup in map
   return fAs->GetQuick(vector<double>{Q})[0];


   // working examples
   // Mind: Call PAR(Alpha_s) somewhere in advance to guarantee correctly updated values
   return QUICK(Alpha_s, ({Q}) )[0];
   //return QUICK_VAR(Alpha_s,1,Q )[0]; // alternative example
   //return QUICK_VEC(Alpha_s, Q )[0]; // alternative example

}


//______________________________________________________________________________


bool fastNLOAlpos::InitPDF() {
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // LHAPDF interface:
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinitialized the set PDF-set.

   // nothing todo
   PAR(Alpha_s);
   PAR(PDF);
   fPDF = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".PDF"));
   fAs  = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".Alpha_s"));
   return true;
}


//______________________________________________________________________________



vector<double> fastNLOAlpos::GetXFX(double x, double Q) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //

   // Example access (0)
   // avoid repeated lookup in map (i.e. string comparisons)
   return fPDF->GetQuick({x,Q});

   // Example access (1)
   // Mind: Call PAR(PDF) somewhere in advance to guarantee correctly updated values
   return QUICK(PDF, ({x,Q}) );

   // Example access (2)
   //return QUICK_VAR(PDF, 2, x,Q );

   // Example access (3)
   // vector<double> xpmuf = {x,Q};
   // return QUICK_VEC(PDF, xpmuf );

   // Example access (4)
   // SET(PDF.xp,x,0);
   // SET(PDF.muf,Q,0);
   // vector<double> vals = VALUES(PDF);
   // return vals;

}


//______________________________________________________________________________


