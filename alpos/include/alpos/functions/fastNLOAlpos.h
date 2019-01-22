// Author: Daniel Britzger
// DESY, 15/01/2015

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO interface for alpos                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef FASTNLOALPOS
#define FASTNLOALPOS

#include "fastnlotk/fastNLOReader.h"
#include "alpos/ATheoryHandler.h" // AFuncD*

class fastNLOAlpos : public fastNLOReader {

public:
   fastNLOAlpos(std::string name);
   fastNLOAlpos(std::string name, std::string LHAPDFfile, int PDFSet = 0);
   
   void SetAlposName(const std::string& name) {fAlposName = name;}
   const std::string& GetAlposName() const { return fAlposName;}

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   virtual bool InitPDF();
   std::vector<double> GetXFX(double xp, double muf) const ;

private:
   std::string fAlposName;						//! Name of AlposFunction that contains this instance. An alpos name with getter 'GetName()' is required to obtain parameter and function values
   AFuncD* fPDF = NULL;
   AFuncD* fAs  = NULL;
};

#endif
