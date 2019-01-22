// DB 16.01.2015

#include "alpos/functions/AfastNLOalt.h"
#include <iostream>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLOalt::fRequirements = {"Filename","PDF","Alpha_s","ScaleFacMuR","ScaleFacMuF","Units","iOrd","MuRFuncForm","MuFFuncForm"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOalt::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLOalt::fFunctionName = "fastNLOalt"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLOalt::AfastNLOalt(const std::string& name) : AParmFuncBase<double>(name), fastNLOReader() {
   // fastNLO initializations
   fUnits               = fastNLO::kPublicationUnits;
   fMuRFunc             = fastNLO::kScale1;
   fMuFFunc             = fastNLO::kScale1;
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   fUseHoppet           = false;
}


// __________________________________________________________________________________________ //
AfastNLOalt::~AfastNLOalt() {
}


// ___________________________________________________________________________________________ //
bool AfastNLOalt::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   cout<<" Init AfastNLOalt!"<<endl;
   string filename = PAR_S(Filename);
   // redo fastNLO initialisation
   //   fastNLOBase::SetFilename(filename); // just 'set' it
   ReadTable();
   fastNLOReader::SetFilename(filename); // do some remaining initialization
   // Things which do  not make sense to change later
   SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
   if (GetIsFlexibleScaleTable()) {
      SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
      SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
   }

   CONST(Filename);
   CONST(MuRFuncForm);
   CONST(MuFFuncForm);

   return true;
}


// __________________________________________________________________________________________ //
bool AfastNLOalt::Update() {  //alpos
   if ( CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF) )
      SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
   if ( CHECK(iOrd) ) {
      int iOrd = PAR(iOrd);
      //SetContributionON(fastNLO::kFixedOrder , 0 , true);  // LO is always ON
      if ( iOrd==0 ) {
         SetContributionON(fastNLO::kFixedOrder , 1, false);
         SetContributionON(fastNLO::kFixedOrder , 2, false);
      }
      else if ( iOrd==1 ){
         SetContributionON(fastNLO::kFixedOrder , 1,  true);
         SetContributionON(fastNLO::kFixedOrder , 2, false);
      }
      else if ( iOrd==2 ){
         SetContributionON(fastNLO::kFixedOrder , 1,  true);
         SetContributionON(fastNLO::kFixedOrder , 2,  true);
      }
   }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(PDF);
   UPDATE(Alpha_s);

   // get cross sections
   CalcCrossSection();
   fValue = GetCrossSection();
   fError.resize(fValue.size());

   return true;
}

//______________________________________________________________________________


double AfastNLOalt::EvolveAlphas(double Q) const { //fastNLO
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return QUICK(Alpha_s, ({Q}) )[0];
}


//______________________________________________________________________________


bool AfastNLOalt::InitPDF() { //fastNLO
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // nothing todo
   PAR(Alpha_s);
   PAR(PDF);
   return true;
}


//______________________________________________________________________________



vector<double> AfastNLOalt::GetXFX(double x, double Q) const { //fastNLO
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return QUICK(PDF, ({x,Q}) );

}


//______________________________________________________________________________
