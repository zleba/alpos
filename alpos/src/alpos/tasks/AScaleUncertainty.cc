#include "alpos/tasks/AScaleUncertainty.h"
#include "alpos/AFactory.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AFitter.h"

/**
 AScaleUncertainty

 */

using namespace std;

const string AScaleUncertainty::fTaskType = "ScaleUncertainty";

//____________________________________________________________________________________ //
AScaleUncertainty::AScaleUncertainty(const string& aname ) : ATask(aname) {
   //! constructor
   SetClassName("AScaleUncertainty");

   //! Important: create always new result-object here!
   fResult = new AScaleUncertaintyResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AScaleUncertainty::~AScaleUncertainty(){
   //! destructor.
   //! do not delete the AResult object
}


//____________________________________________________________________________________ //
bool AScaleUncertainty::Init(){
   // --- setting up task
   return true;
}


//____________________________________________________________________________________ //
bool AScaleUncertainty::Execute(){
   // --- execute task
   const string c_mur = STRING_NS(ParMuRFac,NS());
   const string c_muf = STRING_NS(ParMuFFac,NS());
   const string fittername = STRING_NS(fittername,NS());
   const string fittertype = STRING_NS(fittertype,NS());
   vector<string> fitpar;
   if (EXISTARRAY_NS(VaryParameters, NS())) {
      // get parameters to vary, if given
      fitpar = STRING_ARR_NS(VaryParameters, NS());
   }
   else {
      // otherwise, take the fitter's fit parameters
      fitpar = STRING_ARR_NS(FitParameters, fittername);
   }

   map<string, double> v0; // <parname,value>
   for ( auto ip : fitpar ) v0[ip] = PAR_ANY(ip);

   // --- do new fit
   vector<double> cmur = { 0.5, 1.0, 1.0, 2.0, 0.5, 2.0, 1.0 };
   vector<double> cmuf = { 1.0, 0.5, 2.0, 1.0, 0.5, 2.0, 1.0 };
   vector<map<string,double> > vV(cmur.size());
   vector<map<string,double> > vE(cmur.size());
   for ( unsigned int ic = 0 ; ic<cmur.size() ; ic++ ) {
      SET_ANY(c_mur,cmur[ic],0);
      SET_ANY(c_muf,cmuf[ic],0);
      const auto& super = TheoryHandler::Handler()->GetSuperPair();
      ATask* fitter = AFactory::TaskFactory(fittername,fittertype);
      fitter->Init();
      fitter->Execute();
      // --- get value & difference
      for ( auto ip : fitpar ) vV[ic][ip] = PAR_ANY(ip);
      for ( auto ip : fitpar ) vE[ic][ip] = PAR_ANY(ip)-v0[ip];
      //delete fitter->GetResult();
      delete fitter;
   }

   // --- printout
   info["Execute"]<<"Result:"<<endl;

   printf("################################################################################\n");
   printf("# Task \"Scale Uncertainty\": Print scale variation results\n");
   printf("################################################################################\n");

   printf("(cmur,cmuf)");
   for ( auto ip : fitpar ) printf("     %-17s          ",ip.c_str());
   printf("\n");

   // --- individual results
   double vVmin, vVmax;
   for ( auto ip : fitpar ) vVmin = v0[ip];
   for ( auto ip : fitpar ) vVmax = v0[ip];
   for ( unsigned int ic = 0 ; ic<cmur.size() ; ic++ ) {
      for ( auto ip : fitpar ) vVmin = min(vV[ic][ip],vVmin);
      for ( auto ip : fitpar ) vVmax = max(vV[ic][ip],vVmax);
      printf("(%4.2f,%4.2f)",cmur[ic],cmuf[ic]);
      for ( auto ip : fitpar )
         printf("\t%7.5f %+7.5f [%+8.2e] (%+5.2f%%)", vV[ic][ip], vE[ic][ip], vE[ic][ip], vE[ic][ip]/v0[ip]*100. );
      printf("\n");
   }

   // --- 2-points (symmetric mur,muf variation)
   printf("( 2-pts + )");
   for ( auto ip : fitpar ) {
      double vmax = max(max(vV[4][ip],vV[5][ip]),vV[6][ip]);
      printf("\t%7.5f %+7.5f [%+8.2e] (%+5.2f%%)",vmax,vmax-v0[ip],vmax-v0[ip],(vmax/v0[ip]-1.)*100.);
   }
   printf("\n");

   printf("( 2-pts - )");
   for ( auto ip : fitpar ) {
      double vmin = min(min(vV[4][ip],vV[5][ip]),vV[6][ip]);
      printf("\t%7.5f %+7.5f [%+8.2e] (%+5.2f%%)",vmin,vmin-v0[ip],vmin-v0[ip],(vmin/v0[ip]-1.)*100.);
   }
   printf("\n");

   // --- 6-points (includes asymmetric mur,muf variation)
   printf("( 6-pts + )");
   for ( auto ip : fitpar ) {
      printf("\t%7.5f %+7.5f [%+8.2e] (%+5.2f%%)",vVmax,vVmax-v0[ip],vVmax-v0[ip],(vVmax/v0[ip]-1.)*100.);
   }
   printf("\n");

   printf("( 6-pts - )");
   for ( auto ip : fitpar ) {
      printf("\t%7.5f %+7.5f [%+8.2e] (%+5.2f%%)",vVmin,vVmin-v0[ip],vVmin-v0[ip],(vVmin/v0[ip]-1.)*100.);
   }
   printf("\n");

   printf("################################################################################\n");

   // reset theory to inital values
   TheoryHandler::Handler()->SetTheorySet(this->GetResult()->GetInitialTheorySet());

   return true;

}


//____________________________________________________________________________________ //
