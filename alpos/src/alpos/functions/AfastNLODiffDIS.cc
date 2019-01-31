// DB 15.01.2015

#include "alpos/functions/AfastNLODiffDIS.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include <alpos/AError.h>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLODiffDIS::fRequirements = {//"Filename",                      // fastNLO table filenmae
                                                          "DPDF",                           // a 'PDF' function with Quick-function. Moslty LHAPDF6
                                                          "Alpha_s",                       // a alpha_s(mu_r) function.
                                                          "ScaleFacMuR","ScaleFacMuF",     // Scale factors for ren. and fact. scale
                                                          // "Units",                         // publication or absolute units (to obtain same units as your data table)
                                                          "iOrd",                          // order
							  "MuRFuncForm","MuFFuncForm",      // mu_r and mu_f functional form for fastNLO flexible scale tables
							  "xpom_min","xpom_max","nxpom","logxpom" // xpom slicing
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLODiffDIS::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLODiffDIS::fFunctionName = "fastNLODiffDIS"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLODiffDIS::AfastNLODiffDIS(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AfastNLODiffDIS");
}


// __________________________________________________________________________________________ //
AfastNLODiffDIS::~AfastNLODiffDIS() {
   //if ( fnlo ) delete fnlo;
   for ( auto i : fnlos ) delete i;
}


// ___________________________________________________________________________________________ //
bool AfastNLODiffDIS::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Requires in datafile defintions  of
   //!  'Filename' (one entry), 'Filenames' (array) or 'Tables' (table)
   //! and further 'Units' and if applicable: 'MuRFuncForm' and 'MuFFuncForm'
   //!

   using namespace AlposTools;

   double Units = DOUBLE_NS(Units,GetAlposName());

   vector<string> filenames;
   vector<int> firstbins, lastbins;
   if ( EXIST_NS(Filename,GetAlposName() ))
      filenames.push_back(STRING_NS(Filename,GetAlposName()));
   else if ( EXIST_NS(Filenames,GetAlposName() ) )
      filenames = STRING_ARR_NS(Filenames,GetAlposName()); // direct access to array
   else if ( EXIST_NS(Tables,GetAlposName()) ) {
      filenames = STRING_COL_NS(FnloTables,Filenames,GetAlposName());
      firstbins = INT_COL_NS(FnloTables,FirstBin,GetAlposName());
      lastbins  = INT_COL_NS(FnloTables,LastBin,GetAlposName());
   }


   CONST(xpom_min);
   CONST(xpom_max);
   CONST(nxpom);
   CONST(logxpom);

   //string filename = PAR_S(Filename);
   //cout<<"fastnlo: Filename: "<<filename<<endl;
   for (auto& fn : filenames) {
      AlposTools::CheckFileExit(fn);
      info["Init"] << "Reading table file: " << fn << endl;
      //fastNLODiffH12006FitB* f = new fastNLODiffH12006FitB(fn);
      fastNLOAlposDPDF* f = new fastNLOAlposDPDF(fn);
      fnlos.push_back(f);
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
      //fnloreaders[i]->SetXPomLinSlicing( INT(NumberOfXPomSlices), 0. ,  .03 ); //12.  
      if ( PAR(logxpom) == 0 ) 
	 f->SetXPomLinSlicing( PAR(nxpom),PAR(xpom_min),PAR(xpom_max) ); //12.  
      else
	 f->SetXPomLogSlicing( PAR(nxpom),PAR(xpom_min),PAR(xpom_max) ); //12.  
      //int FitID = BOOL(DoFitB) ? 2 : 1;
      //f->SetFitID(2);
      if (f->GetIsFlexibleScaleTable()) {
         // f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         // f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
	 f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
	 f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
      }
   }

   // init binmap if specified
   if (!firstbins.empty()) {
      for (unsigned int ig = 0; ig < filenames.size(); ig++) {
         vector<bool> bmap(fnlos[ig]->GetNObsBin());
         for (unsigned int ib = 0; ib < bmap.size(); ib++) {
            bmap[ib] = (ib >= firstbins[ig] && ib <= lastbins[ig]);
         }
         fBinmap += bmap;
      }
   }

   /*
   //if ( fnlo ) delete fnlo;
   fnlo = new fastNLOAlpos(filenames[0]);
   fnlo->SetAlposName(this->GetAlposName());
   //Things which do not make sense to change later
   //fnlo->SetFilename(PAR_S(Filename)); // filename should not be changed
   fnlo->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
   if (fnlo->GetIsFlexibleScaleTable()) {
      fnlo->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
      fnlo->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
   }
   CONST(MuRFuncForm);
   CONST(MuFFuncForm);
   CONST(Units);
   CONST(Filename);
   */
   CONST(iOrd);

   return true;
}

bool AfastNLODiffDIS::CalcCrossSections() {

   using namespace AlposTools;

   // get cross sections
   fValue.clear();
   for ( auto fnlo : fnlos ) {
      // fnlo->CalcCrossSection();
      // fValue += fnlo->GetCrossSection();
      fValue += fnlo->GetDiffCrossSection();
   }
   // apply binmap if applicable
   if ( !fBinmap.empty() ) {
      int ii=0;
      for ( int ib = 0 ; ib<fValue.size() ; ib++ ){
         if ( fBinmap[ib] ) fValue[ii++]=fValue[ib];
      }
      fValue.resize(ii);
   }

   fError.resize(fValue.size());

   return true;
}

// __________________________________________________________________________________________ //
bool AfastNLODiffDIS::Update() {

   using namespace AlposTools;

   if ( fValue.empty() || (fValue.size()==1 && fValue[0]==0 )) SetOrder(); // must be initialized in 'update' only

   if ( CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF) )
      for ( auto fnlo : fnlos ) {
         debug["Update"] << "Setting ScaleFacMuR, ScaleFacMuF to (" << PAR(ScaleFacMuR) << ", " << PAR(ScaleFacMuF) << ")..." << std::endl;
         fnlo->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(DPDF);
   UPDATE(Alpha_s);

   CalcCrossSections();

   fHasMultErrors = false;

   return true;
}


// __________________________________________________________________________________________ //
void AfastNLODiffDIS::SetOrder() {
   //! Set correct order of fastNLO calculation
   for ( auto fnlo : fnlos ) {
   // if ( CHECK(iOrd) ) {
      //! Check on existence of various pQCD contributions in table (Id = -1 if not existing)
      //! Check on existence of LO (Id = -1 if not existing)
      int ilo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
      if (ilo < 0) {
         error["SetOrder"] << "LO not found, nothing to be done!" << endl;
         exit(1);
      } else {
         info["SetOrder"] << "The LO contribution has Id: " << ilo << endl;
      }
      //! Check on existence of NLO (Id = -1 if not existing)
      int inlo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
      if (inlo < 0) {
         info["SetOrder"] << "No NLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NLO contribution has Id: " << inlo << endl;
      }
      //! Check on existence of NNLO (Id = -1 if not existing)
      int innlo = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToNextToLeading);
      if (innlo < 0) {
         info["SetOrder"] << "No NNLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NNLO contribution has Id: " << innlo << endl;
      }
      //! Switch selected contributions ON, if possible
      //! LO & NLO are ON by default
      //! Fixed-order
      int iOrd = PAR(iOrd);
      if ( iOrd==0 ) {
         if ( !(inlo<0) )  {fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, false);}
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==1 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==2 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, true);
         } else {
            error["SetOrder"] << "NNLO requested, but not found. Ignoring call!" << endl;
            //exit(1);
         }
      }
//      }
//   }
   }
}
