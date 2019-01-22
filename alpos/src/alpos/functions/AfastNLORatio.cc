// DB 15.01.2015

#include "alpos/functions/AfastNLORatio.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include <alpos/functions/ALhapdf6.h>
#include <alpos/AError.h>

/*!
 AfastNLORatio

 This function computes the ratio of two fastNLO tables.

 The fastNLO tables are given as the steering parameters
 "FilenameNumerator" and "FilenameDenominator".

 Alternatively, a series of fastNLO tables may be specified for both the
 numerator and the denominator. In this case, the relevant steering parameters
 are lists of filenames and are called "FilenamesNumerator" and "FilenamesDenominator"
 (note the plural forms).

 Note: in case of several fastNLO tables, the number of files provided for the
       numerator and denominator must be the same. This is also required of the number
       of observable bins for each pair of table files.

 */


using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLORatio::fRequirements = {//"Filename",                      // fastNLO table filenmae
      "PDF",                           // a 'PDF' function with Quick-function. Moslty LHAPDF6
      "Alpha_s",                       // a alpha_s(mu_r) function.
      "ScaleFacMuR", "ScaleFacMuF",     // Scale factors for ren. and fact. scale
      // "Units",                         // publication or absolute units (to obtain same units as your data table)
      "iOrd",                          // order
      "iThr",                          // Use threshold corrections if available in fastNLO table
//							  "MuRFuncForm","MuFFuncForm"      // mu_r and mu_f functional form for fastNLO flexible scale tables
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLORatio::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLORatio::fFunctionName = "fastNLORatio"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLORatio::AfastNLORatio(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AfastNLORatio");
}


// __________________________________________________________________________________________ //
AfastNLORatio::~AfastNLORatio() {
   //if ( fnlo ) delete fnlo;
   for (auto i : fnlosNumerator) delete i;
   for (auto i : fnlosDenominator) delete i;
}


// ___________________________________________________________________________________________ //
bool AfastNLORatio::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Requires in datafile defintions  of
   //!  'Filename' (one entry), 'Filenames' (array) or 'Tables' (table)
   //! and further 'Units' and if applicable: 'MuRFuncForm' and 'MuFFuncForm'
   //!

   using namespace AlposTools;

   double Units = DOUBLE_NS(Units, GetAlposName());

   unsigned int nFiles = 0;
   unsigned int nObsBins = 0;
   vector<string> filenamesNumerator;
   vector<string> filenamesDenominator;
   vector<int> firstbins, lastbins;

   // Note: general assumption: if more than one table is given for numerator and denominator,
   // the number of observable bins must match up *for each table pair*

   if (EXIST_NS(Tables, GetAlposName())) {
      // -- fastNLO interpolation tables are given in a single table
      filenamesNumerator = STRING_COL_NS(FnloTables, FilenameNumerator, GetAlposName());
      filenamesDenominator = STRING_COL_NS(FnloTables, FilenameDenominator, GetAlposName());
      firstbins = INT_COL_NS(FnloTables, FirstBin, GetAlposName());
      lastbins = INT_COL_NS(FnloTables, LastBin, GetAlposName());
   }
   else {
      // -- fastNLO interpolation tables are given separately
      if (EXIST_NS(FilenameNumerator, GetAlposName()))
         filenamesNumerator.push_back(STRING_NS(FilenameNumerator, GetAlposName()));
      else if (EXIST_NS(FilenamesNumerator, GetAlposName()))
         filenamesNumerator = STRING_ARR_NS(FilenamesNumerator, GetAlposName());
      else {
         error["Init"] << "Required steering parameter 'Filename(s)Numerator' not given. Exiting..." << endl;
         exit(1);
      }

      if (EXIST_NS(FilenameDenominator, GetAlposName()))
         filenamesDenominator.push_back(STRING_NS(FilenameDenominator, GetAlposName()));
      else if (EXIST_NS(FilenamesDenominator, GetAlposName()))
         filenamesDenominator = STRING_ARR_NS(FilenamesDenominator, GetAlposName());
      else {
         error["Init"] << "Required steering parameter 'Filename(s)Numerator' not given. Exiting..." << endl;
         exit(1);
      }

      if (filenamesNumerator.empty()) {
         error["Init"] << "Couldn't read steering parameter 'Filename(s)Numerator' or no files given. Exiting..." <<
         endl;
         exit(1);
      }
      if (filenamesDenominator.empty()) {
         error["Init"] << "Couldn't read steering parameter 'Filename(s)Denominator' or no files given. Exiting..." <<
         endl;
         exit(1);
      }
      if (filenamesNumerator.size() != filenamesDenominator.size()) {
         error["Init"] << "Number of files in 'FilenamesNumerator' (" << filenamesNumerator.size()
         << ") and 'FilenamesDenominator' (" << filenamesNumerator.size()
         << ") do not match.. Exiting..." << endl;
         exit(1);
      }
      nFiles = filenamesNumerator.size();
   }

   for (auto fn : filenamesNumerator) {
      info["Init"] << "Reading numerator table file: " << fn << endl;
      fastNLOAlpos* f = new fastNLOAlpos(fn);
      fnlosNumerator.push_back(f);
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }
   }

   for (auto fn : filenamesDenominator) {
      info["Init"] << "Reading denominator table file: " << fn << endl;
      fastNLOAlpos* f = new fastNLOAlpos(fn);
      fnlosDenominator.push_back(f);
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }
   }

   unsigned int nBinsNumerator = 0;
   unsigned int nBinsDenominator = 0;

   // check if total bin numbers are compatible
   for (unsigned int iFile = 0; iFile < nFiles; iFile++) {
      nBinsNumerator += fnlosNumerator[iFile]->GetNObsBin();
      nBinsDenominator += fnlosDenominator[iFile]->GetNObsBin();
   }
   if (nBinsNumerator != nBinsDenominator) {
      error["Init"] << "Total number of observable bins in the numerator tables (" << nBinsNumerator
      << ") does not match that of the denominator tables (" << nBinsDenominator
      << ")! Exiting..." << endl;
      exit(1);
   }

   nObsBins = nBinsNumerator;


   // init binmap if specified
   if (!firstbins.empty()) {
      for (unsigned int iFile = 0; iFile < nFiles; iFile++) {
         vector<bool> bmap(fnlosNumerator[iFile]->GetNObsBin());
         for (unsigned int iObsBin = 0; iObsBin < bmap.size(); iObsBin++) {
            bmap[iObsBin] = (iObsBin >= firstbins[iFile] && iObsBin <= lastbins[iFile]);
         }
         fBinmap += bmap;
      }
   }

   CONST(iOrd);
   CONST(iThr);

   SetOrder();

   return true;
}

bool AfastNLORatio::CalcCrossSections() {

   using namespace AlposTools;

   std::vector<fastNLOAlpos*> allfnlos = fnlosNumerator;
   allfnlos += fnlosDenominator;

   // get cross sections
   fValue.clear();
   for (unsigned int iFile = 0; iFile < fnlosNumerator.size(); iFile++) {
      fnlosNumerator[iFile]->CalcCrossSection();
      fnlosDenominator[iFile]->CalcCrossSection();
      std::vector<double> values = fnlosNumerator[iFile]->GetCrossSection();
      values /= fnlosDenominator[iFile]->GetCrossSection();
      fValue += values;
   }

   // apply binmap if applicable
   if (!fBinmap.empty()) {
      int iObsBinFiltered = 0;
      for (int iObsBin = 0; iObsBin < fValue.size(); iObsBin++) {
         if (fBinmap[iObsBin]) fValue[iObsBinFiltered++] = fValue[iObsBin];
      }
      fValue.resize(iObsBinFiltered);
   }

   fError.resize(fValue.size());

   return true;
}

// __________________________________________________________________________________________ //
bool AfastNLORatio::Update() {

   using namespace AlposTools;

   std::vector<fastNLOAlpos*> allfnlos = fnlosNumerator;
   allfnlos += fnlosDenominator;

   if (CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF))
      for (auto fnlo : allfnlos)
         fnlo->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR), PAR(ScaleFacMuF));

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(PDF);
   UPDATE(Alpha_s);

   CalcCrossSections();

   // rescale errors (if any) to new predictions
   if (HasErrors()) {
      fHasMultErrors = true;  // theory errors are always multiplicative
      debug["Update"] << "Rescaling all uncertainties (" << fAllErrors.size() << ")..." << std::endl;
      RescaleErrors();
   }
   else {
      fHasMultErrors = false;
   }

   return true;
}

// __________________________________________________________________________________________ //
void AfastNLORatio::SetOrder() {
   //! Set correct order of fastNLO calculation

   using namespace AlposTools;

   std::vector<fastNLOAlpos*> allfnlos = fnlosNumerator;
   allfnlos += fnlosDenominator;

   for (auto fnlo : allfnlos) {
      // if ( CHECK(iOrd) ) {
      //! Check on existence of various pQCD contributions in table (Id = -1 if not existing)
      //! Check on existence of LO (Id = -1 if not existing)
      int ilo = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
      if (ilo < 0) {
         error["SetOrder"] << "LO not found, nothing to be done!" << endl;
         exit(1);
      } else {
         info["SetOrder"] << "The LO contribution has Id: " << ilo << endl;
      }
      //! Check on existence of NLO (Id = -1 if not existing)
      int inlo = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
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
      //! Check on existence of threshold corrections
      int ithc1 = fnlo->ContrId(fastNLO::kThresholdCorrection, fastNLO::kLeading);
      int ithc2 = fnlo->ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
      if (ithc1 < 0) {
         info["SetOrder"] << "1-loop threshold corrections not found!" << endl;
      } else {
         info["SetOrder"] << "1-loop threshold corrections have Id: " << ithc1 << endl;
      }
      if (ithc2 < 0) {
         info["SetOrder"] << "2-loop threshold corrections not found!" << endl;
      } else {
         info["SetOrder"] << "2-loop threshold corrections have Id: " << ithc2 << endl;
      }

      //! Switch selected contributions ON, if possible
      //! LO & NLO are ON by default
      //! Fixed-order
      int iOrd = PAR(iOrd);
      if (iOrd == 0) {
         if (!(inlo < 0)) { fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, false); }
         if (!(innlo < 0)) { fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false); }
      } else if (iOrd == 1) {
         if (!(inlo < 0)) {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if (!(innlo < 0)) { fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false); }
      } else if (iOrd == 2) {
         if (!(inlo < 0)) {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if (!(innlo < 0)) {
            fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, true);
         } else {
            error["SetOrder"] << "NNLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
      }
      // Threshold corrections
//      if ( CHECK(iThr) ) {
      int iThr = PAR(iThr);
      if (!(iThr < 0)) {
         if (iThr == 0 && iOrd == 0) {
            if (!(ithc1 < 0)) {
               fnlo->SetContributionON(fastNLO::kThresholdCorrection, ithc1, true);
            } else {
               error["SetOrder"] << "LO threshold corrections requested, but not found. Nothing to be done!" << endl;
               exit(1);
            }
         } else if (iThr == 1 && iOrd == 1) {
            if (!(ithc2 < 0)) {
               fnlo->SetContributionON(fastNLO::kThresholdCorrection, ithc2, true);
            } else {
               error["SetOrder"] << "NLO threshold corrections requested, but not found. Nothing to be done!" << endl;
               exit(1);
            }
         } else {
            error["SetOrder"] << "Inconsistent request for threshold corrections. Nothing to be done!" << endl;
            exit(1);
         }
      }
//      }
//   }
   }
}

