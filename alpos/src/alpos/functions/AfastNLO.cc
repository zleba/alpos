// DB 15.01.2015

#include "alpos/functions/AfastNLO.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include <alpos/functions/ALhapdf6.h>
#include <alpos/AError.h>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLO::fRequirements = {//"Filename",                      // fastNLO table filenmae
                                                          "PDF",                           // a 'PDF' function with Quick-function. Moslty LHAPDF6
                                                          "Alpha_s",                       // a alpha_s(mu_r) function.
                                                          "ScaleFacMuR","ScaleFacMuF",     // Scale factors for ren. and fact. scale
                                                          // "Units",                         // publication or absolute units (to obtain same units as your data table)
                                                          "iOrd",                          // order
                                                          "iThr",                          // Use threshold corrections if available in fastNLO table
                                                        "MuRFuncForm","MuFFuncForm"      // mu_r and mu_f functional form for fastNLO flexible scale tables
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLO::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLO::fFunctionName = "fastNLO"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLO::AfastNLO(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AfastNLO");
}


// __________________________________________________________________________________________ //
AfastNLO::~AfastNLO() {
   //if ( fnlo ) delete fnlo;
   for ( auto i : fnlos ) delete i;
}


// ___________________________________________________________________________________________ //
bool AfastNLO::Init() {
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


   //string filename = PAR_S(Filename);
   //cout<<"fastnlo: Filename: "<<filename<<endl;
   for (auto& fn : filenames) {
      AlposTools::CheckFileExit(fn);
      info["Init"] << "Reading table file: " << fn << endl;
      fastNLOAlpos* f = new fastNLOAlpos(fn);
      fnlos.push_back(f);
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
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
   CONST(iThr);

   // is PDF uncertainty calculation requested? (default: false)
   const string pdfUncParamName = this->GetFunctionName()+".PDFUncertainties";
   if (TheoryHandler::Handler()->FindParameter(pdfUncParamName) != "") {
      fNeedPDFUncertainties = bool(TheoryHandler::Handler()->GetParmD(pdfUncParamName)->GetValue());
      // PDF uncertainties are multiplicative, so set flag
      fHasMultErrors = fNeedPDFUncertainties;
   }
   return true;
}

bool AfastNLO::CalcCrossSections() {

   using namespace AlposTools;

   // get cross sections
   fValue.clear();
   for ( auto fnlo : fnlos ) {
      fnlo->CalcCrossSection();
      fValue += fnlo->GetCrossSection();
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
bool AfastNLO::Update() {

   using namespace AlposTools;

   if ( fValue.empty() || (fValue.size()==1 && fValue[0]==0 )) SetOrder(); // must be initialized in 'update' only

   if ( CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF) )
      for ( auto fnlo : fnlos ) {
         debug["Update"] << "Setting ScaleFacMuR, ScaleFacMuF to (" << PAR(ScaleFacMuR) << ", " << PAR(ScaleFacMuF) << ")..." << std::endl;
         fnlo->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(PDF);
   UPDATE(Alpha_s);

   CalcCrossSections();

   // if user requested calculation of PDF uncertainty and these have not yet been calculated
   if (fNeedPDFUncertainties) {
      if (!fCalculatedPDFUncertainties)
         CalcPDFUncertainties();
   }

   // rescale errors (if any) to new predictions
   if (HasErrors()) {
      fHasMultErrors = true;  // theory errors are always multiplicative
      if (fRescalePDFUncertainties) {
         debug["Update"] << "Rescaling all uncertainties (" << fAllErrors.size() << ")..." << std::endl;
         RescaleErrors();
      }
   }
   else {
      fHasMultErrors = false;
   }

   return true;
}

// __________________________________________________________________________________________ //
void AfastNLO::CalcPDFUncertainties(double cl, bool alternative) {
   //! Calculate PDF uncertainties
   /*!
    *  This calculates the uncertainty of the theoretical
    *  cross section predictions due to the PDFs. The uncertainty
    *  is calculated a a covariance matrix and stored.
    *
    *  The calculation routines are implemented in LHAPDF6 and are
    *  only called from here. The error model is determined by
    *  LHAPDF6 depending on the nature of the PDF used (usually
    *  Hessian for parameterized PDFs).
    *
    *  Errors are stored as matrix-type AError objects
    *  containing the computed covariance matrix.
    *
    *  \param cl Confidence level to assume for uncertainties, in percent.
    *            The default value of 68.2689... corresponds to a one-sigma uncertainty.
    *
    *  \param alternative Whether to use an alternative calculation method (see LHAPDF
    *                     documentation).
    *
    *  \note This method only works for PDFs accessed using LHAPDF6,
    *        since the calculation routine itself in implemented there.
    *
    *  \note The calculation takes a relatively long time, which is why
    *        it should only be done once. Should the cross section predictions
    *        change after calculation, the PDF uncertainty is not recalculated,
    *        but rescaled accordingly (see AError::RescaleError()). The underlying
    *        assumption is that the relative errors do not change significantly.
    *
    *  \see AError::RescaleError()
    */
   debug["CalcPDFUncertainties"]<<"Calculating cross section uncertainty due to PDFs."<<endl;
   debug["CalcPDFUncertainties"]<<"CL is "<<cl<<","<<(alternative ? " not " : " ")<<"using alternative method."<<endl;

   vector<LHAPDF::PDFUncertainty> PDFUnc;
   const string pdfAlposName = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".PDF"))->GetAlposName();
   const string pdfFuncName = TheoryHandler::Handler()->GetFuncD(this->GetAlposName() + std::string(".PDF"))->GetFunctionName();

   //bool LHAType=TheoryHandler::Handler()->GetFuncD(pdfFuncName)->GetFunctionName().find(ALhapdf6::fFunctionName) !=std::string::npos;
   bool LHAType = pdfFuncName.find(ALhapdf6::fFunctionName) != std::string::npos;
   if ( LHAType ) {
      // if using LHAPDF6, do the uncertainty calculation

      // code adapted from fastNLOLHAPDF::GetPDFUncertaintyLHAPDF

      // pointer to function providing PDF
      ALhapdf6 *lha = (ALhapdf6 *) TheoryHandler::Handler()->GetFuncD(pdfAlposName);
      const unsigned int nMem = lha->GetPDFSet()->size();
      const unsigned int nObsBins = fValue.size();

      // Alpos name of PDF interface
      //string lhaName = TheoryHandler::Handler()->GetFuncD(pdfFuncName)->GetAlposName();

      // member ID to return to after calculation
      int iOrigMem = TheoryHandler::Handler()->GetParmD(pdfAlposName+string(".")+string("PDFSet"))->GetValue();
      //int iOrigMem = lha->GetPDF()->memberID();  // get directly from LHAPDF (?)

      // calculate and store CSs for each member
      vector<vector<double>> CSs(nObsBins);
      for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ )
         CSs[iObs].resize(nMem);

      for ( unsigned int iMem = 0 ; iMem<nMem ; iMem++ ) {
         TheoryHandler::Handler()->GetParmD(pdfAlposName+string(".")+string("PDFSet"))->SetValue(iMem,0,false);
         this->CalcCrossSections();
         for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ )
            CSs[iObs][iMem] = fValue[iObs];
      }

      // calculate uncertainty
      for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) {
         LHAPDF::PDFUncertainty err_struct = lha->GetPDFSet()->uncertainty(CSs[iObs],cl,alternative);
         double err_val = err_struct.errsymm;
         PDFUnc.push_back(err_struct);
         fError[iObs] = err_val; // add symmetric error to fError
      }

      // calculate covariance matrix
      TMatrixDSym covMat(nObsBins);
      for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) {
         covMat[iObs][iObs] = pow(PDFUnc[iObs].errsymm,2);
         for ( unsigned int jObs = iObs+1 ; jObs<nObsBins ; jObs++ ) {
            double corr = lha->GetPDFSet()->correlation(CSs[iObs], CSs[jObs]);  // correlation of PDFs between CS bins A and B
            double entry = PDFUnc[iObs].errsymm * PDFUnc[jObs].errsymm * corr;
            covMat[iObs][jObs] = entry;
            covMat[jObs][iObs] = entry; //TODO: Check if necessary for symmetric matrix
         }
      }

      // update fAllErrors
      string errorName = "PDFUncer_" + GetFunctionName();
      AError tmp(errorName, GetAlposName()); // TODO: Check if ErrorSet=GetAlposName() produces expected behavior
      debug["CalcPDFUncertainties"]<<"Updating fAllErrors. Error name is '"<<tmp.GetErrorName()<<"'"<<endl;
      if ( fAllErrors.find(tmp.GetErrorName()) == fAllErrors.end() ) {
         fAllErrors[tmp.GetErrorName()] = tmp;
      }
      AError& pdfError = fAllErrors[tmp.GetErrorName()];

      pdfError.SetIsTheo();
      //                      err     mat     val     rel?   corr?  nature
      pdfError.SetMatrixError(fError, covMat, fValue, false, false, "");

      // set the name of the 'up' and 'down' data table columns for this error to
      // dummy values. This was originally required because AError objects were always
      // associated with an AData object
      pdfError.SetColUpDn(errorName, errorName);  // should have no effect otherwise (?)

      // reset "active" PDF member and re-calculate CS
      TheoryHandler::Handler()->GetParmD(pdfAlposName+string(".")+string("PDFSet"))->SetValue(iOrigMem,0,false);
      this->CalcCrossSections();

      debug["CalcPDFUncertainties"]<<"Covariance matrix calculated."<<endl;
      fCalculatedPDFUncertainties = true;
   }
   else {
      // if not using LHAPDF6, warn and do nothing
      warn["CalcPDFUncertainties"]<<"Only LHAPDF-style PDFs are currently implemented."<<endl;
   }
}

// __________________________________________________________________________________________ //
void AfastNLO::SetOrder() {
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
      // Threshold corrections
//      if ( CHECK(iThr) ) {
         int iThr = PAR(iThr);
         if ( !(iThr<0) ) {
            if ( iThr==0 && iOrd==0 ) {
               if ( !(ithc1<0) )  {
                  fnlo->SetContributionON(fastNLO::kThresholdCorrection, ithc1, true);
               } else {
                  error["SetOrder"] << "LO threshold corrections requested, but not found. Nothing to be done!" << endl;
                  exit(1);
               }
            } else if ( iThr==1 && iOrd==1 ) {
               if ( !(ithc2<0) )  {
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
