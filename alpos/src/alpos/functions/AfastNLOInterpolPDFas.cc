// DB 15.01.2015

#include "alpos/functions/AfastNLOInterpolPDFas.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include <iostream>
#include <fastnlotk/read_steer.h>
#include <alpos/functions/ALhapdf6.h>
//#include <boost/math/distributions/chi_squared.hpp>
#include <TSpline.h>
#include "TF1.h"
#include "TList.h"

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLOInterpolPDFas::fRequirements = {//"Filename",                    // fastNLO table name
                        "LHAPDFasSeries",              // LHAPDF set with alpha_s series
                        "AlphasMz",                    // alpha_s(mz) of current cross sections
                        "ScaleFacMuR","ScaleFacMuF",   // scale factors
                        "Units",                       // units consistent wiht input data
                        "iOrd",                        // order of calculation
                        "iThr",                        // order of threshold corrections
                        "MuRFuncForm","MuFFuncForm",   // mur and muf functional form for fastNLO flexible scale tables
                        "FitFunc"                      // ROOT::TF1 fit function for interpolation ('pol2')
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOInterpolPDFas::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLOInterpolPDFas::fFunctionName = "fastNLOInterpolPDFas"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLOInterpolPDFas::AfastNLOInterpolPDFas(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
}


// __________________________________________________________________________________________ //
AfastNLOInterpolPDFas::~AfastNLOInterpolPDFas() {
   // clean up 'Splines' objects (if they exist)
   for (unsigned int iSpline = 0; iSpline < fSplines.size(); iSpline++) {
      if (fSplines[iSpline]) delete fSplines[iSpline];
   }
}


// ___________________________________________________________________________________________ //
bool AfastNLOInterpolPDFas::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.

   // is PDF uncertainty calculation requested? (default: false)
   const string pdfUncParamName = this->GetFunctionName()+".PDFUncertainties";
   if (TheoryHandler::Handler()->FindParameter(pdfUncParamName) != "") {
      fNeedPDFUncertaintiesAs = TheoryHandler::Handler()->GetParmS(pdfUncParamName)->GetValue();
   }
   return true;
}


// ___________________________________________________________________________________________ //
bool AfastNLOInterpolPDFas::ReInit() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   info["ReInit"]<<"Re-initializing AfastNLOInterpolPDFas"<<endl;

   using namespace AlposTools;

   //string filename = PAR_S(Filename);
   string lhapdfstring = PAR_S(LHAPDFasSeries);

   vector<string> filenames;
   if (EXIST_NS(Filename, GetAlposName()))
      filenames.push_back(STRING_NS(Filename, GetAlposName()));
   else if (EXIST_NS(Filenames, GetAlposName())) {
      filenames = STRING_ARR_NS(Filenames, GetAlposName());
   }

   fGraphs.clear();
   // split lhapdfstring into vector<string>
   vector<string> lhapdfnames;
   //boost::split(lhapdfnames, lhapdfstring, boost::is_any_of(", "), boost::token_compress_on);

   // some PDFs are published as separate filesets for each value of alpha_s
   // e.g. for CT10..., which uses CT10..._as_0112 etc.
   // others have PDFs for all values of alpha_s in a single fileset
   // e.g. for MSTW2008nlo_asmzrange
   //     -> need to handle cases differently

   int nTotalAlphasPoints = 0;       // # of total alpha_s points for interpolation (from all PDF filesets)
   bool qSinglePDFFileset;

   // check # of pdf filesets and # of members
   if (lhapdfnames.size() == 1) {
      // one PDF fileset with several alpha_s values
      info["ReInit"]<<"Using single multi-member PDF fileset '"<<lhapdfnames[0]<<"' as an alpha_s(Mz)-series PDF."<<endl;

      qSinglePDFFileset = true;
   }
   else {
      // several PDF filesets, each with one alpha_s value
      info["ReInit"]<<"Using several single-member PDF filesets, one for each alpha_s(Mz)."<<endl;

      qSinglePDFFileset = false;
   }

   // preliminary run through all lhapdfnames
   //vector<fastNLOLHAPDF*> fnlos;
   set<double> processedAsmzValues;  // set to register alpha_s(M_Z) values
   for (string lhapdfname : lhapdfnames) {
      debug["ReInit"]<<"Reading LHAPDF Info for PDF fileset '"<<lhapdfname<<"'..."<<endl;
      LHAPDF::PDFSet* lhapdfset = new LHAPDF::PDFSet(lhapdfname);

      int nMembers = lhapdfset->size();

      // check # of pdf filesets and # of members
      if (qSinglePDFFileset) {
         // single-file, multi-member mode
         if (nMembers == 1) {
            error["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has a single member."
            <<"No interpolation possible. Exiting."<<endl;
            exit(1);
         }
         else {
            nTotalAlphasPoints = nMembers;
            debug["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has "<<nMembers<<" members. "<<endl;

         }
      }
      else {
         // multi-file, single-member mode
         if (nMembers != 1) {
            warn["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has more than one member, but only"
            <<"one PDF is needed per alpha_s(M_Z) point. Using member #0..."<<endl;
            nMembers = 1;
         }
         else {
            debug["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has 1 member."<<endl;
         }
         nTotalAlphasPoints += 1;
      }

      // check for alpha_s(M_Z) duplicates
      for ( int iMember = 0 ; iMember<nMembers ; iMember++ ) {
         LHAPDF::PDF *lhapdfmember = lhapdfset->mkPDF(iMember);
         double asmz = lhapdfmember->alphasQ(91.1876);

         if ( processedAsmzValues.find(asmz) != processedAsmzValues.end() ) {
            // already processed a PDF member with the same alpha_s(M_Z) -> exclude from TGraph
            error["ReInit"]<<"Already encountered this value of alpha_s(M_Z)! "
            <<"Multiple chi2 points for the same value of alpha_s(M_Z) not supported. Exiting."<<endl;
            exit(1);
         }
         processedAsmzValues.insert(asmz);   // store processed alpha_s(M_Z)
      }
   }  // end preliminary loop through lhapdfnames
   debug["ReInit"]<<"Preliminary loop done."<<endl;
   info["ReInit"]<<"Interpolating "<<nTotalAlphasPoints<<" alpha_s(M_Z) points."<<endl;

   //debug["ReInit"]<<"Loading fastNLOLHAPDF for PDF fileset '"<<lhapdfname<<"'..."<<endl;
   debug["ReInit"]<<"Loading fastNLOLHAPDF."<<endl;


   std::vector<fastNLOLHAPDF*> fnlos;

   int nBins = 0;
   for (auto fn : filenames) {
      AlposTools::CheckFileExit(fn);
      info["Init"] << "Reading table file: " << fn << endl;
      fastNLOLHAPDF* f = new fastNLOLHAPDF(fn, lhapdfnames[0], 0);
      f->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }

      SetOrder(f);

      f->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));

      nBins += f->GetNObsBin();
      fnlos.push_back(f);
   }

   fGraphs.resize(nBins);

   // run through all lhapdfnames and
   for ( int iLhapdfname = 0 ; iLhapdfname < lhapdfnames.size() ; iLhapdfname++ ) {
      info["ReInit"]<<"Processing PDF fileset '"<<lhapdfnames[iLhapdfname]<<"'..."<<endl;

      for (auto fnlo : fnlos) {
         fnlo->SetLHAPDFFilename(lhapdfnames[iLhapdfname]);
      }

      info["ReInit"]<<"Running over all " << fnlos[0]->GetNPDFMembers() <<" PDF members."<<endl;
      for ( int iFnloMember = 0 ; iFnloMember<fnlos[0]->GetNPDFMembers() ; iFnloMember++ ) {
         info["ReInit"]<<"Processing PDF set '"<<lhapdfnames[iLhapdfname]<<"', member #"<<iFnloMember<<endl;

         vector<double> cs;

         for (auto fnlo : fnlos) {
            fnlo->SetLHAPDFMember(iFnloMember);
            fnlo->CalcCrossSection();
            cs += fnlo->GetCrossSection();
         }

         double asmz = fnlos[0]->GetAlphasMz();
         info["ReInit"]<<"Member #"<<iFnloMember<<" has alpha_s(M_Z) = "<<asmz<<endl;

         // decide how to fill CS Points into fGraphs
         int iGraphPoint;
         if (qSinglePDFFileset) {
            iGraphPoint = iFnloMember;
         }
         else {
            iGraphPoint = iLhapdfname;
         }

         // Fill points into TGraph
         debug["ReInit"]<<"Adding CS points to TGraph"<<endl;
         for ( int iCSBin = 0 ; iCSBin<cs.size() ; iCSBin++ ) {
            fGraphs[iCSBin].SetPoint(iGraphPoint, asmz, cs[iCSBin]);
            debug["ReInit"]<<"fGraphs["<<iCSBin<<"].SetPoint("<<iGraphPoint<<", "<<asmz<<", cs["<<iCSBin<<"]="<<cs[iCSBin]<<");"<<endl;
         }

         // only process first member, if not in single PDF set mode
         if (!qSinglePDFFileset) {
            break;
         }

      } // end loop over PDFMembers

   } // end loop over fnlos

   // sort the graphs by asmz in ascending order
   for (unsigned int iCSBin = 0; iCSBin < fGraphs.size(); iCSBin++) {
      fGraphs[iCSBin].Sort();
   }

   // do interpolation fits
   string FitFunc = PAR_S(FitFunc);
//   for ( auto& fg : fGraphs ) {

   // should not be necessary, but include for safety
   for (unsigned int iSpline = 0; iSpline < fSplines.size(); iSpline++) {
      if (fSplines[iSpline]) delete fSplines[iSpline];
   }

   fSplines.resize(nBins);
   for (unsigned int iCSBin = 0; iCSBin < nBins; iCSBin++) {
      auto& fg = fGraphs[iCSBin];
      if ((FitFunc == "TSpline3") || (FitFunc == "Spline3") || (FitFunc == "spline3")) {
         //if (fSplines.size() != nBins) fSplines.resize(nBins);
         fSplines[iCSBin] = new TSpline3(fg.GetName(), &fg);
         this->fHasSplines = true;
      }
      else if ((FitFunc == "TSpline5") || (FitFunc == "Spline5") || (FitFunc == "spline5")) {
         //if (fSplines.size() != nBins) fSplines.resize(nBins);
         fSplines[iCSBin] = new TSpline5(fg.GetName(), &fg);
         this->fHasSplines = true;
      }
      else {
         fg.Fit(FitFunc.c_str(), "Q");
         this->fHasSplines = false;
      }
   }

   // cleanup
   for (auto fnlo : fnlos) {
      delete fnlo;
   }
   
   fValue.resize(nBins);
   fError.resize(nBins);
   
   return true;
}


// __________________________________________________________________________________________ //
bool AfastNLOInterpolPDFas::Update() {

   // Other parameters will typically not change during the program run
   if ( CHECK(ScaleFacMuR) ||
        CHECK(ScaleFacMuF) ||
        CHECK(iOrd) ||
        CHECK(MuRFuncForm) ||
        CHECK(MuFFuncForm) )
      ReInit();

   double asMz = PAR(AlphasMz);
   string FitFunc = PAR_S(FitFunc);
   for ( unsigned int ib = 0 ; ib<fGraphs.size() ; ib++ ) {
      if (this->fHasSplines) {
         fValue[ib] = fSplines[ib]->Eval(asMz);
      }
      else {
         fValue[ib] = fGraphs[ib].GetFunction(FitFunc.c_str())->Eval(asMz);
      }
   }

   // if user requested calculation of PDF uncertainty and these have not yet been calculated
   if (!fCalculatedPDFUncertainties) {
      if ((fNeedPDFUncertaintiesAs == "Eigenvectors") || (fNeedPDFUncertaintiesAs == "eigenvectors"))
         CalcHessianPDFUncertaintiesAsEigenvectors();
      else if ((fNeedPDFUncertaintiesAs == "Covariance") || (fNeedPDFUncertaintiesAs == "covariance"))
         CalcPDFUncertaintiesAsCovariance();
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

   //fnlo->PrintCrossSections();
   //exit(1);

   return true;
}

// __________________________________________________________________________________________ //
void AfastNLOInterpolPDFas::SetOrder(fastNLOLHAPDF* fnlo) {

   //! Set correct order of fastNLO calculation
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
         error["SetOrder"] << "NNLO requested, but not found. Nothing to be done!" << endl;
         exit(1);
      }
   }

  // Threshold corrections
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
}

// __________________________________________________________________________________________ //
bool AfastNLOInterpolPDFas::PrintValues() {
  /*
   * Print interpolated cross sections at current value of alpha_s(M_Z)
   */

   double asMz = PAR(AlphasMz);
   string FitFunc = PAR_S(FitFunc);
   info["PrintValues"]<<"alpha_s(M_Z) is "<<asMz<<endl;
   for ( unsigned int ib = 0 ; ib<fGraphs.size() ; ib++ ) {
      info["PrintValues"]<<"fValue[iBin="<<ib<<"] = "<<fValue[ib]<<endl;
   }

   return true;
}

// __________________________________________________________________________________________ //
void AfastNLOInterpolPDFas::CalcPDFUncertaintiesAsCovariance(double cl, bool alternative) {
   //! Calculate PDF uncertainties as a single covariance matrix
   /*!
    *  This calculates the uncertainty of the theoretical
    *  cross section predictions due to the PDFs. The uncertainty
    *  is calculated as a covariance matrix and stored.
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
   debug["CalcPDFUncertaintiesAsCovariance"]<<"Calculating cross section uncertainty due to PDFs."<<endl;
   debug["CalcPDFUncertaintiesAsCovariance"]<<"CL is "<<cl<<","<<(alternative ? " not " : " ")<<"using alternative method."<<endl;

   using namespace AlposTools;

   vector<LHAPDF::PDFUncertainty> PDFUnc;
   // TODO: make sure this works for subsets...
   const string lhapdfEigenvectorsName = TheoryHandler::Handler()->GetParmS(this->GetFunctionName() + std::string(".PDFEigenvectorSet"))->GetValue();

   vector<string> filenames;
   if (EXIST_NS(Filename, GetAlposName()))
      filenames.push_back(STRING_NS(Filename, GetAlposName()));
   else if (EXIST_NS(Filenames, GetAlposName())) {
      filenames = STRING_ARR_NS(Filenames, GetAlposName());
   }

   // Load and initialize PDF sets for uncertainty calculation

   std::vector<fastNLOLHAPDF*> fnlos;

   int nBins = 0;
   for (auto fn : filenames) {
      info["Init"] << "Reading table file: " << fn << endl;
      fastNLOLHAPDF* f = new fastNLOLHAPDF(fn, lhapdfEigenvectorsName, 0);
      f->CalcCrossSection();
      f->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }
      SetOrder(f);
      f->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      nBins += f->GetNObsBin();
      fnlos.push_back(f);
   }

   // code adapted from fastNLOLHAPDF::GetPDFUncertaintyLHAPDF

   // pointer to function providing PDF

   const unsigned int nObsBins = fValue.size();

   // FIXME: no getter for PDFSet in fastNLOLHAPDF. workaround: initialize own PDFSet
   LHAPDF::PDFSet pdfSet(lhapdfEigenvectorsName);
   double setCL = pdfSet.errorConfLevel() / 100.;
   const unsigned int nMem = pdfSet.size();

   // calculate and store CSs for each member
   vector<vector<double>> CSs(nObsBins, vector<double>(nMem));

   // Calc cross sections
   for ( unsigned int iMem = 0 ; iMem<nMem ; iMem++ ) {
      std::vector<double> CS_mem;
      for (auto fnlo : fnlos) {
         fnlo->SetLHAPDFMember(iMem);
         fnlo->CalcCrossSection();
         CS_mem += fnlo->GetCrossSection();
      }
      for (unsigned int iObs = 0 ; iObs < nObsBins ; iObs++ ) {
         CSs[iObs][iMem] = CS_mem[iObs];
      }
   }

   // calculate uncertainty
   for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) {
      LHAPDF::PDFUncertainty err_struct = pdfSet.uncertainty(CSs[iObs],cl,alternative);
      double err_val = err_struct.errsymm;
      PDFUnc.push_back(err_struct);
      fError[iObs] = err_val; // add symmetric error to fError
   }

   // calculate covariance matrix
   TMatrixDSym covMat(nObsBins);
   for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) {
      covMat[iObs][iObs] = pow(PDFUnc[iObs].errsymm,2);
      for ( unsigned int jObs = iObs+1 ; jObs<nObsBins ; jObs++ ) {
         double corr = pdfSet.correlation(CSs[iObs], CSs[jObs]);  // correlation of PDFs between CS bins A and B
         double entry = PDFUnc[iObs].errsymm * PDFUnc[jObs].errsymm * corr;
         covMat[iObs][jObs] = entry;
         covMat[jObs][iObs] = entry;
      }
   }

   // update fAllErrors
   string errorName = "PDFUncer_" + GetFunctionName();
   AError tmp(errorName, GetAlposName()); // TODO: Check if ErrorSet=GetAlposName() produces expected behavior
   debug["CalcPDFUncertaintiesAsCovariance"]<<"Updating fAllErrors. Error name is '"<<tmp.GetErrorName()<<"'"<<endl;
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

   debug["CalcPDFUncertaintiesAsCovariance"]<<"Covariance matrix calculated."<<endl;
   fCalculatedPDFUncertainties = true;
}

// __________________________________________________________________________________________ //
void AfastNLOInterpolPDFas::CalcHessianPDFUncertaintiesAsEigenvectors(double cl) {
   //! Calculate 'hessian'-type PDF uncertainties as one error per eigenvector
   /*!
    *  This calculates the uncertainty of the theoretical
    *  cross section predictions due to the PDFs. The uncertainty
    *  is calculated as a number of separate, fully correlated,
    *  errors.
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
    */
   debug["CalcHessianPDFUncertaintiesAsEigenvectors"]<<"Calculating cross section uncertainty due to PDFs."<<endl;
   debug["CalcHessianPDFUncertaintiesAsEigenvectors"]<<"CL is: "<<cl<<"%"<<endl;

   using namespace AlposTools;

   vector<LHAPDF::PDFUncertainty> PDFUnc;
   // TODO: make sure this works for subsets...
   const string lhapdfEigenvectorsName = TheoryHandler::Handler()->GetParmS(this->GetFunctionName() + std::string(".PDFEigenvectorSet"))->GetValue();

   vector<string> filenames;
   if (EXIST_NS(Filename, GetAlposName()))
      filenames.push_back(STRING_NS(Filename, GetAlposName()));
   else if (EXIST_NS(Filenames, GetAlposName())) {
      filenames = STRING_ARR_NS(Filenames, GetAlposName());
   }

   // Load and initialize PDF sets for uncertainty calculation

   std::vector<fastNLOLHAPDF*> fnlos;

   int nBins = 0;
   for (auto fn : filenames) {
      AlposTools::CheckFileExit(fn);
      info["Init"] << "Reading table file: " << fn << endl;
      fastNLOLHAPDF* f = new fastNLOLHAPDF(fn, lhapdfEigenvectorsName, 0);
      f->CalcCrossSection();
      f->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }
      SetOrder(f);
      f->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      nBins += f->GetNObsBin();
      fnlos.push_back(f);
   }

   double errScale = 1.0;

   const unsigned int nObsBins = fValue.size();

   // FIXME: no getter for PDFSet in fastNLOLHAPDF. workaround: initialize own PDFSet
   LHAPDF::PDFSet pdfSet(lhapdfEigenvectorsName);
   double setCL = pdfSet.errorConfLevel() / 100.;
   const unsigned int nMem = pdfSet.size();
   string errType = pdfSet.errorType();

   if (errType != "hessian") {
      error["CalcHessianPDFUncertaintiesAsEigenvectors"] << "Requested 'hessian' errors, but PDF set '"
                                                         << lhapdfEigenvectorsName << "' doesn't seem to support "
                                                         << "these. Not calculating PDF uncertainties..." << endl;
      return;
   }

   // calculate and store CSs for each member
   vector<vector<double>> CSs(nObsBins, vector<double>(nMem));

   // Calc cross sections
   for ( unsigned int iMem = 0 ; iMem<nMem ; iMem++ ) {
      std::vector<double> CS_mem;
      for (auto fnlo : fnlos) {
         fnlo->SetLHAPDFMember(iMem);
         fnlo->CalcCrossSection();
         CS_mem += fnlo->GetCrossSection();
      }
      for (unsigned int iObs = 0 ; iObs < nObsBins ; iObs++ ) {
         CSs[iObs][iMem] = CS_mem[iObs];
      }
   }

   // central cross section values
   std::vector<double> CSs_central(nObsBins);
   for (int iObs = 0; iObs < nObsBins; iObs++) {
      CSs_central[iObs] = CSs[iObs][0];
   }

   // --- calculate uncertainties from per-member observables

   bool success = true;

   // handle error scaling, if necessary
   if (setCL * 100 != cl) {
      //boost::math::chi_squared chiSquareDistribution(1);
      double qsetCL = TMath::ChisquareQuantile(setCL,1);// RADEK change   boost::math::quantile(chiSquareDistribution, setCL);
      double qreqCL = TMath::ChisquareQuantile(cl/100,1); // RADEK boost::math::quantile(chiSquareDistribution, cl/100.);
      errScale = sqrt(qreqCL / qsetCL);
      debug["CalcSymHessianAndAddAsEigenvectors"] << "PDFSet CL is " << setCL * 100 << " %, but " << cl
      << " % was requested. " << "Errors will be scaled by a factor of "
      << errScale << "." << std::endl;
   }

   // go through all eigenvectors (hence all pdf members) and calculate observables
   unsigned int nEV = (nMem - 1) / 2;
   string errorPrefix = "PDFUncer_" + GetFunctionName();
   for (int iEV = 1; iEV <= nEV; iEV++) {
      std::vector<double> errUp(nObsBins);
      std::vector<double> errDn(nObsBins);
      for (int iObs = 0; iObs < nObsBins; iObs++) {
         errUp[iObs] = (CSs[iObs][2 * iEV - 1] - CSs_central[iObs]) * errScale;
         errDn[iObs] = (CSs[iObs][2 * iEV] - CSs_central[iObs]) * errScale;
      }
      std::vector<double> errSymm(nObsBins);
      for (int iObs = 0; iObs < nObsBins; iObs++) {
         errSymm[iObs] = (CSs[iObs][iEV] - CSs_central[iObs]) * errScale;
      }
      // determine error name
      std::string errName = errorPrefix + std::string("_EV_") + std::to_string(iEV);
      // add error and log success
      // Note: use linear averaging, since this is used in LHAPDF
      success &= this->AddAsymError(errName, errUp, errDn, CSs_central, false, 1., "",
                                    AError::kLinear);

      // set error flags: theoretical, multiplicative
      std::map<std::string, AError>& errs = const_cast<std::map<std::string, AError>&>(this->GetAllErrors());
      errs[errName].SetIsTheo(true);
      errs[errName].SetIsMult(true);
   }

      debug["CalcHessianPDFUncertaintiesAsEigenvectors"]<<"Eigenvector errors calculated."<<endl;
      fCalculatedPDFUncertainties = true;
}

// __________________________________________________________________________________________ //
void AfastNLOInterpolPDFas::PrintInterpolatedPoints(unsigned int nPoints) {
   //! Sample interpolation function at `nPoints` points and print cross section values
   /*!
    * This method is provided for use in tasks or for debugging/plotting purposes.
    * It samples the chosen interpolation function at approximately `nPoints`
    * inside an interval covering all alpha_s(M_Z) values used for interpolation.
    *
    * The sampling interval is chosen to be slightly larger than the actual
    * alpha_s(M_Z) range, in order to show extrapolation behavior as well. This means that
    * the outermost points of the sampling interval lie outside the actual alpha_s(M_Z)
    * range.
    *
    * Note: owing to adaptive sampling, the interpolation points are not evenly spaced
    * over the entire range.
    *
    * \param  nPoints  Approximate number of total points to use for sampling. Owing to
    *                  the adaptive sampling technique used, the actual number of points
    *                  used will vary depending on the number of support points.
    */

   cout << endl << endl;
   printf("#%12s ", "Bin:");
   for (unsigned int ib = 0; ib < this->N(); ib++) {
      printf(" %12d ", ib);
   }

   string FitFunc = PAR_S(FitFunc);

   std::vector<double> asmzSupportPoints;

   for (unsigned int iGraphPoint = 0; iGraphPoint < fGraphs[0].GetN(); iGraphPoint++) {
      double asmz;
      double cs;
      fGraphs[0].GetPoint(iGraphPoint, asmz, cs);
      asmzSupportPoints.push_back(asmz);
   }

   // set 'plot' data range to be a little bigger than the actual data range
   double scale_factor = 1.2;  // 20% increase
   double asmzLo = (asmzSupportPoints.back()*(1+scale_factor) + asmzSupportPoints.front()*(1-scale_factor))/2.;
   double asmzHi = (asmzSupportPoints.back()*(1-scale_factor) + asmzSupportPoints.front()*(1+scale_factor))/2.;

   asmzSupportPoints.insert(asmzSupportPoints.begin(), asmzLo);
   asmzSupportPoints.push_back(asmzHi);


   int nPlotPointsInSegment = ceil((double) nPoints / (double)(asmzSupportPoints.size() - 1));
   cout << endl << "# using " << nPoints << " in " << asmzSupportPoints.size() - 1 << " per segment for a total of " << nPlotPointsInSegment << " per segment.";
   for (unsigned int iSegment = 0; iSegment < asmzSupportPoints.size() - 1; iSegment++) {
      double asmzFrom = asmzSupportPoints[iSegment];
      double asmzTo = asmzSupportPoints[iSegment+1];
      cout << endl << "# new segment (" << iSegment << ") from asmz = " << asmzFrom << " to " << asmzTo;

      for (unsigned int iPointInSegment = 0; iPointInSegment < nPlotPointsInSegment; iPointInSegment++) {
         double asmz = asmzFrom + iPointInSegment * (asmzTo - asmzFrom) / nPlotPointsInSegment;
         cout << endl;
         printf(" %12f ", asmz);
         for (unsigned int ib = 0; ib < this->N(); ib++) {
            if (this->fHasSplines) {
               printf(" %12f ", fSplines[ib]->Eval(asmz));
            }
            else {
               printf(" %12f ", fGraphs[ib].GetFunction(FitFunc.c_str())->Eval(asmz));
            }
         }
      }
   }
   cout << endl;
}

// __________________________________________________________________________________________ //
void AfastNLOInterpolPDFas::PrintSupportPoints() {
   //! Print cross sections for each of the requested alpha_s values/PDFs
   /*!
    * This method is provided for use in tasks or for debugging/plotting purposes.
    * It prints the cross section theory predictions obtained by using each PDF given
    * in ``LHAPDFasSeries``, acompanied by the respective value of alpha_s(M_Z).
    *
    * These are the 'support points' used for the cross section interpolation.
    *
    *    \sa  AfastNLOInterpolPDFas::PrintInterpolatedPoints()
    */

   cout << endl << endl;
   printf("#%12s ", "Bin:");
   for (unsigned int ib = 0; ib < this->N(); ib++) {
      printf(" %12d ", ib);
   }

   string FitFunc = PAR_S(FitFunc);

   for (unsigned int ip = 0; ip < fGraphs[0].GetN(); ip++) {
      double asmz0;
      double cs0;
      cout << endl;
      fGraphs[0].GetPoint(ip, asmz0, cs0);
      printf(" %12f ", asmz0);
      for (unsigned int ib = 0; ib < this->N(); ib++) {
         double asmz;
         double cs;
         fGraphs[ib].GetPoint(ip, asmz, cs);
         if (AlposTools::IdenticalToNSignificantDigits(asmz0, asmz, 2))
            printf(" %12f ", cs);
         else
            printf(" (%10f) ", cs);
      }
   }
   cout << endl;
}
