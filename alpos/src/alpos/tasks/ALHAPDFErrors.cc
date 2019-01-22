//
// Created by DS on 4/18/16.
//

#include <alpos/functions/ALhapdf6.h>
#include "alpos/tasks/ALHAPDFErrors.h"

#include <boost/math/distributions/chi_squared.hpp>
#include <alpos/ASubsetFunction.h>

/*!
 ALHAPDFErrors

 */

const std::string ALHAPDFErrors::fTaskType = "LHAPDFErrors";

//____________________________________________________________________________________ //
//ALHAPDFErrors::ALHAPDFErrors(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
ALHAPDFErrors::ALHAPDFErrors(const std::string& aname) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName(ALHAPDFErrors::fTaskType);
   //! Important: create always new result-object here!
   fResult = new ALHAPDFErrorsResult(aname, GetTaskType());
}


//____________________________________________________________________________________ //
ALHAPDFErrors::~ALHAPDFErrors() {
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool ALHAPDFErrors::Init() {
   info["Init"] << "Hello. ALHAPDFErrors::Init()." << std::endl;

   fErrorPrefix = EXIST_NS(ErrorPrefix, NS()) ? STRING_NS(ErrorPrefix, NS()) : "LHAPDFUncer";
   fLhapdfFunctionName = STRING_NS(LHAPDFfunction, NS());
   fCalculateAs = EXIST_NS(CalculateAs, NS()) ? STRING_NS(CalculateAs, NS()) : "";

   if ((fCalculateAs != "Eigenvectors") && (fCalculateAs != "eigenvectors") &&
       (fCalculateAs != "Covariance") && (fCalculateAs != "covariance") &&
       (fCalculateAs == "")) {
      error["Execute"] << "Unknown PDF error calculation '" << fCalculateAs << "' requested. Must be one of "
      << "'Eigenvectors', 'Covariance' or omitted. Treating as omitted..." << std::endl;
      fCalculateAs = "";
   }

   fOverridePDFSet = STRING_NS(OverridePDFSet, NS());

   if ((fOverridePDFSet == "")) {
      error["Init"] << "Required steering parameter 'OverridePDFSet' is not set. Please specify a PDF set for PDF "
                    << "error calculation, or specify 'no' to use the PDF set already loaded in LHAPDF." << std::endl;
      exit(1);
   }
   else if ((fOverridePDFSet == "No") || (fOverridePDFSet == "no")) {
      info["Init"] << "Using PDF set loaded in Alpos function '" << fLhapdfFunctionName << "' for error calculation." << std::endl;
   }
   else {
      info["Init"] << "Attempting to use PDF set '" << fOverridePDFSet << "' for error calculation." << std::endl;
   }

   return true;
}


//____________________________________________________________________________________ //
bool ALHAPDFErrors::Execute() {
   //! Calculate PDF uncertainties
   /*!
    *  This task calculates the uncertainty of theoretical
    *  cross section predictions due to PDFs. The resulting uncertainty
    *  is stored by this task in each theory object of the 'super'
    *  data-theory pair and is stored in each theory.
    *
    *  The uncertainties calculated depend on the error type supported
    *  by the current PDF set loaded in the 'ALhapdf6' function.
    *
    *  For 'replica'-type uncertainties, only one global covariance matrix
    *  is created. The calculation of this global error matrix is done
    *  by calling the relevant LHAPDF6 routines.
    *
    *  For 'hessian'-type uncertainties, one error source is generated
    *  per eigenvector. These are asymmetric for 'hessian' errors and
    *  symmetric for 'symmhessian' error sources. These are individually
    *  calculated for each eigenvector using a reimplementation of the
    *  same algorithm that LHAPDF6 uses internally.
    *
    *  \note The error matrix resulting from 'adding up' these individual
    *        error sources is equal to the global covariance matrix provided
    *        by LHAPDF6 for 'hessian' errors. Splitting the PDF uncertainties
    *        into individual eigenvector contributions allows for more flexible
    *        manipulation of uncertainties, e.g. during fitting.
    *
    */


   if ( Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
      info["Execute"] << "Skipping calculation of PDF errors, since 'IgnoreTheoryErrors' is set." << std::endl;
      return true;
   }
   else
      info["Execute"] << "Calculating PDF uncertainties..." << std::endl;

   bool success = false;  // return success flag

   // LHAPDF instance information
   ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD(fLhapdfFunctionName);

   std::string initialPdfSet = PAR_ANY_S(fLhapdfFunctionName + std::string(".LHAPDFFile"));

   // PDF set meta-information (e.g. Confidence Level, error type)
   if ((fOverridePDFSet != "No") && (fOverridePDFSet != "no")) {
      // override PDF set with the specified one
      SET_ANY_S(fLhapdfFunctionName + std::string(".LHAPDFFile"), fOverridePDFSet, "");
      lhapdfInstance->Update();
   }

   std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
   std::string errorType = lhapdfInstance->GetPDFSet()->errorType();

   // check error type
   if ((errorType != "replicas") && (errorType != "hessian") && (errorType != "symmhessian")) {
      error["Execute"] << "Error type '" << errorType << "' of PDF Set '" << pdfSetName <<
      "' is unknown or not supported." << std::endl;
      exit(1);
   }

   // get theory array from all data-theory pairs
   const auto& dtps = TheoryHandler::Handler()->GetDataTheoryPairs();
   std::vector<AFuncD*> th_funcs;
   for (auto& dtp : dtps) {
      // check if function is actually a subset
      if (dtp.second.second->GetFunctionName().find(ASubsetFunction::fFunctionName) != std::string::npos) {
         std::string theoAlposName = dtp.second.second->GetRequirements()[0];
         AFuncD* funcToBeAdded = TheoryHandler::Handler()->GetFuncD(theoAlposName);
         if (std::find(th_funcs.begin(), th_funcs.end(), funcToBeAdded) == th_funcs.end()) {
            // if theory function not yet added -> add and inform
            debug["Init"] << "Function '" << dtp.second.second->GetAlposName() << "' is a subset of '" << funcToBeAdded->GetAlposName()
                          << "'. Adding to uncertainty calculation..." << std::endl;
            th_funcs.push_back(funcToBeAdded);
         }
         else {
            // else, skip and inform
            debug["Init"] << "Function '" << theoAlposName << "' is a subset of '" << funcToBeAdded->GetAlposName()
                          << "', but this function is already considered. Skipping..." << std::endl;
         }
      }
      else {
         debug["Init"] << "Function '" << dtp.second.second->GetAlposName() << "' is a not subset. "
                       << "Adding to uncertainty calculation..." << std::endl;
         th_funcs.push_back(dtp.second.second);
      }
   }
   if ((fCalculateAs == "Covariance") || (fCalculateAs == "covariance") || (fCalculateAs == "")) {
      // the super-theory must know about the bin-to-bin correlations
      th_funcs.push_back(TheoryHandler::Handler()->GetFuncD("SuperTheory"));  // requires non-const ptr
      //th_funcs.push_back(TheoryHandler::Handler()->GetSuperPair().second); // this is a const ptr
   }


   // 'hessian', 'symmhessian' or 'replicas'
   if (errorType == "replicas") {
      // no eigenvectors -> only full covariance matrix possible
      info["Execute"] << "Error type is: " << errorType << std::endl;

      // 'replicas' error cannot be calculated as 'eigenvectors'
      if ((fCalculateAs == "Eigenvector") || (fCalculateAs == "eigenvector")) {
         error["Execute"] << "PDF error type '" << errorType << "' incompatible with requested error calculation '"
         << fCalculateAs << "': Must be 'Covariance' or omitted." << std::endl;
         exit(1);
      }
      else if ((fCalculateAs == "Covariance") || (fCalculateAs == "covariance") || (fCalculateAs == "")) {
         // default for 'replicas': one error source with full covariance matrix
         success = this->CalcAndAddAsCovariance(th_funcs);
      }
   }
   else if (errorType == "hessian") {
      // has eigenvectors -> one asymmetric error source per eigenvector
      info["Execute"] << "Error type is: " << errorType << std::endl;

      // check requested calculation
      if ((fCalculateAs == "Eigenvectors") || (fCalculateAs == "eigenvectors") || (fCalculateAs == "")) {
         // one error source per PDF eigenvector
         success = this->CalcHessianAndAddAsEigenvectors(th_funcs);
      }
      else if ((fCalculateAs == "Covariance") || (fCalculateAs == "covariance")) {
         // alternatively, calculate one error source with full covariance matrix
         success = this->CalcAndAddAsCovariance(th_funcs);
      }
   }
   else if (errorType == "symmhessian") {
      // has eigenvectors (symmetric) -> one symmetric error source per eigenvector
      info["Execute"] << "Error type is: " << errorType << std::endl;

      // check requested calculation
      if ((fCalculateAs == "Eigenvectors") || (fCalculateAs == "eigenvectors") || (fCalculateAs == "")) {
         // one error source per PDF eigenvector
         success = this->CalcSymHessianAndAddAsEigenvectors(th_funcs);
      }
      else if ((fCalculateAs == "Covariance") || (fCalculateAs == "covariance")) {
         // alternatively, calculate one error source with full covariance matrix
         success = this->CalcAndAddAsCovariance(th_funcs);
      }
   }

   // restore original PDF set
   SET_ANY_S(fLhapdfFunctionName + std::string(".LHAPDFFile"), initialPdfSet, "");

   return success;
}

//______________________________________________________________________________________________________//
bool ALHAPDFErrors::CalcAndAddAsCovariance(std::vector<AFuncD*> th_funcs) {
   //! Get full covariance matrix from LHAPDF6 and add as Alpos error to theory object.
   /*!
    *  LHAPDF6 provides routines for calculating the full PDF covariance matrix.
    *
    *  This method relies on those routines to construct this matrix
    *  and add it to the theory function as a single Alpos error.
    *
    *  \param  th_func  Pointer to the theory function. This must depend
    *                   on the LHAPDF6 function instance referred to in
    *                   the task parameter 'LHAPDFfunction' in the steering
    *                   file.
    *
    *  \note  For PDF sets using an 'eigenvector' approach to providing their
    *         uncertainties, one of the corresponding methods
    *         `CalcHessianAndAddAsEigenvectors` and `CalcHessianAndAddAsEigenvectors`
    *         should be used.
    */

   bool success = true;

   // LHAPDF instance information
   ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD(fLhapdfFunctionName);
   int initialPdfMember = PAR_ANY(fLhapdfFunctionName + std::string(".PDFSet"));

   // PDF set meta-information (e.g. Confidence Level, error type)
   std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
   unsigned int nMembers = lhapdfInstance->GetPDFSet()->size();
   double setCL = lhapdfInstance->GetPDFSet()->errorConfLevel() / 100.;
   double reqCL = 0.6826894 ;//boost::math::erf(1 / sqrt(2));  // required confidence level: 68.2689... %
   double errScale = 1.0;  // scaling due to set CL != required CL

   // cross sections for each member
   std::vector<std::vector<std::vector<double>>> observableArrays = this->CalcObservableArrays(th_funcs);

   // loop over theory functions
   for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) {
      AFuncD* th_func = th_funcs[iTheoryFunction];
      std::vector<std::vector<double>>& observableArray = observableArrays[iTheoryFunction];
      const unsigned int nObsBins = th_func->N();

      // central cross section values
      std::vector<double> observablesCentral(nObsBins);
      for (int iObs = 0; iObs < nObsBins; iObs++) {
         observablesCentral[iObs] = observableArray[iObs][0];
      }

      // --- calculate uncertainties from per-member observables

      // handle error scaling, if necessary
      if (lhapdfInstance->GetPDFSet()->errorConfLevel() < 0) {
         debug["CalcAndAddAsCovariance"] << "PDF set CL is negative (this is the case for 'replicas'-type PDF sets). "
                                         << "Uncertainties will be estimated to correspond to the requested "
                                         << 100*reqCL << " % CL." << std::endl;
      }
      else if (setCL != reqCL) {
         boost::math::chi_squared chiSquareDistribution(1);
         double qsetCL = boost::math::quantile(chiSquareDistribution, setCL);
         double qreqCL = boost::math::quantile(chiSquareDistribution, reqCL);
         errScale = sqrt(qreqCL / qsetCL);
         debug["CalcAndAddAsCovariance"] << "PDFSet CL is " << setCL * 100 << " %, but " << reqCL * 100
                                         << " % was requested. " << "Errors will be scaled by a factor of "
                                         << errScale << "." << std::endl;
      }

      // calculate uncertainties (routines from LHAPDF6)
      std::vector<double> errVals(nObsBins);
      for (unsigned int iObs = 0; iObs < nObsBins; iObs++) {
         // 'alternative=false' : replicas central value is *mean*, not median
         LHAPDF::PDFUncertainty err_struct = lhapdfInstance->GetPDFSet()->uncertainty(observableArray[iObs], 100*reqCL,
                                                                                      false);
         // TODO: use Alpos averaging, or get symmetrized from LHAPDF?
         errVals[iObs] = err_struct.errsymm;
         //errVals[iObs] = insert_Alpos_averaging_function_here(err_struct.errplus, err_struct.errminus);
      }

      // calculate covariance matrix, errors
      TMatrixDSym covMat(nObsBins);
      for (unsigned int iObs = 0; iObs < nObsBins; iObs++) {
         covMat[iObs][iObs] = pow(errVals[iObs], 2);
         for (unsigned int jObs = iObs + 1; jObs < nObsBins; jObs++) {
            double corr = lhapdfInstance->GetPDFSet()->correlation(observableArray[iObs],
                                                                   observableArray[jObs]);  // correlation of PDFs between obs bins
            double entry = errVals[iObs] * errVals[jObs] * corr;
            covMat[iObs][jObs] = entry;
            covMat[jObs][iObs] = entry;
         }
      }

      // determine error name
      std::string errName = fErrorPrefix + std::string("_Cov");

      // add error source
      success &= th_func->AddMatrixError(errName, errVals, covMat, observablesCentral, false, false, "");

      // set error flags: theoretical, multiplicative
      std::map<std::string, AError>& errs = const_cast<std::map<std::string, AError>&>(th_func->GetAllErrors());
      if ( errs.count(errName) ) errs[errName].SetIsTheo(true);
      if ( errs.count(errName) ) errs[errName].SetIsMult(true);
   }

   // reset initial PDF member
   SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), initialPdfMember, 0);

   return success;
}

//______________________________________________________________________________________________________//
bool ALHAPDFErrors::CalcHessianAndAddAsEigenvectors(std::vector<AFuncD*> th_funcs) {
   //! Assuming PDF set has 'hessian' uncertainties, calculate Alpos errors and add them to theory object.
   /*!
    *  Apart from the central value, 'hessian' PDF sets contain two
    *  additional member per eigenvector, indicating the asymmetric
    *  PDF uncertainty.
    *
    *  This method adds one Alpos error source per eigenvector
    *  to the theory object.
    *
    *  \param  th_func  Pointer to the theory function. This must depend
    *                   on the LHAPDF6 function instance referred to in
    *                   the task parameter 'LHAPDFfunction' in the steering
    *                   file.
    */
   info["Execute"] << "Calculating PDF errors as '" << fCalculateAs << "'." << std::endl;

   bool success = true;  // return success flag

   // LHAPDF instance information
   ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD(fLhapdfFunctionName);
   int initialPdfMember = PAR_ANY(fLhapdfFunctionName + std::string(".PDFSet"));

   // PDF set meta-information (e.g. Confidence Level, error type)
   std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
   unsigned int nMembers = lhapdfInstance->GetPDFSet()->size();
   double setCL = lhapdfInstance->GetPDFSet()->errorConfLevel() / 100.;
   double reqCL = 0.6826894 ;//boost::math::erf(1 / sqrt(2));  // required confidence level: 68.2689... %
   double errScale = 1.0;  // scaling due to set CL != required CL

   // cross sections for each member
   std::vector<std::vector<std::vector<double>>> observableArrays = this->CalcObservableArrays(th_funcs);

   // loop over theory functions
   for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) {
      AFuncD* th_func = th_funcs[iTheoryFunction];
      std::vector<std::vector<double>>& observableArray = observableArrays[iTheoryFunction];
      const unsigned int nObsBins = th_func->N();

      // central cross section values
      std::vector<double> observablesCentral(nObsBins);
      for (int iObs = 0; iObs < nObsBins; iObs++) {
         observablesCentral[iObs] = observableArray[iObs][0];
      }

      // --- calculate uncertainties from per-member observables

      // handle error scaling, if necessary
      if (setCL != reqCL) {
         boost::math::chi_squared chiSquareDistribution(1);
         double qsetCL = boost::math::quantile(chiSquareDistribution, setCL);
         double qreqCL = boost::math::quantile(chiSquareDistribution, reqCL);
         errScale = sqrt(qreqCL / qsetCL);
         debug["CalcHessianAndAddAsEigenvectors"] << "PDFSet CL is " << setCL * 100 << " %, but " << reqCL * 100
                                                  << " % was requested. " << "Errors will be scaled by a factor of "
                                                  << errScale << "." << std::endl;
      }

      // go through all eigenvectors (hence all pdf members) and calculate observables
      unsigned int nEV = (nMembers - 1) / 2;
      for (int iEV = 1; iEV <= nEV; iEV++) {
         std::vector<double> errUp(nObsBins);
         std::vector<double> errDn(nObsBins);
         for (int iObs = 0; iObs < nObsBins; iObs++) {
            // TODO: review these
            errUp[iObs] = (observableArray[iObs][2 * iEV - 1] - observablesCentral[iObs]) * errScale;
            errDn[iObs] = (observableArray[iObs][2 * iEV] - observablesCentral[iObs]) * errScale;
            //errUp[iObs] =  std::max(std::max(observableArray[iObs][2*iEV-1] - observablesCentral[iObs], observableArray[iObs][2*iEV] - observablesCentral[iObs]), 0.0) * errScale;
            //errDn[iObs] = -std::max(std::max(observablesCentral[iObs] - observableArray[iObs][2*iEV-1], observablesCentral[iObs] - observableArray[iObs][2*iEV]), 0.0) * errScale;
         }
         // determine error name
         std::string errName = fErrorPrefix + std::string("_EV_") + std::to_string(iEV);
         // add error and log success
         // Note: use linear averaging, since this is used in LHAPDF
         success &= th_func->AddAsymError(errName, errUp, errDn, observablesCentral, false, 1., "",
                                          AError::kLinear);

         // set error flags: theoretical, multiplicative
         std::map<std::string, AError>& errs = const_cast<std::map<std::string, AError>&>(th_func->GetAllErrors());
         if ( errs.count(errName) ) errs[errName].SetIsTheo(true);
         if ( errs.count(errName) ) errs[errName].SetIsMult(true);

      }
   }

   // reset initial PDF member
   SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), initialPdfMember, 0);

   return success;
}

//______________________________________________________________________________________________________//
bool ALHAPDFErrors::CalcSymHessianAndAddAsEigenvectors(std::vector<AFuncD*> th_funcs) {
   //! Assuming PDF set has 'symmhessian' uncertainties, calculate Alpos errors and add them to theory object.
   /*!
    *  Apart from the central value, 'symmhessian' PDF sets contain one
    *  additional member per eigenvector, indicating the symmetric
    *  PDF uncertainty.
    *
    *  This method adds one Alpos error source per eigenvector
    *  to the theory object.
    *
    *  \param  th_func  Pointer to the theory function. This must depend
    *                   on the LHAPDF6 function instance referred to in
    *                   the task parameter 'LHAPDFfunction' in the steering
    *                   file.
    */
   info["Execute"] << "Calculating PDF errors as '" << fCalculateAs << "'." << std::endl;

   bool success = true;  // return success flag

   // LHAPDF instance information
   ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD(fLhapdfFunctionName);
   int initialPdfMember = PAR_ANY(fLhapdfFunctionName + std::string(".PDFSet"));

   // PDF set meta-information (e.g. Confidence Level, error type)
   std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
   unsigned int nMembers = lhapdfInstance->GetPDFSet()->size();
   double setCL = lhapdfInstance->GetPDFSet()->errorConfLevel() / 100.;
   double reqCL = 0.6826894 ;//boost::math::erf(1 / sqrt(2));  // required confidence level: 68.2689... %
   double errScale = 1.0;  // scaling due to set CL != required CL

   // cross sections for each member
   std::vector<std::vector<std::vector<double>>> observableArrays = this->CalcObservableArrays(th_funcs);

   // loop over theory functions
   for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) {
      AFuncD* th_func = th_funcs[iTheoryFunction];
      std::vector<std::vector<double>>& observableArray = observableArrays[iTheoryFunction];
      const unsigned int nObsBins = th_func->N();

      // central cross section values
      std::vector<double> observablesCentral(nObsBins);
      for (int iObs = 0; iObs < nObsBins; iObs++) {
         observablesCentral[iObs] = observableArray[iObs][0];
      }

      // --- calculate uncertainties from per-member observables

      // handle error scaling, if necessary
      if (setCL != reqCL) {
         boost::math::chi_squared chiSquareDistribution(1);
         double qsetCL = boost::math::quantile(chiSquareDistribution, setCL);
         double qreqCL = boost::math::quantile(chiSquareDistribution, reqCL);
         errScale = sqrt(qreqCL / qsetCL);
         debug["CalcSymHessianAndAddAsEigenvectors"] << "PDFSet CL is " << setCL * 100 << " %, but " << reqCL * 100
                                                     << " % was requested. " << "Errors will be scaled by a factor of "
                                                     << errScale << "." << std::endl;
      }

      // go through all eigenvectors (hence all pdf members as well) and calculate uncertainties
      unsigned int nEV = nMembers;

      for (int iEV = 1; iEV < nEV; iEV++) {
         std::vector<double> errSymm(nObsBins);
         for (int iObs = 0; iObs < nObsBins; iObs++) {
            errSymm[iObs] = (observableArray[iObs][iEV] - observablesCentral[iObs]) * errScale;
         }
         // determine error name
         std::string errName = fErrorPrefix + std::string("_EV_") + std::to_string(iEV);
         // add error and log success
         success &= th_func->AddSymError(errName, errSymm, observablesCentral, false, 1., "");

         // set error flags: theoretical, multiplicative
         std::map<std::string, AError>& errs = const_cast<std::map<std::string, AError>&>(th_func->GetAllErrors());
         if ( errs.count(errName) ) errs[errName].SetIsTheo(true);
         if ( errs.count(errName) ) errs[errName].SetIsMult(true);
      }
   }

   // reset initial PDF member
   SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), initialPdfMember, 0);

   return success;
}

//______________________________________________________________________________________________________//
std::vector<std::vector<double>> ALHAPDFErrors::CalcObservableArray(AFuncD* th_func) {
   //! Get theory predictions calculated with each member of a PDF set (one theory function only).
   /*!
    *  \sa  ALHAPDFErrors::CalcObservableArrays
    */
   return this->CalcObservableArrays({th_func})[0];
}

//______________________________________________________________________________________________________//
std::vector<std::vector<std::vector<double>>> ALHAPDFErrors::CalcObservableArrays(std::vector<AFuncD*> th_funcs) {
   //! Get theory predictions calculated with each member of a PDF set (all theory functions simultaneously).
   /*!
    *  This method goes through all PDF members of an LHAPDF6 PDF set
    *  and calculates theory predictions for each theory function in `th_funcs`
    *  simultaneaously.
    *
    *  Returns a list of lists organized as:
    *
    *       observableArrays[iTheoryFunction][iObservableBin][iPDFMember]
    *
    *  \param  th_func  Pointer to the theory function object which provides the
    *                   predictions. This must depend on the LHAPDF6 function
    *                   instance referred to in the task parameter 'LHAPDFfunction'
    *                   in the steering file.
    *
    */

   ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD(fLhapdfFunctionName);
   int initialPdfMember = PAR_ANY(fLhapdfFunctionName + std::string(".PDFSet"));

   std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
   unsigned int nMembers = lhapdfInstance->GetPDFSet()->size();

   // containers for predictions

   // loop over theory functions
   /*
   std::vector<std::vector<std::vector<double>>> observableArrays;
   for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) {
      AFuncD* th_func = th_funcs[iTheoryFunction];
      if ( th_func->GetAlposName().find("Super") != std::string::npos) 
	 info["CalcAndAddAsCovariance"]<<"Calculating bin-to-bin correlations between data sets in "<<th_func->GetAlposName()<<"."<<std::endl;
      const unsigned int nObsBins = th_func->N();
      std::vector<std::vector<double>> observableArray(nObsBins, std::vector<double>(nMembers));

      // go through all pdf members and calcuate predictions
      for (int iMember = 0; iMember < nMembers; iMember++) {
         SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), iMember, 0);
         const std::vector<double>& memberObservables = th_func->GetValues();
         for (int iObs = 0; iObs < nObsBins; iObs++) {
            observableArray[iObs][iMember] = memberObservables[iObs];
         }
      }

      // add array for current function to large return array
      observableArrays.push_back(observableArray);
   }
   */

   std::vector<std::vector<std::vector<double>>> observableArrays(th_funcs.size());
   for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) observableArrays[iTheoryFunction].resize(th_funcs[iTheoryFunction]->N());

   // go through all pdf members and calcuate predictions
   for (int iMember = 0; iMember < nMembers; iMember++) {
      SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), iMember, 0);

      for (unsigned int iTheoryFunction = 0; iTheoryFunction < th_funcs.size(); iTheoryFunction++) {
	 AFuncD* th_func = th_funcs[iTheoryFunction];
	 if ( th_func->GetAlposName().find("Super") != std::string::npos) 
	    info["CalcAndAddAsCovariance"]<<"Calculating bin-to-bin correlations between data sets in "<<th_func->GetAlposName()<<"."<<std::endl;
	 const unsigned int nObsBins = th_func->N();
         const std::vector<double>& memberObservables = th_func->GetValues();
         for (int iObs = 0; iObs < nObsBins; iObs++) {
	    observableArrays[iTheoryFunction][iObs].push_back(memberObservables[iObs]);
         }
      }
   }


   // reset initial PDF member
   SET_ANY(fLhapdfFunctionName + std::string(".PDFSet"), initialPdfMember, 0);
      lhapdfInstance->Update();

   return observableArrays;
}
