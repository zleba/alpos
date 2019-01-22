// DB 15.01.2015

#include "alpos/ASuperData.h"
#include "alpos/Alpos.h"
#include <iostream>

using namespace std;


// __________________________________________________________________________________________ //
//const std::vector<std::string> ASuperData::fRequirements = {}; //< Set requirements
const std::vector<std::string> ASuperData::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ASuperData::fFunctionName = "SuperData"; //< The function's name


// __________________________________________________________________________________________ //
ASuperData::ASuperData(const std::string& name) :  AData(name,true) { //AParmFuncBase<double>(name) { 
   SetName(name); // revert the "_Data" statement
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
	SetClassName("ASuperData");
}


// __________________________________________________________________________________________ //
ASuperData::~ASuperData() {
}


// ___________________________________________________________________________________________ //
bool ASuperData::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   
   return true;
   //return Update();
}


// __________________________________________________________________________________________ //
bool ASuperData::Update() {
   using namespace AlposTools;
   fValue.clear();
   fError.clear();
   for ( const auto& i : GetRequirements()) {
      //#define VALUES(X)        TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetValues()
      fValue += VALUES_ANY(GetAlposName()+std::string(".")+i);
   }
   fError.resize(fValue.size()); // which one should go here?

	// update pointers to data tables of "child" AData objects
	fChildren.clear();
	for ( const auto dtp : TheoryHandler::Handler()->GetDataTheoryPairs() ) {
		fChildren.push_back(dtp.second.first);
	}

   // --- calculate 'new' errors
   CalculateSuperErrors();

   return true;
}


// __________________________________________________________________________________________ //
void ASuperData::CalculateSuperErrors() {
   //! Calculate errors of the super array
   //! These are:
   using namespace AlposTools;
   
   //vector<double> empty(fValue.size());
   // --- calculate fAllErrors
   // --- sym and asym errors:
   map<string,vector<double> > eup, edn; // all new error sources
   map<string,vector<vector<double>>> emat; //correlation matrix
   map<string,double> ecorrfrac; // all new corrfrac's
   map<string,string> enature; // all new corrfrac's
   vector<double> allsetszero;
   map<string,set<bool> > eIsStat; // stst
   map<string,set<bool> > eIsTheo; // theo
   map<string,set<bool> > eIsMult; // mult

   for ( const auto& i : GetRequirements() ) { // loop over all datasets
      AData* dat = (AData*)TheoryHandler::Handler()->GetFuncD(GetAlposName()+std::string(".")+i);
      const std::map<std::string, AError>& errset = dat->GetAllErrors();
		// if a dataset has multiplicative errors, set the flag for 'SuperData' as well
		if (dat->HasMultErrors()) {
			this->SetHasMultErrors();
		}
      // if new error source add it to eup,edn
      for ( const auto& ies : errset ) {
	 		if ( eup.count(ies.first) == 0 ) {
	    		//cout<<"Info. New SuperData Error souce: "<<ies.first<<endl;
	    		eup[ies.first] += allsetszero;
	    		edn[ies.first] += allsetszero;
	    		ecorrfrac[ies.first] =  ies.second.GetCorrelatedFraction();
	    		enature[ies.first] = ies.second.GetNature();
	    		eIsStat[ies.first].insert(ies.second.GetIsStat());
	         eIsTheo[ies.first].insert(ies.second.GetIsTheo());
	         eIsMult[ies.first].insert(ies.second.GetIsMult());
	    		// build an empty matrix for matrix type uncertainties
	    		if ( ies.second.GetIsMatType() ) {
	       		emat[ies.first].resize(fValue.size());
	       		for ( auto& v : emat[ies.first] )
		  				v.resize(fValue.size());
	       		for ( unsigned int x = 0 ; x<emat[ies.first].size() ; x++ ) // initialize diagonal elements with 1
		  				emat[ies.first][x][x] = 1;
	    		}
	 		}
      }

      // todo. If one uncertainty is of 'matrix' type in one of the datasets,
      //       but not in all, then 'transform' to matrix type uncertainty.

      // build error vectors and matrices
      vector<double> zero(dat->GetValues().size());
      for ( const auto& iea : eup ) {
	 		// if error is present in current dataset, add it
	 		if ( errset.count(iea.first) > 0 ) {//&& !errset.at(iea.first).GetIsMatType()) {
	    		// check is corrfrac and nature is identical
	    		if ( ecorrfrac[iea.first] != errset.at(iea.first).GetCorrelatedFraction() )
	       		warn["CalculateSuperErrors"]<<"Correlated fraction of error source '"<<iea.first<<"' is not identical within datasets."<<endl;
	    
	    		if ( enature[iea.first] != errset.at(iea.first).GetNature() )
	       		warn["CalculateSuperErrors"]<<"'Nature' of error source '"<<iea.first<<"' is not identical within datasets."<<endl;
	    
	    		//eup[iea.first] += errset.at(iea.first).GetErrorRelUp();
	    		//edn[iea.first] += errset.at(iea.first).GetErrorRelDn();
				eup[iea.first] += errset.at(iea.first).GetError("RelUpTotAll");
				edn[iea.first] += errset.at(iea.first).GetError("RelDnTotAll");


			}
	 		// otherwise add zero's
	 		else {
			   eup[iea.first] += zero;
			   edn[iea.first] += zero;
	 		}
	 		// if this source is of matrix type in one data set, then also copy matrix elements
	 		//if (errset.count(iea.first) > 0 && errset.at(iea.first).GetIsMatType() ) {
	 		if (errset.count(iea.first) > 0 && emat.count(iea.first) > 0 ) {
	    		ecorrfrac[iea.first] = -1;
	    		const TMatrixDSym& corr = errset.at(iea.first).GetCorrelations();
	    		vector<vector<double> >& mat = emat[iea.first];
			const TMatrixDSym* matinit =
			   fAllErrors.count(iea.first) ?
			   &fAllErrors[iea.first].GetCorrelations() :
			   NULL;
			debug["CalculateSuperErrors"]<<"Found matrix type uncertainty. Keeping bin-to-bin correlations between measurements."<<endl;
	    		unsigned int noff = allsetszero.size();
	    		for ( int x = 0 ; x<corr.GetNrows() ;x++ ) {
			   for ( int y = 0 ; y<corr.GetNrows() ;y++ ) {
			      mat[x+noff][y+noff] = corr[x][y];
			      //cout<<"x="<<x<<", y="<<y<<"\tcorr="<<corr[x][y]<<"\tdiag="<<corr[y][x]<<endl;
			   }
			   // // keep bin-to-bin correlations among data sets for matrix type uncertainites
			   if ( matinit ) {
			      for ( int Y = 0 ; Y<mat.size() ;Y++ ) {
				 if ( Y>= noff && Y < noff+corr.GetNrows() )
				    continue;
				 mat[x+noff][Y] = (*matinit)[x+noff][Y];
				 //cout<<iea.first<<"\tx="<<x<<",   "<<" "<<"\txSuper="<<x+noff<<", ySuper="<<Y<<"\tcorr="<<mat[x+noff][Y]<<endl;
			      }
			   }
			}
			}
      }
      
      // resize also 'zero's
      allsetszero += zero;
   }
   
   // todo: propagate 'error averaging (AError::kLinear) or get from steering...'

   // --- create AErrors and fAllErrors
   fAllErrors.clear();
   for ( const auto& iea : eup ){
      AError err(iea.first,"");
      fAllErrors[iea.first] = err;
      if ( ecorrfrac[iea.first] >= 0  )
	 //fAllErrors[iea.first].SetAsymError(iea.second,edn[iea.first],fValue,true,ecorrfrac[iea.first],enature[iea.first],AError::kSignImprovedQuadratic);
	 fAllErrors[iea.first].SetAsymError(iea.second,edn[iea.first],fValue,true,ecorrfrac[iea.first],enature[iea.first],Alpos::Current()->Settings()->ErrorSymmetrization);
      else {
	 fAllErrors[iea.first].SetMatrixError(iea.second ,emat[iea.first],fValue,true,true,enature[iea.first]);
      }
      fAllErrors[iea.first].SetIsStat(*eIsStat[iea.first].begin());
      fAllErrors[iea.first].SetIsTheo(*eIsTheo[iea.first].begin());
      fAllErrors[iea.first].SetIsMult(*eIsMult[iea.first].begin());
      if ( eIsStat[iea.first].size() > 1 )
	      warn["CalculateSuperErrors"]<<"Error '"<<iea.first<<"' has specified different IsStat flags!"<<endl;
      if ( eIsTheo[iea.first].size() > 1 )
	      warn["CalculateSuperErrors"]<<"Error '"<<iea.first<<"' has specified different IsTeho flags!"<<endl;
      if ( eIsMult[iea.first].size() > 1 )
	      warn["CalculateSuperErrors"]<<"Error '"<<iea.first<<"' has specified different IsMult flags!"<<endl;
   }

	// --- clear the error cache
	fSumErrors.clear();
	fSumErrMats.clear();

	// not needed under new interface -> calculation happens on first call to GetSumErrorMatrix(...)
   //CalculateMatrices();

	// TODO: remnants of old interface -> remove once these become accessible through new interface
	// fErrMat, fErrMatRel (matrix-type uncertainties) excluding stat
	fErrMat.clear();
	fErrMatRel.clear();
	fErrMat.resize(fValue.size());
	fErrMatRel.resize(fValue.size());
	for ( const auto& ierr : fAllErrors ) {
		if ( !ierr.second.GetIsStat() && ierr.second.GetIsMatType()) {
			for ( unsigned int iea = 0 ; iea<fValue.size() ; iea++ ) {
				fErrMat[iea]    += pow(ierr.second.GetMatError()[iea],2) ;
				fErrMatRel[iea] += pow(ierr.second.GetMatErrorRel()[iea],2) ;
			}
		}
	}
	for ( unsigned int iea = 0 ; iea<fValue.size() ; iea++ ) {
		fErrMat[iea]    = sqrt(fErrMat[iea]);
		fErrMatRel[iea] = sqrt(fErrMatRel[iea]);
	}


}


// __________________________________________________________________________________________ //
