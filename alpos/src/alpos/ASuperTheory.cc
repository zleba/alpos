// DB 15.01.2015

#include "alpos/ASuperTheory.h"
#include <iostream>
#include "alpos/AlposTools.h"
#include "alpos/Alpos.h"

using namespace std;


// __________________________________________________________________________________________ //
//const std::vector<std::string> ASuperTheory::fRequirements = {}; //< Set requirements
const std::vector<std::string> ASuperTheory::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ASuperTheory::fFunctionName = "SuperTheory"; //< The function's name


// __________________________________________________________________________________________ //
ASuperTheory::ASuperTheory(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("ASuperTheory");
}


// __________________________________________________________________________________________ //
ASuperTheory::~ASuperTheory() {
}


// ___________________________________________________________________________________________ //
bool ASuperTheory::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //Update();
   return true;
}


// __________________________________________________________________________________________ //
bool ASuperTheory::Update() {
   using namespace AlposTools;
   fValue.clear();
   fError.clear();

   /*
   for(auto f : fValue)
      cout <<"RADEK-f: "<<  f << endl;
   for ( const auto& i : GetRequirements()) {
      cout << "RADEK-i: "<< i << endl;
   }
   */

   vector<pair<string,vector<double>>> vals;
   for ( const auto& i : GetRequirements()) {
      vals.push_back(make_pair(i, vector<double>({})));
   }

//#pragma omp parallel for
   for(int k = 0; k < vals.size(); ++k) {
      string i = vals[k].first;
      vals[k].second = VALUES_ANY(GetAlposName()+std::string(".")+i);
      //cout << "Radek inside super-i: "<< i << endl;
   }

   //Merge them
   for(int k = 0; k < vals.size(); ++k) {
      fValue += vals[k].second;
   }

   /* RADEK comment
   for ( const auto& i : GetRequirements()) {
      //#define VALUES(X)        TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetValues()
      fValue += VALUES_ANY(GetAlposName()+std::string(".")+i);
      cout << "Radek inside super-i: "<< i << endl;
   }
   */

   fError.resize(fValue.size()); // not yet implemented

   // update pointers to data tables of "child" theory function objects
   fChildren.clear();
   for ( const auto dtp : TheoryHandler::Handler()->GetDataTheoryPairs() ) {
      fChildren.push_back(dtp.second.second);
   }

   // --- calculate 'new' errors
   CalculateSuperErrors();

   return true;
}


// __________________________________________________________________________________________ //
void ASuperTheory::CalculateSuperErrors() {
   //! Calculate errors of the super array
   //! These are:
   using namespace AlposTools;

   //vector<double> empty(fValue.size());
   // --- calculate fAllErrors
   // --- sym and asym errors:
   map<string,vector<double>> eup, edn; // all new error sources
   map<string,vector<vector<double>>> emat; //correlation matrix
   map<string,double> ecorrfrac; // all new corrfrac's
   map<string,string> enature; // all new corrfrac's
   vector<double> allsetszero;
   map<string,set<bool> > eIsStat; // stst
   map<string,set<bool> > eIsTheo; // theo
   map<string,set<bool> > eIsMult; // mult


   // first loop: extract meta-information from 'child' Theory functions and set up empty error lists/matrices
   for ( const auto& i : GetRequirements() ) { // loop over all datasets
      //AData* dat = (AData*)TheoryHandler::Handler()->GetFuncD(GetAlposName()+std::string(".")+i);
      AParmFuncBase* th = (AParmFuncBase*)TheoryHandler::Handler()->GetFuncD(GetAlposName() + std::string(".") + i);
      const std::map<std::string, AError>& errset = th->GetAllErrors();

      for ( const auto& ies : errset ) {
         // if new error source add it to eup,edn
         if ( eup.count(ies.first) == 0 ) {
            //cout<<"Info. New SuperTheory Error souce: "<<ies.first<<endl;
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

      // second loop: build error vectors and matrices
      vector<double> zero(th->GetValues().size());
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
	    info["CalculateSuperErrors"]<<"Found matrix type uncertainty. Keeping bin-to-bin correlations between measurements."<<endl;
            unsigned int noff = allsetszero.size();
            for ( int x = 0 ; x<corr.GetNrows() ;x++ ) {
               for ( int y = 0 ; y<corr.GetNrows() ;y++ ) {
                  mat[x+noff][y+noff] = corr[x][y];
                  //cout<<iea.first<<"\tx="<<x<<", y="<<y<<"\txSuper="<<x+noff<<", ySuper="<<y+noff<<"\tcorr="<<corr[x][y]<<"\tdiag="<<corr[y][x]<<endl;
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


   // --- third loop: create 'super' AErrors (in fAllErrors)
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

   // // --- clear the error cache
   // fSumErrors.clear();
   // fSumErrMats.clear();

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
