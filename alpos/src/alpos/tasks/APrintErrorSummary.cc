#include "alpos/tasks/APrintErrorSummary.h"
#include "alpos/AFactory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/AData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
/** 
 APrintErrorSummary

 Print information on uncertainties for all datasets.

 */

using namespace std;

const string APrintErrorSummary::fTaskType = "PrintErrorSummary";

//____________________________________________________________________________________ //
//APrintErrorSummary::APrintErrorSummary(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
APrintErrorSummary::APrintErrorSummary(const string& aname ) : ATask(aname) {
   //! constructor
   //! Important: create always new result-object here!
   fResult = new APrintErrorSummaryResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
APrintErrorSummary::~APrintErrorSummary() {
   //! destructor.
   //! do not delete the AResult object
}


//____________________________________________________________________________________ //
bool APrintErrorSummary::Init() {
   // --- setting up minimzer
   return true;
}


//____________________________________________________________________________________ //
void APrintErrorSummary::FillSummary(const string& dat, const string& th) {
   FillSummary((AData*)TheoryHandler::Handler()->GetFuncD(dat),TheoryHandler::Handler()->GetFuncD(th));
}


//____________________________________________________________________________________ //
void APrintErrorSummary::FillSummary(AData* dat, AFuncD* th) {
}


//____________________________________________________________________________________ //
bool APrintErrorSummary::Execute() {

   cout<<endl;
   cout<<" +--------------------  ErrorSummary  ---------------------- "<<endl;

   // super
   const auto& supdat  = TheoryHandler::Handler()->GetSuperPair().first;
   const auto& suptheo = TheoryHandler::Handler()->GetSuperPair().second;
   fdatsets.push_back(supdat->GetAlposName());
   for ( int iThDt = 0 ; iThDt<2 ; iThDt++ ) {
      const std::map<std::string, AError>& errset = iThDt ? 
	 suptheo->GetAllErrors() :
	 supdat->GetAllErrors();
      for ( const auto& ierr : errset ) {
	 for ( int iThDt = 0 ; iThDt<2 ; iThDt++ ) {
	    cout<<"error name: "<<ierr.first<<endl;
	    fsummary[ierr.first];//second.GetErrorName()]; // add error
	    char buf[10];
	    double cf = ierr.second.GetCorrelatedFraction();
	    if ( int(cf) - cf == 0 ) 
	       sprintf(buf,"%1d",int(cf));
	    else sprintf(buf,"%3.1f",cf);
	    string txt = buf;
	    txt+=",";
	    //txt+=ierr.second.GetNature();
	    txt+=ierr.second.GetType();
	    fsummary[ierr.first][supdat->GetAlposName()] = txt; // add dataset
	 }
      }
   }

   for ( const auto& id : TheoryHandler::Handler()->GetDataTheoryPairs() )  {
      for ( int iThDt = 0 ; iThDt<2 ; iThDt++ ) {
	 const std::map<std::string, AError>& errset = iThDt ? 
	    id.second.first->GetAllErrors() :
	    id.second.second->GetAllErrors();
	 fdatsets.push_back(id.second.first->GetAlposName());
	 for ( const auto& ierr : errset ) {
	    fsummary[ierr.first]; // add error
	    char buf[10];
	    double cf = ierr.second.GetCorrelatedFraction();
	    if ( int(cf) - cf == 0 ) 
	       sprintf(buf,"%1d",int(cf));
	    else sprintf(buf,"%3.1f",cf);
	    string txt = buf;
	    txt+=",";
	    //txt+=ierr.second.GetNature();
	    txt+=ierr.second.GetType();
	    fsummary[ierr.first][id.second.first->GetAlposName()] = txt; // add dataset
	 }
      }
   }

   cout<<endl;
   // print
   string forme = "%-38s ";
   string formt = "%35s ";
   string formt2 = "%8s ";
   // header 
   // column names
   int nmax = 0;
   for ( auto id : fdatsets ) nmax = max(nmax,int(id.size()));
   vector<string> slds(fdatsets.size());
   for ( unsigned int i = 0 ; i< slds.size() ; i++ ) {
      slds[i].resize(nmax-fdatsets[i].size(),' ');
      slds[i]+=fdatsets[i];
   }

   for ( int i = 0 ; i< nmax ; i++ ) {
      printf(forme.c_str(),"");
      for ( auto id : slds ) {
	 printf(formt2.c_str(),id.substr(i,1).c_str());
      }
      cout<<endl;
   }
   // printf(forme.c_str(),"");
   // for ( auto id : fdatsets ) {
   //    printf(formt.c_str(),id.c_str());
   // }
   cout<<endl;
   // entries
   for ( auto is : fsummary) {
      printf(forme.c_str(),is.first.c_str());
      for ( auto id : fdatsets ) {
	 printf(formt2.c_str(),is.second[id].c_str());
      }
      cout<<endl;
   }
   cout<<endl;

   cout<<" +---------------------------------------------------------- "<<endl;
   cout<<endl;

   // --- job done successfully
   return true;
   
}


//____________________________________________________________________________________ //
