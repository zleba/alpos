#include "alpos/tasks/AReplaceDataWithTheoryValues.h"
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
#include <numeric>
#include "TH1D.h"

/** 
 AReplaceDataWithTheoryValues

 Replace cross section values of data (column 'Sigma') by 
 current theoretical predictions.
 Reset all (absolute) uncertainties.
 All uncertainties are re-calculated, as these are by
 default stored as relative uncertainties.

 */

using namespace std;

const string AReplaceDataWithTheoryValues::fTaskType = "ReplaceDataWithTheoryValues";

//____________________________________________________________________________________ //
AReplaceDataWithTheoryValues::AReplaceDataWithTheoryValues(const string& aname ) : ATask(aname) {
   //! constructor
   //! Important: create always new result-object here!
   fResult = new AReplaceDataWithTheoryValuesResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AReplaceDataWithTheoryValues::~AReplaceDataWithTheoryValues() {
   //! destructor.
   //! do not delete the AResult object
}


//____________________________________________________________________________________ //
bool AReplaceDataWithTheoryValues::Init() {
   // --- set up task
   return true;
}

//____________________________________________________________________________________ //
bool AReplaceDataWithTheoryValues::Execute() {
   
   info["Execute"]<<"Replacing all data cross sections with current theory predictions."<<endl;

   // pointers to superdata and supertheory
   /*
   const auto& supdata = TheoryHandler::Handler()->GetSuperPair().first;
   const auto& suptheo = TheoryHandler::Handler()->GetSuperPair().second;

   // get pointers to DataTables
   std::vector<AData*> dataChildren = supdata->GetChildren();           //TODO: handle subsets correctly
   std::vector<AParmFuncBase<double>*> theoChildren = suptheo->GetChildren();   //TODO: handle subsets correctly

   // go through each DataTable in SuperData (one per AData object)
   for (int iChild = 0 ; iChild < dataChildren.size() ; iChild++ ) {

      const std::vector<double>& dataPts = dataChildren[iChild]->GetValues();
      const std::vector<double>& theoPts = theoChildren[iChild]->GetValues();

      cout<<"Name: "<<dataChildren[iChild]->GetAlposName()<<endl;
      cout<<"before null: data="<<dataPts[0]<<"\ttheo="<<theoPts[0]<<endl;
      cout<<"before eins: data="<<dataPts[1]<<"\ttheo="<<theoPts[1]<<endl;
      bool IsConst = dataChildren[iChild]->GetIsConstant();
      if ( IsConst ) { 
	 // set theory as data values
	 dataChildren[iChild]->SetIsConstant(false);
	 dataChildren[iChild]->ChangeValues(theoPts,vector<double>(theoPts.size()));
      }
      else {
	 // revert changes!
	 // set data to data again
	 const std::map<std::string, std::vector<double> >& OrigData = dataChildren[iChild]->GetDataTable();
	 dataChildren[iChild]->ChangeValues(OrigData.at("Sigma"),vector<double>(theoPts.size()));
	 dataChildren[iChild]->SetIsConstant(true);
      }
      dataChildren[iChild]->ClearErrorCache();

      cout<<"after  null: data="<<dataPts[0]<<"\ttheo="<<theoPts[0]<<endl;
      cout<<"after  eins: data="<<dataPts[1]<<"\ttheo="<<theoPts[1]<<endl;
   }
   */

   const auto& DataTheoryPair = TheoryHandler::Handler()->GetDataTheoryPairs();
   const auto& SubsetPair = TheoryHandler::Handler()->GetSubsetPairs();

   // run over all data-functions, and below choose the right one...
   set<string> dtnames;
   for ( auto dt :  DataTheoryPair ) {
      dtnames.insert(dt.first) ;
   }
   for ( auto dt :  SubsetPair ) {
      dtnames.insert(dt.first) ;
   }
   cout<<endl;

   for ( auto dt :  SubsetPair ) {
      for ( auto dt2 :  dt.second ) {
	 dtnames.insert(dt2.first) ;
      }
   }

   // const std::map<std::string, std::pair<AData*, AFuncD*> >& GetDataTheoryPairs() const { return fDataTheoryPairs; }; 
   // const std::map<std::string, std::map<std::string, std::pair<ASubsetData*, ASubsetFunction*> > >& GetSubsetPairs() const { return fSubsetPairs; };
   for ( auto ipair : dtnames ) {


      AData* dt  = (AData* )TheoryHandler::Handler()->GetFuncD(ipair+"_Data");
      AFuncD* th = (AData* )TheoryHandler::Handler()->GetFuncD(ipair);
      
      //const AData* dt = ipair.second.first;
      const std::vector<double>& dataPts = dt->GetValues(); // dt->GetValues()
      const std::vector<double>& theoPts = th->GetValues();

      // string req = GetRequirements()[0];
      // AData* dat = (AData*)TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+req);

      bool IsConst = dt->GetIsConstant();
      vector<string> req = dt->GetRequirements();
      // cout<<endl;
      // cout<<"Name:  "<<dt->GetAlposName()<<endl;
      // cout<<"Name:  "<<dt->GetFunctionName()<<endl;
      // cout<<"Req:   "<<req.size()<<endl;
      // cout<<"Const: "<<IsConst<<endl;
      // for ( auto s : req ) cout<<"Req: "<<s<<endl;
      // cout<<"before null: data="<<dataPts[0]<<"\ttheo="<<theoPts[0]<<endl;
      // cout<<"before eins: data="<<dataPts[1]<<"\ttheo="<<theoPts[1]<<endl;
      
      if ( req.size() == 0 ) {  /// the 'actual' data
	 if ( IsConst ) { 
	    // set theory as data values
	    dt->SetIsConstant(false);
	    dt->ChangeValues(theoPts,vector<double>(theoPts.size()));
	 }
	 else {
	    cout<<"  **  REVERT CHANGES ** "<<endl;
	    // revert changes!
	    // set data to data again
	    const std::map<std::string, std::vector<double> >& OrigData = dt->GetDataTable();
	    dt->ChangeValues(OrigData.at("Sigma"),vector<double>(theoPts.size()));
	    dt->SetIsConstant(true);
	 }
      }
      dt->ClearErrorCache();
      // cout<<"after  null: data="<<dataPts[0]<<"\ttheo="<<theoPts[0]<<endl;
      // cout<<"after  eins: data="<<dataPts[1]<<"\ttheo="<<theoPts[1]<<endl;

   }
   
   // --- job done successfully
   return true;

}


//____________________________________________________________________________________ //
