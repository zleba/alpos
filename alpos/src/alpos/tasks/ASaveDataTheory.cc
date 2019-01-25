#include "alpos/tasks/ASaveDataTheory.h"
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
#include "TNtuple.h"





/* 
 ATask

 */

using namespace std;

const string ASaveDataTheory::fTaskType = "ASaveDataTheory";

//____________________________________________________________________________________ //
//ASaveDataTheory::ASaveDataTheory(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
ASaveDataTheory::ASaveDataTheory(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("ASaveDataTheory");
   //! Important: create always new result-object here!
   fResult = new ASaveDataTheoryResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
ASaveDataTheory::~ASaveDataTheory(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool ASaveDataTheory::Init(){
   info<<"Hello. ASaveDataTheory::Init()."<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool ASaveDataTheory::Execute(){
   debug["Execute"]<<"Now getting 'WelcomeString' from steering and printing it:"<<endl;
   cout<<endl;
   cout<<"  "<<STRING_NS(WelcomeString,NS())<<endl;
   cout<<endl;

   TDirectory* file = Alpos::Current()->Settings()->rootoutput;
   file->mkdir(GetTaskName().c_str())->cd();


   const auto& supdata = TheoryHandler::Handler()->GetSuperPair().first;
   const auto& suptheo = TheoryHandler::Handler()->GetSuperPair().second;

   std::vector<AData*> dataChildren = supdata->GetChildren();           //TODO: handle subsets correctly
   std::vector<AParmFuncBase<double>*> theoChildren = suptheo->GetChildren();   //TODO: handle subsets correctly

   int iChild = 0;

   const std::vector<double>* dataPts = &dataChildren[iChild]->GetValues();
   const std::vector<double>* theoPts = &theoChildren[iChild]->GetValues();

   auto dataTable = dataChildren[iChild]->GetDataTable();

   vector<double> q2      = dataTable.at("Q2");
   vector<double> xpom    = dataTable.at("xp");
   vector<double> beta    = dataTable.at("beta");
   vector<double> dataUnc = dataTable.at("tot");

   TTree *ThDataTab = new TTree("ThDataTab","table with data and theory");

   double xp_, q2_, beta_, xpSigData_, xpSigDataErr_, xpSigTh_, xpSigThErr_;
   ThDataTab->Branch("xp",&xp_,"xp/D");
   ThDataTab->Branch("Q2",&q2_,"Q2/D");
   ThDataTab->Branch("beta",&beta_,"beta/D");
   ThDataTab->Branch("xpSigData",&xpSigData_,"xpSigData/D");
   ThDataTab->Branch("xpSigDataErr",&xpSigDataErr_,"xpSigDataErr/D");
   ThDataTab->Branch("xpSigTh",&xpSigTh_,"xpSigTh/D");
   ThDataTab->Branch("xpSigThErr",&xpSigThErr_,"xpSigThErr/D");

   for(int i = 0; i < dataPts->size(); ++i) {
       //cout << " "<<beta[i] <<" "<< q2[i] <<" "<<  xpom[i] <<" : "<< dataPts->at(i) << endl;
       xp_ = xpom[i];
       q2_ = q2[i];
       beta_ = beta[i];
       xpSigData_ = dataPts->at(i);
       xpSigDataErr_ = dataUnc[i] * 0.01; //from % to relErr
       xpSigTh_ = theoPts->at(i);
       xpSigThErr_ = 0;
       ThDataTab->Fill();
   }
   //cout << "Helenka " << endl;
   //exit(0);

   const std::vector<double>* errsData = &dataChildren[iChild]->GetSumError("AA", "AbsAvTot");
   const std::vector<double>* errsTheo = &theoChildren[iChild]->GetSumError("AA", "AbsAvTot");

   file->Write();

   return true;
}


//____________________________________________________________________________________ //
