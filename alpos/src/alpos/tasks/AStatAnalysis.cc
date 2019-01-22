#include "alpos/tasks/AStatAnalysis.h"
#include "alpos/AFactory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/AData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include <iostream>
#include "TMath.h"

/** 
 AStatAnalysis

 Perform a statistical analysis of the theory.

 */

using namespace std;

const string AStatAnalysis::fTaskType = "StatAnalysis";

//____________________________________________________________________________________ //
//AStatAnalysis::AStatAnalysis(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AStatAnalysis::AStatAnalysis(const string& aname ) : ATask(aname) {
   //! constructor

   //! Important: create always new result-object here!
   fResult = new AStatAnalysisResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AStatAnalysis::~AStatAnalysis() {
   //! destructor.
   //! do not delete the AResult object
   if ( fChisq ) delete fChisq;
}


//____________________________________________________________________________________ //
bool AStatAnalysis::Init() {
   // --- setting up minimzer
   vector<string> parnames;
   string chisqdef = STRING_NS(Chisq,NS());
   fChisq = AFactory::ChisqFactory(chisqdef,parnames,
				   (AData*)TheoryHandler::Handler()->GetFuncD("SuperData"),
				   TheoryHandler::Handler()->GetFuncD("SuperTheory"));
   if ( !fChisq ) {
      cout<<"ERROR. AStatAnalysis:Init(). Failed to initialize chisq with name: "<<chisqdef<<endl;
      exit(1);
   }
   
   return true;
}


//____________________________________________________________________________________ //
void AStatAnalysis::Print(const string& dat, const string& th) {
   Print((AData*)TheoryHandler::Handler()->GetFuncD(dat),TheoryHandler::Handler()->GetFuncD(th));
}


//____________________________________________________________________________________ //
void AStatAnalysis::Print(AData* dat, AFuncD* th) {

   if ( dat->N() == 0 ) return; // nothing todo
   vector<string> parnames;
   const string& ftyp = th->GetFunctionName();
   // chisq
   string chisqdef = STRING_NS(Chisq,NS());
   AChisqBase* chi = AFactory::ChisqFactory(chisqdef,parnames,dat,th);

   double chisq = chi->DoEval(NULL);
   std::map<std::string, double> nui = chi->GetNuisanceParameters();

   AChisqBase* chiSU = NULL;
   double chisqSU = 0;
   std::map<std::string, double> nuiSU;
   if ( EXIST_NS(Chisq2,NS()) && STRING_NS(Chisq2,NS()) != "" ) {
      string chisqdef = STRING_NS(Chisq2,NS());
      chiSU = AFactory::ChisqFactory(chisqdef,parnames,dat,th);
      chisqSU = chiSU->DoEval(NULL);
      nuiSU = chiSU->GetNuisanceParameters();

   }


   int ndf = chi->Data()->N();
   // pull
   APull* pull = (APull*)AFactory::ChisqFactory("Pull",parnames,dat,th);

   double pMean = pull->CalcMeanStatUncorr();
   double pRMS  = pull->CalcRMSStatUncorr();
   double pMedi = pull->CalcMedianStatUncorr();
   double pMin = pull->CalcMinStatUncorr();
   double pMax = pull->CalcMaxStatUncorr();
   
   double pMeanTot = pull->CalcMeanTotErr();
   double pRMSTot  = pull->CalcRMSTotErr();
   double pMediTot = pull->CalcMedianTotErr();
   double pMinTot = pull->CalcMinTotErr();
   double pMaxTot = pull->CalcMaxTotErr();

   // printout
   cout<<" + --- "<<dat->GetAlposName()<<" | "<<th->GetAlposName()<<" ("<<ftyp<<")"<<endl;
   cout<<" |"<<endl;
   cout<<" | Summary: ";
   printf("%-30s\t%8.2f\t%6d\n",th->GetAlposName().c_str(),chisq,ndf);
   cout<<" |  ndf                    "<<ndf<<endl;
   cout<<" |  chisq-def              "<< STRING_NS(Chisq,NS()) <<endl;
   cout<<" |     chisq               "<<chisq<<endl;
   cout<<" |     chisq/ndf           "<<chisq/ndf<<endl;
   cout<<" |     Prob                "<<TMath::Prob(chisq,ndf)<<endl;
   chi->PrintNuisanceParameters(true,&nui);
   if ( chisqSU ) {
      cout<<" |  chisq-def StatUncor    "<< STRING_NS(Chisq2,NS()) <<endl;
      cout<<" |     chisq               "<<chisqSU<<endl;
      cout<<" |     chisq/ndf           "<<chisqSU/ndf<<endl;
      cout<<" |     Prob                "<<TMath::Prob(chisqSU,ndf)<<endl;
      chiSU->PrintNuisanceParameters(true,&nuiSU);
   }
   cout<<" |  Pull (stat+unc, tot):  "<<endl;
   cout<<" |     mean                "<<pMean<<", \t"<<pMeanTot<<endl;           
   cout<<" |     median:             "<<pMedi<<", \t"<<pMediTot<<endl;           
   cout<<" |     RMS:                "<<pRMS<<", \t"<<pRMSTot<<endl;           
   cout<<" |     max:                "<<pMax<<", \t"<<pMaxTot<<endl;           
   cout<<" |     min:                "<<pMin<<", \t"<<pMinTot<<endl;           
   cout<<" |  "<<endl;
   delete chi;
   delete pull;
   if ( chisqSU ) delete chiSU;
}


//____________________________________________________________________________________ //
bool AStatAnalysis::Execute() {

   cout<<endl;
   cout<<" +--------------------  StatAnalysis  ---------------------- "<<endl;
   // super
   const auto& super = TheoryHandler::Handler()->GetSuperPair();
   Print(super.first,super.second);
   // single datasets
   for ( const auto& id : TheoryHandler::Handler()->GetDataTheoryPairs() ) 
      Print(id.second.first,id.second.second);
   // subsets
   for ( const auto& is: TheoryHandler::Handler()->GetAllSubsetPairs() ) 
 	 Print(is.second.first,is.second.second);
   cout<<" +---------------------------------------------------------- "<<endl;
   cout<<endl;
   


   // --- job done successfully
   return true;
   
}


//____________________________________________________________________________________ //
