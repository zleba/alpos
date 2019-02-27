#include "alpos/tasks/ASaveDataTheory.h"
#include "alpos/AFactory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/AData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include "alpos/ATheory.h"
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include "TH1D.h"
#include "TNtuple.h"


extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}

//Multiplied by xPom
double getSigRed2006(int ifit, double xpom, double q2, double z, double sSqrt = 319) //1==FitA, 2==FitB
{
    double xPq[13], f2[2], fl[2], c2[2], cl[2];
    qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);

    double t = -1;
    int Int = 1;
    int ipom;

    double fluxIP, fluxIR;
    ipom = 1;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &fluxIP);
    ipom = 2;
    h12006flux_(&xpom, &t, &Int, &ifit, &ipom, &fluxIR);


    //double Ep = (q2 < 120) ? 820 : 920;
    double Ee = 27.5;
    const double s = sSqrt;//  4*Ep * Ee;

    const double mp2 = pow(0.92, 2);
    double x = z*xpom;
    double y = q2/(s-mp2)/x;
    double yplus  = 1+pow(1-y,2);


    double xpSigRed_IP =  fluxIP*xpom * (f2[0]  - y*y/yplus*fl[0]);
    double xpSigRed_IR =  fluxIR*xpom * (f2[1]  - y*y/yplus*fl[1]);

    return xpSigRed_IP  + xpSigRed_IR;
}

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

   std::vector<AData*> dataChildren = supdata->GetChildren();           
   std::vector<AParmFuncBase<double>*> theoChildren = suptheo->GetChildren();   

   cout << "dataSize theorSize " << dataChildren.size() << " "<< theoChildren.size() << endl;


   for(int iChild = 0; iChild < dataChildren.size(); ++iChild) {

      TString name = "sample_";
      name +=  theoChildren[iChild]->GetAlposName();
      if(name.Contains(':'))
         name = name(0, name.First(':'));
      name = name.ReplaceAll('-','_');

      const std::vector<double>* dataPts = &dataChildren[iChild]->GetValues();
      const std::vector<double>* theoPts = &theoChildren[iChild]->GetValues();

      auto dataTable = dataChildren[iChild]->GetDataTable();

      if(!name.Contains("jets")) {

         vector<double> q2      = dataTable.at("Q2");
         vector<double> xpom    = dataTable.at("xp");
         vector<double> beta    = dataTable.at("beta");
         vector<double> dataUnc = dataTable.at("tot");

         assert(q2.size() == dataPts->size());


         TTree *ThDataTab = new TTree(name,"table with data and theory");

         double xp_, q2_, beta_, xpSigData_, xpSigDataErr_, xpSigTh_, xpSigThErr_;
         double xpSigThOrgA_, xpSigThOrgB_;

         ThDataTab->Branch("xp",&xp_,"xp/D");
         ThDataTab->Branch("Q2",&q2_,"Q2/D");
         ThDataTab->Branch("beta",&beta_,"beta/D");
         ThDataTab->Branch("xpSigData",&xpSigData_,"xpSigData/D");
         ThDataTab->Branch("xpSigDataErr",&xpSigDataErr_,"xpSigDataErr/D");
         ThDataTab->Branch("xpSigTh",&xpSigTh_,"xpSigTh/D");
         ThDataTab->Branch("xpSigThErr",&xpSigThErr_,"xpSigThErr/D");

         ThDataTab->Branch("xpSigThOrgA",&xpSigThOrgA_,"xpSigThOrgA/D");
         ThDataTab->Branch("xpSigThOrgB",&xpSigThOrgB_,"xpSigThOrgB/D");

         double sSqrt = 319;
         if(name.Contains("225"))
            sSqrt = 225;
         else if(name.Contains("252"))
            sSqrt = 252;
         else if(name.Contains("H1incDDIS-HERA-I-SpacMB") || name.Contains("H1incDDIS-HERA-I-SpacTrg"))
            sSqrt = 301;


         for(int i = 0; i < dataPts->size(); ++i) {
            //cout << " "<<beta[i] <<" "<< q2[i] <<" "<<  xpom[i] <<" : "<< dataPts->at(i) << endl;
            xp_ = xpom[i];
            q2_ = q2[i];
            beta_ = beta[i];
            xpSigData_ = dataPts->at(i);
            xpSigDataErr_ = dataUnc[i] * 0.01; //from % to relErr
            xpSigTh_ = theoPts->at(i);
            xpSigThErr_ = 0;

            xpSigThOrgA_ = getSigRed2006(1, xp_, q2_, beta_, sSqrt);
            xpSigThOrgB_ = getSigRed2006(2, xp_, q2_, beta_, sSqrt);

            ThDataTab->Fill();
         }
      }
      else { //for jets
         vector<double> q2_min      = dataTable.at("Q2_min");
         vector<double> q2_max      = dataTable.at("Q2_max");
         vector<double> pt_min      = dataTable.at("Pt_min");
         vector<double> pt_max      = dataTable.at("Pt_max");

         vector<double> dataUnc = dataTable.at("tot.(%)");
         assert(q2_min.size() == dataPts->size());

         TTree *ThDataTab = new TTree(name,"table with data and theory");


         double xp_, q2_, beta_, SigData_, SigDataErr_, SigTh_, SigThErr_;


         double q2_min_, q2_max_, pt_min_, pt_max_; 

         ThDataTab->Branch("q2_min",&q2_min_,"q2_min/D");
         ThDataTab->Branch("q2_max",&q2_max_,"q2_max/D");

         ThDataTab->Branch("pt_min",&pt_min_,"pt_min/D");
         ThDataTab->Branch("pt_max",&pt_max_,"pt_max/D");


         ThDataTab->Branch("SigData",&SigData_,"SigData/D");
         ThDataTab->Branch("SigDataErr",&SigDataErr_,"SigDataErr/D");
         ThDataTab->Branch("SigTh",&SigTh_,"SigTh/D");
         ThDataTab->Branch("SigThErr",&SigThErr_,"SigThErr/D");


         for(int i = 0; i < dataPts->size(); ++i) {
            //cout << " "<<beta[i] <<" "<< q2[i] <<" "<<  xpom[i] <<" : "<< dataPts->at(i) << endl;
            q2_min_ = q2_min[i];
            q2_max_ = q2_max[i];
            pt_min_ = pt_min[i];
            pt_max_ = pt_max[i];

            SigData_ = dataPts->at(i);
            SigDataErr_ = dataUnc[i] * 0.01; //from % to relErr
            SigTh_ = theoPts->at(i);
            SigThErr_ = 0;

            ThDataTab->Fill();
         }
      }

   }

   //cout << "Helenka " << endl;
   //exit(0);

   //const std::vector<double>* errsData = &dataChildren[iChild]->GetSumError("AA", "AbsAvTot");
   //const std::vector<double>* errsTheo = &theoChildren[iChild]->GetSumError("AA", "AbsAvTot");

   //string asfunc = STRING_NS(Alphas,NS());
   //cout << asfunc << endl;

   vector<string> parameters = {
           "QcdnumInit.AlphasMz"  ,
           "QcdnumInit.ScaleFacMuR",
           "QcdnumInit.ScaleFacMuF",
           "QcdnumInit.mcharm",
           "QcdnumInit.mbottom",
           "QcdnumInit.Q0",
           "DPDF.Flux_pom1_ap",
           "DPDF.Flux_pom1_b0",
           "DPDF.Flux_reg1_a0",
           "DPDF.Flux_reg1_ap",
           "DPDF.Flux_reg1_b0"};


   TH1D *hFixedPars = new TH1D("hFixedPars", "Fixed Parameters varied by theor-unc.", parameters.size(), -0.5, parameters.size() - 0.5);
   for(int i = 0; i < parameters.size(); ++i) {
       hFixedPars->SetBinContent(i+1, PAR_ANY(parameters[i]));
       hFixedPars->GetXaxis()->SetBinLabel(i+1, parameters[i].c_str());
   }

   /*
   cout << "Radek : " << PAR_ANY("QcdnumInit.AlphasMz")  << endl;
   cout << "Radek : " << PAR_ANY("QcdnumInit.ScaleFacMuR")  << endl;
   cout << "Radek : " << PAR_ANY("QcdnumInit.ScaleFacMuF")  << endl;
   cout << "Radek : " << PAR_ANY("QcdnumInit.mcharm")  << endl;
   cout << "Radek : " << PAR_ANY("QcdnumInit.mbottom")  << endl;
   cout << "Radek : " << PAR_ANY("QcdnumInit.Q0")  << endl;
   cout << "Radek : " << PAR_ANY("DPDF.Flux_pom1_ap")  << endl;
   cout << "Radek : " << PAR_ANY("DPDF.Flux_pom1_b0")  << endl;
   cout << "Radek : " << PAR_ANY("DPDF.Flux_reg1_a0")  << endl;
   cout << "Radek : " << PAR_ANY("DPDF.Flux_reg1_ap")  << endl;
   cout << "Radek : " << PAR_ANY("DPDF.Flux_reg1_b0")  << endl;
   exit(0);
   */

   file->Write();

   return true;
}


//____________________________________________________________________________________ //
