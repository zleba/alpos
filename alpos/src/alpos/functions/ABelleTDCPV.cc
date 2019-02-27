#include "alpos/functions/ABelleTDCPV.h"
#include <iostream>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

using namespace std;

// __________________________________________________________________________________________ //
const std::vector<std::string> ABelleTDCPV::fRequirements = {
   "tau","dm","A","S",
   "tau0","dm0","A0","S0",
   "W","fsig",
   "MCFile","MCFileTreename"
}; //< List of all AParm's which this function depends on
const std::vector<std::string> ABelleTDCPV::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ABelleTDCPV::fFunctionName = "BelleTDCPV"; //< The function's name


// __________________________________________________________________________________________ //
ABelleTDCPV::ABelleTDCPV(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("ABelleTDCPV");
}


// __________________________________________________________________________________________ //
ABelleTDCPV::~ABelleTDCPV() {
}


// ___________________________________________________________________________________________ //
bool ABelleTDCPV::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   // 'real' QCDNUM init is called in constructor

   debug["Init"]<<"This is ABelleTDCPV::Init(). Please implement the Likelihood calculation here"<<endl;
   fValue.resize(1);
   fError.resize(1);

   // generated values should not be changed
   CONST(tau0);
   CONST(dm0);
   CONST(A0);
   CONST(S0);

   // --- get info from 'data file'
   vector<string> datafiles = STRING_COL_NS(Data,datafile,GetAlposName());
   vector<double> sigmaData = DOUBLE_COL_NS(Data,Sigma,GetAlposName());
   string datadir;
   if ( EXIST_NS(DataDir,GetAlposName()) ) datadir = STRING_NS(DataDir,GetAlposName())+"/";
   
   // --- read data
   if ( sigmaData.size()!=1 || sigmaData[0]!=0 ) {
      error["Update"]<<"Currently, input .dat-file must specify one file in Data table. Like:"<<endl;
      error["Update"]<<"   Sigma         datafile       stat"<<endl;
      error["Update"]<<"     0           test.root        0"<<endl;
      exit(1);
   }

   for ( string datafile : datafiles ) {
      datafile=datadir+datafile;
      AlposTools::CheckFileExit(datafile);
      TFile* DtFile = TFile::Open(datafile.c_str(),"READ");
      if (!DtFile || DtFile->IsZombie()) {
         error["Init"]<<"Cannot read data root file: "<<datafile<<endl;
         exit(1);
      }
      else
         info["Init"]<<"Opened data file: "<<datafile<<endl;
      fDataFiles.push_back(DtFile);
   }
   

   // --- reading MC file
   string MCFilename = PAR_S(MCFile);
   AlposTools::CheckFileExit(MCFilename);

   // --- read MC
   TFile* MCFile = TFile::Open(MCFilename.c_str(),"READ");
   if (!MCFile || MCFile->IsZombie()) {
      error["Init"]<<"Cannot read MC root file: "<<MCFilename<<endl;
      exit(1);
   }
   else
      info["Init"]<<"Opened MC file: "<<MCFilename<<endl;

   fMcFiles.push_back(MCFile);

   return true;
}


// __________________________________________________________________________________________ //
bool ABelleTDCPV::Update() {
   //! Update! This method must include the calculation of the likelihood function for any data point

   debug["Update"]<<"This is ABelleTDCPV::Init(). Please implement the Likelihood calculation here"<<endl;
   
   // ---- parameters, which are input to the
   //      calculation and can be fitted
   double tau  = PAR(tau);
   double dm   = PAR(dm);
   double A    = PAR(A);
   double S    = PAR(S);
   double tau0 = PAR(tau0);
   double dm0  = PAR(dm0);
   double A0   = PAR(A0);
   double S0   = PAR(S0);
   double W    = PAR(W);
   double fsig = PAR(fsig);



   const string DataTreename = STRING_NS(TTreename,GetAlposName());
   if ( DataTreename.empty() ) {
      error["Update"]<<"No TTreename specified in data file: "<<GetAlposName()<<endl;
         exit(1);
   }
   for ( TFile* fDataFile : fDataFiles ) {
      TTreeReader DataReader(DataTreename.c_str(), fDataFile);
      TTreeReaderValue<Int_t>   dt_nB0_mcTagPDG(DataReader, "B0_mcTagPDG");
      TTreeReaderValue<Float_t> dt_tB0_DeltaT(DataReader, "B0_DeltaT");
      TTreeReaderValue<Float_t> dt_tB0_DeltaTErr(DataReader, "B0_DeltaTErr");
      TTreeReaderValue<Int_t>   dt_nB0_mcPDG(DataReader, "B0_mcPDG");
      TTreeReaderValue<Int_t>   dt_nB0_mcErrors(DataReader, "B0_mcErrors");
      TTreeReaderValue<Int_t>   dt_nB0_Jpsi_mu0_nPXDHits(DataReader, "B0_Jpsi_mu0_nPXDHits");
      TTreeReaderValue<Int_t>   dt_nB0_Jpsi_mu1_nPXDHits(DataReader, "B0_Jpsi_mu1_nPXDHits");
      TTreeReaderValue<Float_t> dt_tB0_TagVPvalue(DataReader, "B0_TagVPvalue");
      TTreeReaderValue<Float_t> dt_tB0_TagVz(DataReader, "B0_TagVz");
      TTreeReaderValue<Float_t> dt_tB0_TagVez(DataReader, "B0_TagVez");
   }

   
   if ( fMcFiles.size()!=1 ) warn["Update"]<<"Currently, only the evaluation of a single MC file is implemented."<<endl;
   TFile* MCFile = fMcFiles[0];
   const string MCTreename = PAR_S(MCFileTreename);
   TTreeReader MCReader(MCTreename.c_str(), MCFile);
   TTreeReaderValue<Float_t> tB0_DeltaT(MCReader, "B0_DeltaT"); // dt
   TTreeReaderValue<Float_t> tB0_DeltaTErr(MCReader, "B0_DeltaTErr"); //ddt
   TTreeReaderValue<Float_t> tB0_TruthDeltaT(MCReader, "B0_TruthDeltaT"); // dtmc
   TTreeReaderValue<Float_t> tB0_TagVPvalue(MCReader, "B0_TagVPvalue");
   TTreeReaderValue<Float_t> tB0_TagVz(MCReader, "B0_TagVz");
   TTreeReaderValue<Float_t> tB0_TruthTagVz(MCReader, "B0_TruthTagVz");
   TTreeReaderValue<Float_t> tB0_TagVez(MCReader, "B0_TagVez");
   TTreeReaderValue<Int_t>   nB0_mcTagPDG(MCReader, "B0_mcTagPDG"); // charge: PDG ~ qmc
   TTreeReaderValue<Int_t>   nB0_mcPDG(MCReader, "B0_mcPDG");
   TTreeReaderValue<Int_t>   nB0_mcErrors(MCReader, "B0_mcErrors");
   TTreeReaderValue<Int_t>   nB0_Jpsi_mu0_nPXDHits(MCReader, "B0_Jpsi_mu0_nPXDHits");
   TTreeReaderValue<Int_t>   nB0_Jpsi_mu1_nPXDHits(MCReader, "B0_Jpsi_mu1_nPXDHits");

   // --- read MC
   int nev=0; // counter
   int ngoodev=0; // counter
   info["Update"]<<"Reading "<<MCReader.GetEntries(false)<<" MC events..."<<endl;
   while (MCReader.Next()) {
      nev+=1;
      if( (float)*tB0_TagVPvalue  < 0.0  ) continue; // crazy event
      if( !(abs((int)*nB0_mcPDG)==511 && (int)*nB0_mcErrors<2 )) continue;
      if( !((int)*nB0_Jpsi_mu0_nPXDHits>0.9 || (int)*nB0_Jpsi_mu1_nPXDHits>0.9)) continue;
      bool badVtx = (float)*tB0_DeltaTErr > 2.75 && abs((float)*tB0_TagVz)<0.01 && (float)*tB0_TagVPvalue < 0.066  ;
      if( badVtx == kTRUE ) continue;
      
      int    qmc    = ((int)*nB0_mcTagPDG<0) ? -1 : 1;   // tag B0    =  511, or B0bar = -511 
      double dt     = (float)*tB0_DeltaT;
      int    dtsign = dt < 0. ? -1 : 1;
      
      double adt    = abs(dt);
      int    idt    = 0;
      if      ( adt <  1. ) idt =            adt/0.05;  //  0 - 19
      else if ( adt <  2. ) idt = 20 + (adt-1.0)/0.10;  // 20 - 29
      else if ( adt <  3. ) idt = 30 + (adt-2.0)/0.20;  // 30 - 34
      else if ( adt <  4. ) idt = 35 + (adt-3.0)/0.25;  // 35 - 38
      else if ( adt <  5. ) idt = 39 + (adt-4.0)/0.50;  // 39 - 40
      else if ( adt <  6. ) idt = 41 + (adt- 5.)/1.0;   // 41
      else if ( adt <  8. ) idt = 42 + (adt- 6.)/2.0;   // 42
      else if ( adt < 10. ) idt = 43 + (adt- 8.)/4.0;   // 43
      else if ( adt < 14. ) idt = 44 + (adt-10.)/8.0;   // 44
      else if ( adt < 22. ) idt = 45 + (adt-14.)/8.0;   // 45
      else if ( adt < 38. ) idt = 46 + (adt-22.)/8.0;   // 46
      else if ( adt < 70. ) idt = 47 + (adt-38.)/8.0;   // 47
      else                  idt = 48;                   // 48

      double pdtFuncGen = Pdt(dt, qmc, tau0, dm0, A0, S0);
      double pdtFuncNew = Pdt(dt, qmc, tau , dm , A,  S );
      double wdt    = 0.;
      if ( pdtFuncGen > 0) wdt = pdtFuncNew / pdtFuncGen ;
      double invwdt = 1.;
      if ( pdtFuncGen > 0.01) invwdt = 1. / pdtFuncGen ;

      Int_t ibin = idt*dtsign;
      if( dtsign==-1 ) ibin = ibin -1;
      ngoodev++;      
   }
   info["Update"]<<"Done. Read "<<nev<<" MC events, with "<<ngoodev<<" good events."<<endl;

   // fValue must contain a number of probabilites,
   // which are then input to the -2log(L) function.
   fValue[0] = exp(-0.5*pow((A0-S0)/dm,2)) / (sqrt(2*M_PI)*dm);
   fError[0] = 0;

   return true;
}


// ______________________________________________________________________________________ //
double ABelleTDCPV::Pdt(double dt, double q, double tau, double dm, double A, double S){
   double Pdt = exp(-abs(dt)/tau)/(4*tau); 
   Pdt *= ( 1 + q*(A*cos(dm*(dt))+S*sin(dm*(dt))) );
   return Pdt;
}



// ______________________________________________________________________________________ //
std::vector<double> ABelleTDCPV::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);
   
   std::vector<double> ret(1);
   //if ( mur.size() != 1 ) {
   {
      cout<<"Error in ABelleTDCPV::GetQuick(vector<double>). Quick acces is not implemented."<<endl;
      exit(1);
      return ret;
   }
   return ret;
}


