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
   "MCFile","MCFileTreename",
   "DropDataOverflowBins", // drop overflow bins of data 0=false, 1=true;
   "nDataMax" // maximum number of data events
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

   vector<TFile*> DataFiles;
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
      DataFiles.push_back(DtFile);
   }
   


   // --- fill data histogram
   int nDataMax = PAR(nDataMax);
   const string DataTreename = STRING_NS(TTreename,GetAlposName());
   if ( DataTreename.empty() ) {
      error["Update"]<<"No TTreename specified in data file: "<<GetAlposName()<<endl;
      exit(1);
   }
   for ( TFile* fDataFile : DataFiles ) {
      TTreeReader DataReader(DataTreename.c_str(), fDataFile);
      TTreeReaderValue<Int_t>   nB0_mcTagPDG(DataReader, "B0_mcTagPDG");
      TTreeReaderValue<Float_t> tB0_DeltaT(DataReader, "B0_DeltaT");
      TTreeReaderValue<Float_t> tB0_DeltaTErr(DataReader, "B0_DeltaTErr");
      TTreeReaderValue<Int_t>   nB0_mcPDG(DataReader, "B0_mcPDG");
      TTreeReaderValue<Int_t>   nB0_mcErrors(DataReader, "B0_mcErrors");
      TTreeReaderValue<Int_t>   nB0_Jpsi_mu0_nPXDHits(DataReader, "B0_Jpsi_mu0_nPXDHits");
      TTreeReaderValue<Int_t>   nB0_Jpsi_mu1_nPXDHits(DataReader, "B0_Jpsi_mu1_nPXDHits");
      TTreeReaderValue<Float_t> tB0_TagVPvalue(DataReader, "B0_TagVPvalue");
      TTreeReaderValue<Float_t> tB0_TagVz(DataReader, "B0_TagVz");
      TTreeReaderValue<Float_t> tB0_TagVez(DataReader, "B0_TagVez");

      info["Init"]<<"Reading "<<DataReader.GetEntries(false)<<" data events..."<<endl;
      int nev=0; // counter
      int ngoodev=0; // counter
      float minddt =  10000;
      float maxddt = -10000;
      float mindt  =  10000;
      float maxdt  = -10000;
      while (DataReader.Next()) {
         nev++;
         float  dt     = (float)*tB0_DeltaT;
         float  ddt    = (float)*tB0_DeltaTErr;
         // todo: this is still qmc !
         int    qdt    = ((int)*nB0_mcTagPDG<0) ? -1 : 1;   // tag B0    =  511, or B0bar = -511 
         
         if( (float)*tB0_TagVPvalue  < 0.0  ) continue; // crazy event
         if( !(abs((int)*nB0_mcPDG)==511 && (int)*nB0_mcErrors<2 )) continue;
         if( !((int)*nB0_Jpsi_mu0_nPXDHits>0.9 || (int)*nB0_Jpsi_mu1_nPXDHits>0.9)) continue;
         bool badVtx = (float)*tB0_DeltaTErr > 2.75 && abs((float)*tB0_TagVz)<0.01 && (float)*tB0_TagVPvalue < 0.066  ;
         if( badVtx == kTRUE ) continue;

         // todo: read variables for 'classes' from data root file
         double dtpar = 0;
         int rqi      = 0;
         int ftag     = 0;
         int PDGtag   = (int)*nB0_mcTagPDG;
         const int ic = GetClassId(dtpar,rqi,ftag,PDGtag);
         
         // in case: init histograms
         if ( fH2Dt_dtddt.count(ic) == 0 ) {
            info["Init"]<<"Initializing histograms for class ID: "<<ic<<endl;
            fH2Dt_dtddt[ic]   = MakedtddtTH2D(Form("fH2Dt_dtddt_%03d",ic));
            fH2MC_dtddt_p[ic] = MakedtddtTH2D(Form("fH2MC_dtddt_p_%03d",ic));
            fH2MC_dtddt_m[ic] = MakedtddtTH2D(Form("fH2MC_dtddt_m_%03d",ic));
         }
         
         double evwgt = 1;
         fH2Dt_dtddt[ic].Fill(dt,ddt,evwgt);
         maxdt  = max(maxdt,dt);
         mindt  = min(maxdt,dt);
         maxddt = max(maxddt,ddt);
         minddt = min(maxddt,ddt);
         ngoodev++;
         if ( nDataMax>0 && ngoodev>=nDataMax ) break;
      }
      info["Update"]<<"Done. Read "<<nev<<" MC events, with "<<ngoodev<<" good events."<<endl;
      info["Update"]<<"      min(dt)="<<mindt<<"\tmax(dt)="<<maxdt<<"\tmin(ddt)="<<minddt<<"\tmax(ddt)="<<maxddt<<endl;
      fDataFile->Close();
   }
   // ---- sanity checks, and warning messages...
   for ( auto hh : fH2Dt_dtddt ) {
      TH1D* pYDt   = hh.second.ProjectionY();
      if ( pYDt->GetBinContent(0) != 0 )
         warn["Init"]<<"Non-zero bin content ("<<pYDt->GetBinContent(0)
                     <<")in ddt underflow bin for data and class="<<hh.first<<endl;
      if ( pYDt->GetBinContent(pYDt->GetNbinsX()+1) != 0 )
         warn["Init"]<<"Non-zero bin content ("<<pYDt->GetBinContent(pYDt->GetNbinsX()+1)
                     <<") in ddt overflow bin for data and class="<<hh.first<<endl;
      TH1D* pXDt   = hh.second.ProjectionX();
      if ( pXDt->GetBinContent(0) != 0 )
         warn["Init"]<<"Non-zero bin content ("<<pXDt->GetBinContent(0)
                     <<") in dt underflow bin for data and class="<<hh.first<<endl;
      if ( pXDt->GetBinContent(pXDt->GetNbinsX()+1) != 0 )
         warn["Init"]<<"Non-zero bin content ("<<pXDt->GetBinContent(pXDt->GetNbinsX()+1)
                     <<") in dt overflow bin for data and class="<<hh.first<<endl;
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
TH2D ABelleTDCPV::MakedtddtTH2D(std::string name){
   static       int    dt_nBins    = INT_NS(dt_nBins    ,GetAlposName());
   static const double dt_min      = DOUBLE_NS(dt_min   ,GetAlposName());
   static const double dt_max      = DOUBLE_NS(dt_max   ,GetAlposName());
   static       int    ddt_nBins   = INT_NS(ddt_nBins   ,GetAlposName());
   static const double ddt_min     = DOUBLE_NS(ddt_min  ,GetAlposName());
   static const double ddt_max     = DOUBLE_NS(ddt_max  ,GetAlposName());
   static vector<double> dt_binning = EXIST_NS(dt_binning,GetAlposName()) ?
      DOUBLE_ARR_NS(dt_binning,GetAlposName()) :
      vector<double>();
   static vector<double> ddt_binning = EXIST_NS(ddt_binning,GetAlposName()) ?
      DOUBLE_ARR_NS(ddt_binning,GetAlposName()) :
      vector<double>();
   if ( dt_binning.size()  ) dt_nBins  = dt_binning.size() - 1;
   if ( ddt_binning.size() ) ddt_nBins = ddt_binning.size() - 1;
   info["MakedtddtTH2D"]<<"Making histogram "<<name<<" with "
                        <<dt_nBins<< " dt bins, and "
                        <<ddt_nBins<< " ddt bins."<<endl;
   if ( dt_binning.size() && ddt_binning.empty() )
      return TH2D(name.c_str(), name.c_str(), dt_nBins, &dt_binning[0], ddt_nBins, ddt_min, ddt_max);
   else if ( dt_binning.empty() && ddt_binning.size() )
      return TH2D(name.c_str(), name.c_str(), dt_nBins, dt_min, dt_max, ddt_nBins, &ddt_binning[0]);
   else if ( dt_binning.size() && ddt_binning.size() )
      return TH2D(name.c_str(), name.c_str(),  dt_nBins, &dt_binning[0], ddt_nBins,  &ddt_binning[0]);
   else
      return TH2D(name.c_str(), name.c_str(), dt_nBins, dt_min, dt_max, ddt_nBins, ddt_min, ddt_max);
}


// __________________________________________________________________________________________ //
int  ABelleTDCPV::GetClassId(double dtpar, int rqi, int ftag, int PDGtag) {
   // todo: return some ID for any class.
   // remark: the ID's do not need to be subsequent.
   int clid = 0;
   if ( PDGtag<0) clid=0;
   else           clid=1;
   return clid;
}



// __________________________________________________________________________________________ //
bool ABelleTDCPV::Update() {
   //! Update! This method must include the calculation of the likelihood function for any data point

   debug["Update"]<<"This is ABelleTDCPV::Init(). Please implement the Likelihood calculation here"<<endl;
   
   // ---- parameters, which are input to the
   //      calculation and can be fitted
   const double tau  = PAR(tau);
   const double dm   = PAR(dm);
   const double A    = PAR(A);
   const double S    = PAR(S);
   const double tau0 = PAR(tau0);
   const double dm0  = PAR(dm0);
   const double A0   = PAR(A0);
   const double S0   = PAR(S0);
   const double W    = PAR(W);
   const double fsig = PAR(fsig);

   // ---- clearing MC histos
   for ( auto& hist : fH2MC_dtddt_p ) hist.second.Reset();
   for ( auto& hist : fH2MC_dtddt_m ) hist.second.Reset();


   // ---- calculate histograms with probabilies
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
      double ddt    = (float)*tB0_DeltaTErr;
      double dtmc   = (float)*tB0_TruthDeltaT;

      double dtpar = 0;
      int    rqi = 0;
      int    ftag = 0;
      int    PDGtag = (int)*nB0_mcTagPDG;
      const int ic = GetClassId(dtpar,rqi,ftag,PDGtag);
      
      double pdtFuncGen = Pdt(dtmc, qmc, tau0, dm0, A0, S0);
      double pdtFuncNew = Pdt(dtmc, qmc, tau , dm , A,  S );

      double mcwgt    = 1.;
      if ( pdtFuncGen != 0) mcwgt *= pdtFuncNew / pdtFuncGen ;
      else warn["Update"]<<"pdtFuncGen is zero. Skipping reweighting."<<endl;

      if (qmc>0) fH2MC_dtddt_p[ic].Fill(dt,ddt,mcwgt);
      else       fH2MC_dtddt_m[ic].Fill(dt,ddt,mcwgt);
      
      ngoodev++;      
   }
   info["Update"]<<"Done. Read "<<nev<<" MC events, with "<<ngoodev<<" good events."<<endl;


   // calculate probabilities
   double logL_p = 0;
   double logL_m = 0;
   int xymin = PAR(DropDataOverflowBins)==1 ? 1 : 0;
   for ( auto pH2 : fH2Dt_dtddt ) {
      int ic = pH2.first;
      const TH2D& h2dt = pH2.second;
      int xmax = h2dt.GetNbinsX()+1;
      int ymax = h2dt.GetNbinsY()+1;
      if ( PAR(DropDataOverflowBins)==1 ) {
         ymax-=1;
         xmax-=1;
      }

      TH1D* pYMC_p = fH2MC_dtddt_p[ic].ProjectionY();
      TH1D* pYMC_m = fH2MC_dtddt_m[ic].ProjectionY();
      TH1D* pYDt   = h2dt.ProjectionY();
      //cout<<"nBins pY: "<<pYMC_p->GetNbinsX()<<endl;
      for ( int iy = xymin; iy<=ymax; iy++ ) {
         if ( pYDt->GetBinContent(iy) == 0 ) continue;// no data
         double sumDtMC_p = pYMC_p->GetBinContent(iy);
         double sumDtMC_m = pYMC_m->GetBinContent(iy);
         // if ( sumDtMC_p == 0 ) ...
         // if ( sumDtMC_m == 0 ) ...
         for ( int ix = xymin; ix<=xmax; ix++ ) {
            const double nDt = h2dt.GetBinContent(ix,iy);
            if ( nDt == 0 ) continue;// no data
            if ( sumDtMC_p !=0 ) {
               double prob_p = fH2MC_dtddt_p[ic].GetBinContent(ix,iy) / sumDtMC_p;
               if ( prob_p == 0 ) warn["Update"]<<"MC probability is zero for non-zero data event."<<endl;
               else logL_p += nDt * log(prob_p);
            }
            if ( sumDtMC_m !=0 ) {
               double prob_m = fH2MC_dtddt_m[ic].GetBinContent(ix,iy) / sumDtMC_m;
               if ( prob_m == 0 ) warn["Update"]<<"MC probability is zero for non-zero data event."<<endl;
               else logL_m += nDt * log(prob_m);
            }
            // info["Update"]<<"ix="<<ix<<"\tiy="<<iy
            //                <<"\tprob_p="<<prob_p<<"\tprob_m="<<prob_m
            //                <<"\tsum_p="<<pYMC_p->GetBinContent(iy)
            //                <<"\tsum_m="<<pYMC_m->GetBinContent(iy)
            //                <<endl;
         }
      }
   }
   
   // double L_p = exp(logL_p);
   // double L_m = exp(logL_m);
   // double pdfSig = (1-W)*L_p + W*L_m;
   // double model = fsig*pdfSig + (1-fsig)*0;
   // info["Update"]<<"logL_p="<<logL_p
   //               <<"\tlogL_m="<<logL_m
   //               <<"\tL_m="<<L_m
   //               <<"\tL_m="<<L_p
   //               <<"\tpdfSig="<<pdfSig
   //               <<"\tmodel="<<model
   //               <<"\t-2log(L)="<<-2*log(model)
   //               <<endl;

   if ( W>=1 || W <=0 ) {
      error["Update"]<<"W exceeds limits ]0;1[. W="<<W<<endl;
      exit(1);
   }
   double pdfSig_p = logL_p + log(1-W);
   double pdfSig_m = logL_m + log(W);
   double pdfSig   = pdfSig_p + pdfSig_m; // ?? correct ??
   // if ( fsig>=1 || fsig <=0 ) {
   //    error["Update"]<<"fsig exceeds limits ]0;1[. W="<<W<<endl;
   //    exit(1);
   // }
   //double model = (log(fsig) + pdfSig) * (log(1-fsig) + 0);
   double model = -2*pdfSig;

   info["Update"]<<"pdfSig_p="<<pdfSig_p
                 <<"\tpdfSig_m="<<pdfSig_m
                 <<"\tpdfSig="<<pdfSig
                 <<"\tmodel="<<model
                 <<endl;

      
   
   // fValue must contain a number of probabilites,
   // which are then input to the -2log(L) function.
   // double dS=0.01;
   // fValue[0] = exp(-0.5*pow((A-A0)/dm,2)) / (sqrt(2*M_PI)*dm) * exp(-0.5*pow((S-S0)/dS,2)) / (sqrt(2*M_PI)*dS);
   // fError[0] = 0;
   
   fValue[0]= model;

   return true;
}


// ______________________________________________________________________________________ //
double ABelleTDCPV::Pdt(double dt, double q, double tau, double dm, double A, double S){
   double Pdt = exp(-abs(dt)/tau)/(4*tau); 
   Pdt *= ( 1 + q*(A*cos(dm*(dt))+S*sin(dm*(dt))) );
   if ( Pdt <= 0 )warn["Pdt"]<<"Weigthing function less than zero: "<<Pdt<<" Pdt(dt="<<dt
                             <<",q="<<q
                             <<",tau="<<tau
                             <<",dm="<<dm
                             <<",A="<<A
                             <<",S="<<S<<endl;
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


