// DB 15.01.2015

#include "alpos/functions/AfastNLODiffDIS.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include <alpos/AError.h>

#include <functional>
#include <numeric>
#include "TH1D.h"
#include "TVectorD.h"

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLODiffDIS::fRequirements = {//"Filename",                      // fastNLO table filenmae
   "DPDF",                           // a 'PDF' function with Quick-function. Moslty LHAPDF6
   "Alpha_s",                       // a alpha_s(mu_r) function.
   "ScaleFacMuR","ScaleFacMuF",     // Scale factors for ren. and fact. scale
   // "Units",                         // publication or absolute units (to obtain same units as your data table)
   "iOrd",                          // order
   "MuRFuncForm","MuFFuncForm",      // mu_r and mu_f functional form for fastNLO flexible scale tables
   "xpom_min","xpom_max","nxpom","logxpom" // xpom slicing
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLODiffDIS::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLODiffDIS::fFunctionName = "fastNLODiffDIS"; //< The function's name


vector<double> CalculateSpecialVar(fastNLOAlposDPDF  &fnlodiff,
                                   std::function<double(double,double,double)>  func);//xpom,q2,tableVar

vector<double> CalculateSpecialXpom(fastNLOAlposDPDF  &fnlodiff);

// __________________________________________________________________________________________ //
AfastNLODiffDIS::AfastNLODiffDIS(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AfastNLODiffDIS");
}


// __________________________________________________________________________________________ //
AfastNLODiffDIS::~AfastNLODiffDIS() {
   //if ( fnlo ) delete fnlo;
   for ( auto i : fnlos ) delete i;
}


// ___________________________________________________________________________________________ //
bool AfastNLODiffDIS::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Requires in datafile defintions  of
   //!  'Filename' (one entry), 'Filenames' (array) or 'Tables' (table)
   //! and further 'Units' and if applicable: 'MuRFuncForm' and 'MuFFuncForm'
   //!

   using namespace AlposTools;

   double Units = DOUBLE_NS(Units,GetAlposName());

   vector<string> filenames;
   vector<int> firstbins, lastbins;
   if ( EXIST_NS(Filename,GetAlposName() ))
      filenames.push_back(STRING_NS(Filename,GetAlposName()));
   else if ( EXIST_NS(Filenames,GetAlposName() ) )
      filenames = STRING_ARR_NS(Filenames,GetAlposName()); // direct access to array
   else if ( EXIST_NS(Tables,GetAlposName()) ) {
      filenames = STRING_COL_NS(FnloTables,Filenames,GetAlposName());
      firstbins = INT_COL_NS(FnloTables,FirstBin,GetAlposName());
      lastbins  = INT_COL_NS(FnloTables,LastBin,GetAlposName());
   }


   CONST(xpom_min);
   CONST(xpom_max);
   CONST(nxpom);
   CONST(logxpom);

   //string filename = PAR_S(Filename);
   //cout<<"fastnlo: Filename: "<<filename<<endl;
   for (auto& fn : filenames) {
      AlposTools::CheckFileExit(fn);
      info["Init"] << "Reading table file: " << fn << endl;
      //fastNLODiffH12006FitB* f = new fastNLODiffH12006FitB(fn);
      fastNLOAlposDPDF* f = new fastNLOAlposDPDF(fn);
      fnlos.push_back(f);
      cout << "Radek " << this->GetAlposName() << endl;
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
      //fnloreaders[i]->SetXPomLinSlicing( INT(NumberOfXPomSlices), 0. ,  .03 ); //12.  
      if ( PAR(logxpom) == 0 ) 
         f->SetXPomLinSlicing( PAR(nxpom),PAR(xpom_min),PAR(xpom_max) ); //12.  
      else
         f->SetXPomLogSlicing( PAR(nxpom),PAR(xpom_min),PAR(xpom_max) ); //12.  
      f->SetTIntegratedRange(-1.);
      f->SetProtonE(920.);

      //int FitID = BOOL(DoFitB) ? 2 : 1;
      //f->SetFitID(2);
      if (f->GetIsFlexibleScaleTable()) {
         // f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         // f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
      }
   }

   // init binmap if specified
   if (!firstbins.empty()) {
      for (unsigned int ig = 0; ig < filenames.size(); ig++) {
         vector<bool> bmap(fnlos[ig]->GetNObsBin());
         for (unsigned int ib = 0; ib < bmap.size(); ib++) {
            bmap[ib] = (ib >= firstbins[ig] && ib <= lastbins[ig]);
         }
         fBinmap += bmap;
      }
   }

   /*
   //if ( fnlo ) delete fnlo;
   fnlo = new fastNLOAlpos(filenames[0]);
   fnlo->SetAlposName(this->GetAlposName());
   //Things which do not make sense to change later
   //fnlo->SetFilename(PAR_S(Filename)); // filename should not be changed
   fnlo->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
   if (fnlo->GetIsFlexibleScaleTable()) {
   fnlo->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
   fnlo->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
   }
   CONST(MuRFuncForm);
   CONST(MuFFuncForm);
   CONST(Units);
   CONST(Filename);
   */
   CONST(iOrd);

   return true;
}

bool AfastNLODiffDIS::CalcCrossSections() {

   using namespace AlposTools;

   const auto MX = [] (double xpom, double Q2, double y) {
      const double s = 4*27.6*920.;
      return  sqrt( max(0.0,y * s * xpom - Q2)); };

   const auto ZP = [] (double xpom, double Q2, double xi12) {
      return xi12/xpom; };

   const auto BETA= [] (double xpom, double Q2, double xBjor) {
      return xBjor/xpom; };


   // get cross sections
   fValue.clear();
   cout << "start" << endl;


   vector<vector<double>> vals(fnlos.size());

//#pragma omp parallel for
   for(int i = 0; i < fnlos.size(); ++i) {
      auto fnlo = fnlos[i];
   //for ( auto fnlo : fnlos ) {
      // fnlo->CalcCrossSection();
      // fValue += fnlo->GetCrossSection();

      string aname = fnlo->GetAlposName();
      string fname = fnlo->GetFilename();
      cout << "Radek " << fname << endl;
      if(fname.find("xi2zIP") != string::npos) { //isZpom
         vals[i] =  CalculateSpecialVar(*fnlo, ZP);
      }
      else if(fname.find("yMx") != string::npos) { //isMx
         vals[i] =  CalculateSpecialVar(*fnlo, MX); 
      }
      else if(fname.find("xbjbeta") != string::npos) { //isBeta
         vals[i] =  CalculateSpecialVar(*fnlo, BETA); 
      }
      else if(aname.find("xpom") != string::npos) {
         vals[i] = CalculateSpecialXpom(*fnlo);
      }

      else {
         vals[i] = fnlo->GetDiffCrossSection();
      }
   }
   //Merge outputs
   for(auto v : vals)
      fValue += v;


   //cout << endl;
   // apply binmap if applicable
   if ( !fBinmap.empty() ) {
      int ii=0;
      for ( int ib = 0 ; ib<fValue.size() ; ib++ ){
         if ( fBinmap[ib] ) fValue[ii++]=fValue[ib];
      }
      fValue.resize(ii);
   }

   fError.resize(fValue.size());

   return true;
}

// __________________________________________________________________________________________ //
bool AfastNLODiffDIS::Update() {

   using namespace AlposTools;

   if ( fValue.empty() || (fValue.size()==1 && fValue[0]==0 )) SetOrder(); // must be initialized in 'update' only

   if ( CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF) )
      for ( auto fnlo : fnlos ) {
         debug["Update"] << "Setting ScaleFacMuR, ScaleFacMuF to (" << PAR(ScaleFacMuR) << ", " << PAR(ScaleFacMuF) << ")..." << std::endl;
         fnlo->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(DPDF);
   UPDATE(Alpha_s);

   CalcCrossSections();

   fHasMultErrors = false;

   return true;
}


// __________________________________________________________________________________________ //
void AfastNLODiffDIS::SetOrder() {
   //! Set correct order of fastNLO calculation
   for ( auto fnlo : fnlos ) {
      // if ( CHECK(iOrd) ) 
      //! Check on existence of various pQCD contributions in table (Id = -1 if not existing)
      //! Check on existence of LO (Id = -1 if not existing)
      int ilo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
      if (ilo < 0) {
         error["SetOrder"] << "LO not found, nothing to be done!" << endl;
         exit(1);
      } else {
         info["SetOrder"] << "The LO contribution has Id: " << ilo << endl;
      }
      //! Check on existence of NLO (Id = -1 if not existing)
      int inlo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
      if (inlo < 0) {
         info["SetOrder"] << "No NLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NLO contribution has Id: " << inlo << endl;
      }
      //! Check on existence of NNLO (Id = -1 if not existing)
      int innlo = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToNextToLeading);
      if (innlo < 0) {
         info["SetOrder"] << "No NNLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NNLO contribution has Id: " << innlo << endl;
      }
      //! Switch selected contributions ON, if possible
      //! LO & NLO are ON by default
      //! Fixed-order
      int iOrd = PAR(iOrd);
      if ( iOrd==0 ) {
         if ( !(inlo<0) )  {fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, false);}
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==1 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==2 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, true);
         } else {
            error["SetOrder"] << "NNLO requested, but not found. Ignoring call!" << endl;
            //exit(1);
         }
      }
      //      
      //   
   }
}


TMatrixD CreateInterpolationMatrix(vector<double> &fastBinsLo, vector<double> &fastBinsHi)
{
    const int Nb = fastBinsLo.size();

    assert(Nb >2);

    vector<double> wb(Nb);
    for(int i = 0; i < Nb; ++i)
        wb[i] = fastBinsHi[i] - fastBinsLo[i];

    TMatrixD mat(Nb,Nb);

    mat(0,0)     = 1;
    mat(Nb-1,Nb-1) = 1;

    for(int i = 1; i < Nb-1; ++i) {
        mat(i,i-1) = 1./4 * wb[i] / ( wb[i-1] + wb[i] );
        mat(i,i+1) = 1./4 * wb[i] / ( wb[i+1] + wb[i] );
        mat(i,i  ) = 1./4 * ( 2 + wb[i-1] / (wb[i-1]+wb[i]) + wb[i+1] / (wb[i+1]+wb[i]) );
    }

    return mat.Invert();

}






vector<double> CalculateSpecialVar(fastNLOAlposDPDF  &fnlodiff,
                                   std::function<double(double,double,double)>  func) //xpom,q2,tableVar
{

    //  If you want to receive your cross section in
    //   pb/GeV or in pb. Here we choose pb/GeV
    fnlodiff.SetUnits(fastNLO::kPublicationUnits);


    //fnlodiff.SetExternalFuncForMuR (&Function_Mu);
    //fnlodiff.SetExternalFuncForMuF (&Function_Mu);



    vector<double> fastBinsLo = fnlodiff.GetObsBinsLoBounds(0);
    vector<double> fastBinsHi = fnlodiff.GetObsBinsUpBounds(0);


    // calculate and access the cross section
    typedef std::map<double,  std::vector < std::map< double, double > > >  array3D;

    array3D xsQ2;
    //std::map<double,  std::vector < std::map< double, double > > > xsQ2 = new array3D[3+npdfall];


   string aname = fnlodiff.GetAlposName();
   string fname = fnlodiff.GetFilename();

    //TODO Need to get data binning to xMin, xMax
   vector<double> xMin; //= {0., 0.3, 0.6, 0.9};
   vector<double> xMax; //= {0.3, 0.6, 0.9, 1.0};

   if(fname.find("xi2zIP") != string::npos) { //isZpom
      xMin = DOUBLE_COL_NS(Data,zpom_min,aname);
      xMax = DOUBLE_COL_NS(Data,zpom_max,aname);
   }
   else if(fname.find("yMx") != string::npos) { //isMx
      xMin = DOUBLE_COL_NS(Data,mx_min,aname);
      xMax = DOUBLE_COL_NS(Data,mx_max,aname);
   }
   else if(fname.find("xbjbeta") != string::npos) { //isBeta
      xMin = DOUBLE_COL_NS(Data,beta_min,aname);
      xMax = DOUBLE_COL_NS(Data,beta_max,aname);
   }
   else {
      cout << "Unknow type" << endl;
      exit(1);
   }


   int nBins = xMin.size();

    vector<double> bins(nBins+1);
    for(int i = 0; i < nBins; ++i)
        bins[i] = xMin[i];
    bins[nBins] = xMax[nBins-1];

    int hash = rand();

    TH1::AddDirectory(false);
    TH1D hist("histogram", "histogram", nBins, bins.data());

    //fnlodiff.SetLHAPDFMember(0);

    //fnlodiff.SetScaleFactorsMuRMuF(1.0, 1.0);
    vector<double> NloXs=fnlodiff.GetDiffCrossSection();
    xsQ2 = fnlodiff.Get3DCrossSection();

    /*
    double Sum=0;
    for(int i = 0; i <NloXs.size(); ++i) {
        Sum+= NloXs[i] * (fastBinsHi[i]-fastBinsLo[i]);
        cout << i <<" : "<<  fastBinsLo[i]<<" "<<fastBinsHi[i]<< "  <>  "<<   NloXs[i] << endl;
    }
    cout << "Sum is " << Sum << endl;
    exit(0);
    */


    //For interpretation
    TMatrixD mat = CreateInterpolationMatrix(fastBinsLo, fastBinsHi);
    TVectorD nloVec(fastBinsLo.size() );


    //Fill array of Q2
    vector<double> Q2arr;
    for(const auto &a :   (xsQ2.begin()->second)[0])
        Q2arr.push_back(a.first);


    const int Niter = 400;
    srand(1);


    for( const auto &v : xsQ2 ) {
        double xpom = v.first;
        auto xs2D   = v.second;

        for(unsigned i = 0; i < fastBinsLo.size(); ++i) {
            for(const auto &item : xs2D[i]) {
                double Q2 = item.first;
                double xs = item.second*(fastBinsHi[i]-fastBinsLo[i]) / Niter;
                for(int k = 0; k < Niter; ++k) {
                    double varHist = fastBinsLo[i] + rand()/(RAND_MAX+0.) * (fastBinsHi[i]-fastBinsLo[i]);
                    double VarCalc = func(xpom, Q2, varHist);
                    //cout << "Zpom " << VarCalc <<  endl; 
                    hist.Fill(VarCalc, xs);
                }
            }
        }
    }

    vector<double> xsc;
    xsc.resize(nBins);

    //double Units = DOUBLE_NS(Units,GetAlposName());

    for(int i = 0; i < nBins; ++i) {
        double corr = 1; // 1./ (xMax[i]-xMin[i]); Todiff
        xsc[i]=hist.GetBinContent(i+1) * corr;
    }


    return xsc;
}


vector<double> CalculateSpecialXpom(fastNLOAlposDPDF  &fnlodiff)
{
    string name = "xpom";

    vector<double> logxpomLRG = {-2.30, -2.10, -1.90, -1.70, -1.52};
    vector<double> logxpomFPS = {-2.3, -1.9, -1.6, -1.4, -1.2, -1.0};
    vector<double> xpomVFPS   = {0.010,  0.014, 0.019, 0.024};
    vector<double> xpomZEUS   = {0.0025, 0.0050, 0.0079, 0.0126, 0.0199, 0.0300};
   //xpom binning
    vector<double> xMin = logxpomLRG;
    xMin.erase(xMin.end()-1); //remove last
    vector<double> xMax = logxpomLRG;
    xMax.erase(xMax.begin()); //remove first

    int nBins = xMin.size();
    vector<double>  nloXpom(nBins);

    //cout << "RADEKHERE " << __LINE__ << endl;


    for(int i = 0; i < nBins; ++i) {
    //cout << "RADEKHERE " << __LINE__ << endl;

        if(name.find("xpom") != string::npos) {
            fnlodiff.SetXPomLinSlicing( 5, xMin[i], xMax[i]); // VFPS range
            fnlodiff.SetTIntegratedRange(-0.6);
        }
        else if(name.find("logxpom") != string::npos) {
            fnlodiff.SetXPomLogSlicing( 5, pow(10,xMin[i]) ,  pow(10,xMax[i]) ); 
            fnlodiff.SetTIntegratedRange(-1.);
    //cout << "RADEKHERE " << __LINE__ << endl;
            
        }
        else {
            cout << "error in xpom " << name <<  endl;
        }

        vector<double> Xs = fnlodiff.GetDiffCrossSection();
        //assume absolute units 
        double thTot = accumulate(Xs.begin(),Xs.end(), 0.);

        double bw = 1;//xMax[i] - xMin[i]; Normalization

        nloXpom[i]   = thTot / bw;

    }

    return nloXpom;

}
