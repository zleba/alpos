#include "alpos/functions/ABelleTDCPV.h"
#include <iostream>

using namespace std;

// __________________________________________________________________________________________ //
const std::vector<std::string> ABelleTDCPV::fRequirements = {"tau","dm","A","S","W","fsig"}; //< List of all AParm's which this function depends on
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

   info["Init"]<<"This is ABelleTDCPV::Init(). Please implement the Likelihood calculation here"<<endl;
   fValue.resize(1);
   fError.resize(1);
 
   return true;
}


// __________________________________________________________________________________________ //
bool ABelleTDCPV::Update() {
   //! Update! This method must include the calculation of the likelihood function for any data point

   info["Update"]<<"This is ABelleTDCPV::Init(). Please implement the Likelihood calculation here"<<endl;

   // ---- parameters, which are input to the
   //      calculation and can be fitted
   double tau = PAR(tau);
   double dm  = PAR(dm);
   double A   = PAR(A);
   double S   = PAR(S);
   double W   = PAR(W);
   double fsig= PAR(fsig);

   // --- get info from 'data file'
   vector<string> datafiles = STRING_COL_NS(Data,datafile,GetAlposName());
   vector<double> sigmaData = DOUBLE_COL_NS(Data,Sigma,GetAlposName());
   string datadir = STRING_NS(DataDir,GetAlposName());

   // --- read data
   if ( sigmaData.size() && sigmaData[0]==0 ) {
      // read file from colum datafiles
      for ( const string& datafile : datafiles ) {
         cout<<"Opening data file: "<<datafile<<endl;
         // --- open data file
         TFile* f = NULL;

         // --- calculate likelihood for every event
         double L = 0;
      }
   }

   // fValue must contain a number of probabilites,
   // which are then input to the -2log(L) function.
   fValue[0] = exp(-0.5*pow((A-S)/dm,2)) / (sqrt(2*M_PI)*dm);
   fError[0] = 0;

   return true;
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


