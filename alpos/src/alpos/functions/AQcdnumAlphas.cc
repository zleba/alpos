// DB 15.01.2015

#include "alpos/functions/AQcdnumAlphas.h"
#include <iostream>

using namespace std;

extern "C" {
   // call qcinit(6,' ');//        !initialize 
   // call setord(I_FIT_ORDER);//         !LO, NLO, NNLO       
   //void qcdnum_ini_();
   // void getintegrateddisxsection_(double*);
   // void setabr_(double*, double*);
   // void zmdefq2_(double*, double*);

   // void qcinit_(int*,char*);
   // void setalf_(double*,double*);
   // void getalf_(double*,double*);
   // void setcbt_(int* nfin,int* mc,int* mb,int* mt);
   double asfunc_(double* q2, int* nf, int* ierr);
   // double setord_(int* iord);
   // void gxmake_(double* xmin,int *iwt,int* nx,int* nin,int* nxout,int* ispline);//        !x-grid 
   // //void gqmake_(double* qarr,int* wgt,int* n,int* nqin,int* nqout);//          !mu2-grid  
   // void gqmake_(double* qarr,double* wgt,int* n,int* nqin,int* nqout);//          !mu2-grid  
   // int iqfrmq_(double* imc);// grid indices
}

// __________________________________________________________________________________________ //
const std::vector<std::string> AQcdnumAlphas::fRequirements = {"QcdnumInitializer","mur"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumAlphas::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumAlphas::fFunctionName = "QcdnumAlphas"; //< The function's name


// __________________________________________________________________________________________ //
AQcdnumAlphas::AQcdnumAlphas(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   fValue.resize(1);
   fError.resize(1);
   /*
   int lun = 6;// qcdnum output on standard screen
   char nix[] = " ";
   qcinit_(&lun,nix);
   // we must init x and muf grids to set nf in alpha_s 
   // --- these are just dummy values;
   vector<double> xmin = {0.00001,0.0001,0.001,0.01,0.1};
   vector<int>    iwgt = {1,1,1,1,1,1};
   vector<double> dwgt = {1,1,1,1,1,1};
   int nx = 3;
   int nin = 5;
   int nxout;
   int ispline = 2 ;
   gxmake_(&xmin[0],&iwgt[0],&nx,&nin,&nxout,&ispline);
   
   vector<double> qarr = {1,10,100,1000,10000,100000};
   gqmake_(&qarr[0],&dwgt[0],&nx,&nin,&nxout);//          !mu2-grid  
   */
}


// __________________________________________________________________________________________ //
AQcdnumAlphas::~AQcdnumAlphas() {
}


// ___________________________________________________________________________________________ //
bool AQcdnumAlphas::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   // 'real' QCDNUM init is called in constructor

   // update QCDNUM
   // it is not possible to call another function during the init step
   //int init = PAR(QcdnumInitializer);

   debug["Init"]<<endl;

   /*
   double mc2 = PAR(mcharm)*PAR(mcharm);
   double mb2 = PAR(mbottom)*PAR(mbottom);
   double mt2 = PAR(mtop)*PAR(mtop);
   int imc = iqfrmq_(&mc2);
   int imb = iqfrmq_(&mb2);
   int imt = iqfrmq_(&mt2);
   int nfin  = PAR(nfFix);
   setcbt_(&nfin,&imc,&imb,&imt);

   double Mz2 = PAR(Mz);
   int iOrd = PAR(iOrd)+1;
   Mz2*=Mz2;
   double AsMz = PAR(AlphasMz);
   setord_(&iOrd);
   setalf_(&AsMz,&Mz2);
   */   
   return true;
}


// __________________________________________________________________________________________ //
bool AQcdnumAlphas::Update() {

   debug["Update"]<<endl;

   //cout<<"AQcdnumAlphas.Update() Check(QcdnumInit): " <<CHECK(QcdnumInitializer)<<endl;
   if ( CHECK(QcdnumInitializer) ) 
      PAR(QcdnumInitializer);

   double Q2 = PAR(mur);
   Q2*=Q2;
   int ierr =0;
   int nf = 0;
   fValue[0] = asfunc_(&Q2,&nf,&ierr);
   fError[0] = 0;

   return true;
}

// ______________________________________________________________________________________ //
std::vector<double> AQcdnumAlphas::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);
   std::vector<double> ret(1);
   if ( mur.size() != 1 ) {
      cout<<"Error in AQcdnumAlphas::GetQuick(vector<double>). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   double Q2 = mur[0]*mur[0];
   int ierr =0;
   int nf = 0;
   ret[0] = asfunc_(&Q2,&nf,&ierr);
   return ret;
}


