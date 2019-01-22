// DB 22.09.2015

#include "alpos/functions/AEprc.h"
#include "alpos/ATheory.h"
#include <cmath>
#include <iostream>

using namespace std;

extern "C" {
   void eprc_init_alpos_full_(
			 int* LP1,int* LP2,int* LP3,int* LP4,int* LP5,
			 int* LP6,int* LP7,int* LP8,int* LP9,int* LP10,
			 int* LP11,int* LP12,int* LP13,int* LP14,int* LP15,
			 int* LP16,int* LP17,int* LP18,int* LP19,int* LP20,
			 double* me_i,double* mmu_i,double* mtau_i,
			 double* mu_i,double* md_i,double* ms_i,double* mc_i,double* mb_i,double* mt_i,
			 double* mZ_i,double* mW_i,double* mH_i,
			 double* alpha_i, double* gf_i,
			 double* au_i, double* vu_i, double* ad_i, double* vd_i,
			 int* doprint);
   void eprc_init_alpos_(
			 int* LP1,int* LP2,int* LP3,int* LP4,int* LP5,
			 int* LP6,int* LP7,int* LP8,int* LP9,int* LP10,
			 int* LP11,int* LP12,int* LP13,int* LP14,int* LP15,
			 int* LP16,int* LP17,int* LP18,int* LP19,int* LP20,
			 int* doprint);
   void setpar_alpos_(
		      double* me_i,double* mmu_i,double* mtau_i,
		      double* mu_i,double* md_i,double* ms_i,double* mc_i,double* mb_i,double* mt_i,
		      double* mZ_i,double* mW_i,double* mH_i,
		      double* alpha_i, double* gf_i,
		      double* au_i, double* vu_i, double* ad_i, double* vd_i,
		      int* doprint);
   double aemrun_(double* q2);
   void getpar_alpos_(double* q2_in, double* DR_o, double* SWEFF2_o, double* AKAPPA_o,
		      double* SW2_o, double* MW_o, double* MZ_o);


}

// __________________________________________________________________________________________ //
const std::vector<std::string> AEprc::fRequirements = {"mur",     // renormalization scale
						       "LPAR1",   // BORN CROSS SECTION
						       "LPAR2",   // 1-LOOP CORRECTIONS
						       "LPAR3",   // EXPONENTIATION:
						       "LPAR4",   // MW OR GMU FIXED
						       "LPAR5",   // QCD CORR IN DELTA-R
						       "LPAR6",   // PARTON DISTRIB.
						       "LPAR7",   // SIGMA-GAMMA
						       "LPAR8",   // SIGMA-GAMMA-Z
						       "LPAR9",   // SIGMA-Z
						       "LPAR10",  // SIGMA-W
						       "LPAR11",  // QED CORRECTIONS
						       "LPAR12",  // LEPTONIC QED
						       "LPAR13",  // HADRONIC QED
						       "LPAR14",  // LEPTON-QUARK-INTRF
						       "LPAR15",  // WEAK CORRECTIONS
						       "LPAR16",  // WEAK BOXES
						       "LPAR17",  // GAMMA OR/AND Z
						       "LPAR18",  // LEPTONIC LLA INCL.
						       "LPAR19",  // X/Y DEFINITION
						       "LPAR20",  // FIXED Q2 IN HARD BS
						       "me","mmu","mtau", // lepton masses
						       "mu","md","ms","mc","mb","mt", // quark masses
						       "mZ","mW","mH", // boson masses
						       "alpha", // alpha_em(Q=0) ~ 1./137
						       "gf", // gf
						       "au","ad","vu","vd", // quark couplings
//						       "convfac", // conversion 
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AEprc::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AEprc::fFunctionName = "EPRC"; //< The function's name
int AEprc::fNinstances = 0;

// __________________________________________________________________________________________ //
AEprc::AEprc(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("AEprc");
   if ( ++fNinstances > 1 ) {
      cout<<"Error. AEprc::AEprc(). Only one AlphaEmRun function is allowed."<<endl;
      exit(1);
   }
   fValue.resize(1);
   fError.resize(1);
}


// __________________________________________________________________________________________ //
AEprc::~AEprc() {
}


// ___________________________________________________________________________________________ //
bool AEprc::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   // 'real' QCDNUM init is called in constructor

   int P1 = PAR(LPAR1);   // BORN CROSS SECTION   
   int P2 = PAR(LPAR2);   // 1-LOOP CORRECTIONS   
   int P3 = PAR(LPAR3);   // EXPONENTIATION:      
   int P4 = PAR(LPAR4);   // MW OR GMU FIXED      
   int P5 = PAR(LPAR5);   // QCD CORR IN DELTA-R  
   int P6 = PAR(LPAR6);   // PARTON DISTRIB.      
   int P7 = PAR(LPAR7);   // SIGMA-GAMMA          
   int P8 = PAR(LPAR8);   // SIGMA-GAMMA-Z        
   int P9 = PAR(LPAR9);   // SIGMA-Z              
   int P10 = PAR(LPAR10);  // SIGMA-W              
   int P11 = PAR(LPAR11);  // QED CORRECTIONS      
   int P12 = PAR(LPAR12);  // LEPTONIC QED         
   int P13 = PAR(LPAR13);  // HADRONIC QED         
   int P14 = PAR(LPAR14);  // LEPTON-QUARK-INTRF   
   int P15 = PAR(LPAR15);  // WEAK CORRECTIONS     
   int P16 = PAR(LPAR16);  // WEAK BOXES           
   int P17 = PAR(LPAR17);  // GAMMA OR/AND Z       
   int P18 = PAR(LPAR18);  // LEPTONIC LLA INCL.   
   int P19 = PAR(LPAR19);  // X/Y DEFINITION       
   int P20 = PAR(LPAR20);  // FIXED Q2 IN HARD BS  
   double me = PAR(me);
   double mmu = PAR(mmu);
   double mtau = PAR(mtau); // lepton masses         
   double mu = PAR(mu);
   double md = PAR(md);
   double ms = PAR(ms);
   double mc = PAR(mc);
   double mb = PAR(mb);
   double mt = PAR(mt); // quark masses            
   double mZ = PAR(mZ);
   double mW = PAR(mW);
   double mH = PAR(mH); // boson masses                                                                             
   double alpha = PAR(alpha); // alpha_em(Q=0) ~ 1./13                                                                 
   double gf = PAR(gf); // gf
   double au = PAR(au); //
   double ad = PAR(ad); //
   double vu = PAR(vu); //
   double vd = PAR(vd); //
   //double conv = PAR(convfac);
   int doprint = 0;

    // eprc_init_alpos_full_(
    // 		     &P1,&P2,&P3,&P4,&P5,
    // 		     &P6,&P7,&P8,&P9,&P10,
    // 		     &P11,&P12,&P13,&P14,&P15,
    // 		     &P16,&P17,&P18,&P19,&P20,
    // 		     &me,&mmu,&mtau, // lepton masses
    // 		     &mu,&md,&ms,&mc,&mb,&mt, // quark masses
    // 		     &mZ,&mW,&mH, // boson masses
    // 		     &alpha, // alpha_em(Q=0) ~ 1./137
    // 		     &gf,
    // 		     &au,&vu,&ad,&vd,
    // 		     &doprint);

   eprc_init_alpos_(
    		     &P1,&P2,&P3,&P4,&P5,
    		     &P6,&P7,&P8,&P9,&P10,
    		     &P11,&P12,&P13,&P14,&P15,
    		     &P16,&P17,&P18,&P19,&P20,
    		     &doprint);

   return true;
}


// __________________________________________________________________________________________ //
bool AEprc::Update() {

   double me = PAR(me);
   double mmu = PAR(mmu);
   double mtau = PAR(mtau); // lepton masses         
   double mu = PAR(mu);
   double md = PAR(md);
   double ms = PAR(ms);
   double mc = PAR(mc);
   double mb = PAR(mb);
   double mt = PAR(mt); // quark masses            
   double mZ = PAR(mZ);
   double mW = PAR(mW);
   double mH = PAR(mH); // boson masses                                                                             
   double alpha = PAR(alpha); // alpha_em(Q=0) ~ 1./13                                                                 
   double gf = PAR(gf); // gf                      
   //double conv = PAR(convfac);
   double au = PAR(au); //
   double ad = PAR(ad); //
   double vu = PAR(vu); //
   double vd = PAR(vd); //
   static int doprint = 1;

   setpar_alpos_(&me,&mmu,&mtau, // lepton masses
		 &mu,&md,&ms,&mc,&mb,&mt, // quark masses
		 &mZ,&mW,&mH, // boson masses
		 &alpha,&gf,
		 &au,&vu,&ad,&vd,
		 &doprint);
   doprint=0;

   double q2 = PAR(mur);
   q2*=q2;
   //double arun = aemrun_(&q2);
   fValue.resize(4);
   fError.resize(4);

   fValue[0] = aemrun_(&q2);
   fError[0] = 0;

   // fValue[0]:  aem(q2)
   // fValue[1]:  aem(deltaR)
   // fValue[2]:  aem(sin2thweff(q2))
   // fValue[3]:  kappa
   // fValue[3]:  sin2thw
   double _mz,_mw,_sw2;
   getpar_alpos_(&q2,&fValue[1],&fValue[2],&fValue[3],&_sw2,&_mw,&_mz);

 
   return true;
}

// ______________________________________________________________________________________ // 
std::vector<double> AEprc::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);
   std::vector<double> ret(fValue.size());
   if ( n != 1 ) {
      cout<<"Error in AEprc::GetQuick(...). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   // for a simple example see: http://en.wikipedia.org/wiki/Variadic_function
   va_list ap;
   va_start(ap, n);
   double mur = va_arg(ap, double);
   va_end (ap);
   double q2 = mur*mur;
   ret[0] = aemrun_(&q2);
   return ret;
}


// ______________________________________________________________________________________ // 
std::vector<double> AEprc::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using: 
   //!   ::GetQuick(vector<double> mur);   
   std::vector<double> ret(fValue.size());
   if ( mur.size() != 1 ) {
      cout<<"Error in AEprc::GetQuick(vector<double>). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   double q2 = mur[0]*mur[0];
   ret[0] = aemrun_(&q2);
   return ret;
}
