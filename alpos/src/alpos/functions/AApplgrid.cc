// DB 07.12.2015

#include "alpos/functions/AApplgrid.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include "alpos/AlposTools.h"
#include "alpos/ATheory.h"

using namespace std;


// __________________________________________________________________________________________ //
// ---  global functions for Applgrid
string _g_appl_pdf_name = ""; //!<pdf function for upcoming calls
string _g_appl_as_name = ""; //!<alpha_s function for upcoming calls
AFuncD* _g_appl_pdf_ = NULL;
AFuncD* _g_appl_as_ = NULL;
void appl_pdf(const double& x, const double& Q, double* f) {
   if ( _g_appl_pdf_ == NULL ) {
      say::error["appl_pdf"]<<"No PDF function specified."<<endl;
      return;
   }
   double xp = x > 1 ? 1 : x;
   const vector<double>& xfx = _g_appl_pdf_->GetQuick(2,xp,Q);
   for ( unsigned int i = 0 ; i<13 ; i++ ) 
      f[i]=xfx[i];
   return;
}

void appl_pdf_bar(const double& x, const double& Q, double* f) {
   if ( _g_appl_pdf_ == NULL ) {
      say::error["appl_pdf_bar"]<<"No PDF function specified."<<endl;
      return;
   }
   double xp = x > 1 ? 1 : x;
   const vector<double>& xfx = _g_appl_pdf_->GetQuick(2,xp,Q);
   for ( unsigned int i = 0 ; i<13 ; i++ ) 
      f[12-i]=xfx[i];
   return;
}


double appl_alphas(const double& Q) {
   if ( _g_appl_as_ == NULL ) {
      say::error["appl_alphas"]<<"No Alpha_s function specified."<<endl;
      return 0;
   }
   vector<double> muf = {Q};
   return _g_appl_as_->GetQuick( muf )[0];
}


// __________________________________________________________________________________________ //
const std::vector<std::string> AApplgrid::fRequirements = {
// "Filename",                      // fastNLO table filenmae
 							  "PDF",                           // a 'PDF' function with Quick-function. Moslty LHAPDF6
 							  "Alpha_s",                       // a alpha_s(mu_r) function.
 							  "ScaleFacMuR","ScaleFacMuF",     // Scale factors for ren. and fact. scale
// 							  "Units",                         // publication or absolute units (to obtain same units as your data table)
 							  "iOrd",                          // order 
// 							  "iThr",                          // Use threshold corrections if available in fastNLO table
// 							  "MuRFuncForm","MuFFuncForm"      // mu_r and mu_f functional form for fastNLO flexible scale tables
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AApplgrid::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApplgrid::fFunctionName = "Applgrid"; //< The function's name


// __________________________________________________________________________________________ //
AApplgrid::AApplgrid(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AApplgrid");
}


// __________________________________________________________________________________________ //
AApplgrid::~AApplgrid() {
   //if ( fnlo ) delete fnlo;
   for ( auto g : fgrids ) delete g;
}


// ___________________________________________________________________________________________ //
bool AApplgrid::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;

   using namespace AlposTools;
   vector<string> gridfiles;
   vector<int> firstbins, lastbins;
   if ( EXIST_NS(Gridfile,GetAlposName() )) 
      gridfiles.push_back(STRING_NS(Gridfile,GetAlposName()));
   else if ( EXIST_NS(Gridfiles,GetAlposName() ) )
      gridfiles = STRING_ARR_NS(Gridfiles,GetAlposName()); // direct access to array
   else if ( EXIST_NS(Grids,GetAlposName()) ) {
      gridfiles = STRING_COL_NS(Grids,Gridfiles,GetAlposName());
      firstbins = INT_COL_NS(Grids,FirstBin,GetAlposName());
      lastbins  = INT_COL_NS(Grids,LastBin,GetAlposName());
   }

   for ( auto f : gridfiles ) {
      info["Init"]<<"Rading grid file: "<<f<<endl;
      fgrids.push_back(new appl::grid(f));
   }
   
   // init binmap if specified
   if ( !firstbins.empty() ) {
      //for ( auto grid : fgrids ) {
      for ( unsigned int ig = 0; ig<gridfiles.size() ; ig++ ){
	 vector<bool> bmap(fgrids[ig]->Nobs());
	 for ( unsigned int ib=0 ; ib<bmap.size() ; ib++ ) {
	    bmap[ib] = (ib>=firstbins[ig] && ib <= lastbins[ig] );
	 }
	 fBinmap += bmap;
      }
   }
   return true;
}


// __________________________________________________________________________________________ //
bool AApplgrid::Update() {

   using namespace AlposTools;
   
   // --- 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(PDF);
   UPDATE(Alpha_s);

   // --- set PDF to appl_pdf
   _g_appl_pdf_name = TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string("PDF"))->GetAlposName();
   _g_appl_as_name = TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string("Alpha_s"))->GetAlposName();
   _g_appl_pdf_ = TheoryHandler::Handler()->GetFuncD(_g_appl_pdf_name);
   _g_appl_as_ = TheoryHandler::Handler()->GetFuncD(_g_appl_as_name);


   //void appl_pdf(const double& x, const double& Q, double* f) {

   double mur = PAR(ScaleFacMuR);
   double muf = PAR(ScaleFacMuF);
   int iorder = PAR(iOrd);
   fValue.clear();
   for ( auto g : fgrids ) {
      vector<double> xs = g->vconvolute(appl_pdf, appl_pdf_bar, appl_alphas, iorder, mur, muf );
      fValue += xs;
      //cout<<"xs.size: "<<xs.size()<<"\tfValue.size: "<<fValue.size()<<endl;      
   }

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
