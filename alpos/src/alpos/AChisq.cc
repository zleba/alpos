#include "alpos/AChisq.h"
#include "alpos/ATheory.h"
#include "alpos/ASuperData.h"
#include "alpos/AlposTools.h"
#include "alpos/AInvMatrices.h"
#include <cmath>

// Apccpp
//#include <apccpp/Apccpp.h>
#include <Apccpp.h>

/*
 AChisq

 */

using namespace std;

/**
 * A base class for chisq functions
 */

//____________________________________________________________________________________ //
std::vector<double> AChisqBase::fNuisance = std::vector<double>(0);
std::vector<std::string> AChisqBase::fNuisanceName = std::vector<std::string>(0);

//____________________________________________________________________________________ //
std::map<std::string, double> AChisqBase::GetNuisanceParameters() {
   //! return nuisance parameters:
   //!  < name , value >
   //! 
   //! mind; fNuisance is calulated by 'AChisqBase', while
   //! names of nuisance parameters are provided by the chisq-definition itself.
   //!
   std::map<std::string, double> ret;
   if ( fNuisance.size() != fNuisanceName.size() ){
      say::error["AChisqBase::GetNuisance"]<<"Names and values of nuisance parameters are not consistent. nNames: "<<fNuisanceName.size()<<"\tnNuisance:" <<fNuisance.size()<<endl;
      exit(303);
   }
   
   for ( unsigned int i = 0 ; i<fNuisance.size() ; i++ ) {
      ret[fNuisanceName[i]] = fNuisance[i];
   }
   return ret;
}


//____________________________________________________________________________________ //
void AChisqBase::PrintNuisanceParameters(bool printnone, std::map<std::string, double>* nuisance) const {
   //! Print nuisance parameters if these have been calculated
   map<string,double> nui = nuisance ? *nuisance : this->GetNuisanceParameters();
   if ( nui.empty() ) {
      if ( printnone ) 
	 say::info["AChisqBase::PrintNuisanceParameters"]<<"No Nuisance parameters have been calculated."<<endl;// by chisq-definition '"<<GetChisqName()<<"'."<<endl;
   }
   else {
      //info["AChisqBase::PrintNuisanceParameters"]<<"Printing nuisance parameters calculated by '"<<GetChisqName()<<"'."<<endl;
      for ( auto ie : nui ) {
         printf("   %-22s\t% 5.3f\n",ie.first.c_str(),ie.second);
      }
   }
}


//____________________________________________________________________________________ //
std::vector<double> AChisqBase::CalcNuisance(const std::vector<const std::vector<double>* >& gk, const TMatrixD& VInv, const std::vector<double>& dmt)  {
   //! determine nuisance parameters
   //! Input:
   //!   gk:   Correlated errors for each source (gk[error][bin])
   //!   VInv: Inverse covariance matrix
   //!   dmt:  (data-theory)
   //!
   //! The nuiscane parameters are calculated by solving following eqation:
   //!     (gVg+E) * b = a
   //! with:
   //!     gVg_ik = G^T_i V^-1_ik G_k [an e*e matrix with e=number of correlated errors]
   //!     E   = e*e unit matrix
   //!     G_i = Correlated Errors
   //!     a   = G^T_i V^-1_ik (d-t)_i
   //!     b   = nuisance parameters
   //!
   //! Covariance matrix holds all 'matrix type and uncorrelated uncertainties'


  // --- get data and uncertainties
   const unsigned int NE = gk.size();
   if ( NE == 0 ) return vector<double>();
   const unsigned int nb = dmt.size();
   const TMatrixD& A = VInv;
   // --- GkA
   vector<vector<double> > GkA(NE); // GkA = GkT_j * V^_ij
   for ( unsigned int ie = 0 ; ie<NE ; ie++) {
      GkA[ie].resize(nb);
      const vector<double>& vGk = *gk[ie];
      for ( unsigned int x = 0 ; x<nb ; x++) {
         for ( unsigned int y = 0 ; y<nb ; y++) {
            GkA[ie][x] +=vGk[y]*A[x][y];
         }
      }
   }
   // --- gVg
   TMatrixDSym gVg(NE); // gVg-ij = GkA_i*G_j
   for ( unsigned int e1 = 0 ; e1<NE ; e1++) {
      const vector<double>& vGk = *gk[e1];
      for ( unsigned int e2 = 0 ; e2<=e1 ; e2++) {
         double val = 0;
         for ( unsigned int x = 0 ; x<nb ; x++) {
            val += GkA[e2][x] * vGk[x];
         }
         if (e1==e2) val+=1; //  unit matrix: B = gVg-1 (B->gVg) (from \sum b^2)
         gVg(e1,e2)=val;
         gVg(e2,e1)=val;
      }
   }
   //gVg.Print();
   // --- BI
   TMatrixDSym BI = AlposTools::InvertChol(gVg); // BI*b=a
   //TMatrixD BI = AlposTools::InvertLU(gVg); // BI*b=a
   //BI.Print();
   // --- a
   vector<double> a(NE); // a_k =  GkT_i * V^_ij (d-t)_j
   for ( unsigned int ei =0 ;ei<NE; ei++ ) {
      for ( unsigned int x = 0 ; x<nb ; x++) {
         a[ei] += GkA[ei][x] * dmt[x];
      }
   }
   // --- b
   return SolveStoreNuisance(BI,a);
}


//____________________________________________________________________________________ //
std::vector<double> AChisqBase::CalcNuisance(const std::vector<const std::vector<double>* >& gk, const vector<double>& VDiag, const std::vector<double>& dmt)  {
   //! determine nuisance parameters
   //! Input:
   //!   gk:   Correlated errors for each source (gk[error][bin])
   //!   VDiag: Covariances (Squared uncorrelated uncertainties. no inverse!)
   //!   dmt:  (data-theory)
   //!
   //! The nuiscane parameters are calculated by solving following eqation:
   //!     (gVg+E) * b = a
   //! with:
   //!     gVg_ik = G^T_i V^-1_ik G_k [an e*e matrix with e=number of correlated errors]
   //!     E   = e*e unit matrix
   //!     G_i = Correlated Errors
   //!     a   = G^T_i V^-1_ik (d-t)_i
   //!     b   = nuisance parameters
   //!
   //! Covariance matrix holds all 'matrix type and uncorrelated uncertainties'

  // --- get data and uncertainties
   const unsigned int NE = gk.size();
   if ( NE == 0 ) return vector<double>();
   const unsigned int nb = dmt.size();
   const vector<double>& A = VDiag;

   // --- GkA
   vector<vector<double> > GkA(NE); // GkA = GkT_j * V^_ij
   for ( unsigned int ie = 0 ; ie<NE ; ie++) {
      GkA[ie].resize(nb);
      const vector<double>& vGk = *gk[ie];
      for ( unsigned int x = 0 ; x<nb ; x++) {
         GkA[ie][x] +=vGk[x]/A[x];
      }
   }
   // --- gVg
   TMatrixDSym gVg(NE); // gVg-ij = GkA_i*G_j
   for ( unsigned int e1 = 0 ; e1<NE ; e1++) {
      const vector<double>& vGk = *gk[e1];
      for ( unsigned int e2 = 0 ; e2<=e1 ; e2++) {
         double val = 0;
         for ( unsigned int x = 0 ; x<nb ; x++) {
            val += GkA[e2][x] * vGk[x];
         }
         if (e1==e2) val+=1; //  unit matrix: B = gVg-1 (B->gVg) (from \sum b^2)
         gVg(e1,e2)=val;
         gVg(e2,e1)=val;
      }
   }
   //gVg.Print();
   // --- BI
   TMatrixDSym BI = AlposTools::InvertChol(gVg); // BI*b=a
   //TMatrixD BI = AlposTools::InvertLU(gVg); // BI*b=a
   //BI.Print();
   // --- a
   vector<double> a(NE); // a_k =  GkT_i * V^_ij (d-t)_j
   for ( unsigned int ei =0 ;ei<NE; ei++ ) {
      for ( unsigned int x = 0 ; x<nb ; x++) {
         a[ei] += GkA[ei][x] * dmt[x];
      }
   }
   // --- b
   return SolveStoreNuisance(BI,a);
}


//____________________________________________________________________________________ //
std::vector<double> AChisqBase::CalcNuisance(const std::vector<std::vector<double> >& gk, const TMatrixD& VInv, const std::vector<double>& dmt)  {
   //! determine nuisance parameters
   //! Input:
   //!   gk:   Correlated errors for each source (gk[error][bin])
   //!   VInv: Inverse covariance matrix
   //!   dmt:  (data-theory)
    vector<const vector<double>* > bla;
    for (unsigned int i=0; i<gk.size() ; i++ )
       bla.push_back(&(gk[i]));
   //for ( auto i : gk ) bla.push_back(&i);
   // cout<<"Hier geht was nicht!"<<endl;
   // exit(1);
   return CalcNuisance(bla,VInv,dmt);
 }


//____________________________________________________________________________________ //
std::vector<double> AChisqBase::CalcNuisance(const std::vector<std::vector<double> >& gk, const vector<double>& VDiag, const std::vector<double>& dmt)  {
   //! determine nuisance parameters
   //! Input:
   //!   gk:   Correlated errors for each source (gk[error][bin])
   //!   VInv: covariaances (squared uncorrelated uncertainties)
   //!   dmt:  (data-theory)
    vector<const vector<double>* > bla;
    for (unsigned int i=0; i<gk.size() ; i++ )
       bla.push_back(&(gk[i]));
   return CalcNuisance(bla,VDiag,dmt);
 }


//____________________________________________________________________________________ //
std::vector<double> AChisqBase::SolveStoreNuisance(const TMatrixDSym& BI, const std::vector<double> a) {
   //! calculate b = BI*a
   //! and store b as result in fNuisance
   //! return fNuisance
   const unsigned int NE = a.size();
   //vector<double> fNuisance(NE); // b_k = BI_k * a_k
   fNuisance.clear();
   fNuisance.resize(NE);
   double sumsq=0;
   for ( unsigned int e = 0 ; e<NE ; e++) {
      //fNuisance[e]=0;
      for ( unsigned int ex = 0 ; ex<NE ; ex++) {
         fNuisance[e] += BI(e,ex)* a[ex];
      }
      //cout<<"b["<<e<<"]="<<b[e]<<endl;
      sumsq+=pow(fNuisance[e],2);
   }
   //cout<<"sum_bi^2 = "<<sumsq<<endl;
   return fNuisance;

 }


//____________________________________________________________________________________ //
/**
 *
 * AChisqCov
 *
 * A simple chisq using the covariance matrix
 *  chisq  = ( d - t ) V^-1  ( d-t )
 */

double AChisqCov::DoEval(const double *p) const {
   //! Calculate chisq
   //!    chisq = sum_ij [m_i - t_i] * V^-1_ij * [m_i - t_i]
   //! V is the covariance matrix containing all uncertainties
   //!
   //! \note any uncertainties marked as multiplicative
   //!       will be treated as additive by this chi-square
   //!       definition.

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )
      SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get data and theory values
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get uncertainties from data and theory objects

   // pointers to 'sum' matrices entering into the 'final' covariance matrix
   std::set<const TMatrixDSym*> matCollection;

   // additive errors can be considered directly
   const TMatrixDSym* InvCov = &AInvMatrices::Instance()->GetInvMatrix(
      {
         &fData->GetSumErrorMatrix("AA", "AbsAvTotAll"),
         &fTheo->GetSumErrorMatrix("AA", "AbsAvTotAll")
      }
   );

   // --- loop over all super-vector data points and calculate chisq
   double chisq = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += pow((da[x]-th[x]),2) * (*InvCov)[x][x] ;
      //cout<<"dt= "<<da[x]<<"\tth= "<<th[x]<<"\terr= "<<InvCov[x][x] <<endl;
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * (da[x]-th[x]) * (*InvCov)[x][y] * (da[y]-th[y]);
      }
   }

   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqCMS
 *
 * A simple chisq using the covariance matrix
 *  chisq  = ( d - t ) V^-1  ( d-t )
 *  where V is composed from multiplicative and additive uncertainties
 *  HOWEVER! uncorrelated and statistical uncertainties are always treated as additive !  
 *  Theoretical uncertainties are treated as 'additive' as well, i.e. rel*theo. 
 */

double AChisqCMS::DoEval(const double *p) const {
   //! Calculate chisq
   //!    chisq = sum_ij [m_i - t_i] * V^-1_ij * [m_i - t_i]
   //! V is the covariance matrix containing all uncertainties
   //!   where V is composed from multiplicative and additive uncertainties
   //!   HOWEVER! uncorrelated and statistical uncertainties are always treated as additive !  
   //!   Theoretical uncertainties are treated as 'additive' as well, i.e. rel*theo. 
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )
      SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get data and theory values
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get uncertainties from data and theory objects
   // additive errors can be considered directly
   const TMatrixDSym& CovStat      = fData->GetSumErrorMatrix("AS", "AbsAvTotAll");
   const TMatrixDSym& CovUnc       = fData->GetSumErrorMatrix("EY", "AbsAvUncAll");
   const TMatrixDSym& CovMultSysRel= fData->GetSumErrorMatrix("EY", "RelAvCorAll");
   //const TMatrixDSym* CovAddSys  = ; // this is not forseen here !
   const TMatrixDSym& CovTheo      = fData->GetSumErrorMatrix("TA", "AbsAvTotAll"); // this is not forseen here !
   const TMatrixDSym& CovFTheo     = fTheo->GetSumErrorMatrix("AA", "AbsAvTotAll"); // this is not forseen here !

   TMatrixDSym CovMultSys(CovMultSysRel);
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      for ( unsigned int y = 0 ; y<th.size() ; y++ ) {
	 CovMultSys[x][y] = CovMultSysRel[x][y] * th[x] * th[y];
      }
   }

   const TMatrixDSym* InvCov = &AInvMatrices::Instance()->GetInvMatrix(
      {&CovStat, &CovUnc, &CovMultSys, &CovTheo, &CovFTheo}, true );


   // --- loop over all super-vector data points and calculate chisq
   double chisq = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += pow((da[x]-th[x]),2) * (*InvCov)[x][x] ;
      //cout<<"dt= "<<da[x]<<"\tth= "<<th[x]<<"\terr= "<<InvCov[x][x] <<endl;
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * (da[x]-th[x]) * (*InvCov)[x][y] * (da[y]-th[y]);
      }
   }

   return chisq;
}

//____________________________________________________________________________________ //
/**
 *
 * AChisqCovMult
 *
 * A simple chisq using the covariance matrix
 *  chisq  = ( d - t ) V^-1  ( d-t )
 *
 * The main difference to AChisqCov is the rescaling of all
 * absolute errors marked 'multiplicative' to correspond to the
 * current theory values.
 */

double AChisqCovMult::DoEval(const double *p) const {
   //! Calculate chisq
   //!    chisq = sum_ij [m_i - t_i] * V^-1_ij * [m_i - t_i]
   //! V is the covariance matrix containing all uncertainties
   //!
   //! \note if uncertainties are marked as multiplicative,
   //!       this means they are given relative to the theory
   //!       and must therefore the inverse covariance matrix
   //!       must be calculated for every iteration.
   //!       This can bring about a significant increase in the
   //!       time complexity.

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )
      SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get data and theory values
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get uncertainties from data and theory objects

   bool recalcInvCov=false;  // need to recalculate inverse covariance matrix?
   // pointers to 'sum' matrices entering into the 'final' covariance matrix
   std::set<const TMatrixDSym*> matCollection;

   // additive errors can be considered directly
   matCollection.insert(&fData->GetSumErrorMatrix("AA", "AbsAvTotAdd"));

   // handle multiplicative errors stored in AData (if they exist)
   TMatrixDSym dataCovMatMul;
   if (fData->HasMultErrors()) {
      dataCovMatMul.ResizeTo(fData->N(), fData->N());
      dataCovMatMul = fData->GetSumErrorMatrix("AA", "AbsAvTotMul");
      // rescale error reference values to current theory values
      AlposTools::CovRescale(dataCovMatMul, da, th);
      matCollection.insert(&dataCovMatMul);
      recalcInvCov=true;
   }

   // handle errors stored in Theory objects (if they exist)
   // and assume these are multiplicative
   TMatrixDSym theoCovMatMul;
   if (fTheo->HasErrors()) {
      // get relative covariance matrix and scale to current theory cross sections
      theoCovMatMul.ResizeTo(fData->N(), fData->N());
      theoCovMatMul = fTheo->GetSumErrorMatrix("AA", "RelAvTotAll");
      AlposTools::CovRelToCov(theoCovMatMul, th);
      matCollection.insert(&theoCovMatMul);
      recalcInvCov=true;
   }

   // --- get/construct inverse covariances
   const TMatrixDSym* InvCov;
   if (!recalcInvCov) {
      // if no multiplicative errors so far, use AInvMatrices for fast access
      InvCov = &AInvMatrices::Instance()->GetInvMatrix(matCollection);
   }
   else {
      // otherwise, need to recalculate anyway, so don't use AInvMatrices -> create temporary matrix
      InvCov = new TMatrixDSym(AlposTools::InvertChol(matCollection));
   }


   // --- loop over all super-vector data points and calculate chisq
   double chisq = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += pow((da[x]-th[x]),2) * (*InvCov)[x][x] ;
      //cout<<"dt= "<<da[x]<<"\tth= "<<th[x]<<"\terr= "<<InvCov[x][x] <<endl;
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * (da[x]-th[x]) * (*InvCov)[x][y] * (da[y]-th[y]);
      }
   }

   // destroy temporary matrix
   if (recalcInvCov) {
      delete InvCov;
   }

   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqCovStatUncorr
 *
 * A simple chisq using the covariance matrix
 *  chisq  = ( d - t ) V^-1  ( d-t )
 *  where V is composed only from uncorrelated and statistical uncertainties
 */

double AChisqCovStatUncorr::DoEval(const double *p) const {
   //! Calculate chisq
   //!    chisq = sum_ij [m_i - t_i] * V^-1_ij * [m_i - t_i]
   //! V containing all uncertainties

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )
      SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get data and theory values
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data for uncertainties and inverse covariances
   //const TMatrixD& InvCov = fData->GetInverseCovarianceStatUncorr();  // old interface
   const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix(
      { &fData->GetSumErrorMatrix("ES","AbsAvTot"), &fData->GetSumErrorMatrix("EY","AbsAvUnc")} );

   // --- loop over all super-vector data points and calculate chisq
   double chisq = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      for ( unsigned int y = 0 ; y<th.size() ; y++ ) {
         chisq += (da[x]-th[x]) * InvCov[x][y] * (da[y]-th[y]);
      }
   }
   return chisq;
}
//____________________________________________________________________________________ //




//____________________________________________________________________________________ //
/**
 *
 * AChisqLogNormal
 *
 * A chisq as used in [arxiv:1406.4709]
 *
 * Correlated uncertainties are included in the covariance matrix
 * while being constructed from the relative uncertainties
 *
 */

AChisqLogNormal::AChisqLogNormal(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo) {
}

unsigned int AChisqLogNormal::NDim() const {
   return fFitPar.size();
}

ROOT::Math::IMultiGenFunction* AChisqLogNormal::Clone() const {
   //! You must implement a 'Clone' method
   return new AChisqLogNormal(*this);//fFitPar,fData,fTheo);
};


double AChisqLogNormal::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(log(m_i) - log(t_i)] * V^-1_ij * [(log(m_i) - log(t_i)]
   //! With V containing stat, correlated and uncorrelated uncertainties

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data for uncertainties
   //const map<string,AError>& errors = fData->GetAllErrors();

   // --- get data for uncertainties and inverse covariances
   //const TMatrixD& InvCov = fData->GetInverseCovarianceRel();
   //const TMatrixDSym& InvCov = fData->GetSumErrorMatrix("EA","RelAvTotInv");

   const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix(
         {&fData->GetSumErrorMatrix("AA","RelAvTot"),
          &fTheo->GetSumErrorMatrix("AA","RelAvTot")});

   //const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix({&fData->GetSumErrorMatrix("EA","RelAvTot"),&fData->GetSumErrorMatrix("EA","RelAvTot")});

   double chisq = 0;
   // --- loop over all super-vector data points and calculate chisq
   vector<double> dldxltx(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      //dldxltx[x] = log(da[x]) - log(th[x]);
      dldxltx[x] = log(da[x]/th[x]);
      if ( !isfinite(dldxltx[x]) ){
	 warn["DoEval"]<<"Contribution to chisq is not finite in bin "<<x<<" data="<<da[x]<<", theo="<<th[x]<<"\tlog(d/t)="<<dldxltx[x]<<"\t data name: "<<Theo()->GetAlposName()<<endl;
	 dldxltx[x] = 1./sqrt(InvCov[x][x])*8.; // some non-zero contribution to chisq
      }
      chisq += dldxltx[x]*dldxltx[x] * InvCov[x][x];
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * dldxltx[x] * InvCov[x][y] * dldxltx[y];
         //cout<<"chisq: "<<chisq<<"\tx,y: "<<x<<","<<y<<"\tda[x]="<<da[x]<<"\tth="<<th[x]<<endl;
      }
   }

   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    cout<<"x="<<x<<" da: "<<da[x]<<",\t th: "<<th[x]<<"\t chisq: "<<(da[x]-th[x])*(da[x]-th[x])/InvCov[x][x]<<endl;
   // }

   if ( std::isnan(chisq)  ){
      cerr<<"Error.  AChisqLogNormal::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqNormalLogNormal
 *
 * An updated chisq as used in [arxiv:1406.4709]
 *
 * Correlated uncertainties are included in the covariance matrix
 * while being constructed from the relative uncertainties
 *
 * Data points without statistical uncertainties are treated
 * with normal distributed errors.
 * These are assumed be flagged further as 'additive'.
 *
 * This consistency must be ensured by the user 
 * (no stat -> all syst. errors are flaged as additive)!
 * Furthermore, these syst. uncertainties may not be 
 * correlated to datapoints with log-normal uncertainties
 *
 *
 */
double AChisqNormalLogNormal::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(log(m_i) - log(t_i)] * V^-1_ij * [(log(m_i) - log(t_i)]
   //! With V containing stat, correlated and uncorrelated uncertainties

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   vector<double> th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   vector<double> da = Data()->GetValues();// VALUES_ANY("SuperData");
   //const vector<double>& dstat = fData->GetUncertaintyStat();
   const vector<double>& dstat = fData->GetSumError("AS", "AbsAvTot");
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      if ( dstat[x]==0 ) {
	 //cout<<"NoStat Error. Assume all uncertianties for datapoint "<<x<<" as Normal-distributed."<<endl;
	 th[x]=exp(th[x]);
	 da[x]=exp(da[x]);
      }
   }   

   // --- get data for uncertainties
   //const map<string,AError>& errors = fData->GetAllErrors();

   // --- get data for uncertainties and inverse covariances
   //const TMatrixD& InvCov = fData->GetInverseCovarianceRel();
   //const TMatrixDSym& InvCov = fData->GetSumErrorMatrix("EA","RelAvTotInv");

   const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix(
      {&fData->GetSumErrorMatrix("EA","RelAvTotMul"), &fData->GetSumErrorMatrix("EA","AbsAvTotAdd")} );
   //const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix({&fData->GetSumErrorMatrix("EA","RelAvTot"),&fData->GetSumErrorMatrix("EA","RelAvTot")});

   double chisq = 0;
   // --- loop over all super-vector data points and calculate chisq
   vector<double> dldxltx(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      //dldxltx[x] = log(da[x]) - log(th[x]);
      //dldxltx[x] = log(fabs(da[x])) - log(fabs(th[x]));
      dldxltx[x] = log(da[x]/th[x]);
      if ( std::isnan(dldxltx[x]) ) {
	 cout<<"Warning! Data point "<<x<<" has 'nan' for log(d)-log(t). log(d)="<<log(fabs(da[x]))<<", log(t)="<<log(fabs(th[x]))<<", th="<<th[x]<<endl;
	 TheoryHandler::Handler()->PrintCurrentTheorySet();
	 cout<<endl;
	 cout<<"Exiting, due to nan in chisq."<<endl;
	 cout<<endl;
	 exit(4);
	 dldxltx[x]=1.e11;
      }
      chisq += pow(dldxltx[x],2) * InvCov[x][x];
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * dldxltx[x] * InvCov[x][y] * dldxltx[y];
         //cout<<"chisq: "<<chisq<<"\tx,y: "<<x<<","<<y<<"\tda[x]="<<da[x]<<"\tth="<<th[x]<<endl;
      }
   }

   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    cout<<"x="<<x<<" da: "<<da[x]<<",\t th: "<<th[x]<<"\t chisq: "<<(da[x]-th[x])*(da[x]-th[x])/InvCov[x][x]<<endl;
   // }

   if ( std::isnan(chisq)  ){
      cerr<<"Error.  AChisqLogNormal::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqLogNormalStatUncorr
 *
 * A chisq as used in [arxiv:1406.4709]
 *
 * Only stat and uncorr uncertainties are considered
 */

double AChisqLogNormalStatUncorr::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(log(m_i) - log(t_i)] * V^-1_ij * [(log(m_i) - log(t_i)]
   //! With V containing stat, correlated and uncorrelated uncertainties

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data for uncertainties and inverse covariances
   //const TMatrixD& InvCov = fData->GetInverseCovarianceStatUncorRel();
   const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix(
         {&fData->GetSumErrorMatrix("AS","RelAvTot"),
          &fData->GetSumErrorMatrix("AY","RelAvUnc"),
          &fTheo->GetSumErrorMatrix("AS","RelAvTot"),
          &fTheo->GetSumErrorMatrix("AY","RelAvUnc")});

   double chisq = 0;
   // --- loop over all super-vector data points and calculate chisq
   vector<double> dldxltx(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      dldxltx[x] = log(da[x]/th[x]);
      chisq += pow(dldxltx[x],2) * InvCov[x][x];
      for ( unsigned int y = 0 ; y<x ; y++ ) {
         chisq += 2 * dldxltx[x] * InvCov[x][y] * dldxltx[y];
      }
   }

   if ( std::isnan(chisq)  ){
      cerr<<"Error.  AChisqLogNormal::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqSimple
 *
 * A very simple chisq:
 *
 * chisq = sum_ij [(m_i) - (t_i)]^2 / d^2_ii
 *
 * All uncertainties are treated as uncorrelated and 'additive'
 */

double AChisqSimple::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(m_i) - (t_i)]^2 / d^2_ii
   //! with d containing all uncertainties

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data for uncertainties and inverse covariances
   //const vector<double>& Errors = fData->GetUncertaintyTot();  // old interface
   const vector<double>& dataErrors = fData->GetSumError("AA", "AbsAvTot");
   const vector<double>& theoErrors = fTheo->GetSumError("AA", "AbsAvTot");

   double chisq = 0;
   // --- loop over all super-vector data points and calculate chisq
   vector<double> dldxltx(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += pow(da[x] - th[x], 2) / (pow(dataErrors[x], 2) + pow(theoErrors[x], 2));
   }

   if ( std::isnan(chisq)  ){
      cerr<<"Error.  AChisqSimple::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqLogNormalNuisance
 *
 * A chisq as used in [arxiv:1406.4709]
 *
 * Only alpha_s is a free parameter
 * Corrleated uncertainties are fitted with free nuisance parameter
 *
 */

AChisqLogNormalNuisanceFit::AChisqLogNormalNuisanceFit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo) {
   const map<string, AError>& daErrors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for (auto ie : daErrors) {
      if (ie.second.GetCorrelatedFraction() > 0)
         fFitPar.push_back(ie.first);
   }
   const map<string, AError>& thErrors = fTheo->GetAllErrors();
   for (auto ie : thErrors) {
      if (ie.second.GetCorrelatedFraction() > 0)
         fFitPar.push_back(ie.first);
   }
}

double AChisqLogNormalNuisanceFit::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(log(m_i) - log(t_i) - E_i ] * V^-1_ij * [(log(m_i) - log(t_i) - E_j] + sum_k b_k
   //! With V containing stat and uncorrelated uncertainties
   //! and E_i = sum_k D_i,k, with D_i,k being the corrlated uncertainty of error k

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'data' and 'theory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data for uncertainties and inverse covariances
   //const TMatrixD& InvCov = fData->GetInverseCovarianceStatUncorRel();
   const TMatrixDSym& InvCov = AInvMatrices::Instance()->GetInvMatrix(
      {&fData->GetSumErrorMatrix("AS","RelAvTot"),
       &fData->GetSumErrorMatrix("AY","RelAvUnc"),
       &fTheo->GetSumErrorMatrix("AS","RelAvTot"),
       &fTheo->GetSumErrorMatrix("AY","RelAvUnc")});

   double chisq = 0;
   // --- get nuisance parameters from errors
   map<string, const AError*> errors;
   for (auto& err : fData->GetAllErrors()) {
      errors[err.first] = &err.second;
   }
   for (auto& err : fTheo->GetAllErrors()) {
      // FIXME: errors with same name in theory as in data get overwritten!
      errors[err.first] = &err.second;
   }
   vector<double> Ei(th.size());
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         // const vector<double>& eup = errors.at(fFitPar[k]).GetErrorCorrRelUp();
         // const vector<double>& edn = errors.at(fFitPar[k]).GetErrorCorrRelDn();
         // Ei[i] += (eup[i]-edn[i])/2.*p[k] /* + (eup[i]+edn[i])/2.*p[k]*p[k]*/;
         //const vector<double>& e = errors.at(fFitPar[k]).GetErrorCorrRelAvg();  // old interface
         const vector<double>& e = errors.at(fFitPar[k])->GetError("RelAvCor");

         Ei[i] += e[i]*p[k] /* + (eup[i]+edn[i])/2.*p[k]*p[k]*/;
      }
   }
   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double dldxltxEix = log(da[x]/th[x]) - Ei[x];
      for ( unsigned int y = 0 ; y<th.size() ; y++ ) {
         //chisq += (log(da[x])-log(th[x])-Ei[x]) * InvCov[x][y] * (log(da[y])-log(th[y])-Ei[y]);
         chisq += dldxltxEix * InvCov[x][y] * (log(da[y]/th[y])-Ei[y]);
      }
   }

   if ( std::isnan(chisq)  ){
      error["DoEval"]<<"Chisq is 'nan'."<<endl;
      //exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqLogNormalNuisance
 *
 * A chisq as used in [arxiv:1406.4709]
 *
 * Only alpha_s is a free parameter
 * Nuisance parameters of corrleated uncertainties are calculated analytically
 *
 */

double AChisqLogNormalNuisance::DoEval(const double *p) const {
   //! Calculate chisq
   //! chisq = sum_ij [(log(m_i) - log(t_i) - E_i ] * V^-1_ij * [(log(m_i) - log(t_i) - E_j] + sum_k b_k
   //! With V containing stat and uncorrelated uncertainties
   //! and E_i = sum_k D_i,k, with D_i,k being the corrlated uncertainty of error k
   //! b_k are calculated analytically from 'relative uncertainties'!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) SET_ANY(fFitPar[ipar],p[ipar],0);
   // --- get 'data' and 'theory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");
   //const TMatrixD& VInv = fData->GetInverseCovarianceStatUncorRel();
   const TMatrixDSym& VInv = AInvMatrices::Instance()->GetInvMatrix(
         {&fData->GetSumErrorMatrix("AS","RelAvTot"),
          &fData->GetSumErrorMatrix("AY","RelAvUnc"),
          &fTheo->GetSumErrorMatrix("AS","RelAvTot"),
          &fTheo->GetSumErrorMatrix("AY","RelAvUnc")});

   // --- nuisance parameters
   // --- get nuisance parameters from errors
   map<string, const AError*> errors;
   for (auto& err : fData->GetAllErrors()) {
      errors[err.first] = &err.second;
   }
   for (auto& err : fTheo->GetAllErrors()) {
      // FIXME: errors with same name in theory as in data get overwritten!
      errors[err.first] = &err.second;
   }
   vector<const vector<double>* > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second->GetCorrelatedFraction() > 0 ) {
	 //gerr.push_back(&ie.second.GetErrorCorrRelAvg());  // old interface
    gerr.push_back(&ie.second->GetError("RelAvCor"));
	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
      }
   }
   vector<double> dmt(da.size());
   //for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i] = log(fabs(da[i])) - log(fabs(th[i]));
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i] = log(da[i]/th[i]);
   vector<double> b = CalcNuisance(gerr,VInv,dmt);

   // --- calc chisq
   double chisq = 0;
   for ( unsigned int e = 0 ; e<gerr.size() ; e++)  chisq += b[e]*b[e]; // sum b_i
   vector<double> Ei(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) { // shifts
      for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
         Ei[x] += b[e] * gerr[e]->at(x);
      }
   }
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) { // rT V r
      double cx = dmt[x] - Ei[x];
      for ( unsigned int y = 0 ; y<=x ; y++ ) {
         double cy = cx * VInv(x,y) * (dmt[y] - Ei[y]);;
         chisq += cy;
         if (x!=y) chisq += cy;
      }
   }

   if ( std::isnan(chisq)  ){
      error["DoEval"]<<"Chisq is 'nan'."<<endl;
      AChisqBase::PrintNuisanceParameters(true);
      //exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqHERAFitterDefaultFit
 *
 * A chisq as used HERAFitter by default, following arXiv:1410:4412
 *
 * This chisq ignores correlation matrices.
 * Nuisance parameters are determined by minuit
 *
 */

AChisqHERAFitterDefaultFit::AChisqHERAFitterDefaultFit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& errors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqHERAFitterDefaultFit::DoEval(const double *p) const {
   //! Calculate chisq following arXiv:1410:4412
   //!
   //!                       [d_i - t_i(1- sum_k(g_ik b_k) ) ]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                  D_i,unc^2 t_i^2 + D_i,stat^2 t_i m_i (1-sum_k g_ik b_k)
   //!
   //! Correlations are ignored: only fully correlated or uncorrelated uncertainties are considered.
   //! If a correlation matrix is specified, the uncertainty is considered as uncorrelated.

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();

   double chisq = 0;
   // --- nuisance parameters
   vector<double> Ei(th.size()); // E_i = sum_j g_ij b_j
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         //const vector<double>& gk = errors.at(fFitPar[k]).GetErrorCorrRelAvg();  // old interface
         const vector<double>& gk = errors.at(fFitPar[k]).GetError("RelAvCor");
         Ei[i] += gk[i]*p[k] /* + (eup[i]+edn[i])/2.*p[k]*p[k]*/;
      }
   }

   //const vector<double>& dunc  = fData->GetUncertaintyUncorrRel();  // old interface
   const vector<double>& dunc  = fData->GetSumError("AY", "RelAvUnc");

   //const vector<double>& dmat  = fData->GetUncertaintyMatRel();  // old interface (left in, for now)
   //const vector<double>& dmat  = fData->GetSumError("??", "???????????");  // FIXME: this is not accessible via new interface

   //const vector<double>& dstat = fData->GetUncertaintyStatRel();  // old interface
   const vector<double>& dstat = fData->GetSumError("AS", "RelAvTot");

   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x] - th[x]*(1-Ei[x]);
      //cout<<"d,t\t"<<da[x]<<"\t"<<th[x]<<endl;
      double denom = (dunc[x]*dunc[x])*th[x]*th[x] + dstat[x]*dstat[x]*th[x]*da[x]*(1-Ei[x]);
      //double denom = (dunc[x]*dunc[x]+dmat[x]*dmat[x])*th[x]*th[x] + dstat[x]*dstat[x]*th[x]*da[x];
      chisq += num*num/denom;
      //cout<<"A: "<<0<<"\tM: "<< dunc[x]*dunc[x]<<"\tP: "<<dstat[x]*dstat[x]<<"\tda="<<da[x]<<"\tth="<<th[x]<<endl;
  }

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterDefaultFit::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqHERAFitterDefault
 * A chisq as used HERAFitter by default, following arXiv:1410:4412
 * This chisq ignores correlation matrices.
 * Nuisance parameter are calculated
 */
double AChisqHERAFitterDefault::DoEval(const double *p) const {
   //! Calculate chisq following arXiv:1410:4412
   //!
   //!                       [d_i - t_i(1- sum_k(g_ik b_k) ) ]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                  D_i,unc^2 t_i^2 + D_i,stat^2 t_i m_i (1-sum_k g_ik b_k)
   //!
   //! Matrix-type correlations are considered only for stat-uncertainties, and otherwise ignored.
   //!
   using namespace AlposTools;

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();

   // --- calculate nuisance
   // --- 'constants' for all iterations
   vector<vector<double> > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 ) {
         //gerr.push_back(ie.second.GetErrorCorrRelAvg());  // old interface
         gerr.push_back(ie.second.GetError("RelAvCor"));
         if ( ie.second.GetIsMult() ) gerr.back()*=th;// multiplicative
	 else {  // additive
	    gerr.back()*=da;
	    for ( unsigned int x = 0 ; x<th.size() ; x++ ) gerr.back()[x]*=-1.;
	 }
	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
      }
   }

   //vector<double> Eunc = fData->GetUncertaintyUncorrRel();  // old interface
   vector<double> Eunc = fData->GetSumError("AY", "RelAvUnc");

   //vector<double> Emat = fData->GetUncertaintyMatRel();  // old interface (left in, for now)
   //vector<double> Emat = fData->GetSumError("??", "???????????");  // FIXME: this is not accessible via new interface

   //const vector<double>& dstat = fData->GetUncertaintyStatRel();  // old interface
   const vector<double>& dstat = fData->GetSumError("AS", "RelAvTot");

   Eunc *= th; // ('linear' rescale)
   //Emat *= th; // ('linear' rescale)
   // --- things to keep (b, Ei, Vinv, d-t)
   vector<double> b(gerr.size());
   vector<double> Ei0(th.size()); // E_i = sum_j g_ij b_j
   vector<double> VDiag(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      // cout<<x<<"\tdt="<<da[x]<<"\tth="<<th[x]<<"\tstat="<<dstat[x]<<"\tEunc="<<Eunc[x]<<"\tEUncRel="<<fData->GetUncertaintyUncorrRel()[x]<<"\tsize: "<< fData->GetUncertaintyUncorrRel().size()<<"\tname: "<<fData->GetAlposName()<<endl;
      VDiag[x] = Eunc[x]*Eunc[x] + dstat[x]*dstat[x]*th[x]*da[x];
   }
   vector<double> dmt(da);
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]-=th[i];
   // --- log term
   vector<double> logt(th.size());
   const bool bUseLogTerm = false;//true;
   if ( bUseLogTerm ){
      //const vector<double> Dunc = fData->GetUncertaintyUncorr();  // old interface
      const vector<double> Dunc = fData->GetSumError("AY", "AbsAvUnc");

      //const vector<double> Dmat = fData->GetUncertaintyMat();  // old interface (left in, for now)
      //const vector<double> Dmat = fData->GetSumError("??", "???????????");  // FIXME: this is not accessible via new interface

      //const vector<double>& Dstat = fData->GetUncertaintyStat();  // old interface
      const vector<double>& Dstat = fData->GetSumError("AS", "AbsAvTot");

      for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
	 if ( VDiag[x] > 0 )
	    logt[x] = log( VDiag[x] / (Dunc[x]*Dunc[x]+Dstat[x]*Dstat[x]) );
      }
   }
   // --- loop
   const int nIter = 1; // >= 2
   for ( int nn = 0 ; nn<nIter ;nn++ ) {
      b = CalcNuisance(gerr,VDiag,dmt);
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) { // keep Ei
         Ei0[x]=0;
         for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
            Ei0[x] += b[e] * gerr[e][x];
         }
      }
      // for ( unsigned int x = 0 ; x<th.size() ; x++ ) // recalculate
      //    VDiag[x] = Eunc[x]*Eunc[x] + Emat[x]*Emat[x] + dstat[x]*dstat[x]*(th[x]+Ei0[x])*da[x];
   }


   // --0 calculate actual chisq
   double chisq = 0;
   // --- nuisance parameters
   for ( unsigned int e = 0 ; e<gerr.size() ; e++)  chisq += b[e]*b[e];
   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += pow((dmt[x] - Ei0[x]),2) / VDiag[x];
      // --- log term
      if ( bUseLogTerm ) chisq += logt[x];
   }
   // printout
   if ( bUseLogTerm ) {
      double sumlog = 0;
      double sumb = 0;
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) sumlog += logt[x];
      for ( unsigned int e = 0 ; e<gerr.size() ; e++)  sumb += b[e]*b[e];
      cout<<" Chisq:  "<<chisq<<"\t\tNuisance parameters:  "<<sumb<<"\t\t Log term:  "<<sumlog<<endl;
   }

   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    cout<<"x="<<x<<"\tda="<<da[x]<<"\tth="<<th[x]<<"\tEi0="<<Ei0[x]<<"\t1/V="<<1/VDiag[x]<<endl;
   // }
   // exit(1);

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterDefault::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}


//____________________________________________________________________________________ //
/**
 *
 * AChisqHERAFitterLogDefault
 * A chisq as used HERAFitter by default, following arXiv:1410:4412
 * This chisq ignores correlation matrices.
 * Nuisance parameter are calculated
 */
double AChisqHERAFitterLogDefault::DoEval(const double *p) const {
   //! Calculate chisq following arXiv:1410:4412
   //!
   //!                       [d_i - t_i(1- sum_k(g_ik b_k) ) ]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                  D_i,unc^2 t_i^2 + D_i,stat^2 t_i m_i (1-sum_k g_ik b_k)
   //!
   //! Matrix-type correlations are considered only for stat-uncertainties, and otherwise ignored.
   //!
   using namespace AlposTools;


   double chisq = AChisqHERAFitterDefault::DoEval(p);

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");


   // --- calculate nuisance
   // --- 'constants' for all iterations
// <<<<<<< .mine
// =======
//    vector<vector<double> > gerr;
//    for ( const auto& ie : errors ) {
//       if ( ie.second.GetCorrelatedFraction() > 0 ) {
//          //gerr.push_back(ie.second.GetErrorCorrRelAvg());
//          gerr.push_back(ie.second.GetError("RelAvCor"));
//          gerr.back()*=th;
// 	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
//       }
//    }
// >>>>>>> .r136
   //vector<double> Eunc = fData->GetUncertaintyUncorrRel();  // old interface
   vector<double> Eunc = fData->GetSumError("AY", "RelAvUnc");
   //vector<double> Emat = fData->GetUncertaintyMatRel();
   //const vector<double>& dstat = fData->GetUncertaintyStatRel();  // old interface
   const vector<double>& dstat = fData->GetSumError("AS", "RelAvUnc");
   Eunc *= th; // ('linear' rescale)
   //Emat *= th; // ('linear' rescale)
   // --- things to keep (b, Ei, Vinv, d-t)
   vector<double> VDiag(th.size());
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      // cout<<x<<"\tdt="<<da[x]<<"\tth="<<th[x]<<"\tstat="<<dstat[x]<<"\tEunc="<<Eunc[x]<<"\tEUncRel="<<fData->GetUncertaintyUncorrRel()[x]<<"\tsize: "<< fData->GetUncertaintyUncorrRel().size()<<"\tname: "<<fData->GetAlposName()<<endl;
      VDiag[x] = Eunc[x]*Eunc[x] + dstat[x]*dstat[x]*th[x]*da[x];
   }
   vector<double> dmt(da);
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]-=th[i];
   // --- log term
   vector<double> logt(th.size());
   const bool bUseLogTerm = true;
   if ( bUseLogTerm ){
      //const vector<double>& Dunc  = fData->GetUncertaintyUncorr();  // old interface
      const vector<double>& Dunc  = fData->GetSumError("AY", "AbsAvUnc");
      //const vector<double>& Dmat  = fData->GetUncertaintyMat();  // old interface (left in, for now)
      //const vector<double>& Dmat = fData->GetSumError("??", "???????????");  // FIXME: this is not accessible via new interface
      //const vector<double>& Dstat  = fData->GetUncertaintyStat();  // old interface
      const vector<double>& Dstat  = fData->GetSumError("AS", "AbsAvTot");
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
	 if ( VDiag[x] > 0 )
	    logt[x] = log( VDiag[x] / (Dunc[x]*Dunc[x]+Dstat[x]*Dstat[x]) );
      }
   }

   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      chisq += logt[x]; // --- log term
   }
   // printout
   if ( bUseLogTerm ) {
      double sumlog = 0;
      double sumb = 0;
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) sumlog += logt[x];
      for ( unsigned int e = 0 ; e<fNuisance.size() ; e++)  sumb += fNuisance[e]*fNuisance[e];
      cout<<" Chisq:  "<<chisq<<"\t\tNuisance parameters:  "<<sumb<<"\t\t Log term:  "<<sumlog<<endl;
   }

   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    cout<<"x="<<x<<"\tda="<<da[x]<<"\tth="<<th[x]<<"\tEi0="<<Ei0[x]<<"\t1/V="<<1/VDiag[x]<<endl;
   // }
   // exit(1);

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterLogDefault::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqHERAFitterDefaultMatrix
 * A chisq as used HERAFitter by default, following arXiv:1410:4412
 * This chisq consideres statistical correlations
 * Nuisance parameter are calculated analytically in an iterative method (nIter)
 */
double AChisqHERAFitterDefaultMatrix::DoEval(const double *p) const {
   //! Calculate chisq following arXiv:1410:4412, similar to
   //!
   //!                       [d_i - t_i(1- sum_k(g_ik b_k) ) ]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                      D_i,unc^2 t_i^2 + D_i,stat^2 t_i m_i
   //!
   //! but matrix-type correlations are considered only for stat-uncertainties.
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();
   // const vector<double>& dunc  = fData->GetUncertaintyUncorrRel();
   // const vector<double>& dmat  = fData->GetUncertaintyMatRel();
   // const vector<double>& dstat = fData->GetUncertaintyStatRel();

   // --- calculate nuisance
   // --- first iteration
   /*
   TMatrixDSym Vstat = Data()->GetCovarianceStatRel();
   vector<double> Eunc = fData->GetUncertaintyUncorrRel();
   vector<double> Emat = fData->GetUncertaintyMatRel();
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      Eunc[x] *= th[x] ;
      Emat[x] *= th[x] ;
   }
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      for ( unsigned int y = 0 ; y<=x ; y++ ) {
         double sq = sqrt(th[x]*da[y])*sqrt(th[y]*da[x]);
         double Vxy = Vstat(x,y);
         Vstat(x,y) = Vxy * sq;
         Vstat(y,x) = Vstat(x,y);
      }
      Vstat(x,x) += Eunc[x]*Eunc[x];
      Vstat(x,x) += Emat[x]*Emat[x];
   }
   TMatrixD VInv = AlposTools::InvertChol(Vstat);

   // -- using 'no matrix _instead_'
   // TMatrixDSym VInv(th.size());
   // vector<double> Estat = fData->GetUncertaintyStatRel();
   // vector<double> Eunc = fData->GetUncertaintyUncorrRel();
   // vector<double> Emat = fData->GetUncertaintyMatRel();
   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    Estat[x] *= sqrt(th[x]*da[x]) ;
   //    Eunc[x] *= th[x] ;
   //    Emat[x] *= th[x] ;
   //    VInv(x,x) += 1./(Emat[x]*Emat[x]+Eunc[x]*Eunc[x]+Estat[x]*Estat[x]);
   // }
   // -- --

   vector<vector<double> > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 ) {
         gerr.push_back(ie.second.GetErrorCorrRelAvg());
         for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
            gerr.back()[x]*=th[x];
         }
      }
   }
   vector<double> dmt(th.size());
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]=da[i]-th[i];
   vector<double> b = CalcNuisance(gerr,VInv,dmt);
   */


   // --- 'constants' for all iterations
   vector<vector<double> > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 ) {
         //gerr.push_back(ie.second.GetErrorCorrRelAvg());  // old interface
         gerr.push_back(ie.second.GetError("RelAvCor"));
	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
         for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
            gerr.back()[x]*=th[x];
         }
      }
   }
   //vector<double> Eunc = fData->GetUncertaintyUncorrRel();  // old-interface
   vector<double> Eunc = fData->GetSumError("AY", "RelAvUnc");
   //vector<double> Emat = fData->GetUncertaintyMatRel();
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      Eunc[x] *= th[x] ; // ('linear' rescale)
      //Emat[x] *= th[x] ; // ('linear' rescale)
   }
   // --- things to keep for chisq calculation (b, Ei, Vinv, d-t)
   vector<double> b(gerr.size());
   vector<double> Ei0(th.size()); // E_i = sum_j g_ij b_j
   TMatrixD VInv(th.size(),th.size());
   vector<double> dmt(th.size());
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]=da[i]-th[i];
   // loop
   const int nIter = 2; // >= 2
   for ( int nn = 0 ; nn<nIter ;nn++ ) {
      //TMatrixDSym Vstat = Data()->GetCovarianceStatRel();  // old-interface
      TMatrixDSym Vstat = Data()->GetSumErrorMatrix("AS", "RelAvTot");
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
         for ( unsigned int y = 0 ; y<=x ; y++ ) {
            //double sq = sqrt((th[x]+Ei0[x])*da[y])*sqrt((th[y]+Ei0[x])*da[x]);
            double Vxy = Vstat(x,y);
            double sq = sqrt((th[x]+Ei0[x])*da[y])*sqrt((th[y]+Ei0[y])*da[x]); // 'poisson' rescale
            Vstat(x,y) = Vxy * sq; // 'poisson' rescale
            Vstat(y,x) = Vstat(x,y);
         }
         Vstat(x,x) += Eunc[x]*Eunc[x];
         //Vstat(x,x) += Emat[x]*Emat[x];
      }
      TMatrixDSym VInvI = AlposTools::InvertChol(Vstat);
      b = CalcNuisance(gerr,VInvI,dmt);
      if ( nn==nIter-1) { // keep inverse matrix
         for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
            for ( unsigned int y = 0 ; y<th.size() ; y++ ) {
               VInv(x,y)=VInvI(x,y);
            }
         }
      }
      for ( unsigned int x = 0 ; x<th.size() ; x++ ) { // keep Ei
         Ei0[x]=0;
         for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
            Ei0[x] += b[e] * gerr[e][x];
         }
      }
   }


   // // --- second iteration with nuisance parameters
   // vector<double> Ei0(th.size()); // E_i = sum_j g_ij b_j
   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
   //    Ei0[x] += b[e] * gerr[e][x];
   //    }
   // }
   // //
   // TMatrixD Vstatb = Data()->GetCovarianceStatRel(); // A=V^-1
   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    for ( unsigned int y = 0 ; y<=x ; y++ ) {
   //    double sq = sqrt((th[x]+Ei0[x])*da[y])*sqrt((th[y]+Ei0[y])*da[x]);
   //    double Vxy = Vstat(x,y);
   //    Vstatb(x,y) = Vxy * sq;
   //    Vstatb(y,x) = Vstat(y,x)
   //    }
   //    Vstat(x,x) += Eunc[x]*Eunc[x];
   //    Vstat(x,x) += Emat[x]*Emat[x];
   // }
   // VInv = AlposTools::InvertChol(Vstat);
   // b = CalcNuisance(gerr,VInv,dmt);


   // --0 calculate actual chisq
   double chisq = 0;
   // --- nuisance parameters
   for ( unsigned int e = 0 ; e<gerr.size() ; e++)  chisq += b[e]*b[e];

   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double cx = (dmt[x] - Ei0[x]);
      for ( unsigned int y = 0 ; y<=x ; y++ ) {
         double cy = cx * VInv(x,y) * (dmt[y] - Ei0[y]);;
         chisq += cy;
         if (x!=y) chisq += cy;
      }
   }


   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterDefaultMatrix::DoEval(). Chisq is 'nan'. Returning 1e7"<<endl;
      return 1.e7;
      exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqHERAFitterFull
 *
 * A chisq as used in HERAFitter, following arXiv:1410:4412 and HERAFitter Manual
 *
 * This chisq ignores correlation matrices.
 * This chisq requires the 'Nature' of the error as 'M', 'A' or 'P', denoting
 * the multiplicative, additive, or poissonian nature.
 *
 */

AChisqHERAFitterFull::AChisqHERAFitterFull(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& errors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqHERAFitterFull::DoEval(const double *p) const {
   //! Calculate chisq following the prescription from HERAFitter:
   //!
   //! Consider following options for Nature for each source of uncertainty:
   //! Let D be the 'relative' uncertainty
   //!
   //! Multiplicative ('M'):
   //!    UncorErr  =  D * t
   //!    CorrErr Scaling rule:  t -> t (1 - Db )
   //! Additive ('A'):
   //!    UncorErrAbsErr  =  D * m
   //!    CorrErr Scaling rule:  m -> m (1 - Db )
   //! Poisson ('P'):
   //!    UncorErrAbsErr  =  D * sqrt(mt)
   //!    CorrErr Scaling rule:  t -> t sqrt(1 - Db )
   //!                           m -> m sqrt(1 - Db )
   //!
   //! The uncorrelated uncertainties are not rescaled !!
   //!
   //!                                    [~d_i - ~t_i ) ]^2
   //! chisq ~~  sum_i ------------------------------------------------------- + sum_k b_k^2
   //!                  D_P,unc^2 d_i ~t_i + D_M,unc^2 t_i^2 + D_A,unc^2 d_i^2
   //!
   //! where ~t and ~d are the bias-corrected/scaled values of t and d respectively.
   //!
   //! Correlations are ignored: only fully correlated or uncorrelated fractions of the uncertainties are considered.
   //! If a correlation matrix is specified, the uncertainty is considered as uncorrelated.
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   vector<double> th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   vector<double> da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();

   double chisq = 0;
   // --- correlated uncertainites: nuisance parameters
   vector<double> Edi(th.size()); // E_i = sum_j g_ij b_j
   vector<double> Eti(th.size()); // E_i = sum_j g_ij b_j
   //const vector<double>& dstat = fData->GetUncertaintyStatRel();
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         //const vector<double>& gk = errors.at(fFitPar[k]).GetErrorCorrRelAvg();  // old interface
         const vector<double>& gk = errors.at(fFitPar[k]).GetError("RelAvCor");
         string nat = errors.at(fFitPar[k]).GetNature();
         if ( nat != "A" && nat != "M" && nat != "P" ) {
            cout<<"AChisqHERAFitterFull::DoEval. Cannot understand 'nature' of error '"<<errors.at(fFitPar[k]).GetErrorName()<<"', which is: "<<nat<<endl;
            cout<<"    Treat it as multiplicative ('M')."<<endl;
            nat = "M";
            //errors.at(fFitPar[k]).SetNature(nat); // keep it !
         }
         if ( nat == "A" )
            Edi[i] += gk[i]*p[k];
         else if ( nat == "M" )
            Eti[i] += gk[i]*p[k];
         else if ( nat == "P" ) {
            Edi[i] += sqrt(gk[i]*p[k]);
            Eti[i] += sqrt(gk[i]*p[k]);
         }
      }
   }

   // --- uncorrelated errors
   vector<double> uncSqA(th.size());
   vector<double> uncSqM(th.size());
   vector<double> uncSqP(th.size());
   using namespace AlposTools;
   for ( const auto& ierr : errors ) {
      const string& nat = ierr.second.GetNature();
      vector<double> uncsq;
      if ( ierr.second.GetCorrelatedFraction() >= 0 && ierr.second.GetCorrelatedFraction() != 1 )
         //uncsq = ierr.second.GetErrorUncorrRelAvg();  // old interface
         uncsq = ierr.second.GetError("RelAvUnc");
      else if ( ierr.second.GetIsMatType() )
         //uncsq = ierr.second.GetErrorRelAvg();  // old interface
         uncsq = ierr.second.GetError("RelAvTot");
      else
         continue; // it's fully correlated
      uncsq *= uncsq;
      for ( unsigned int i = 0 ; i< th.size() ; i++ ) {
         if ( nat == "A" )
            uncSqA[i] += uncsq[i];
         else if ( nat == "M" )
            uncSqM[i] += uncsq[i];
         else if ( nat == "P" )
            uncSqP[i] += uncsq[i];
      }
   }

   // // // --- correct for 'bias' from correlated uncertainties
   // for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
   //    da[x] = da[x]*(1-Edi[x]);
   //    th[x] = th[x]*(1-Eti[x]);
   // }

   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x] - th[x]*(1-Eti[x]);
      double denom = uncSqM[x]*th[x]*th[x];
      denom += uncSqA[x]*da[x]*da[x];
      denom += uncSqP[x]*da[x]*th[x]*(1-Eti[x]);
      //cout<<"A: "<<uncSqA[x]<<"\tM: "<< uncSqM[x]<<"\tP: "<<uncSqP[x]<<"\tda="<<da[x]<<"\tth="<<th[x]<<endl;
      chisq += num*num/denom;
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterFull::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 * A chisq as used in HERAFitter, following arXiv:1410:4412 and HERAFitter Manual
 *
 * This chisq ignores correlation matrices.
 * This chisq requires the 'Nature' of the error as 'M', 'A' or 'P', denoting
 * the multiplicative, additive, or poissonian nature.
 * This improved version rescales all uncertainties accordingly.
 *
 */

AChisqHERAFitterFullImproved::AChisqHERAFitterFullImproved(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& errors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqHERAFitterFullImproved::DoEval(const double *p) const {
   //! Calculate chisq following the prescription from HERAFitter:
   //!
   //! Consider following options for Nature for each source of uncertainty:
   //! Let D be the 'relative' uncertainty
   //!
   //! Multiplicative ('M'):
   //!    UncorErr  =  D * t
   //!    CorrErr Scaling rule:  t -> t (1 - Db )
   //! Additive ('A'):
   //!    UncorErrAbsErr  =  D * m
   //!    CorrErr Scaling rule:  m -> m (1 - Db )
   //! Poisson ('P'):
   //!    UncorErrAbsErr  =  D * sqrt(mt)
   //!    CorrErr Scaling rule:  t -> t sqrt(1 - Db )
   //!                           m -> m sqrt(1 - Db )
   //!
   //!                       [~d_i - ~t_i ) ]^2
   //! chisq ~~  sum_i -------------------------------------------- + sum_k b_k^2
   //!                  D_i,unc^2 x_i^2 + sum_l D_i,stat^2 x_i^2
   //!
   //! where ~t and ~d are the bias-corrected/scaled values of t and d respectively.
   //!
   //! Correlations are ignored: only fully correlated or uncorrelated fractions of the uncertainties are considered.
   //! If a correlation matrix is specified, the uncertainty is considered as uncorrelated.
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   vector<double> th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   vector<double> da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();

   double chisq = 0;
   // --- correlated uncertainites: nuisance parameters
   vector<double> Edi(th.size()); // E_i = sum_j g_ij b_j
   vector<double> Eti(th.size()); // E_i = sum_j g_ij b_j
   //const vector<double>& dstat = fData->GetUncertaintyStatRel();
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         //const vector<double>& gk = errors.at(fFitPar[k]).GetErrorCorrRelAvg();
         const vector<double>& gk = errors.at(fFitPar[k]).GetError("RelAvCor");
         string nat = errors.at(fFitPar[k]).GetNature();
         if ( nat != "A" && nat != "M" && nat != "P" ) {
            cout<<"AChisqHERAFitterFullImproved::DoEval. Cannot understand 'nature' of error '"<<errors.at(fFitPar[k]).GetErrorName()<<"', which is: "<<nat<<endl;
            cout<<"    Treat it as multiplicative ('M')."<<endl;
            nat = "M";
            //errors.at(fFitPar[k]).SetNature(nat); // keep it !
         }
         if ( nat == "A" )
            Edi[i] += gk[i]*p[k];
         else if ( nat == "M" )
            Eti[i] += gk[i]*p[k];
         else if ( nat == "P" ) {
            Edi[i] += sqrt(gk[i]*p[k]);
            Eti[i] += sqrt(gk[i]*p[k]);
         }
      }
   }

   // --- uncorrelated errors
   vector<double> uncSqA(th.size());
   vector<double> uncSqM(th.size());
   vector<double> uncSqP(th.size());
   using namespace AlposTools;
   for ( const auto& ierr : errors ) {
      const string& nat = ierr.second.GetNature();
      vector<double> uncsq;
      if ( ierr.second.GetCorrelatedFraction() >= 0 && ierr.second.GetCorrelatedFraction() != 1 )
         //uncsq = ierr.second.GetErrorUncorrRelAvg();  // old interface
         uncsq = ierr.second.GetError("RelAvUnc");
      else if ( ierr.second.GetIsMatType() )
         //uncsq = ierr.second.GetErrorRelAvg();  // old interface
         uncsq = ierr.second.GetError("RelAvTot");
      else
         continue; // it's fully correlated
      uncsq *= uncsq;
      for ( unsigned int i = 0 ; i< th.size() ; i++ ) {
         if ( nat == "A" )
            uncSqA[i] += uncsq[i];
         else if ( nat == "M" )
            uncSqM[i] += uncsq[i];
         else if ( nat == "P" )
            uncSqP[i] += uncsq[i];
      }
   }

   // // --- correct for 'bias' from correlated uncertainties
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      da[x] = da[x]*(1-Edi[x]);
      th[x] = th[x]*(1-Eti[x]);
   }

   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x] - th[x];
      double denom = uncSqM[x]*th[x]*th[x];
      denom += uncSqA[x]*da[x]*da[x];
      denom += uncSqP[x]*da[x]*th[x];
      //cout<<"A: "<<uncSqA[x]<<"\tM: "<< uncSqM[x]<<"\tP: "<<uncSqP[x]<<"\tda="<<da[x]<<"\tth="<<th[x]<<endl;
      chisq += num*num/denom;
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqHERAFitterFullImproved::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqNuisanceRelFit
 *
 * A simple chisq using nuisance parameters
 * The chisq is less biased than the covariance method, since
 * shifts of correlated uncertainties are applied to the theory preiction
 *
 * This chisq ignores correlation matrices.
 *
 */

AChisqSimpleNuisanceMultFit::AChisqSimpleNuisanceMultFit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& errors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqSimpleNuisanceMultFit::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //!                       [d_i - t_i (1- sum_k(g_ik b_k) )]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                        D_i,unc^2 d_i^2 + D_i,stat^2 d_i^2
   //!
   //! Correlations are ignored: only fully correlated or uncorrelated uncertainties are considered.
   //! If a correlation matrix is specified, the uncertainty is considered as uncorrelated.
   //!
   //! Mind that uncorrelated and correlated uncertainties are treated differently (mult[rel] vs. add[abs]).

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();
   //const vector<double>& dunc  = fData->GetUncertaintyUncorr();  // old-interface
   const vector<double>& dunc  = fData->GetSumError("AY", "AbsAvUnc");
   //const vector<double>& dmat  = fData->GetUncertaintyMat();
   //const vector<double>& dstat = fData->GetUncertaintyStat();  // old-interface
   const vector<double>& dstat = fData->GetSumError("AS", "AbsAvTot");

   double chisq = 0;
   // --- nuisance parameters
   vector<double> Ei(th.size()); // E_i = sum_j g_ij b_j
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         //const vector<double>& gk = errors.at(fFitPar[k]).GetErrorCorrRelAvg();  // old-interface
         const vector<double>& gk = errors.at(fFitPar[k]).GetError("RelAvCor");
         Ei[i] += gk[i]*p[k] /* + (eup[i]+edn[i])/2.*p[k]*p[k]*/;
      }
   }
   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x] - th[x]*(1-Ei[x]);
      double denom = dunc[x]*dunc[x] + dstat[x]*dstat[x];
      chisq += num*num/denom;
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error. AChisqSimpleNuisanceMultFit::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}


//____________________________________________________________________________________ //
/**
 *
 * AChisqNuisanceMult
 *
 * A simple chisq using nuisance parameters
 * The chisq is less biased than the covariance method, since
 * shifts of correlated uncertainties are applied to the theory preiction
 * The nuisance parameters are calculated analytically
 *
 * This chisq ignores correlation matrices.
 *
 */

AChisqNuisanceMult::AChisqNuisanceMult(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   // const map<string,AError>& errors = fData->GetAllErrors();
   // fTheoPar = FitPar;
   // for ( auto ie : errors ) {
   //    if ( ie.second.GetCorrelatedFraction() > 0 )
   //    fFitPar.push_back(ie.first);
   // }
}


double AChisqNuisanceMult::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //! chisq =  [d_i - t_i(1- sum_k(g_ik b_k) ) ]^T V^-1 [d_i - t_i(1- sum_k(g_ik b_k) ) ] + sum_k b_k^2
   //!
   //! Chisq is equivalent to 'simple covariance' chisq with:  chisq = xT V^-1 x
   //!
   //! The nuiscane parameters are calculated by solving following eqation:
   //!     (gVg+E) * b = a
   //! with:
   //!     gVg_ik = G^T_i V^-1_ik G_k [an e*e matrix with e=number of correlated errors]
   //!     E   = e*e unit matrix
   //!     G_i = Correlated Errors
   //!     a   = G^T_i V^-1_ik (d-t)_i
   //!     b   = nuisance parameters
   //!
   //! Covariance matrix holds all 'matrix type and uncorrelated uncertainties'
   //!
   using namespace AlposTools;
   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )  SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- calcute nuisance
   const map<string,AError>& errors = fData->GetAllErrors();
   //const TMatrixD& A = Data()->GetInverseCovarianceStatUncorr(); // A=V^-1
   const TMatrixDSym& A = AInvMatrices::Instance()->GetInvMatrix(
      { &fData->GetSumErrorMatrix("ES","AbsAvTot"),
        &fData->GetSumErrorMatrix("EY","AbsAvUnc")} );

   vector<vector<double> > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 ) {
         //gerr.push_back(ie.second.GetErrorCorrRelAvg());  // old-interface
         gerr.push_back(ie.second.GetError("RelAvCor"));
         gerr.back()*=th;
	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
      }
   }
   // --- calculate nuisance
   vector<double> dmt = th;
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]-=da[i];
   vector<double> b = CalcNuisance(gerr,A,dmt);

   // --0 calculate actual chisq
   double chisq = 0;
   // --- nuisance parameters
   for ( unsigned int e = 0 ; e<gerr.size() ; e++)  chisq += b[e]*b[e];

   // --- loop over all super-vector data points and calculate chisq
   vector<double> Ei(th.size()); // E_i = sum_j g_ij b_j
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
         Ei[x] += b[e] * gerr[e][x];
      }
   }
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double cx = (dmt[x] - Ei[x]);
      for ( unsigned int y = 0 ; y<=x ; y++ ) {
         double cy = cx * A(x,y) * (dmt[y] - Ei[y]);;
         chisq += cy;
         if (x!=y) chisq += cy;
      }
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error. AChisqNuisanceMult::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqNuisanceAdd
 *
 * A chisq using nuisance parameters, which is consistent to the covariance expression.
 * The nuisance parameters are calculated analytically
 *
 */

AChisqNuisanceAdd::AChisqNuisanceAdd(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   // const map<string,AError>& errors = fData->GetAllErrors();
   // fTheoPar = FitPar;
   // for ( auto ie : errors ) {
   //    if ( ie.second.GetCorrelatedFraction() > 0 )
   //    fFitPar.push_back(ie.first);
   // }
}


double AChisqNuisanceAdd::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //! chisq =  [d_i(1- sum_k(g_ik b_k) ) - t_i ]^T V^-1 [d_i(1- sum_k(g_ik b_k) ) - t_i ] + sum_k b_k^2
   //!
   //! Chisq is equivalent to 'simple covariance' chisq with:  chisq = xT V^-1 x
   //!
   //! The nuiscane parameters are calculated by solving following eqation:
   //!     (gVg+E) * b = a
   //! with:
   //!     gVg_ik = G^T_i V^-1_ik G_k [an e*e matrix with e=number of correlated errors]
   //!     E   = e*e unit matrix
   //!     G_i = Correlated Errors
   //!     a   = G^T_i V^-1_ik (d-t)_i
   //!     b   = nuisance parameters
   //!
   //! Covariance matrix holds all 'matrix type and uncorrelated uncertainties'
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ )  SET_ANY(fFitPar[ipar],p[ipar],0);

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- calcute nuisance
   const map<string,AError>& errors = fData->GetAllErrors();
   //const TMatrixD& A = Data()->GetInverseCovarianceStatUncorr(); // A=V^-1
    const TMatrixDSym& A = AInvMatrices::Instance()->GetInvMatrix(
       { &fData->GetSumErrorMatrix("ES","AbsAvTot"),
         &fData->GetSumErrorMatrix("EY","AbsAvUnc")} );

   vector<const vector<double>* > gerr;
   for ( const auto& ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 ) {
	 //gerr.push_back(&ie.second.GetErrorCorrAbsAvg());  // old-interface
	 gerr.push_back(&ie.second.GetError("AbsAvCor"));
	 if ( fNuisance.empty() ) fNuisanceName.push_back(ie.first);
      }
   }
   vector<double> dmt = da;
   for ( unsigned int i = 0 ; i<dmt.size() ; i++ ) dmt[i]-=th[i];
   vector<double> b = CalcNuisance(gerr,A,dmt);

   // --0 calculate actual chisq
   double chisq = 0;
   // --- nuisance parameters
   for ( unsigned int e = 0 ; e<gerr.size() ; e++)  chisq += b[e]*b[e];

   // --- loop over all super-vector data points and calculate chisq
   vector<double> Ei(th.size()); // E_i = sum_j g_ij b_j
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      for ( unsigned int e = 0 ; e<gerr.size() ; e++) {
         Ei[x] += b[e] * gerr[e]->at(x);
      }
   }
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double cx = (dmt[x] - Ei[x]);
      for ( unsigned int y = 0 ; y<=x ; y++ ) {
         double cy = cx * A(x,y) * (dmt[y] - Ei[y]);;
         chisq += cy;
         if (x!=y) chisq += cy;
      }
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error. AChisqNuisanceAdd::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}



//____________________________________________________________________________________ //
/**
 *
 * AChisqSimpleNuisanceAddFit
 *
 * A simple chisq using nuisance parameters
 * The chisq give identical results as the covariance method, if only
 * correlated and uncorrelated uncertainties are specified.
 *
 * This nicely illustrates the bias of the 'simple covariance chisq'
 *
 * This chisq ignores correlation matrices.
 *
 */

AChisqSimpleNuisanceAddFit::AChisqSimpleNuisanceAddFit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& errors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : errors ) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqSimpleNuisanceAddFit::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //!                       [d_i(1- sum_k(g_ik b_k) ) - t_i ]^2
   //! chisq =  sum_i ----------------------------------------------------------- + sum_k b_k^2
   //!                        D_i,unc^2 d_i^2 + D_i,stat^2 d_i^2
   //!
   //! Chisq is equivalent to 'simple' chisq with covariance matrix  chisq = xT V^-1 x
   //!
   //! Correlations are ignored: only fully correlated or uncorrelated uncertainties are considered.
   //! If a correlation matrix is specified, the uncertainty is considered as uncorrelated.
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();
   //const vector<double>& dunc  = fData->GetUncertaintyUncorr();  // old-interface
   const vector<double>& dunc  = fData->GetSumError("AY", "AbsAvUnc");
   //const vector<double>& dmat  = fData->GetUncertaintyMat();
   //const vector<double>& dstat = fData->GetUncertaintyStat();  // old-interface
   const vector<double>& dstat = fData->GetSumError("AS", "AbsAvTot");

   double chisq = 0;
   // --- nuisance parameters
   vector<double> Ei(th.size()); // E_i = sum_j g_ij b_j
   for ( unsigned int k = fTheoPar.size() ; k<fFitPar.size() ;k++ ) {
      chisq += p[k]*p[k];
      for ( unsigned int i = 0 ; i<th.size() ; i++ ) {
         //const vector<double>& gk = errors.at(fFitPar[k]).GetErrorCorrRelAvg();  // old-interface
         const vector<double>& gk = errors.at(fFitPar[k]).GetError("RelAvCor");
         Ei[i] += gk[i]*p[k] /* + (eup[i]+edn[i])/2.*p[k]*p[k]*/;
      }
   }
   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x]*(1-Ei[x]) - th[x];
      double denom = dunc[x]*dunc[x] + dstat[x]*dstat[x];
      chisq += num*num/denom;
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error. AChisqSimpleNuisanceAddFit::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}




//____________________________________________________________________________________ //
/**
 *
 * AChisqD0Fit
 *
 *
 * A chisq as used by D0 alpha_s fits from jets and in ATLAS dijet azimuthal decorrelations
 *
 * This chisq takes exp. and theo. uncertainties into account.
 *
 * The nuisance parameter (for both: theoretical and exp. uncertainties) are fitted
 *
 *
 * More details in CDS: ATL-COM-PHYS-2012-1647, appendix I (restriced ATLAS internal)
 *
 */

AChisqD0Fit::AChisqD0Fit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   const map<string,AError>& daErrors = fData->GetAllErrors();
   fTheoPar = FitPar;
   for ( auto ie : daErrors) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
   const map<string,AError>& thErrors = fTheo->GetAllErrors();
   for ( auto ie : thErrors) {
      if ( ie.second.GetCorrelatedFraction() > 0 )
         fFitPar.push_back(ie.first);
   }
}


double AChisqD0Fit::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //! Simplified expression:
   //!
   //!                 [d_i - t_i (1 + sum_k(gth_ik a_k) )/(1 + sum_k(gda_ik b_l) ) ]^2
   //! chisq =  sum_i ----------------------------------------------------------------- + sum_k a_k^2 + sum_l b_l^2  
   //!                            D_i,stat^2 d_i^2 + D_i,unc^2 d_i^2 
   //!
   //! gth are theoretical correlated uncertainties
   //! tda are experimental systematic correlated uncertainties
   //!
   //! In addition, this prescription takes asymmetric uncertainties into account
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   map<string, const AError*> errors;
   for (auto& err : fData->GetAllErrors()) {
      errors[err.first] = &err.second;
   }
   for (auto& err : fTheo->GetAllErrors()) {
      // FIXME: errors with same name in theory as in data get overwritten!
      errors[err.first] = &err.second;
   }

   // const TMatrixDSym& Vs = fData->GetSumErrorMatrix("ES","RelAvTot");
   // const TMatrixDSym& Vu = fData->GetSumErrorMatrix("EY","RelAvUnc");

   const vector<double>& dunc  = fData->GetSumError("EY","AbsAvUnc");
   const vector<double>& dstat = fData->GetSumError("AS","AbsAv"); // consider also 'TS' here.

   double chisq = 0;
   // --- nuisance parameters
   vector<double> Eai(th.size()); // E_i = sum_j gth_ij a_k
   vector<double> Ebi(th.size()); // E_i = sum_j gda_ij b_l
   for (unsigned int k = fTheoPar.size(); k < fFitPar.size(); k++) {
      chisq += p[k] * p[k];
      for (unsigned int i = 0; i < th.size(); i++) {
         //const vector<double>& gkav = daErrors.at(fFitPar[k]).GetError("RelAvCor");
         const vector<double>& gkup = errors.at(fFitPar[k])->GetError("RelUpCor");
         const vector<double>& gkdn = errors.at(fFitPar[k])->GetError("RelDnCor");
         double dij = p[k] * (gkup[i] - gkdn[i]) / 2. + (gkup[i] + gkdn[i]) / 2. * p[k] * p[k];

         if (errors.at(fFitPar[k])->GetIsTheo())  // Theo uncertainties
            Eai[i] += dij;
         else //data uncertainties
            Ebi[i] += dij;

      }
   }

   // --- loop over all super-vector data points and calculate chisq
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double num = da[x] - th[x]*(1+Eai[x])/(1.+Ebi[x]);
      double denom = dunc[x]*dunc[x] + dstat[x]*dstat[x];
      chisq += num*num/denom;
   }

   if ( std::isnan(chisq)  ){
      cout<<"Error. AChisqD0Fit::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
}

/**
 *
 * AChisqD0StatCorrFit
 *
 *
 * A chisq as used by D0 alpha_s fits from jets and in ATLAS dijet azimuthal decorrelations
 *
 * This chisq takes exp. and theo. uncertainties into account.
 *
 * The nuisance parameter (for both: theoretical and exp. uncertainties) are fitted
 *
 *
 * More details in CDS: ATL-COM-PHYS-2012-1647, appendix I (restriced ATLAS internal)
 *
 */

AChisqD0StatCorrFit::AChisqD0StatCorrFit(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo)
{
   fTheoPar = FitPar;

   // -- create nuisance parameters

   // errors stored in AData
   const map<string, AError>& daErrors = fData->GetAllErrors();
   for (auto ie : daErrors) {
      // one nuisance parameter for each correlated systematic uncertainty
      if ((ie.second.GetCorrelatedFraction() > 0) && (!ie.second.GetIsStat())) {
         fFitPar.push_back(ie.first);
      }
   }

   // errors stored in AFuncD:
   const map<string, AError>& thErrors = fTheo->GetAllErrors();
   for (auto ie : thErrors) {
      // one nuisance parameter for each correlated systematic uncertainty
      if ((ie.second.GetCorrelatedFraction() > 0) && (!ie.second.GetIsStat())) {
         fFitPar.push_back(ie.first);
      }
   }
}


double AChisqD0StatCorrFit::DoEval(const double *p) const {
   //! Calculate simple chisq with nuisance parameters
   //!
   //! Simplified expression:
   //!
   //!
   //! chisq =  sum_i sum_l  A_i * (V^-1)_il * A_l + sum_j e_j^2 + sum_k a_k^2
   //!
   //!   A_i = d_i - t_i * (1 + sum_k(gth_ik a_k)) / (1 + sum_k(gda_ij e_j))
   //!
   //! gth are theoretical correlated uncertainties
   //! gda are experimental systematic correlated uncertainties
   //! V is the covariance matrix containing the statistical and uncorrelated
   //! uncertainties.
   //!
   //! In addition, this prescription takes asymmetric uncertainties into account
   //!

   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fTheoPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar], p[ipar], 0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- get data and uncertainties
   map<string, const AError*> errors;
   for (auto& err : fData->GetAllErrors()) {
      errors[err.first] = &err.second;
   }
   for (auto& err : fTheo->GetAllErrors()) {
      // FIXME: errors with same name in theory as in data get overwritten!
      errors[err.first] = &err.second;
   }

   // --- get uncertainties from data and theory objects
   // following errors treated as additive:
   const TMatrixDSym& Vdstat = fData->GetSumErrorMatrix("AS", "AbsAvTotAll");  // statistical errors
   const TMatrixDSym& Vdunc  = fData->GetSumErrorMatrix("AY", "AbsAvUncAll");  // systematic errors (uncorrelated part)
   const TMatrixDSym& Vtstat = fTheo->GetSumErrorMatrix("AS", "AbsAvTotAll");  // statistical errors
   const TMatrixDSym& Vtunc  = fTheo->GetSumErrorMatrix("AY", "AbsAvUncAll");  // systematic errors (uncorrelated part)
   // other errors treated as multiplicative via the nuisance parameters...

   const TMatrixDSym* InvCov = &AInvMatrices::Instance()->GetInvMatrix({&Vdstat, &Vdunc, &Vtstat, &Vtunc});

   double chisq = 0;

   // --- nuisance parameters
   vector<double> Eai(th.size()); // E_i = sum_j gth_ij a_k
   vector<double> Ebi(th.size()); // E_i = sum_j gda_ij b_l
   for (unsigned int k = fTheoPar.size(); k < fFitPar.size(); k++) {
      chisq += p[k] * p[k];
      for (unsigned int i = 0; i < th.size(); i++) {
         //const vector<double>& gkav = daErrors.at(fFitPar[k]).GetError("RelAvCor");
         const vector<double>& gkup = errors.at(fFitPar[k])->GetError("RelUpCor");
         const vector<double>& gkdn = errors.at(fFitPar[k])->GetError("RelDnCor");
         double dij = p[k] * (gkup[i] - gkdn[i]) / 2. + (gkup[i] + gkdn[i]) / 2. * p[k] * p[k];

         if (errors.at(fFitPar[k])->GetIsTheo()) {
            // theory uncertainties
            Eai[i] += dij;
         }
         else {
            // data uncertainties
            Ebi[i] += dij;
         }
      }
   }

   // --- loop over all super-vector data points and calculate chisq
   for (unsigned int i = 0; i<th.size(); i++) {
      double A_i = da[i] - th[i] * (1.+Eai[i]) / (1.+Ebi[i]);
      chisq += A_i * A_i * (*InvCov)[i][i];
      for (unsigned int j = 0; j<i; j++) {
         double A_j = da[j] - th[j] * (1.+Eai[j]) / (1.+Ebi[j]);
         chisq += 2 * A_i * A_j * (*InvCov)[i][j];
      }
   }

   if (std::isnan(chisq)) {
      error["DoEval"] << "Chisq is 'nan'." << std::endl;
      exit(626);
   }

   return chisq;
}



//____________________________________________________________________________________ //
/** 
 *
 * AChisqApc
 *
 * 
 *
 * This chisq ignores correlation matrices.
 *
 */

AChisqApc::AChisqApc(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo) 
{
   const map<string,AError>& mErrors = fData->GetAllErrors();
   int iNsystErr = 0;
   int iNda =  data->GetValues().size();
   for(auto ie: mErrors)
     {
       if(ie.second.GetCorrelatedFraction() > 0)
         {
         ++iNsystErr;
         }
     }
   emNuisParam = Eigen::MatrixXd::Zero(iNda, iNsystErr);
   int iCount = 0;
   for(auto ie : mErrors)
    {
      if(ie.second.GetCorrelatedFraction() > 0)
        {
          emNuisParam.col(iCount) = 
            //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0],iNda);  // old-interface
            Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0],iNda);
          ++iCount;
        }
    }
}


double AChisqApc::DoEval(const double *p) const {
  
   
   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      //cout<<" ** AChisq(). SetValue "<<fFitPar[ipar]<<": "<<p[ipar]<<endl;
      SET_ANY(fFitPar[ipar],p[ipar],0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");
   
   Eigen::ArrayXd eaTh = Eigen::ArrayXd::Map(&th[0], th.size());
   Eigen::ArrayXd eaDa = Eigen::ArrayXd::Map(&da[0], th.size());

   // --- get data and uncertainties
   const map<string,AError>& errors = fData->GetAllErrors();
   //const vector<double>& dunc  = fData->GetUncertaintyUncorrRel();  // old-interface
   const vector<double>& dunc  = fData->GetSumError("AA", "RelAvUnc");
   //const vector<double>& dstat = fData->GetUncertaintyStatRel();  // old-interface
   const vector<double>& dstat = fData->GetSumError("AS", "RelAvTot");
   
   Eigen::ArrayXd eaDunc = Eigen::ArrayXd::Map(&dunc[0], dunc.size());
   Eigen::ArrayXd eaDstat = Eigen::ArrayXd::Map(&dstat[0], dstat.size());   
   Eigen::ArrayXd eaVxDiag = eaDunc*eaDunc*eaTh*eaTh + eaDstat*eaDstat*eaTh*eaDa;

   Eigen::VectorXd evX(eaDa.size()+emNuisParam.cols());
   evX << eaDa.matrix(), Eigen::VectorXd::Zero(emNuisParam.cols());
   
   Eigen::MatrixXd emVx = Eigen::MatrixXd::Identity(eaDa.size()+emNuisParam.cols(),eaDa.size()+emNuisParam.cols());
   emVx.topLeftCorner(eaDa.size(),eaDa.size()) = eaVxDiag.matrix().asDiagonal();

   Eigen::VectorXd evF(eaDa.size());
   Apccpp ChisqNuis(evX);
   ChisqNuis.SetPrintFlag(0);
   do { 
     evX = ChisqNuis.evGetX();
     evF = evX.head(eaDa.size()).array() - eaTh*(1 - (emNuisParam*evX.segment(eaDa.size(), emNuisParam.cols())).array());
   } while(ChisqNuis.bApcpp(evF, emVx));

   evX = ChisqNuis.evGetX();

   const map<string,AError>& mErrors = fData->GetAllErrors();
   int iCount = eaDa.size();
   for(auto ie : mErrors) {
     if(ie.second.GetCorrelatedFraction() > 0)
       {
         cout << std::left << std::setw(20) << std::setfill(' ') << ie.first;
         cout << std::left << std::setw(20) << std::setfill(' ') << " = ";
         cout << std::left << std::setw(20) << std::setfill(' ') << evX[iCount];
         cout << endl;
         ++iCount;
       }
   }
   cout << "-------------------------------------------------" << endl;
   double chisq = ChisqNuis.dGetChi2();

   if ( std::isnan(chisq)  ){
      cout<<"Error.  AChisqApc::DoEval(). Chisq is 'nan'."<<endl;
      exit(1);
   }
   return chisq;
  
}


//____________________________________________________________________________________ //
/**
 *
 * ALogLikelihood
 *
 */
//ALogLikelihood::ALogLikelihood(const std::vector<std::string>& FitPar ,AData* data, AFuncD* theo ) : AChisqBase(FitPar,data,theo) {
//   fTheoPar = FitPar;
//}
double ALogLikelihood::DoEval(const double *p) const {
   //! Calculate Log likelihood
   // --- set new theory parameters
   for ( unsigned int ipar = 0 ; ipar < fFitPar.size() ; ipar++ ) {
      SET_ANY(fFitPar[ipar], p[ipar], 0);
   }

   // --- get 'superdata' and 'supertheory' arrays
   const vector<double>& th = Theo()->GetValues();// VALUES_ANY("SuperTheory");
   const vector<double>& da = Data()->GetValues();// VALUES_ANY("SuperData");

   // --- loop over all super-vector data points and calculate chisq
   double prob = 1;
   for (unsigned int i = 0; i<th.size(); i++) {
      prob *= th[i];
   }
   return -2*log(prob);
}





//____________________________________________________________________________________ //
/**
    Calculate Pull values
*/
double APull::CalcMeanStatUncorr() const {
   //! Calculate  pull value, using the stat+uncorr uncertainties
   // --- get 'data' and 'theory' arrays
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovarianceStatUncorr();  // old interface
   const vector<double>& estat = Data()->GetSumError("AS", "AbsAv");
   const vector<double>& eunc  = Data()->GetSumError("AY", "AbsAvUnc");
   double ret = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ret += (da[x]-th[x])/sqrt(estat[x]*estat[x]+eunc[x]*eunc[x]);
   }
   return ret/th.size();
}

double APull::CalcMeanTotErr() const {
//! Calculate pull value, using the stat+uncorr+corr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovariance();  // old interface
   const TMatrixDSym& err = Data()->GetSumErrorMatrix("AA", "AbsAvTot");
   double ret = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ret += (da[x]-th[x])/sqrt(err[x][x]);
   }
   return ret/th.size();
}

double APull::CalcMedianStatUncorr() const {
   //! Calculate meadian of pull distribution, using the stat+uncorr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovarianceStatUncorr();  // old interface
   const vector<double>& estat = Data()->GetSumError("AS", "AbsAv");
   const vector<double>& eunc  = Data()->GetSumError("AY", "AbsAvUnc");
   multiset<double> ps;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ps.insert((da[x]-th[x])/sqrt(estat[x]*estat[x]+eunc[x]*eunc[x]));
   }
   std::multiset<double>::iterator it=ps.begin();
   if (  th.size() % 2 == 1 ) {
      std::advance(it,(th.size()-1)/2);
      return *it;
   }
   else {
      std::advance(it,th.size()/2-1);
      double lo = *it;
      std::advance(it,1);
      double up = *it;
      return (lo+up)/2.;
   }
}

double APull::CalcMedianTotErr() const {
   //! Calculate meadian of pull distribution, using the stat+uncorr+corr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovariance();  // old interface
   const TMatrixDSym& err = Data()->GetSumErrorMatrix("AA", "AbsAvTot");
   multiset<double> ps;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ps.insert((da[x]-th[x])/sqrt(err[x][x]));
   }
   std::multiset<double>::iterator it=ps.begin();
   if (  th.size() % 2 == 1 ) {
      std::advance(it,(th.size()-1)/2);
      return *it;
   }
   else {
      std::advance(it,th.size()/2-1);
      double lo = *it;
      std::advance(it,1);
      double up = *it;
      return (lo+up)/2.;
   }
}

double APull::CalcRMSStatUncorr() const {
   //! Calculate RMS, using the stat+uncorr uncertainties

   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovarianceStatUncorr();  // old interface
   const vector<double>& estat = Data()->GetSumError("AS", "AbsAv");
   const vector<double>& eunc  = Data()->GetSumError("AY", "AbsAvUnc");
   double ret = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ret += pow((da[x]-th[x])/sqrt(estat[x]*estat[x]+eunc[x]*eunc[x]),2);
   }
   return sqrt(ret/th.size());
}
double APull::CalcRMSTotErr() const {
   //! Calculate RMS , using the stat+uncorr+corr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovariance();  // old interface
   const vector<double>& err = Data()->GetSumError("AA", "AbsAvTot");
   double ret = 0;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      ret += pow((da[x]-th[x])/sqrt(err[x]*err[x]),2);
   }
   return sqrt(ret/th.size());
}

double APull::CalcMaxStatUncorr() const {
   //! Calculate RMS, using the stat+uncorr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovarianceStatUncorr();  // old interface
   const vector<double>& estat = Data()->GetSumError("AS", "AbsAv");
   const vector<double>& eunc  = Data()->GetSumError("AY", "AbsAvUnc");
   double ret = -1e7;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double pull = (da[x]-th[x])/sqrt(estat[x]*estat[x]+eunc[x]*eunc[x]);
      if ( pull > ret ) ret = pull;
   }
   return ret;
}

double APull::CalcMaxTotErr() const {
   //! Calculate RMS, using the stat+uncorr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovariance();  // old interface
   //const TMatrixDSym& err = Data()->GetSumErrorMatrix("AA", "AbsAvTot");
   const vector<double>& err = Data()->GetSumError("AA", "AbsAvTot");
   double ret = -1e7;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double pull = (da[x]-th[x])/sqrt(err[x]*err[x]);
      if ( pull > ret ) ret = pull;
   }
   return ret;
}

double APull::CalcMinStatUncorr() const {
   //! Calculate RMS, using the stat+uncorr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovarianceStatUncorr();  // old interface
   const vector<double>& estat = Data()->GetSumError("AS", "AbsAv");
   const vector<double>& eunc  = Data()->GetSumError("AY", "AbsAvUnc");
   double ret = 1e7;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double pull = (da[x]-th[x])/sqrt(estat[x]*estat[x]+eunc[x]*eunc[x]);
      if ( pull < ret ) ret = pull;
   }
   return ret;
}

double APull::CalcMinTotErr() const {
   //! Calculate RMS, using the stat+uncorr uncertainties
   const vector<double>& th = Theo()->GetValues();
   const vector<double>& da = Data()->GetValues();
   //const TMatrixDSym& err = Data()->GetCovariance();  // old interface
   //const TMatrixDSym& err = Data()->GetSumErrorMatrix("AA", "AbsAvTot");
   const vector<double>& err = Data()->GetSumError("AA", "AbsAvTot");
   double ret = 1e7;
   for ( unsigned int x = 0 ; x<th.size() ; x++ ) {
      double pull = (da[x]-th[x])/err[x];
      if ( pull < ret ) ret = pull;
   }
   return ret;
}
