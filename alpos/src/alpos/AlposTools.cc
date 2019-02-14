// DB


/**
 * 
 *  AlposTools
 * 
 *  Useful tools for Alpos
 * 
 * 
 **/

#include "alpos/AlposTools.h"
#include "alpos/Alpos.h"
#include <iostream>
#include <cmath>
#include <TDecompLU.h>
#include <TDecompChol.h>
#include "fastnlotk/speaker.h"
//#include <sys/access.h>
#include <unistd.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixDEigen.h>
//#include <set>


//____________________________________________________________________________________ //
bool AlposTools::IdenticalToNSignificantDigits(double a, double b, unsigned int n) {
   return ((a==b) || (std::fabs(a-b) / std::max(std::fabs(a), std::fabs(b)) < 1.0 / std::pow(10, n)));
}

//____________________________________________________________________________________ //
std::vector<double > AlposTools::CovToCorr(TMatrixDSym& mat) {
   //! Transform a covariance matrix into a correlation matrix.
   //! That is, mat[i][j] = mat[i][j] / sqrt(mat[i][i]) / sqrt(mat[j][j])
   //! Obtain a vector of the error (i.e. sqrt(diagonal))
   using namespace std;
   vector<double> err(mat.GetNcols());
   for ( int x = 0 ; x<mat.GetNcols() ; x++ ) {
      err[x] = sqrt(mat[x][x]);
   }

   for ( int x = 0 ; x<mat.GetNcols() ; x++ ) {
      mat[x][x] = 1; // set diagonal elements exactly to one
      for ( int y = 0 ; y<x ; y++ ) {
	      mat[x][y] = mat[x][y]/(err[x]*err[y]);
	      mat[y][x] = mat[x][y];
      }
   }

   return err;
}


//____________________________________________________________________________________ //
void AlposTools::CovToCovRel(TMatrixDSym& mat, const std::vector<double>& sigma) {
   //! Transform a covariance matrix into a 'relative' covariance matrix.
   //! That is, mat[i][j] = mat[i][j] / sigma[i] / sigma [j]

   using namespace std;

   // check compatibility
   if (mat.GetNrows() != sigma.size()) {
      say::debug["AlposTools::CovToCovRel"]<<"Matrix (size "<<mat.GetNrows()<<" x "<<mat.GetNcols()<<") and vector (size "<<sigma.size()<<") are incompatible."<<endl;
      exit(1);
   }
   for ( int x = 0 ; x<mat.GetNrows() ; x++ ) {
      for ( int y = 0 ; y<mat.GetNcols() ; y++ ) {
         mat[x][y] = mat[x][y]/(sigma[x]*sigma[y]);
      }
   }
}

//____________________________________________________________________________________ //
void AlposTools::CovRelToCov(TMatrixDSym& mat, const std::vector<double>& sigma) {
   //! Transform a 'relative' covariance matrix into a covariance matrix.
   /*!
       That is, mat[i][j] = mat[i][j] * sigma[i] * sigma [j]
    */

   using namespace std;

   // check compatibility
   if (mat.GetNrows() != sigma.size()) {
      say::debug["AlposTools::CovRelToCov"]<<"Matrix (size "<<mat.GetNrows()<<" x "<<mat.GetNcols()<<") and vector (size "<<sigma.size()<<") are incompatible."<<endl;
      exit(1);
   }
   for ( int x = 0 ; x<mat.GetNrows() ; x++ ) {
      for ( int y = 0 ; y<mat.GetNcols() ; y++ ) {
         mat[x][y] = mat[x][y]*(sigma[x]*sigma[y]);
      }
   }
}

//____________________________________________________________________________________ //
void AlposTools::CovRescale(TMatrixDSym& mat, const std::vector<double>& sigma_current, const std::vector<double>& sigma_rescale) {
   //! Rescale a covariance matrix in-place.
   /*!
       This means dividing by the current reference cross sections to obtain a
       'relative' covariance matrix and multiplying that by the new, resacled cross
       sections. That is:

         mat[i][j] = mat[i][j] / (sigma_current[i]*sigma_current[j]) * (sigma_rescale[i]*sigma_rescale[j])
    */

   using namespace std;

   // check compatibility
   if (sigma_rescale.size() != sigma_current.size()) {
      say::error["AlposTools::CovRescale"]<<"The reference cross section vector (size "<<sigma_current.size()<<
                                            ") and the rescaled cross section vector (size "<<sigma_rescale.size()<<
                                            ") are incompatible."<<endl;
      exit(1);
   }
   if (mat.GetNrows() != sigma_current.size()) {
      say::error["AlposTools::CovRescale"]<<"Matrix (size "<<mat.GetNrows()<<" x "<<mat.GetNcols()<<") and vector (size "<<sigma_current.size()<<") are incompatible."<<endl;
      exit(1);
   }

   // do rescaling
   for ( int x = 0 ; x<mat.GetNrows() ; x++ ) {
      for ( int y = 0 ; y<mat.GetNcols() ; y++ ) {
         mat[x][y] = mat[x][y]/(sigma_current[x]*sigma_current[y])*(sigma_rescale[x]*sigma_rescale[y]);
      }
   }
}

//____________________________________________________________________________________ //
void AlposTools::CorrToCov(TMatrixDSym& mat, const std::vector<double>& error){
   //! Transform a correlation matrix into its covariance matrix
   //! That is, mat[i][j] = mat[i][j] * error[i] * error [j]
   //! Error must be given in absolute units

   using namespace std;

   // check compatibility
   using namespace std;
   if (mat.GetNrows() != error.size()) {
      say::error["AlposTools::CorrToCov"]<<"Matrix (size "<<mat.GetNrows()<<" x "<<mat.GetNcols()<<") and error vector (size "<<error.size()<<") are imcompatible."<<endl;
      exit(1);
   }
   for ( int x = 0 ; x<mat.GetNcols() ; x++ ) {
      for ( int y = 0 ; y<x+1 ; y++ ) {
	      if ( x==y && mat[x][y]!=1 )
	         say::warn["AlposTools::CorrToCov"]<<"Diagonal element of correlation matrix should be 1, but is "<<mat[x][y]<<endl;
         mat[x][y] = mat[x][y]*error[x]*error[y];
	      mat[y][x] = mat[x][y];
      }
   }
}


//____________________________________________________________________________________ //
TMatrixD AlposTools::InvertLU(const TMatrixD& mat){
   using namespace std;
   // --- invert matrix
   TVectorD eigVec;
   mat.EigenVectors(eigVec);
   cout << "RADEK condition "  << eigVec.Min() <<" "<< eigVec.Max() << endl;
   TDecompLU lu(mat);
   TMatrixD Inv(mat);
   lu.Invert(Inv);

   // --- verify the inversion, compare T = A-1 * A with unit matrix (within global tolerance)
   TMatrixD test(mat, TMatrixD::kMult, Inv);
   TMatrixD unit(TMatrixD::kUnit, mat);

   double tolerance = Alpos::Current()->Settings()->matInvTol;
   Bool_t ok = VerifyMatrixIdentity(unit, test, false, tolerance);
   if (!ok){
      say::error["AlposTools::InvertLU"]<<"Failed to invert the matrix ["<<mat.GetNcols()<<","<<mat.GetNrows()<<"]" << endl;

      TMatrixD diff(test-unit);
      double minEl=0, maxEl=0;
      for (int i=0; i<diff.GetNrows(); i++) {
         for (int j=0; j<diff.GetNcols(); j++) {
            if (diff[i][j] < minEl) minEl=diff[i][j];
            if (diff[i][j] > maxEl) maxEl=diff[i][j];
         }
      }

      say::info["AlposTools::InvertLU"]<<"Value range of (A^-1.A)-1 is ["<<minEl<<", "<<maxEl<<"]"<<endl;
      if ( (-tolerance > minEl) || (tolerance < maxEl) ) {
         say::info["AlposTools::InvertLU"]<<"Inversion failed because tolerance level (+/-"<<tolerance<<") is within this range."<<endl;
      }

      exit(3);
   }
   return Inv;
}


//____________________________________________________________________________________ //
TMatrixDSym AlposTools::InvertChol(const TMatrixDSym& mat){
   using namespace std;
   // --- invert matrix
   say::info["AlposTools::InvertChol"]<<"Inverting matrix."<<endl;
   TDecompChol lu(mat);
   TMatrixDSym Inv(mat);
   lu.Invert(Inv);

   // --- verify the inversion, compare T = A-1 * A with unit matrix                                                        
   TMatrixD test(mat, TMatrixD::kMult, Inv);
   TMatrixD unit(TMatrixD::kUnit, mat);
   Bool_t ok = VerifyMatrixIdentity(unit, test, false, 1e-8);
   if (!ok){
      say::error["AlposTools::InvertChol"]<<"Failed to invert the matrix." << endl;
      say::info["AlposTools::InvertChol"]<<"Trying LU decomposition."<<endl;
      TMatrixD invLU = AlposTools::InvertLU(mat);
      for ( int x = 0 ; x<invLU.GetNcols() ; x++ ) {
	 for ( int y = 0 ; y<invLU.GetNcols() ; y++ ) {
	    Inv(x,y) = invLU(x,y);
	 }
      }
      //exit(3);
   }
   return Inv;
}

TMatrixDSym AlposTools::InvertChol(std::set<const TMatrixDSym*> mat_collection) {
   //! Invert sum of matrices present in collection

   TMatrixDSym Inv(*(*mat_collection.begin()));
   int i = 0;
   for ( auto mat : mat_collection ) {
      if ( i++ == 0 )
         continue;
      Inv += *(mat);
   }

   return AlposTools::InvertChol(Inv);
}

//____________________________________________________________________________________ //
void AlposTools::CheckFileExit(std::string& file){
   //! check if a file exists and is readable (under linux)
   if ( !CheckFile(file) ) {
      say::error["AlposTools::CheckFileExit"]<<"File '"<<file<<"' does not exist. Maybe the wrong path is specified. Exiting."<<std::endl;
      exit(1);
   }
}

//____________________________________________________________________________________ //
bool AlposTools::CheckFile(std::string& file){
   //! check if a file exists and is readable (under linux)
   if ( -1 != access(file.c_str(), R_OK) ) {
      return true;// found file
   }
   else {
      if ( file.front() == '/' ) {
	 // absolute path is given, but file was not found before there
	 return false;
      }
      else { 
	 // relative path. Looking into ALPOS_DIR
	 std::cout<<"Modifying path: alposDir: "<<Alpos::Current()->Settings()->Alpos_dir<<std::endl;
	 file = Alpos::Current()->Settings()->Alpos_dir + file;
	 std::cout<<"file="<<file<<std::endl;
	 return ( -1 != access(file.c_str(), R_OK) );
      }
   }
   return false;
}



//____________________________________________________________________________________ //
void AlposTools::Calc13partonsFromLico(std::vector<double>& xfx13, const TMatrixD& Inv, const std::vector<double>& xf6, const std::vector<bool>& found ) {
   //! calculate 13 partons from LiCo 
   xfx13.resize(13);
   for ( int ip = 0 ; ip<13 ; ip++ ) xfx13[ip]=0;
   for ( int ip = 0 ; ip<13 ; ip++ ) {
      for ( int i = 0 ; i<13 ; i++ )  {
         xfx13[ip] += Inv(ip,i)*xf6[i];
      }
   }
   for ( int i = 0 ; i<13 ; i++ ) {
      if ( fabs(xfx13[i]) < 1.e-30 ) xfx13[i]=0;
      if ( !found[i] ) xfx13[i]=0;
   }
}


//____________________________________________________________________________________ //
void AlposTools::CalcLicoFrom13partons(std::vector<double>& LiCo, const TMatrixD& def, const std::vector<double>& xfx13 ) {
   //! calculate Lico from 13 partons 
   LiCo.resize(13);
   for ( int ip = 0 ; ip<13 ; ip++ ) LiCo[ip]=0;
   for ( int i = 0 ; i<13 ; i++ ) {
      for ( int ip = 0 ; ip<13 ; ip++ )  {
	 LiCo[i] += def(ip,i)*xfx13[ip];
      }
   }
}

//____________________________________________________________________________________ //
void AlposTools::CalcLicoFrom13partons(std::vector<double>& LiCo, const std::vector<double>& def, const std::vector<double>& xfx13 ) {
   //! calculate Lico from 13 partons 
   LiCo.resize(13);
   for ( int ip = 0 ; ip<13 ; ip++ ) LiCo[ip]=0;
   for ( int i = 0 ; i<13 ; i++ ) {
      for ( int ip = 0 ; ip<13 ; ip++ )  {
	 LiCo[i] += def[i*13+ip]*xfx13[ip];
      }
   }
}

//____________________________________________________________________________________ //
std::vector<double> AlposTools::LicoLhaToApfelxx(const std::vector<double>& LiCoLha ) {
   //! calculate PDF linear combination as input to Apfel++
   //! Input vector denots PDFs like in LHAPDF: tb...g d u s c b t (0...6...12)
   
   // q^+ = q + qbar
   // q^- = q - qbar
   static const std::vector<std::vector<double> > M{
      {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}, // gluon
      {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}, // SIGMA = \sum_q q^+
      { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}, // VALENCE = \sum_q q^-
      {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}, // T3  = u^+ - d^+
      {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}, // V3  = u^- - d^-
      {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}, // T8  = u^+ + d^+ - 2 s^+
      {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}, // V8  = u^- + d^- - 2 s^-
      {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}, // T15 = u^+ + d^+ + s^+ - 3 c^+
      {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}, // V15 = u^- + d^- + s^- - 3 c^-
      {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}, // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+
      {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}, // V24 = u^- + d^- + s^- + c^- - 4 b^-
      { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}, // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
      {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,}  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-
   };
   std::vector<double> apfl(13);
   apfl[0] = LiCoLha[6];
   for ( int i = 1 ; i<13; i++ ) {
      for ( int j = 0 ; j<13; j++ ) {
	 apfl[i] += M[i][j]*LiCoLha[j];
      }
   }
   return apfl;
}

//____________________________________________________________________________________ //
std::vector<double> AlposTools::LicoApfelxxToLha(const std::vector<double>& LiCoApfl ) {
   //! calculate PDF linear combination as input to Apfel++
   //! Input vector denots PDFs like in LHAPDF: tb...g d u s c b t (0...6...12)
   
   // q^+ = q + qbar
   // q^- = q - qbar
   static const std::vector<std::vector<double> > M{
      { 0, 1./12, -1./12,     0,     0,     0,     0,     0,   0,  0,  0,   -1./12,  1./12 }, 
      { 0, 1./12, -1./12,     0,     0,     0,     0,     0,   0,   -1./10,  1./10,  1./60, -1./60 },  
      { 0, 1./12, -1./12,     0,     0,     0,     0, -1./8,  1./8,   1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12,     0,     0, -1./6,  1./6, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12,  1./4, -1./4, 1./12,-1./12, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12, -1./4,  1./4, 1./12,-1./12, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60  }, 
      { 1, 0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  }, 
      { 0, 1./12,  1./12, -1./4, -1./4, 1./12, 1./12, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,  1./4,  1./4, 1./12, 1./12, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0, -1./6, -1./6, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0, -1./8,  -1./8,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0,     0,   0,   -1./10, -1./10,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0,     0,   0,  0,  0, -1./12, -1./12  }
   };

   std::vector<double> lha(13);
   for ( int i = 0 ; i<13; i++ ) {
      for ( int j = 0 ; j<13; j++ ) {
	 lha[i] += M[i][j]*LiCoApfl[j];
      }
      if ( fabs(lha[i])<1.e-13 ) lha[i]=0;
   }
   return lha;
}


//____________________________________________________________________________________ //
std::map<int,double> AlposTools::LicoLhaToApfelxxMap(const std::vector<double>& LiCoLha ) {
   //! calculate PDF linear combination as input to Apfel++
   //! Input vector denots PDFs like in LHAPDF: tb...g d u s c b t (0...6...12)
   
   // q^+ = q + qbar
   // q^- = q - qbar
   static const std::vector<std::vector<double> > M{
      {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}, // gluon
      {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}, // SIGMA = \sum_q q^+
      { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}, // VALENCE = \sum_q q^-
      {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}, // T3  = u^+ - d^+
      {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}, // V3  = u^- - d^-
      {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}, // T8  = u^+ + d^+ - 2 s^+
      {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}, // V8  = u^- + d^- - 2 s^-
      {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}, // T15 = u^+ + d^+ + s^+ - 3 c^+
      {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}, // V15 = u^- + d^- + s^- - 3 c^-
      {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}, // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+
      {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}, // V24 = u^- + d^- + s^- + c^- - 4 b^-
      { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}, // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
      {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,}  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-
   };
   std::map<int,double> apfl;
   apfl.insert({0,LiCoLha[6]});
   for ( int i = 1 ; i<13; i++ ) {
      double f = 0;
      for ( int j = 0 ; j<13; j++ ) {
	 f += M[i][j]*LiCoLha[j];
      }
      apfl.insert({i,f});
   }
   return apfl;
}

//____________________________________________________________________________________ //
std::map<int,double> AlposTools::LicoApfelxxToLhaMap(const std::vector<double>& LiCoApfl ) {
   //! calculate PDF linear combination as input to Apfel++
   //! Input vector denots PDFs like in LHAPDF: tb...g d u s c b t (0...6...12)
   
   // q^+ = q + qbar
   // q^- = q - qbar
   static const std::vector<std::vector<double> > M{
      { 0, 1./12, -1./12,     0,     0,     0,     0,     0,   0,  0,  0,   -1./12,  1./12 }, 
      { 0, 1./12, -1./12,     0,     0,     0,     0,     0,   0,   -1./10,  1./10,  1./60, -1./60 },  
      { 0, 1./12, -1./12,     0,     0,     0,     0, -1./8,  1./8,   1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12,     0,     0, -1./6,  1./6, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12,  1./4, -1./4, 1./12,-1./12, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60},   
      { 0, 1./12, -1./12, -1./4,  1./4, 1./12,-1./12, 1./24, -1./24,  1./40, -1./40,  1./60, -1./60  }, 
      { 1, 0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  }, 
      { 0, 1./12,  1./12, -1./4, -1./4, 1./12, 1./12, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,  1./4,  1./4, 1./12, 1./12, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0, -1./6, -1./6, 1./24,  1./24,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0, -1./8,  -1./8,  1./40,  1./40,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0,     0,   0,   -1./10, -1./10,  1./60,  1./60  }, 
      { 0, 1./12,  1./12,     0,     0,     0,     0,     0,   0,  0,  0, -1./12, -1./12  }
   };

   std::map<int,double> lha;
   for ( int i = 0 ; i<13; i++ ) {
      double f = 0;
      for ( int j = 0 ; j<13; j++ ) {
	 f += M[i][j]*LiCoApfl[j];
      }
      if ( fabs(f)<1.e-13 ) f=0;
      lha.insert({i-6,f});
   }
   return lha;
}

