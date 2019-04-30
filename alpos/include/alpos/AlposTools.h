// DB 15.01.2015
#ifndef Alpos_AlposTools
#define Alpos_AlposTools

#include "TMatrixDSymfwd.h"
#include "TMatrixTSym.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include <vector>
#include <set>
#include <iostream>
#include <map>

namespace AlposTools {

   //! --- useful matrix operations/functions
   template <typename T> T sq(const T& v) { return v*v;} // square
   std::vector<double > CovToCorr(TMatrixDSym& mat);
   void CovToCovRel(TMatrixDSym& mat, const std::vector<double>& sigma);
   void CovRelToCov(TMatrixDSym& mat, const std::vector<double>& sigma);
   void CovRescale(TMatrixDSym& mat, const std::vector<double>& sigma_current, const std::vector<double>& sigma_rescale);
   void CorrToCov(TMatrixDSym& mat, const std::vector<double>& error);
   TMatrixD InvertLU(const TMatrixD& mat);
   TMatrixDSym InvertChol(const TMatrixDSym& mat);
   TMatrixDSym InvertChol(std::set<const TMatrixDSym*> mat_collection);

   bool IdenticalToNSignificantDigits(double a, double b, unsigned int n=5);

   //! --- other stuff
   void CheckFileExit(std::string& file);
   bool CheckFile(std::string& file);

   //! --- calculations for PDF LiCo's
   void Calc13partonsFromLico(std::vector<double>& xfx13, const TMatrixD& Inv, const std::vector<double>& xf6, const std::vector<bool>& found );
   void CalcLicoFrom13partons(std::vector<double>& LiCo, const TMatrixD& def, const std::vector<double>& xfx13 );
   void CalcLicoFrom13partons(std::vector<double>& LiCo, const std::vector<double>& def, const std::vector<double>& xfx13 );
   std::vector<double> LicoApfelxxToLha(const std::vector<double>& LiCoApfl );
   std::vector<double> LicoLhaToApfelxx(const std::vector<double>& LiCoLha  );
   std::map<int,double> LicoApfelxxToLhaMap(const std::vector<double>& LiCoApfl );
   std::map<int,double> LicoLhaToApfelxxMap(const std::vector<double>& LiCoLha  );

   double rfluxRawInt(double a0, double ap, double b0,  double x_pom, double tAbsMin, double tAbsMax);
   double rfluxRaw(double a0, double ap, double b0, double x_pom, double tAbs);
   double rfluxInt(double a0, double ap, double b0, double x_pom, double tAbsMin, double tAbsMax);
   double rflux(double a0, double ap, double b0, double x_pom, double tAbs);

   //! --- useful (but ambivalent) vector operations 
   //! Multiply element by element
   template <typename T>
   void operator*=(std::vector<T> &v1, const std::vector<std::string> &v2) {
      std::cout<<"Operator *= is not defined for 'strings'."<<std::endl;
      //for ( unsigned int i = 0 ; i<v1.size() ; i++ ) 
      //v1[i] *= v2[i];
   }

   template <typename T,typename A>
   void operator*=(std::vector<T> &v1, const std::vector<A> &v2) {
      for ( unsigned int i = 0 ; i<v1.size() ; i++ ) 
	 v1[i] *= v2[i];
   }

   //! Divide element by element
   template <typename T, typename A>
   void operator/=(std::vector<T> &v1, const std::vector<A> &v2) {
      for ( unsigned int i = 0 ; i<v1.size() ; i++ ) 
	 v1[i] /= v2[i];
   }

   //! Multiply each element with constant
   template <typename T>
   void operator*=(std::vector<T> &v1, double v2) {
      for ( auto& i : v1 )
	 i*=v2;
   }

   //! Divide each element with constant
   template <typename T>
   void operator/=(std::vector<T> &v1, double v2) {
      for ( auto& i : v1 ) i/=v2;
   }

   //! append one vector to the other
   template <typename T>
   void operator+=(std::vector<T>& v1, const std::vector<T>& v2) {
      size_t sz = v1.size();
      v1.resize(v1.size()+v2.size());
      for ( const auto& i2 : v2 )
	 v1[sz++] = i2;
   }

   //! --- more vector operations
   template <typename T>
   std::vector<T> VectorSubset(const std::vector<T>& vals,const std::vector<bool>& valid) {
      std::vector<T> ret(vals.size());
      unsigned int n=0;
      for ( unsigned int i =0 ; i<vals.size() ; i++ ){
	 if ( valid[i] ) {
	    ret[n++] = vals[i];
	 }
      }
      ret.resize(n);
      return ret;
   }

   
}

#endif
