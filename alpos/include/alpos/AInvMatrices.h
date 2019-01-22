// DB. 2/2016
#ifndef Alpos_AInvMatrices
#define Alpos_AInvMatrices

/** 
 * 
 * 
 * 
 */

#include <map>
#include <set>

#include "alpos/AlposObject.h"
#include "alpos/AlposTools.h"
#include "TMatrixTSym.h"

/**
 AInvMatrices

 on-demand caching for inverse matrices and inverse of sum-of-matrices

*/

class AInvMatrices : public std::map<std::set<const TMatrixDSym* >,TMatrixDSym> , public AlposObject {

private:
   AInvMatrices();
   ~AInvMatrices();

public:
   static AInvMatrices* Instance();
   // const TMatrixDSym& SetInvMatrix(std::set<const TMatrixDSym&> mats);
   // const TMatrixDSym& SetInvMatrix(std::set<const TMatrixDSym* const> mats);
   // const TMatrixDSym& GetInvMatrix(std::set<const TMatrixDSym&> mats);
   const TMatrixDSym& GetInvMatrix(std::set<const TMatrixDSym*> mats, bool recalculate=false);
   const TMatrixDSym& GetInvMatrix(const TMatrixDSym* mat, bool recalculate=false);
   const TMatrixDSym& GetInvMatrix(const TMatrixDSym& mat, bool recalculate=false);

private:
   static AInvMatrices* instance;
   //std::map<std::set<const TMatrixDSym* >,TMatrixDSym> fInvMat;
   
};


#endif
