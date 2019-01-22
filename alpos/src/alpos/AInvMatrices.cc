#include "alpos/AInvMatrices.h"

#include "TMatrixTSym.h"

/** 
 AInvMatrices
 
 calculate sum of matrices. Invert the sum and stores it for further access.

 */

using namespace std;

//____________________________________________________________________________________ //
AInvMatrices* AInvMatrices::instance = NULL;


//____________________________________________________________________________________ //
AInvMatrices::AInvMatrices() : AlposObject("AInvMatrices") {
}


//____________________________________________________________________________________ //
AInvMatrices::~AInvMatrices(){

}


//____________________________________________________________________________________ //
AInvMatrices* AInvMatrices::Instance(){
   if (instance==NULL) instance = new AInvMatrices();
   return instance;
}


//____________________________________________________________________________________ //
const TMatrixDSym& AInvMatrices::GetInvMatrix(std::set<const TMatrixDSym*> mats, bool recalculate){
   //! calculate sum of all given matrices and invert it and cache it.
   //! Return inverse matrix.
   //! Keep inverse for next access.
   //! if recalculate is true, force recalculation of the cached matrix (default: false)

   if ( (this->count(mats) != 1) || (recalculate) ) {
      if ( mats.size()==1 ) {
         const TMatrixDSym& m = *(*mats.begin());
         (*this)[mats].ResizeTo(m);
         (*this)[mats] = AlposTools::InvertChol(m);
      }
      else {
         TMatrixDSym m(*(*mats.begin()));
         int i = 0;
         for ( auto mm : mats ) {
            if ( i++ == 0 ) continue;
            m+= *(mm);
         }
         (*this)[mats].ResizeTo(m);
         (*this)[mats] = AlposTools::InvertChol(m);
      }
   }
   return (*this)[mats]; // return cached inverse
}


//____________________________________________________________________________________ //
const TMatrixDSym& AInvMatrices::GetInvMatrix(const TMatrixDSym* mat, bool recalculate){
   return GetInvMatrix(std::set<const TMatrixDSym*>{mat}, recalculate);
}


//____________________________________________________________________________________ //
const TMatrixDSym& AInvMatrices::GetInvMatrix(const TMatrixDSym& mat, bool recalculate){
   return GetInvMatrix(&mat, recalculate);
}


//____________________________________________________________________________________ //
