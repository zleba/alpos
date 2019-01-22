// DB 15.01.2015
#ifndef Alpos_AData
#define Alpos_AData

/* 
 An implementation of AData

 */

#include "alpos/ATheory.h"
#include "alpos/AError.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixTSym.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"

class AData : public AParmFuncBase<double> {

public:
   AData(const std::string& name);
   AData(const std::string& name, bool klar); //!< klar without meaning, but to specify that "_Data" is not appended
   virtual ~AData();

   // --- AFuncD defaults
   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name
   // ---------------------

   // --- access to errors and data
   virtual const std::map<std::string, std::vector<double> >& GetDataTable() const {return fOrigData;}

   bool AddColumnToDataTable(const std::string& colname, const std::vector<double>&);

   void ChangeValues(const vecT& val, const vecT& err) { SetValues(val,err);}
   
   // --- get subsets
   const std::map<std::string,std::vector<bool> >& GetSubsets() const { return fSubsets;}; //!< Get subsets <Subsetname,ValidPoints>

   // --- statics
   static bool CheckType(std::string& type); //!< check if 'type' is reasonable

private:
   // --- AFuncD defaults
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values
   // ---------------------
   TMatrixDSym ReadMatrix(const std::vector<std::vector<double> >& values, const std::string& format="Matrix", 
			  const std::vector<std::string>& matHeader=std::vector<std::string>(), const std::vector<std::string>& DataHeader=std::vector<std::string>(), 
			  const std::vector<std::vector<double> >& DataVals=std::vector<std::vector<double> >()); //!< Read a (symmetric) matrix in the given 'format' ('Matrix' or 'SingleValues')
   
private:
   std::string fDataName; //!< The 'unique' dataset name
   std::map<std::string, std::vector<double > > fOrigData; //!< the 'Data'-table from the input file

protected:
   const TMatrixD& CheckReturnInverse(TMatrixD& Inv, const TMatrixDSym& Mat) const;
   void InitSubsets(); //!< Init the subsets for different phase space regions
   void InitSubsetWithCuts(); //!< Prepare a subset, where all cuts are applied (this serves as input to SuperData)
   
   void ReadErrors(const std::string& DataName);
   std::map<std::string,std::vector<bool> > fSubsets;

   // std::map<std::vector<double> > fCorrErrors;
   // std::vector<double > fUncorrErrors;
   // std::vector<double > fStatErrors;

};


#endif
