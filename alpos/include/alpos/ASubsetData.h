// DB 15.01.2015
#ifndef Alpos_ASubsetData
#define Alpos_ASubsetData

/* 
 An implementation of the most trivial function:
 It consists only of one value 'theConst'.
 This function may be used for data, which directly
 measure a constant.

 */

#include "alpos/ATheory.h"
#include "alpos/AData.h"

class ASubsetData : public AData {

public:
   ASubsetData(const std::string& name);
   virtual ~ASubsetData();

   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const {return ASubsetData::fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const; //!< Return the function name of the underlying function (type)
   std::string GetRealFunctionName() const { return fFunctionName;}; //!< Return 'SubsetData'
   static const std::string fFunctionName; //!< The function's name
   void InitErrors(); //!< Initialize fAllErrors
   void SetRequirementValidPoints(const std::string& req, const std::vector<bool>& valid);
   // overload getter to return reduced data table
   const std::map<std::string, std::vector<double> >& GetDataTable() const { return fReducedDataTable; };
   
private:
   std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   std::vector<bool> fPointValid;

   std::map<std::string, std::vector<double>> fReducedDataTable; //!< A subset of the original data table

protected:
   

};


#endif
