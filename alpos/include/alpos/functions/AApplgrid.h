// DB 07.12.2015
#ifndef Alpos_AApplgrid
#define Alpos_AApplgrid

/* 
 An implementation of APPLgrid as an Alpos function.
 */

#include "alpos/ATheory.h"
#include "appl_grid/appl_grid.h"

class AApplgrid : public AParmFuncBase<double> {

public:
   AApplgrid(const std::string& name);
   virtual ~AApplgrid();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   const std::vector<bool>& GetBinmap() const { return fBinmap;}; //!< get bin map of valid points obtained from datacard

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   // members
   //fastNLOAlpos* fnlo = NULL;
   std::vector<appl::grid*> fgrids;
   
   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

protected:
   std::vector<bool> fBinmap;
};


#endif
