// DB 22.02.2015
#ifndef Alpos_AEprc
#define Alpos_AEprc

/* 
 * A function to calculate EW parameters with EPRC
 * Only one instance of this function is allowed
 */

#include "alpos/ATheory.h"

class AEprc : public AParmFuncBase<double> {

public:
   AEprc(const std::string& name);
   virtual ~AEprc();

   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name
   virtual std::vector<double> GetQuick(int n,...); //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&); //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values


private:
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   static int fNinstances;
protected:
   

};


#endif
