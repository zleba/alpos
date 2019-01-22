// DB 15.01.2015
#ifndef Alpos_AQcdnumInit
#define Alpos_AQcdnumInit

/* 
 * A function to initialize QCDNUM
 * Only one instance of this function is allowd
 */

#include "alpos/ATheory.h"

class AQcdnumInit : public AParmFuncBase<double> {

public:
   AQcdnumInit(const std::string& name);
   virtual ~AQcdnumInit();

   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name

   double pdfinput(int ipdf, double xp); //!< PDFs at starting scale and pass to QCDNUM

private:
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   static int fNinstances;

   void Setcbt(); //! Set cbt in QCDNUM
   void SetAsMz(); //! Set As(Mz) and Mz in QCDNUM
   int SetXGrid(std::vector<double> xgrid, std::vector<int> wgt, int nNodes);
   int SetQGrid(std::vector<double> qgrid, std::vector<double> wgt, int nNodes);

protected:
   bool fInitEvol; //!< Init full evolution
   double fQ0;
};


#endif
