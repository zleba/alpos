#include "Constraint.h" 

/* Example implementation for function to perform constrained least square fit.

   D. Reichelt 09.2015
 */


// destructor
Constraint::~Constraint() 
{
  if(apcFitter) delete apcFitter;
}

// initialise the fitter (and delete if one already exists)
void Constraint::InitFit()
{
  if(apcFitter) delete apcFitter;
  apcFitter = new Apccpp(evInitialX());
}

// perform the fit (and intialise first, if not done yet)
void Constraint::Fit()
{
  if(!apcFitter) InitFit();
  Eigen::VectorXd evF = evConstrVec(apcFitter->evGetX());
  Eigen::MatrixXd emCov = emCovariance(apcFitter->evGetX());
  do
    {
      do
        {
          evF = evConstrVec(apcFitter->evGetX());
        } while(apcFitter->bApcpp(evF, emCov));
      emCov = emCovariance(apcFitter->evGetX());
    } while(!apcFitter->bIsFinished());
}
