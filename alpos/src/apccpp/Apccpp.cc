#include "Apccpp.h"

/* Implementing the Constraint least square method

   D. Reichelt 09.2015
*/

using std::cout;
using std::endl;


/////////////////////////////////////////////////////////////////////////////////////////
// constructors /////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

/* constructor
   Arguments: evXInit ... vector containing all measured data and parameters
 */
Apccpp::Apccpp(const Eigen::VectorXd& evXInit) 
  : evX(evXInit), evXorig(evXInit), evXunch(evXInit),
    evVar(Eigen::VectorXd::Ones(evXInit.size())),
    evDx(Eigen::VectorXd::Zero(evXInit.size())),
    emInv(Eigen::MatrixXd::Identity(evXInit.size(),evXInit.size()))
{ 
    PrintDebug("Init from Eigen::VectorXd");
}

/////////////////////////////////////////////////////////////////////////////////////////
// public functions /////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


// main functions ///////////////////////////////////////////////////////////////////////

/* steering loop for derivatives
   Arguments: evF ... vector with constraints at current point evX
   Return: true if another evaluation of constraints with current evX is needed
           false if iteration was done successfull
   use e.g. as
     do
       { 
         evF = // calculate values of constraints as function ov this.evX
       } while(this.bApccpp(evF))
 */
bool Apccpp::bApcpp(const Eigen::VectorXd& evF, const Eigen::MatrixXd& emVx) {
  bool bRet = bDerivative(evF);      // calc derivative, bRet determines if 
                                     // another eval is needed
  if(!bRet)                          // if deriv. calc done, solve eqation
    {
      bool sucs = bSolveEq(emVx); 
      if(!sucs)
        {
          PrintError("Failed to solve equation!");
          exit(1);
        }
      ++iNiter;
    }
  return bRet;
}

/* Invert symmetric matrix emRet, where steps are done for the upper left iUppSq*iUppSq
   square matrix
   Arguments: emRet ... some symmetric matrix (only lower triangle is considered) to be 
                        inverted
              iUppSq ... size of the upper left square matrix to which this algorithm 
                         already was applied
              dEps ... precision to which numbers have to be distinct from 0 to be
                       be counted as 0
    Return: Inverse of emRet (as symmetric matrix)
 */
Eigen::MatrixXd Apccpp::emInvertMatrix(Eigen::MatrixXd emFullMat, int iUppSq) 
{
  PrintDebug("Inverting Matrix!");
  const int iDim = emFullMat.rows();                     // dimension of matrix
  int iPivInd;                                           // index of pivot element
  double dInvPivEl;                                      // value of pivot element
  int iNstep = 0;                                        // number of steps done
  Eigen::VectorXd evPivRow;                              // save pivot row 
  Eigen::ArrayXd eaWasPiv = Eigen::ArrayXd::Ones(iDim);  // array to keep track which 
                                                         // elements have been pivot
                                                         // (0 -> already has been, 
                                                         //  1-> still available)
  eaWasPiv.head(iUppSq).setZero();                       // first iUppSq values already 
                                                         // have been pivot

  Eigen::VectorXd evTol = pow(10, -6)*emFullMat.diagonal();
  //  time_t begin = time(0);
  do
    {
      if((emFullMat.diagonal().array()*eaWasPiv).abs().maxCoeff(&iPivInd) 
         > evTol(iPivInd))                               // find next best pivot element
        {        
          dInvPivEl = 1/emFullMat(iPivInd, iPivInd);
          eaWasPiv(iPivInd) = 0;                                           
  
          // do exchange step with pivot elements
          evPivRow = emFullMat.row(iPivInd);
          evPivRow.tail(iDim-iPivInd) = emFullMat.col(iPivInd).tail(iDim-iPivInd);
          emFullMat.selfadjointView<Eigen::Lower>().rankUpdate(evPivRow, -dInvPivEl); 
          emFullMat.row(iPivInd) = evPivRow.transpose()*dInvPivEl;
          emFullMat.col(iPivInd) = evPivRow*dInvPivEl;          
          emFullMat(iPivInd, iPivInd) = -dInvPivEl;
          ++iNstep;
        }
      else 
        {
          
          for(int i=0; i<iDim; ++i)
            {
              if(eaWasPiv(i) > 0)
                {
                  PrintDebug("Matrix (almost) singular!");
                  emFullMat.row(i).setZero();
                  emFullMat.col(i).setZero();
                  eaWasPiv(i) = 0;
                }
            }
        }
    } 
  while((eaWasPiv!=0).any());
  //time_t end = time(0);
  //cout << difftime(end,begin) << ", " << iNstep << endl;
  emFullMat *= -1;
  return emFullMat.selfadjointView<Eigen::Lower>();
}

/* Set fit back to "initial state", e.g. as no iterations had been done 
   (changes in covariance and steering parameters are kept)
 */
void Apccpp::Reset() 
{
  evX = evXorig;
  evXunch = evXorig;
  evDx = evXorig*0;
  evVar = Eigen::VectorXd::Ones(evX.size());
  dChi2Act = 0;
  dChi2Prev = 0;
  dFavAct = 0;
  dFavPrev = 0;
  iNiter = 0;
  iDir = 0;
  iSign = 0;
  PrintDebug("Reset Fit!");
}


// calculating parameters of last iteration /////////////////////////////////////////////

/* Calculate average of constraints
   Return: average of absolute value of all constraints
 */
double Apccpp::dFav() const {
  double dFsum = evFcurr.cwiseAbs().sum();
  return dFsum/evFcurr.size();
}


/* Calculate pull for one fit parameter
   Arguments: iInd ... index of x for which pull is wanted (< size of Vx)
   Return: pull, defined by p_i = \delta x_i / \sqrt(V_ii_unfit-V_ii_fitted)
 */
double Apccpp::dPull(int iInd) const {
  if(iInd < evVar.size()) {
    return evDx(iInd)/sqrt(fabs(evVar(iInd)-emInv(iInd,iInd)));
  }
  else {
    PrintError("Tried to calculate Pull for parameter out of range!");
    exit(1);
  }
}

/* Calculate vector of pulls
   Return: vector containing pulls for all x_i which are covered by the covariance
           matrix
 */
Eigen::VectorXd Apccpp::evPull() const {
  const int iNvx = evVar.size();
  Eigen::VectorXd evPulls(iNvx);
  for(int i=0; i<iNvx; ++i) 
    {
    evPulls(i) = dPull(i);
    }
  return evPulls;
}

/////////////////////////////////////////////////////////////////////////////////////////
// private functions ////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// internal functions //////////////////////////////////////////////////////////////////

/* Get together full linear system to solve after derivatives are calculated and sole it
   Arguments: emVx ... covariance matrix to be used. Will be interpreted as the 
              covariance matrix for the first n = emVx.size() elements in X. Remaining 
              elements are considered as unmeasured -> parameters. If n is larger than
              the size m of X, only the top left m x m  matrix will be used.
   Return: true if job is done
   At the end, the new corrections evDx, the corrected parameters evX and the new 
   covariance matrix emInv are deterined and set accordingly.
   New values for Chi2 and average constraints are calculated and old values are moved
   to the respective -Prev variables.
 */
bool Apccpp::bSolveEq(const Eigen::MatrixXd& emVx) {
  
  PrintDebug("Start solving equation!");
  // get the dimensions of involved objects 
  const int iNf = evFcurr.size();                // No. constraints
  const int iNx = evX.size();                    // No. fit parameters (measured+unmeas.)
  const int iNvx = emVx.cols();                  // size of covar matrix (no. of meas.)
  const int iDim = iNf + iNx;                    // total size of matrix in lhs of eq.

  
  //emVx.fullPivLu().isInvertible() ? 
  //cout << "Matrix invertible!" : cout << "Matrix NOT invertible!";
  //cout << endl;


  // RHS of equation, Nx zeros and linear approximation of constraints
  Eigen::VectorXd evRhs(iDim);
  evRhs.head(iNx).setZero();  
  evRhs.tail(iNf) = -evFcurr + emA*evDx;

  // full matrix (lower triangular) with weights and derivatives
  Eigen::MatrixXd emFullMat(iDim,iDim);
  emFullMat.triangularView<Eigen::Lower>() = Eigen::MatrixXd::Zero(iDim,iDim);
  if(iNvx < iNx) 
    {
      emFullMat.topLeftCorner(iNvx,iNvx)
        .triangularView<Eigen::Lower>() = -emVx;
      // save initial variance
      evVar = emVx.diagonal();
    }
  else 
    {
      emFullMat.topLeftCorner(iNx,iNx)
        .triangularView<Eigen::Lower>() = -emVx.topLeftCorner(iNx,iNx);
      // save initial variance
      evVar = emVx.diagonal().head(iNx);
    }
  emFullMat.bottomLeftCorner(iNf,iNx) = emA;
  emFullMat.bottomLeftCorner(iNf,iNvx) *= emVx.selfadjointView<Eigen::Lower>();
  emFullMat.bottomRightCorner(iNf,iNf).triangularView<Eigen::Lower>()=-emA.leftCols(iNvx)
    *emVx.selfadjointView<Eigen::Lower>()*emA.leftCols(iNvx).transpose();


  /*//define actual matrix to be converted for consistency checks
  Eigen::MatrixXd emTestInv = Eigen::MatrixXd::Zero(iDim,iDim);;
  Eigen::MatrixXd emW = emVx.inverse();
  if(iNvx < iNx) 
    {
      emTestInv.topLeftCorner(iNvx,iNvx)
        .triangularView<Eigen::Lower>() = emW;
    }
  else 
    {
      emTestInv.topLeftCorner(iNx,iNx)
        .triangularView<Eigen::Lower>() = emW.topLeftCorner(iNx,iNx);
    }
  emTestInv.bottomLeftCorner(iNf,iNx) = emA;
  emTestInv = emTestInv.selfadjointView<Eigen::Lower>();
  */

  // solve equation (i.e. invert matrix)
  emInv = emInvertMatrix(emFullMat, iNvx);
  Eigen::VectorXd evSol = emInv*evRhs;
  evDx = evSol.head(iNx);
  //cout << (emTestInv*emInv - Eigen::MatrixXd::Identity(iDim,iDim)).cwiseAbs().maxCoeff() << "  "   
  //     << (emTestInv*emInv - Eigen::MatrixXd::Identity(iDim,iDim)).cwiseAbs().sum() <<  endl;


  // set X to corrected value
   
  evX  = evXorig + evDx;
  evXunch = evX;

  // get chi2 and average constraint
  dChi2Prev = dChi2Act;
  dChi2Act = -evDx.head(iNx).transpose()*emA.leftCols(iNx).transpose()*evSol.tail(iNf);
  dFavPrev = dFavAct;
  dFavAct = dFav();
  const double dChi2Change =  dChi2Act - dChi2Prev;
  const double dFavChange =  dFavAct - dFavPrev;
  PrintStep(dChi2Act, dChi2Change, dFavAct, dFavChange);

  // done succesfull
  return true;
}


/* Do one step in the derivative calc. and return true if another step is needed.
   Arguments: evF ... values of constraints at this.evX as set after the last call
   Return: true if another evaluation is needed, in this case calculate recalculate evF 
           at evX and call again.
           false if no more calculation is needed. Start solving afterwards
 */
bool Apccpp::bDerivative(const Eigen::VectorXd& evF) 
{
  // length of involved vectors
  const int iNx = evX.size();           // number of fit parameters (meas. and unmeas.)
  const int iNf = evF.size();           // number of constraints

  // can be in three states: 0 - first step in iteration, initialize everything 
  //                             and prepare first step
  //                        -1 - prev. calculated F at x+step in direction iDir, 
  //                             add this to deriv. calc and set X to
  //                             x-step in direction iDir
  //                         1 - prev. calculated F ar x-step in direction iDir, 
  //                             subtract this from deriv calc and 
  //                             divide by step size. Then set iDir to next direction,
  //                             check if end of iteration is 
  //                             reched, otherwise set x to x+step in new iDir
  
  switch (iSign) {
  
  case -1:                                      // after F calc in pos. dir. (1. step)
    emA.col(iDir) += evF;
    evX(iDir) = evXunch(iDir)-evStep(iDir);
    iSign = 1;                                  // iSig will be 1 in next call
    return true;                                // signal that another F eval is needed
    break;

  case 1:                                       // after F calc in neg. dir (2. step)
    emA.col(iDir) -= evF;
    emA.col(iDir) /=
      (evXunch(iDir)+evStep(iDir))-(evXunch(iDir)-evStep(iDir));
    evX(iDir) = evXunch(iDir);
    ++iDir;
    if(iDir<iNx)                                 // if not all der. are calculated
      {
        evX(iDir) += evStep(iDir);
        iSign = -1;                              // iSig will be -1 in next call
        return true;                             // signal that another F eval is needed
      }
    else                                         // if end of iteration is reached, 
      {                                          // set everything back to initial values
        iDir = 0;      
        iSign = 0;      
        return false;                            // signal end of iteration
      }
    break;

  case 0:                                       // if this is the first step
    PrintDebug("Start deriavtive calculation now!");
    emA.setZero(iNf,iNx);                       // clear matrix A before filling it
    evStep =  evXunch.cwiseAbs()*pow(10,-1);     // calculate stepsizes for current 
                                                // values of X
    int iMinCoeff;
    while(evStep.array().abs().minCoeff(&iMinCoeff)==0)
      {
        evStep(iMinCoeff) = pow(10,-1);
      }
    evFcurr = evF;                              // save the constraints at the current X

    evX = evXunch;                              // prepare first step
    evX(iDir) += evStep(iDir);                                      
    iSign = -1;                                 // iSign will be -1 at next call
    return true;                                // signal that another F eval is needed
    break;
  }
  
  PrintError("Fatal error in derivative calculation!"); // not possible to continue if 
  exit(1);                                              // none of these cases

} 


// print functions //////////////////////////////////////////////////////////////////////

void Apccpp::PrintStep(double dChi2, double dChi2Change,
                       double dFav, double dFavChange) const
{
  if(iIp==0) {return;}
  const int iSmallWidth = 10;
  const int iLargeWidth = 15;
  if(iNiter==0) 
    {
      cout << std::left << std::setw(iSmallWidth) << std::setfill(' ') << "iter";
      cout << std::left << std::setw(iSmallWidth) << std::setfill(' ') << "ndf"; 
      cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << "chi2";
      cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << "delta";
      cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << "|F|/NF";
      cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << "delta";
      cout << endl;
    }
  cout << std::left << std::setw(iSmallWidth) << std::setfill(' ') << iNiter+1;
  cout << std::left << std::setw(iSmallWidth) << std::setfill(' ') 
       << evFcurr.size()-(evX.size()-evVar.size()); 
  cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << dChi2;
  cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << dChi2Change;
  cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << dFav;
  cout << std::left << std::setw(iLargeWidth) << std::setfill(' ') << dFavChange;
  cout << endl;
  if(iIp==2)
    {
      if(bReachedMaxIt())
        {
          cout << "Maximum number of Iterations reached. Probably not converging." 
               << endl;
        }
    }
}


void Apccpp::PrintResult() const 
{
  const int iWidth = 15;
  const int iNx = evX.size();
  const int iNvx = evVar.size();
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Par name";
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Fitted";
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Std.dev.";
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Initial";
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Std.dev.";
  cout << std::left << std::setw(iWidth) << std::setfill(' ') << "Pull";
  cout << endl;
  for(int i=0; i<iNx; ++i)
    {
      cout << std::left << std::setw(iWidth) << std::setfill(' ') << i;
      cout << std::left << std::setw(iWidth) << std::setfill(' ') << evX(i);
      cout << std::left << std::setw(iWidth) << std::setfill(' ') << sqrt(emInv(i,i));
      cout << std::left << std::setw(iWidth) << std::setfill(' ') << evXorig(i);
      if(i<iNvx) cout << std::left << std::setw(iWidth) 
                      << std::setfill(' ') << sqrt(evVar(i));
      else cout << std::left << std::setw(iWidth) << std::setfill(' ') << ' ';
      if(i<iNvx) cout << std::left << std::setw(iWidth) << std::setfill(' ') << dPull(i);
      else cout << std::left << std::setw(iWidth) << std::setfill(' ') << ' ';
      cout << endl;
    }
}

void Apccpp::PrintDebug(const std::string& sMessage) const
{
  if(iIp<3) {return;}
  cout << "DEBUG: " << sMessage << endl;
}


void Apccpp::PrintError(const std::string& sMessage) const
{
  if(iIp<1) {return;}
  cout << "ERROR: " << sMessage << endl;
}
