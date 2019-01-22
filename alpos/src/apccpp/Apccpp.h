#ifndef APCCPP_H_
#define APCCPP_H_



#include <iostream>
#include <vector>
#include <iomanip>

//#define EIGEN_NO_DEBUG
#include "Eigen/Dense"
#include <ctime>


/*
  class Apccpp
  Implementing fit by the constrained least square method.

  D. Reichelt 09.2015
*/

class Apccpp {

 public:
  // constructors
  Apccpp(const Eigen::VectorXd& evXInit);
  Apccpp(const std::vector<double>& evXInit);
  
  // destructor
  ~Apccpp() { };
  
  // setter
  void SetEpsilon(double dEps) {dEpsilonChi2=dEps; dEpsilonFav=dEps;}
  void SetEpsChi2(double dEps) {dEpsilonChi2=dEps;}
  void SetEpsFav(double dEps) {dEpsilonFav=dEps;}
  void SetMaxIter(int iMaxIt) {iMaxIter=iMaxIt;}
  void SetPrintFlag(int iIpSet) {iIp=iIpSet;}

  // getter
  double dGetX(int iInd) const {return evX(iInd);}
  const Eigen::VectorXd& evGetX() const {return evX;}
  double dGetDx(int iInd) const {return evDx(iInd);}
  const Eigen::VectorXd& evGetDx() const {return evDx;}
  double dGetStdv(int iInd) const {return sqrt(emInv(iInd,iInd));}
  Eigen::VectorXd evGetStdv() const {return emGetCovar().diagonal().cwiseSqrt();}
  Eigen::MatrixXd emGetCovar() const {return emInv.topLeftCorner(evX.size(),
                                                                     evX.size());}
  double dGetChi2() const {return dChi2Act;}

  // status ratings
  bool bIsChi2ChangeSmall() const {return (fabs(dChi2Act-dChi2Prev) < dEpsilonChi2);}
  bool bIsFavSmall() const {return (fabs(dFav()) < dEpsilonFav);}
  bool bIsConverged() const {return bIsChi2ChangeSmall() && bIsFavSmall();}
  bool bReachedMaxIt() const {return (iNiter >= iMaxIter);}
  bool bIsFinished() const {return (bIsConverged() || bReachedMaxIt());}

  // main functions
  bool bApcpp(const Eigen::VectorXd& evF, const Eigen::MatrixXd& emVx);                  
  Eigen::MatrixXd emInvertMatrix(Eigen::MatrixXd emFullMat, int iUppSq);
  void Reset();
  
  void PrintResult() const;       // print out results (fitted x and init. x w/ stdv.)
  
  // functions to calculate parameters for current fit status
  double dFav() const;
  double dPull(int iInd) const;
  Eigen::VectorXd evPull() const;

  

 private:

  // internal functions (made privat because using them is not usefull/possible in 
  // all states)
  bool bSolveEq(const Eigen::MatrixXd& emVx);   // solve the equation for new corrections
  bool bDerivative(const Eigen::VectorXd& evF); // calc. one step in derivative calc.
   
                                         
  // private print functions respecting the value of print flag
  void PrintStep(double dChi2, double dChi2Change,  // print parameters of last iteration
                 double dFav, double dFavChange) const;
  void PrintDebug(const std::string& sMessage) const;     // show a debug message
  void PrintError(const std::string& sMessage) const;
  
  // user provided variables
  Eigen::VectorXd evX;        // vector of variables
  Eigen::VectorXd evXorig;    // initially provided vector X
  Eigen::VectorXd evXunch;    // vectorX without any changes due to derivative calc.

  Eigen::VectorXd evFcurr;    // current constraint values
  Eigen::VectorXd evVar;      // variances according to last provided matrix  

  // limiting parameters for fit
  double dEpsilonChi2=pow(10, -8);  // > change of chi2 for rating "converged"
  double dEpsilonFav=pow(10, -8);   // > av. constr. value for rating "converged"
  int iMaxIter=40;           // maximum number of iteration, if iNiter exceeds rating is
                             // "reached max iter" and "finished"


  //variables obtained during fit
  Eigen::VectorXd evDx;       // vector of corrections
  Eigen::MatrixXd emInv;      // covariance matrix after fit
  Eigen::MatrixXd emA;        // matrix of first derivatives


  // savings of parameters for last iteration
  double dChi2Act=0;
  double dChi2Prev=0;         // chi2 of previous iteration
  double dFavAct=0;
  double dFavPrev=0;          // average constraint value of previous iteration
  int iNiter=0;               // iteration counter

  // steering flags
  int iIp=2;                  // print flag 0=no output/
                              //            1=standard output and terminating errors
                              //            2=warnings on fit behaviour /3=debug

  
  // stuff needed in derivative calculation
  int iDir=0;                 // direction in which derivative is calculated 
  int iSign=0;                // sign for next value of F
  Eigen::VectorXd evStep;     // stepsizes for derivatives
  
};

#endif
