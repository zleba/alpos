#ifndef CONSTRAINT_H
#define CONSTRAINT_H
#include "../Apccpp.h"

/*
  class Constraint
  Abstract class for setup of constraints, calculation of constaints and fitting
  The functions evInitialX, evConstrVec and emCovariance are meant to return the 
  intital vector X from which the Apccpp instance should be constructed, the 
  vector of constraints and the covariance matrix as used in the fit, 
  respectively.
  The functions InitFit and Fit will initialise the fitter and perform the fit.

  D. Reichelt 09.2015
 */

class Constraint {
 public: 
  Constraint() {};
  virtual ~Constraint();
  virtual Eigen::VectorXd evInitialX() const = 0;
  virtual Eigen::VectorXd evConstrVec(const Eigen::VectorXd& X) const = 0;
  virtual Eigen::MatrixXd emCovariance(const Eigen::VectorXd& X) const = 0;

  void InitFit();
  void Fit();
  
  Apccpp* apcFitter=nullptr;

};
#endif
