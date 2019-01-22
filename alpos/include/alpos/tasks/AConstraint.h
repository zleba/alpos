#ifndef Alpos_AConstraint
#define Alpos_AConstraint

/**
 * classes for constraints
 *
 * provide inital data vector, covariance matrix and constraint
 * calculation for use in AApcFitter and AConstLQFitter
 *
 */

#include <Apccpp.h>
#include "alpos/AData.h"

using namespace std;

class AConstraint {
 public:
 AConstraint(const vector<string>& FitPar, AData* Data, AFuncD* Theo) 
   : fFitPar(FitPar), fData(Data), fTheo(Theo) {};
  virtual ~AConstraint() {};

  virtual Eigen::VectorXd InitialX() const = 0;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const = 0;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd& X) const = 0;

 protected:
  std::vector<string> fFitPar;
  AData* fData=nullptr;
  AFuncD* fTheo=nullptr;
};


//____________________________________________________________________________//

class ANormal : public AConstraint {
 public:
  ANormal(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~ANormal() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd& X) const;

 protected:
  Eigen::MatrixXd fCovar;
};

//____________________________________________________________________________//

class ALogNormal : public AConstraint {
 public:
  ALogNormal(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~ALogNormal() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd& X) const;

 protected:
  Eigen::MatrixXd fCovar;
};

//____________________________________________________________________________//

class ALogNormalNuisance : public AConstraint {
 public:
  ALogNormalNuisance(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~ALogNormalNuisance() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd& X) const;

 protected:
  Eigen::MatrixXd fCovar;
  Eigen::MatrixXd fAddErr;
  Eigen::MatrixXd fMultErr;
  Eigen::MatrixXd fMultAddErr;
};


//____________________________________________________________________________//


class ANormalNuisance : public AConstraint {
 public:
  ANormalNuisance(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~ANormalNuisance() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd& X) const;

 protected:
  Eigen::MatrixXd fCovar;
  Eigen::MatrixXd fAddErr;
  Eigen::MatrixXd fMultErr;
  Eigen::MatrixXd fMultAddErr;
};

//____________________________________________________________________________//

class AHeraFitterConstr : public AConstraint {
 public:
  AHeraFitterConstr(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~AHeraFitterConstr() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd&) const;

 protected:
  Eigen::MatrixXd fSystErr;
};

//____________________________________________________________________________//

/*
class ASumRulesHeraFitter : public AConstraint {
 public:
  ASumRulesHeraFitter(const vector<string>& FitPar, AData* Data, AFuncD* Theo);
  virtual ~ASumRulesHeraFitter() {};
  virtual Eigen::VectorXd InitialX() const;
  virtual Eigen::VectorXd Constraint(const Eigen::VectorXd& X) const;
  virtual Eigen::MatrixXd Covariance(const Eigen::VectorXd&) const;

 protected:
  Eigen::MatrixXd fSystErr;
  int fgAind;
  int fuvAind;
  int fdvAind;
  int fUbarAind;
};
*/


#endif
