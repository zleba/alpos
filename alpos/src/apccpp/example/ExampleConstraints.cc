#include "Constraint.h"

/* 
   Some basic explicit examples for the use of Apccpp.
   The actual fit is implemented in the Constraint class, see
   ther general example of using Apccpp.

   D. Reichelt 09.2015
 */

using std::cout;
using std::endl;

class Pyth : public Constraint {
public:
  Pyth(double aIni, double bIni, double cIni, const Eigen::MatrixXd& iniV) 
    : a(aIni), b(bIni), c(cIni), emV(iniV) 
  {
    cout << endl;
    cout << "Pyth: enforcing pythagoras." << endl;
    cout << "---------------------------" << endl;
    cout << endl;
  }
  virtual ~Pyth() {}

  virtual Eigen::VectorXd evInitialX() const
  {
    return (Eigen::VectorXd(3) << a, b, c).finished();
  }
  virtual Eigen::MatrixXd emCovariance(const Eigen::VectorXd& X) const
  {
    return emV;
  }
  virtual Eigen::VectorXd evConstrVec(const Eigen::VectorXd& X) const
  {
    return (Eigen::VectorXd(1)  <<  X(0)*X(0)+X(1)*X(1)-X(2)*X(2)).finished();
  }
protected:
  double a;
  double b;
  double c;
  Eigen::MatrixXd emV;
};



class Avr : public Constraint {
public: 
  Avr(double aIni, double bIni) : a(aIni), b(bIni) 
  {
    cout << endl;
    cout << "Avr: averaging two poisson numbers." << endl;
    cout << "---------------------------" << endl;
    cout << endl;
  }
  virtual ~Avr() {};
   virtual Eigen::VectorXd evInitialX() const
  {
    return (Eigen::VectorXd(2) << a, b).finished();
  }
  virtual Eigen::MatrixXd emCovariance(const Eigen::VectorXd& X) const
  {
    return X.asDiagonal();
  }
  virtual Eigen::VectorXd evConstrVec(const Eigen::VectorXd& X) const
  {
    return (Eigen::VectorXd(1) << X(0) - X(1)).finished();
  }
protected:
  double a;
  double b;
};


class Mass : public Constraint {
public: 
  Mass(double m1Ini, double m2Ini, double mSumIni, Eigen::MatrixXd iniV) 
    : m1(m1Ini), m2(m2Ini), mSum(mSumIni), emV(iniV) 
  {
    cout << endl;
    cout << "Mass: Improvement of 3 mass measurements." << endl;
    cout << "-----------------------------------------" << endl;
    cout << endl;
  }
  virtual ~Mass() {};
   virtual Eigen::VectorXd evInitialX() const
  {
    return (Eigen::VectorXd(3) << m1, m2, mSum).finished();
  }
  virtual Eigen::MatrixXd emCovariance(const Eigen::VectorXd& X) const
  {
    return emV;
  }
  virtual Eigen::VectorXd evConstrVec(const Eigen::VectorXd& X) const
  {
    return (Eigen::VectorXd(1)  <<  X(0) + X(1) - X(2)).finished();
  }
protected:
  double m1;
  double m2;
  double mSum;
  Eigen::MatrixXd emV;
};


int main()
{
  
  // enforcing pythagoras
  Eigen::VectorXd evXpyth(3);
  double a=3.1; 
  double b=4.1; 
  double c=5.1;
  Eigen::MatrixXd emVpyth = (Eigen::VectorXd(3) 
                             << 0.1*0.1, 0.2*0.2, 0.1*0.1).finished().asDiagonal();
  Pyth pyth(a, b, c, emVpyth);
  pyth.Fit();
  pyth.apcFitter->PrintResult();
  cout << endl;
   
  // averaging poisson numbers
  double d=9.0; 
  double e=16.0;
  Avr avr(d, e);
  avr.Fit();
  avr.apcFitter->PrintResult();
  cout << endl;

  // improve mass from two measurements and measurement of the sum
  double m1 = 101;
  double m2 = 99;
  double mSum = 199;
  Eigen::MatrixXd emV = (Eigen::VectorXd(3) << 1.0, 1.0, 1.0).finished().asDiagonal();
  Mass mass(m1, m2, mSum, emV);
  mass.Fit();
  mass.apcFitter->PrintResult();
  cout << endl;


}
