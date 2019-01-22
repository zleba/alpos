// DB. 01/2015
#ifndef Alpos_AChisq
#define Alpos_AChisq

/** 
 * classes for chisq minimization
 * 
 * Various calculations of a chisq used
 * by AFitter
 * 
 */

#include "alpos/AlposObject.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "alpos/AData.h"

// ROOT
#include "Math/IParamFunction.h"

// Eigen (for apccpp)
#include <Eigen/Eigen>

//____________________________________________________________________________________ //
/** 
 * A base class for chisq functions
 */
class AChisqBase : public ROOT::Math::IMultiGenFunction, public AlposObject {
public:
   AChisqBase(const std::vector<std::string>& FitPar, AData* data=NULL, AFuncD* theo=NULL ) : fFitPar(FitPar), fData(data), fTheo(theo) { fNuisance.clear();fNuisanceName.clear();};
   void SetDataTheory(AData* data, AFuncD* theo) {fData = data; fTheo=theo;};
   virtual ~AChisqBase() {};
   virtual unsigned int NDim() const = 0; //!< Return number of fit parameters
   //virtual ROOT::Math::IMultiGenFunction* Clone() const = 0 ;// { return new AChisqBase(data); }; //!< You must implement a 'Clone' method
   const std::vector<std::string>& GetFitPar() { return fFitPar;}
   AData* Data() const { return fData;}
   AFuncD* Theo() const { return fTheo;}
   //private:
   virtual double DoEval(const double *p) const = 0; //!< return chisq for a given set of fit-parameters p
   static std::map<std::string, double> GetNuisanceParameters(); //!< Get recent calculation of nuisance parameters
   void PrintNuisanceParameters(bool printnone=false, std::map<std::string, double>* nuisance=NULL) const; //!< print nuisance parameters

protected:
   std::vector<std::string> fFitPar;
   AData*  fData=NULL;
   AFuncD* fTheo=NULL;
   static std::vector<double> CalcNuisance(const std::vector<std::vector<double> >& gk, const TMatrixD& VInv, const std::vector<double>& dmt) ;
   static std::vector<double> CalcNuisance(const std::vector<const std::vector<double>* >& gk, const TMatrixD& VInv, const std::vector<double>& dmt);
   static std::vector<double> CalcNuisance(const std::vector<std::vector<double> >& gk, const std::vector<double>& VDiag, const std::vector<double>& dmt);
   static std::vector<double> CalcNuisance(const std::vector<const std::vector<double>* >& gk, const std::vector<double>& VDiag, const std::vector<double>& dmt);
   
   static std::vector<double> SolveStoreNuisance(const TMatrixDSym& BI, const std::vector<double> a);

   static std::vector<double> fNuisance; //!< Nuisance parameters 
   static std::vector<std::string> fNuisanceName; //!< Names of nuisance parameters 
};



//____________________________________________________________________________________ //
/** 
 * A chisq using the covariance matrix
 * 
 *  chisq  = ( d - t ) V^-1  ( d-t ) 
 *
 */
class AChisqCov : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "Covariance"; return name;}
   AChisqCov(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) { };
   virtual ~AChisqCov() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqCov(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};

//____________________________________________________________________________________ //
/** 
 * A chisq using the covariance matrix
 * 
 *  chisq  = ( d - t ) V^-1  ( d-t ) 
 *  where V is composed from multiplicative and additive uncertainties
 *  HOWEVER! uncorrelated and statistical uncertainties are always treated as additive !
 *  Theoretical uncertainties are treated as 'additive' as well, i.e. rel*theo.
 *
 */
class AChisqCMS : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "ChisqCMS"; return name;}
   AChisqCMS(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) { };
   virtual ~AChisqCMS() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqCMS(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};

//____________________________________________________________________________________ //
/** 
 * A chisq using the covariance matrix
 * 
 *  chisq  = ( d - t ) V^-1  ( d-t ) 
 *  where V is composed only from the stat and uncorr uncertainties
 *
 */
class AChisqCovStatUncorr : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "CovStatUncorr"; return name;}
   AChisqCovStatUncorr(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) { };
   virtual ~AChisqCovStatUncorr() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqCovStatUncorr(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};

//____________________________________________________________________________________ //
/**
 * A chisq using the covariance matrix
 *
 *  chisq  = ( d - t ) V^-1  ( d-t )
 *
 * where multiplicative contributions to V are rescaled to the theory
 *
 */
class AChisqCovMult : public AChisqBase {
public:
   static const std::string& GetChisqName() { static const std::string name = "CovarianceMult"; return name;}
   AChisqCovMult(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) { };
   virtual ~AChisqCovMult() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqCovMult(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in [arxiv:1406.4709] 
 */
class AChisqLogNormal : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "LogNormal"; return name;}
   AChisqLogNormal(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqLogNormal() {};
   virtual unsigned int NDim() const; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const;// { return new AChisqBase(data); }; //!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};


//____________________________________________________________________________________ //
/** 
 * An updated version as used in [arxiv:1406.4709] 
 * to allow further for normal-distributed uncertainties
 * These are assumed to be datapoints WITHOUT statistical
 * It is further assumed that these datapoints are uncorrelated to other
 * log-normal distributed data points.
 */
class AChisqNormalLogNormal : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "NormalLogNormal"; return name;}
   AChisqNormalLogNormal(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {;};
   virtual ~AChisqNormalLogNormal() {};
   virtual unsigned int NDim() const {return fFitPar.size();}; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const{ return new AChisqNormalLogNormal(*this);};//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in [arxiv:1406.4709] 
 * where only stat. and uncorr uncertainties are considered
 */
class AChisqLogNormalStatUncorr : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "LogNormalStatUncorr"; return name;}
   AChisqLogNormalStatUncorr(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
   virtual ~AChisqLogNormalStatUncorr() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqLogNormalStatUncorr(*this); };// { return new AChisqBase(data); }; //!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in [arxiv:1406.4709]
 */
class AChisqLogNormalNuisance : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "LogNormalNuisance"; return name;}
   AChisqLogNormalNuisance(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
   virtual ~AChisqLogNormalNuisance() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqLogNormalNuisance(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in [arxiv:1406.4709] 
 */
class AChisqLogNormalNuisanceFit : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "LogNormalNuisanceFit"; return name;}
   AChisqLogNormalNuisanceFit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqLogNormalNuisanceFit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqLogNormalNuisanceFit(*this);}
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
protected:
   std::vector<std::string> fTheoPar;
};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in HERAFitter, where nuisance parameters are determined by minuit
 */
class AChisqHERAFitterDefaultFit : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterDefaultFit"; return name;}
   AChisqHERAFitterDefaultFit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqHERAFitterDefaultFit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterDefaultFit(*this);}
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
protected:
   std::vector<std::string> fTheoPar;

};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in HERAFitter, where nuisance parameters are calculated in a two-step
 * analytic approach
 */
class AChisqHERAFitterDefault : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterDefault"; return name;}
   AChisqHERAFitterDefault(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
   virtual ~AChisqHERAFitterDefault() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterDefault(*this);}
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
};
//! Same class as above, but where statiscal uncertainties correlations are considered (if present)
class AChisqHERAFitterDefaultMatrix : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterDefaultMatrix"; return name;}
   AChisqHERAFitterDefaultMatrix(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
   virtual ~AChisqHERAFitterDefaultMatrix() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterDefaultMatrix(*this);}
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in HERAFitter, where nuisance parameters are calculated in a two-step
 * analytic approach
 */
// class AChisqHERAFitterLogDefault : public AChisqBase { 
// public:
//    static const std::string& GetChisqName() { static const std::string name = "HERAFitterLogDefault"; return name;}
//    AChisqHERAFitterLogDefault(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
//    virtual ~AChisqHERAFitterLogDefault() {};
//    virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
//    virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterLogDefault(*this);}
//    virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
// };
class AChisqHERAFitterLogDefault : public AChisqHERAFitterDefault { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterLogDefault"; return name;}
   AChisqHERAFitterLogDefault(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqHERAFitterDefault(FitPar,data,theo) {};
   virtual ~AChisqHERAFitterLogDefault() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterLogDefault(*this);}
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in HERAFitter [arxiv:1406.4709] 
 * This definition requires the specification of the
 * error 'Nature' as 'A,M,P'. 
 */
class AChisqHERAFitterFull : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterFull"; return name;}
   AChisqHERAFitterFull(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqHERAFitterFull() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterFull(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   std::vector<std::string> fTheoPar;

};


//____________________________________________________________________________________ //
/** 
 * A chisq as used in HERAFitter [arxiv:1406.4709] 
 * It rescalaes multiplicative, additive and poisson errors correctly.
 * This definition requires the specification of the
 * error 'Nature' as 'A,M,P'. 
 */
class AChisqHERAFitterFullImproved : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "HERAFitterFullImproved"; return name;}
   AChisqHERAFitterFullImproved(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqHERAFitterFullImproved() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqHERAFitterFullImproved(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   std::vector<std::string> fTheoPar;

};


//____________________________________________________________________________________ //
/** 
 * A chisq where all uncertainties are treated as uncorrelated and 'additive'
 */
class AChisqSimple : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "Simple"; return name;}
   AChisqSimple(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo) : AChisqBase(FitPar,data,theo) {};
   virtual ~AChisqSimple() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqSimple(*this); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:

};


//____________________________________________________________________________________ //
/** 
 * A chisq using nuisance parameters, ignoring correlation matrices
 * and given the same (biased) result as the simple covariance chisq.
 */
class AChisqSimpleNuisanceAddFit : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "SimpleNuisanceAddFit"; return name;}
   AChisqSimpleNuisanceAddFit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqSimpleNuisanceAddFit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqSimpleNuisanceAddFit(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   std::vector<std::string> fTheoPar;

};


//____________________________________________________________________________________ //
/** 
 * A simple chisq using nuisance parameters, ignoring correlation matrices
 * Mind that uncorrelated uncertainties are treated as additive, while
 * correlated uncertainties are treated as multiplicative
 *
 * Nuisance parameters are free parameter for the (minuit) fit
 */
class AChisqSimpleNuisanceMultFit : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "SimpleNuisanceMultFit"; return name;}
   AChisqSimpleNuisanceMultFit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqSimpleNuisanceMultFit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqSimpleNuisanceMultFit(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
protected:
   std::vector<std::string> fTheoPar;
};


//____________________________________________________________________________________ //
/** 
 * A chisq using nuisance parameters
 * and resulting in the same (biased) result as the simple covariance chisq.
 */
class AChisqNuisanceAdd : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "NuisanceAdd"; return name;}
   AChisqNuisanceAdd(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqNuisanceAdd() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqNuisanceAdd(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   //std::vector<std::string> fTheoPar;

};


//____________________________________________________________________________________ //
/** 
 * A simple chisq using nuisance parameters, ignoring correltion matrices
 * and given the same (biased) result as the simple covariance chisq.
 */
class AChisqNuisanceMult : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "NuisanceMult"; return name;}
   AChisqNuisanceMult(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqNuisanceMult() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqNuisanceMult(*this);}//}fFitPar,fData,fTheo); };//!< You must implement a 'Clone' method
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
};


//____________________________________________________________________________________ //
/** 
 *
 * A chisq as used by D0 alpha_s fits from jets and in ATLAS dijet azimuthal decorrelations
 *
 * This chisq takes exp. and theo. uncertainties into account.
 *
 * The nuisance parameter (for both: theoretical and exp. uncertainties) are fitted
 *
 */
class AChisqD0Fit  : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "D0Fit"; return name;}
   AChisqD0Fit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqD0Fit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqD0Fit(*this);}

   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   std::vector<std::string> fTheoPar;

};

//____________________________________________________________________________________ //
/**
 *
 * Analogous to 'D0Fit', but takes statistical correlations into account.
 *
 * The nuisance parameters (for both theoretical and exp. uncertainties) are fitted
 *
 */
class AChisqD0StatCorrFit : public AChisqBase {
public:
   static const std::string& GetChisqName() { static const std::string name = "D0StatCorrFit"; return name;}
   AChisqD0StatCorrFit(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqD0StatCorrFit() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqD0StatCorrFit(*this);}

   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p

protected:
   std::vector<std::string> fTheoPar;
};

//____________________________________________________________________________________ //
/**
 * A chisq based on nuisance parameter calculation by Apcalc
 */
class AChisqApc : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "ApcNuisance"; return name;}
   AChisqApc(const std::vector<std::string>& FitPar, AData* data, AFuncD* theo);
   virtual ~AChisqApc() {};
   virtual unsigned int NDim() const { return fFitPar.size(); }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new AChisqApc(*this);};//!< You must implement a 'Clone' method
   //private:
   virtual double DoEval(const double *p) const; //!< return chisq for a given set of fit-parameters p
protected:
   Eigen::MatrixXd emNuisParam;
};



//____________________________________________________________________________________ //
/** 
 * Calculate the mean value of the pull distribution (ignoring correlations between data points)
 */
class APull : public AChisqBase { 
public:
   static const std::string& GetChisqName() { static const std::string name = "Pull"; return name;}
   APull(AData* data, AFuncD* theo)  : AChisqBase(std::vector<std::string>(),data,theo) {};
   virtual ~APull() {};
   virtual unsigned int NDim() const { return 0; }; //!< Return number of fit parameters
   virtual ROOT::Math::IMultiGenFunction* Clone() const { return new APull(*this);}
   double CalcMeanStatUncorr() const; //!< Calculate pull value, using the stat+uncorr uncertainties
   double CalcMeanTotErr() const; //!< Calculate pull value, using the stat+uncorr+corr uncertainties
   double CalcMedianStatUncorr() const; //!< Calculate meadian of pull distribution, using the stat+uncorr uncertainties
   double CalcMedianTotErr() const; //!< Calculate meadian of pull distribution, using the stat+uncorr+corr uncertainties
   double CalcRMSStatUncorr() const; //!< Calculate RMS, using the stat+uncorr uncertainties
   double CalcRMSTotErr() const; //!< Calculate RMS , using the stat+uncorr+corr uncertainties
   double CalcMaxStatUncorr() const; //!< Calculate max value
   double CalcMaxTotErr() const; //!< Calculate max value
   double CalcMinStatUncorr() const; //!< Calculate min value
   double CalcMinTotErr() const; //!< Calculate min value
   //
   // .. function including theory uncertainties may be reasonable
   //
 
protected:
   virtual double DoEval(const double *p) const { return 0;};  //!< Don't call DoEval
   std::vector<std::string> fTheoPar;

};



#endif
