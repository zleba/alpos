//
// Created by DS on 4/15/16.
//

#ifndef ALPOS_ATHEORYHANDLER_H
#define ALPOS_ATHEORYHANDLER_H

#include <set>
#include <vector>
#include <utility>

#include <alpos/AlposObject.h>

/**********************************************************
 *
 * ATheorySet
 *
 *  A full set of all theory parameters
 * It can be obtained and set to the TheoryHandler.
 *
 ***********************************************************/

struct AThCont {
   std::vector<std::string> Values;
   std::vector<std::string> Errors;
   bool Const;
};

class ATheorySet {
public:
   ATheorySet() { ; };

   ~ATheorySet() { ; };

   bool AddContent(const std::string& aname, const AThCont& cont) {
      if (fCont.count(aname) > 0) return false;
      fCont[aname] = cont;
      return true;
   }

   const std::map<std::string, AThCont>& GetSet() const { return fCont; } //!< Get full set of parameters
   AThCont GetContent(const std::string& aname) { return fCont[aname]; } //!< Get one content

private:
   std::map<std::string, AThCont> fCont;
};

// ------------------------------------------------------------------------------------

// TheoryHandler needs some forward declarations
class AParmNamed;
template<typename T> class AParm;
template<typename T> class AParmFuncBase;

class AData;
class ASuperData;
class ASuperTheory;
class ASubsetData;
class ASubsetFunction;

// shortcuts to templated types
typedef AParm<double> AParmD;
typedef AParm<std::string> AParmS;
typedef AParm<int> AParmI;
typedef AParmFuncBase<double> AFuncD;

/**
 TheoryHandler

 Global singleton class to manage the theory calculations and
 provide functions to initialize theory, data, subset and
 super-functions.

 The theory handler also owns all theory parameters and
 serves as a central "repository" for these. It provides
 an interface to retrieve the parameters by name (as AParm...
 objects).

*/

class TheoryHandler : public AlposObject {

private:
   TheoryHandler();

   ~TheoryHandler();

public:
   static TheoryHandler* Handler();

   bool InitTheory(const std::string& steerfile = std::string());

   bool InitSuperFunctions(); //!< init SuperTheory and SuperData

   bool InitSubsetFunctions(const std::vector<std::string>& requirements, bool DoCuts,
                            bool DoSubsets); //!< init Subset function for data and theory

   std::set<std::string> GetRequirements(
         const std::string& functype); //!< Get requirements of a function (from 'defaults' in steering!)

   bool RegisterParameter(AParmNamed*);

   bool EraseParameter(const std::string& name);

   bool CheckParameter(const std::string& name) const;

   AParmNamed* GetParameter(const std::string& name);

   AParmD* GetParmD(const std::string& name) {
      return (AParmD*)(GetParameter(name));
   }; //!< Get a ParmD parameter (todo: adjust type-casting)

   AParmS* GetParmS(const std::string& name) {
      return (AParmS*)(GetParameter(name));
   }; //!< Get a ParmS parameter (todo: adjust type-casting)

   AFuncD* GetFuncD(const std::string& name) {
      return (AFuncD*)(GetParameter(name));
   }; //!< Get a ParmD parameter (todo: adjust type-casting)

   bool NewAlias(const std::string& alias, const std::string& aliasfor);

   // quick access
   double GetParmDValue(const std::string& name); //!< Get the value of a parameter of type (double)
   double GetParmDError(const std::string& name); //!< Get the error of a parameter of type (double)

   const std::string& GetParmSValue(const std::string& name); //!< Get the value of a parameter of type (string)
   const std::string& GetParmSError(const std::string& name); //!< Get the error of a parameter of type (string)

   std::string FindParameter(const std::string& name) const; //!< Find a parameter's base name
   bool ResolveName(std::string& name) const; //!< Resolve the parameter's name into its base name
   const std::set<std::string>& GetListOfAllParm() const { return fAllParmNames; }; //!< get a set of all (future) existing parameter names
   void PrintListOfAllParms() const; //!< Print the names of all available parameters

   ATheorySet GetTheorySet() const; //!< Get the values and errors of all parameters
   void SetTheorySet(const ATheorySet& set); //!< Set a full set of theory parameters
   void PrintCurrentTheorySet(ATheorySet* set = nullptr) const; //!< Print the names of all available parameters
   const std::pair<ASuperData*, ASuperTheory*>& GetSuperPair() const { return fSuperPair; } //!< Get pointers to Superdata and Supertheory
   const std::map<std::string, std::pair<AData*, AFuncD*> >& GetDataTheoryPairs() const { return fDataTheoryPairs; }; //!< Get pointer to all data-theory pairs (ordered by name)
   const std::map<std::string, std::map<std::string, std::pair<ASubsetData*, ASubsetFunction*> > >& GetSubsetPairs() const { return fSubsetPairs; }; //!< Get pointers to all subsets (ordered by dataset name, and subsetname)

   std::map< std::string, std::pair<ASubsetData*, ASubsetFunction*> > GetAllSubsetPairs() const {
      std::map< std::string, std::pair<ASubsetData*, ASubsetFunction*> > ret;
      for (const auto& id : fSubsetPairs)
         for (const auto& is : id.second)
            ret[is.first] = is.second;
      return ret;
   }; //!< Get pointers to all subsets (but ordered (only) by subsetname)


private:
   static TheoryHandler* instance;
   std::map<std::string, AParmNamed*> fParms;
   std::set<std::string> fAllParmNames;

   bool CheckFunction(
         AFuncD* func); //!< Check if function defaults are available and add function to list of parameters
   //useful map of datasets and their corresponding predictions
   std::pair<ASuperData*, ASuperTheory*> fSuperPair;
   std::map<std::string, std::pair<AData*, AFuncD*> > fDataTheoryPairs;
   std::map<std::string, std::map<std::string, std::pair<ASubsetData*, ASubsetFunction*> > > fSubsetPairs;
};

#endif //ALPOS_ATHEORYHANDLER_H
