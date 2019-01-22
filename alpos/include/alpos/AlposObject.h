// D. Britzger 23.11.2014

#ifndef Alpos_AlposObject
#define Alpos_AlposObject

/**********************************************************
  AlposObject
  Base class of all alpos objects

***********************************************************/
#include <string>
#include "fastnlotk/speaker.h"

class AlposObject : public PrimalScream /*: public TObject, speaker */{
public:
   AlposObject(const std::string& classname = "") : PrimalScream(classname) {}; //! One may specify a classname for the 'speaker' tool.
   ~AlposObject(){};
   private:
};

#endif
