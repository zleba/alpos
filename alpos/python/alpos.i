/* alpos.i */
%module alpos

%{
#define SWIG_FILE_WITH_INIT
#include "alpos/Alpos.h"
#include "alpos/AlposObject.h"
%}

// SWIG wrappers for STL
%include <std_string.i>
%include <std_map.i>
%include <std_vector.i>
%include <std_pair.i>


// Alpos core
%include "alpos/AlposObject.h"
%include "alpos/AError.h"
%include "alpos/Alpos.h"
