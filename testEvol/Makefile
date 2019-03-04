
LIBAPFEL=../apfel/lib/
INCAPFEL=../apfel/include/
LHALIB=$(shell lhapdf-config --libdir)
LHAINC=$(shell lhapdf-config --incdir)

test: test.cpp
	g++ -I${INCAPFEL} -I{LHAINC}  $^ -o $@ -L${LIBAPFEL} -lAPFEL -Wl,-rpath,${LIBAPFEL}  -L${LHALIB}  -lLHAPDF
