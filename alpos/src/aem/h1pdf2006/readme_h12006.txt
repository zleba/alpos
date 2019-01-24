README for H1 2006 DPDF Fit Parameterization
============================================

Version 1.0, 16/06/2006

H1 Collaboration

Contact persons:
- Frank-Peter Schilling (frank-peter.schilling@cern.ch)
- Paul Newman (newmanpr@mail.desy.de)

WWW source:

http://www-h1.desy.de/h1/www/publications/htmlsplit/DESY-06-049.long.html

These files provide a fortran parameterization of the
results of the NLO DGLAP QCD fits labelled 'H1 2006 DPDF Fit A' and
'H1 2006 DPDF Fit B' from the paper

   H1 Collaboration, A. Aktas et al., "Measurement and QCD Analysis of
   the Diffractive Deep-Inelastic Scattering Cross Section at HERA",
   DESY06-049, hep-ex/0606004, subm. to Eur.Phys.J

Please reference to this paper when using this parametrization.

Files provided
--------------

readme_h12006.txt   This readme file
h12006flux.f        Pomeron/Reggeon flux factor
qcd_2006.f          Pomeron DPDFs and Pomeron/Reggeon structure functions
i_2006_fita.f       Grid data for H1 2006 DPDF Fit A
i_2006_fitb.f       Grid data for H1 2006 DPDF Fit B


Usage
=====

For details on the arguments of the provided routines, please
also look at the comments provided in the .f files


a. Flux factor
--------------

Subroutine h12006flux provides the Pomeron or Reggeon flux factor

   f_{IP,IR/p}(xpom,t), 

as defined in Equation 14 of the paper. The routine may provide either
the flux at particular values of xpom and t, or alternatively provide
the flux at xpom, integrated over t in the range t...tmin.  The
convention for the normalization of the t-integrated flux is

   f_{IP,IR/p}(xpom=0.003,t=-1...tmin)*xpom = 1. 

In case of the Reggeon, there is a further normalization factor by
which the flux is multiplied, as determined from the QCD fit.


b. Parton distributions
-----------------------

The Pomeron parton distributions 

   z*f_{i/IP}(z,Q^2) 

are provided by the subroutine qcd_2006. The Reggeon parton
distributions are not provided here but should be taken from PDFLIB
instead (in the fits PDFLIB-Owens(2,1,1) was used).  The pomeron
partons are provided on a (z,Q^2) grid in the range 0.001<z<1.0,
1<Q^2<30000 GeV^2. They should be used with care beyond the fit range
of 0.0043<z<0.8, 8.5<Q^2<1600 GeV^2.  The evolution to Q^2 values
above and below the fit range is performed using the NLO DGLAP
equations.

In order to obtain proton diffractive PDFs, the Pomeron/Reggeon parton
distributions have to be multiplied by the Pomeron/Reggeon flux factor
(see a.), so that

   f_{i/p}(xpom,t,z,Q^2) = f_{IP,IR/p}(xpom,t) * f_{i/IP,IR}(z,Q^2)


c. Structure functions
----------------------

Subroutine qcd_2006 also provides Pomeron and Reggeon
structure functions, calculated at NLO in the MS-bar scheme:

- F_2(beta,Q^2)
- F_L(beta,Q^2)   (longitudinal structure function)
- F^c_2(beta,Q^2) (charm structure function)
- F^c_L(beta,Q^2) (longitudinal charm structure function)

In order to obtain the diffractive proton structure functions,
again these Pomeron/Reggeon structure functions need to be 
multiplied by the corresponding flux factor such that

   F_2-IP(4)(xpom,t,beta,Q^2) = f_{IP/p}(xpom,t) * F_2-IP(beta,Q^2)

and

   F_2-IR(4)(xpom,t,beta,Q^2) = f_{IR/p}(xpom,t) * F_2-IR(beta,Q^2)

Finally, the full diffractive structure function F_2^D(4) is
obtained by adding the Pomeron and Reggeon contributions

   F_2^D(4)(xpom,t,beta,Q^2) = F_2-IP(4) + F_2-IR(4)

Similar relations hold for the t-integrated structure function F_2^D(3).


d. Diffractive Reduced Cross Section
------------------------------------

The diffractive reduced cross section which corresponds
to the measured data points from the paper, is defined as (Equation 5
of the paper):

   sigma_r^D(3) = F_2^D(3) - y^2/Y_+ * F_L^D(3)

where Y_+ = 1+ (1-y)^2.

