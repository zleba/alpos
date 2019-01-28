***********************************************************************
* H1 2006 DPDF Fits parameterization                                  *
* ----------------------------------                                  *
* Reference: H1 Collaboration, A. Aktas et al., "Measurement and QCD  *
* Analysis of the Diffractive Deep-Inelastic Scattering Cross Section *
* at HERA", DESY06-049, hep-ex/0606004, subm. to Eur.Phys.J           *
*                                                                     *
* Contact persons in case of questions:                               *
* Frank-Peter Schilling (frank-peter.schilling@cern.ch)               *
* Paul Newman (newmanpr@mail.desy.de)                                 *
* ------------------------------------------------------------------- *
* This routine returns the pomeron or reggeon flux factors            *
* f_{IP,IR/p}(xpom,t)                                                 *
* Both IP and IR flux factors are normalized such that xpom*flux=1    *
* at xpom=0.003; In addition, the IR flux is multiplied by a further  *
* normalization parameter as determined in the fit to the data        *
* The routine returns either the flux at fixed values (xpom,t), or    *
* the t-integrated flux (integrated from t to tmin)                   *
* ------------------------------------------------------------------- *
* Input:  xpom:  xpomeron value                                       *
*         t:     t value (negative!)                                  *
*         int:   0: flux at (xpom,t) is returned                      *
*                1: t-integrated (t...tmin) flux at xpom is returned  *
*         ifit:  1: Fit A  2: Fit B                                   *
*         ipdf:  LHApdf set 0: central fit                            *
*                           1...32(30) error pdfs for fit A(B)        *
*         ipom:  1: Pomeron flux 2: Reggeon flux                      *
* Output: flux:  flux value                                           *
***********************************************************************

      subroutine h12006flux(xpom,t,int,ifit,ipdf,ipom,flux)

      implicit none

      integer int,ifit,ipdf,ipom
      double precision xpom,t,a0,ap,b0,c,flux,norm,dm,a0pom,nmes

      
cKC      write(*,*)'KAREL flux ipom,ifit ',ipom,' ',ifit

      if ((ipom.lt.1).or.(ipom.gt.2)) then
         print *,'[H12006FLUX] Unknown ipom: ',ipom
         stop
      endif
      if ((ifit.lt.1).or.(ifit.gt.2)) then
         print *,'[H12006FLUX] Unknown ifit: ',ifit
         stop
      endif

      call fluxpars(ifit,ipdf,a0pom,nmes)

      if (ipom.eq.1) then ! pomeron
        a0=a0pom
        ap = 0.06d0 
        b0 = 5.5d0 
        c = 1.0d0 
      else ! meson
        a0 = 0.5d0 
        ap = 0.3d0 
        b0 = 1.6d0 
        c=nmes
      endif

c     normalization
      call rflux2006(0.003d0,a0,ap,b0,-1.d0,1.0d0,1,dm)
      norm=(1./(0.003d0*dm))*c

c     actual flux
      call rflux2006(xpom,a0,ap,b0,t,norm,int,flux)

      return
      end

*******************************************************************************

      subroutine rflux2006(x_pom,a0,ap,b0,tcut,c,int,fl)

      implicit none

      double precision x_pom,a0,tmin,tcut,ap,b0,c,fl,b
      integer int

c     calc min. kinematically  allowed t
      tmin= -(0.93827231D0*x_pom)**2/(1.D0-x_pom)

c     c*xpom**(-(2apom-1))
      fl = c * dexp((2.0d0*a0-1.)*dlog(1.0d0/x_pom))
      b=(b0+2.0d0*ap*dlog(1.0d0/x_pom))

      if (int.eq.0) then
c       at fixed t:  exp(Bt) 
        fl = fl * dexp(b*tcut)
      else
c       t integrated: (1/B)*[exp(-B*tmax)-exp(-B*tmin)]
        fl = fl * (dexp(tmin*b)-dexp(tcut*b))/b
      endif

      return
      end

*******************************************************************************

      subroutine fluxpars(ifit,ipdf,a0pom,nmes)

      double precision a0pom,nmes
      integer ipdf,ifit

      if (ifit.eq.1) then ! fit a

      if(ipdf.ge.0.and.ipdf.lt.17) then 
        a0pom=1.11824
        nmes=0.00169662
      elseif (ipdf.eq.17) then
        a0pom=1.11843
        nmes=0.00169287
      elseif (ipdf.eq.18) then
        a0pom=1.11804
        nmes=0.00170004
      elseif (ipdf.eq.19) then
        a0pom=1.11687
        nmes=0.00169226
      elseif (ipdf.eq.20) then
        a0pom=1.11854
        nmes=0.00168107
      elseif (ipdf.eq.21) then
        a0pom=1.11833
        nmes=0.0016943
      elseif (ipdf.eq.22) then
        a0pom=1.11819
        nmes=0.00169873
      elseif (ipdf.eq.23) then
        a0pom=1.11936
        nmes=0.00169188
      elseif (ipdf.eq.24) then
        a0pom=1.11734
        nmes=0.00170108
      elseif (ipdf.eq.25) then
        a0pom=1.14722
        nmes=0.00159724
      elseif (ipdf.eq.26) then
        a0pom=1.10871
        nmes=0.00170391
      elseif (ipdf.eq.27) then
        a0pom=1.11546
        nmes=0.000976152
      elseif (ipdf.eq.28) then
        a0pom=1.12169
        nmes=0.00299953
      elseif (ipdf.eq.29) then
        a0pom=1.11657
        nmes=0.00131429
      elseif (ipdf.eq.30) then
        a0pom=1.12065
        nmes=0.00245839
      elseif (ipdf.eq.31) then
        a0pom=1.11691
        nmes=0.00179416
      elseif (ipdf.eq.32) then
        a0pom=1.11855
        nmes=0.00165973
      else 
        print*,'[H12006flux] ipdf out of range:',ipdf
        stop
      endif 

      elseif (ifit.eq.2) then ! fit B

      if(ipdf.ge.0.and.ipdf.lt.15) then 
        a0pom=1.11101
        nmes=0.00139764
      elseif (ipdf.eq.15) then
        a0pom=1.11236
        nmes=0.00142534
      elseif (ipdf.eq.16) then
        a0pom=1.10965
        nmes=0.00137153
      elseif (ipdf.eq.17) then
        a0pom=1.10988
        nmes=0.00140738
      elseif (ipdf.eq.18) then
        a0pom=1.11129
        nmes=0.00138136
      elseif (ipdf.eq.19) then
        a0pom=1.11102
        nmes=0.00139576
      elseif (ipdf.eq.20) then
        a0pom=1.11101
        nmes=0.00139999
      elseif (ipdf.eq.21) then
        a0pom=1.11212
        nmes=0.00139274
      elseif (ipdf.eq.22) then
        a0pom=1.11009
        nmes=0.00140232
      elseif (ipdf.eq.23) then
        a0pom=1.13989
        nmes=0.00129396
      elseif (ipdf.eq.24) then
        a0pom=1.10147
        nmes=0.00140518
      elseif (ipdf.eq.25) then
        a0pom=1.10927
        nmes=0.000811837
      elseif (ipdf.eq.26) then
        a0pom=1.11307
        nmes=0.00244171
      elseif (ipdf.eq.27) then
        a0pom=1.10997
        nmes=0.00108895
      elseif (ipdf.eq.28) then
        a0pom=1.11245
        nmes=0.0020084
      elseif (ipdf.eq.29) then
        a0pom=1.10772
        nmes=0.00134241
      elseif (ipdf.eq.30) then
        a0pom=1.1146
        nmes=0.00147194
      else 
        print*,'[H12006flux] ipdf out of range:',ipdf
        stop
      endif 

      endif

      return
      end

