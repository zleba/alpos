      subroutine diffpdferr(xpom,zpom,muf,ifit,ierr,ireg,pdfs)
      implicit none
      integer ierr,ifit,int,ipom,ireg,calling,nul
      double precision beta,muf,pt2max,pt2min,xpom,pdfs(-6:6)
      double precision zpom
      double precision srh
      double precision xz

      double precision t,flux,qq2,x,dxpom,my,xpbin
      double precision XPQ(-6:6),F2(1:2),FL(1:2),C2(1:2),CL(1:2)
      double precision xpmin,xpmax,Pflux,Rflux
      double precision UPV,DNV,SEA,STR,CHM,GL
      external STROWP1


c ------------------------------
c     set consistent names
      beta=zpom
      x=zpom
      xz=xpom*zpom

c     The scale , for example, is set to the minimum accessible 
c     in VFPS analysis. 
         
c     the bq(-nf:nf) array contains beta*distribution 

c     ifit=1                    ! Fit type A
c     ifit=2                    ! Fit type B
      t=-1d0                    ! maximum t
      int=1                     ! t-integrated flux [t..tmin]


c ------------------------------
c     Get Pomeron flux in the proton. It depends on xpom & t but not on beta and Q2

c
c   For the moment only central value (jpdf=0)
c   From pdf-cteq6_err.h
c     h12006flux(&xp, &t_diff,&iint, &ifit, &ipdf, &ipom, &PomFlux);
c     h12006pdf(&z, &q2, &ifit,&ipdf,xpq,f2,fl,c2,cl);

      if(ireg.eq.2) then 
         Pflux=0
      else
         ipom=1
         call h12006flux(xpom,t,int,ifit,ierr,ipom,flux);
         Pflux=flux	    
c     get pomeron PDF's
         qq2=muf*muf
         call h12006pdf(x,qq2,IFIT,ierr,XPQ,F2,FL,C2,CL)
c         print *,"DIFFPDF: flxIP: ",Pflux
c         print *,"DIFFPDF: gluon: ",XPQ(0)
c         print *,"DIFFPDF: up   : ",XPQ(1)
c         print *,"DIFFPDF: dn   : ",XPQ(2)
      endif

      


c ------------------------------
c     Get Reggeon flux in the proton. It depends on xpom & t but not on beta and Q2 

c      ipom=2; //1-pomeron, 2-reggeon
c      int jpdf=0;
c      h12006flux_(&xp, &t_diff,&iint, &ifit, &jpdf, &ipom, &RegFlux);


      if(ireg.eq.1) then
         Rflux=0
      else
         ipom=2
         nul = 0
         call h12006flux(xpom,t,int,ifit,nul,ipom,flux)
         Rflux=flux	    
      endif


c ------------------------------
c     call stand-alone pion structure function by Owen
c     pion pdf's are used for the reggeon term, best we can do 

      CALL STROWP1(x,muf,UPV,DNV,SEA,STR,CHM,GL)


c ------------------------------
c     various contributions are added together 
c     final DPDF(beta) = pomeron_pdf * Pomflux + reggeon_pdf * Regflux

c     There is a scheme inconsistency between the massless nf=4 
c     scheme employed by NLOjet++ and the FFNS nf=3 of FitB DPDF's

c     Charm contribution are obtained from charm structure functions F2c, 
c     here is called C2(1), and returned by H1 fit 2006 
c     ( which is performed in FFNS scheme, nf=3)

c DB: not updated to h12006pdf_ERR_

      pdfs(-6)= 0d0
      pdfs(-5)= 0d0	
      pdfs(-4)= XPQ(-4)*Pflux+(CHM)*RFLUX
      pdfs(-3)= XPQ(-3)*Pflux+(STR)*RFLUX
      pdfs(-2)= XPQ(-2)*Pflux+(SEA)*RFLUX
      pdfs(-1)= XPQ(-1)*Pflux+(SEA)*RFLUX
      pdfs(0) = XPQ(0)*Pflux+(GL)*RFLUX
      pdfs(1) = XPQ(1)*Pflux+(UPV+SEA)*RFLUX
      pdfs(2) = XPQ(2)*Pflux+(DNV+SEA)*RFLUX
      pdfs(3) = XPQ(3)*Pflux+(STR)*RFLUX
      pdfs(4) = XPQ(4)*Pflux+(CHM)*RFLUX
      pdfs(5) = 0d0
      pdfs(6) = 0d0


c      pdfs(-6)= XPQ(-6)*Pflux
c      pdfs(-5)= XPQ(-5)*Pflux
c      pdfs(-4)= (CHM)*RFLUX+XPQ(-4)*Pflux
c      pdfs(-3)= (STR)*RFLUX+XPQ(-3)*Pflux
c      pdfs(-2)= (SEA)*RFLUX+XPQ(-2)*Pflux
c      pdfs(-1)= (SEA)*RFLUX+XPQ(-1)*Pflux
c      pdfs(0) = (GL)*RFLUX+XPQ(0)*Pflux
c      pdfs(1) = (UPV+SEA)*RFLUX+XPQ(1)*Pflux
c      pdfs(2) = (DNV+SEA)*RFLUX+XPQ(2)*Pflux
c      pdfs(3) = (STR)*RFLUX+XPQ(3)*Pflux
c      pdfs(4) = (CHM)*RFLUX+XPQ(4)*Pflux
c      pdfs(5) = XPQ(5)*Pflux
c      pdfs(6) = XPQ(6)*Pflux

      return
      end

      
cc       INCLUDE 'Boris/ddis2jet/qcd_2006.f'
cc       INCLUDE 'Boris/ddis2jet/h12006flux.f'
cc       INCLUDE 'Boris/ddis2jet/i_2006_fita.f'       
cc       INCLUDE 'Boris/ddis2jet/i_2006_fitb.f' 	
cc       INCLUDE 'pion_stand_alone.f'
cc c      INCLUDE 'Boris/ddis2jet/strowp1.f'
       INCLUDE 'Boris/ddis2jet/h12006flux.err.f'
       INCLUDE 'Boris/ddis2jet/lha2006fita.f'
       INCLUDE 'Boris/ddis2jet/lha2006fitb.f'

       INCLUDE 'Boris/ddis2jet/h12006pdf.f'
c===================================================================
      
