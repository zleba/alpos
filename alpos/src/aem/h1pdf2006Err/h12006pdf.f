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
* This routine returns the pomeron dpdf's (z*pdf(z,Q^2)) and the      * 
* pomeron and reggeon structure functions F2_{IP,IR}(beta,Q^2),       *
* FL_{IP,IR}(beta,Q^2), F2c_{IP,IR}(beta,Q^2), FLc_{IP,IR}(beta,Q^2)  *
* Range of validity: 0.0043 < z < 0.8 ; 8.5 < Q^2 < 1600 GeV^2        *
* Outside, the pdf's are in z an extrapolation beyond the measured    *
* range and in Q^2 evolved to higher/lower scales using the NLO DGLAP *
* equations                                                           *
* The grids are provided for 0.001 < z < 1.0 ; 1 < Q^2 < 30000 GeV^2  *
* ------------------------------------------------------------------- *
* Input:  z:     z-pom (parton momentum fraction)                     *
*         q2:    scale^2 where dpdf's are evaluated (in GeV^2)        *
*         ifit:  1: Fit A  2: Fit B                                   *
*         ipdf:  LHApdf set 0: central fit                            *
*                           1...32(30) error pdfs for fit A(B)        *
* Output: xpq(-6:6):  PDG style array of pomeron dpdf's at (z,Q2)     *
*      the following are provided only for ipdf=0                     *
*         f2(2):      Structure function F2  (1=pomeron,2=reggeon)    *
*         fl(2):      Structure function FL  (1=pomeron,2=reggeon)    *
*         c2(2):      Structure function F2c (1=pomeron,2=reggeon)    *
*         cl(2):      Structure function FLc (1=pomeron,2=reggeon)    *
***********************************************************************

      subroutine h12006pdf(z,q2,ifit,ipdf,xpq,f2,fl,c2,cl)

      implicit none

      double precision z,q2,xpq(-6:6),f2(2),fl(2),c2(2),cl(2)
      integer ifit,ipdf,ifit2,maxpdf,i
      character*1000 alposDir,pathA, pathB

      logical ierr
      data ierr/.false./
      save ierr

      data ifit2/99/
      save ifit2

      integer iwarn,maxwarn
      data iwarn /0/
      save iwarn
      parameter ( maxwarn = 10 )

      double precision zmin,zmax,q2min,q2max

      parameter (zmin=0.0043d0 , zmax=0.8d0)
      parameter (q2min=8.5d0 , q2max=1600.0d0)

cKC
cKC      WRITE(*,*)'KAREL ifit,ipdf = ',ifit,ipdf
cKC

      call getenv('ALPOS_DIR', alposDir)
      pathA = trim(alposDir) // trim('/src/aem/h1pdf2006Err/a.data')
      pathB = trim(alposDir) // trim('/src/aem/h1pdf2006Err/b.data')
      

      if (ifit.eq.1) then
         maxpdf=32
         if (ifit.ne.ifit2) then
            WRITE(6,*) '[H12006PDF] Initializing Fit A'
            call i_2006_fita
            ierr=.false.
      open(unit=1,file=pathA
     &           ,status='OLD',
     &           err=444)
            ierr=.true.
            call h1reada
            close(1)
         endif   
      elseif (ifit.eq.2) then
         maxpdf=30
         if (ifit.ne.ifit2) then
            WRITE(6,*) '[H12006PDF] Initializing Fit B'
            call i_2006_fitb
            ierr=.false.
c      open(unit=1,file='/afs/desy.de/user/p/pokorny/h1/'//
c     &'nlo/modules-4.1.0/ddis2jet/b.data'
c     &           ,status='OLD',
c     &           err=444)
c      open(unit=1,file='/h1wgs/h1mpim13/x02/usr/britzger/'//
c     &'alphasfit/JetsAtHighQ2/FitDiffJets/Boris/ddis2jet/b.data'
c     &           ,status='OLD',
c     &           err=444)
      open(unit=1,file=pathB
c      open(unit=1,file='/afs/desy.de/user/b/britzger/'//
c     &'xxl/alpos/svn-neu/Alpos/Alpos/aem/Boris/ddis2jet/b.data'
     &           ,status='OLD',
     &           err=444)
            ierr=.true.
            call h1readb
            close(1)
         endif   
      else
         WRITE(6,*)'[H12006PDF] Error: Unknown fit=',ifit
         stop
      endif

 555  continue

      if (ipdf.lt.0.or.ipdf.gt.maxpdf) then
         WRITE(6,*)'[H12006PDF] Error: pdf out of range [0,',
     &           maxpdf,'] :',ipdf
         stop
      endif

      ifit2=ifit


c     health warning
c DB: health warning commented out (discussed with federico)
c      if (z<zmin.or.z>zmax.or.q2.lt.q2min.or.q2.gt.q2max) then
c         if (iwarn.lt.maxwarn) then
c            WRITE(6,*) '[H12006PDF] Warning: z,q2 outside '//
c     &'range z=[0.0043,0.8], Q2=[8.5,1600]'
c            WRITE(6,*) '            z,q2=',z,q2            
c            iwarn=iwarn+1
c            if (iwarn.eq.maxwarn) then
c              WRITE(6,*)'[H12006PDF] Last warning.'
c            endif
c         endif
c      endif

      do i=-6,6
         xpq(i)=0
      enddo
      do i=1,2
         f2(i)=0
         fl(i)=0
         c2(i)=0
         cl(i)=0
      enddo

      if (ipdf.eq.0) then
        call qcd_2006(z,q2,0,xpq,f2,fl,c2,cl)
      else
        if (ierr) then
         if (ifit.eq.1) then
          call h1pdfa(ipdf)
          call h1evolvea(z,dsqrt(q2),xpq)
         else
          call h1pdfb(ipdf)
          call h1evolveb(z,dsqrt(q2),xpq)
         endif
        else
          WRITE(6,*)'[H12006PDF] Error: Error pdf requested but not av.'
          stop
        endif
      endif

      return

 444  WRITE(6,*)'[H12006PDF] Warning: File with error pdfs not found!'
      goto 555

      end
