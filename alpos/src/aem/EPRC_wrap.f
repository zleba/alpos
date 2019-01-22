      
c-----------------------------------------------------------------------------
c...set mw
      subroutine eprc_set_mw(mw_in) 
      implicit none
      double precision mw_in
      double precision SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      MW=mw_in
      MW2=MW*MW
      end

c...get mw
      subroutine eprc_get_mw(mw_out) 
      implicit none
      double precision mw_out
      double precision SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      mw_out=MW
      end
      
c...set mt
      subroutine eprc_set_mt(mt_in) 
      implicit none
      double precision mt_in
      double precision SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      MT=mt_in
      end

c...get mt
      subroutine eprc_get_mt(mt_out) 
      implicit none
      double precision mt_out
      double precision SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      mt_out=MT
      end
      
c...sw2 following Czarnecki-Marciano
      function eprc_sw2_cm(q2)
      implicit none
      double precision eprc_sw2_cm,q2
      double precision PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      double precision SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      integer LPAR,LPARIN,IPART
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART

      double precision ze,zmu,ztau,zu,zd,zc,zs,zt,zb,zw
      double precision pe,pmu,ptau,pu,pd,pc,ps,pt,pb,pw
      double precision kappaf,kappab,kappa
      double precision sw2cm,cw2cm,swcm,cwcm
      double precision sw2os,cw2os
      COMPLEX*16 SIGFWS,SIGFZS
      
      ze=me2/q2
      zmu=mmy2/q2
      ztau=mtau2/q2
      zu=mu2/q2
      zd=md2/q2
      zc=mc2/q2
      zs=ms2/q2
      zt=mt2/q2
      zb=mb2/q2
      zw=mw2/q2
      pe=dlog((dsqrt(1d0+4d0*ze)+1d0)**2/(4d0*ze))
      pmu=dlog((dsqrt(1d0+4d0*zmu)+1d0)**2/(4d0*zmu))
      ptau=dlog((dsqrt(1d0+4d0*ztau)+1d0)**2/(4d0*ztau))
      pu=dlog((dsqrt(1d0+4d0*zu)+1d0)**2/(4d0*zu))
      pd=dlog((dsqrt(1d0+4d0*zd)+1d0)**2/(4d0*zd))
      pc=dlog((dsqrt(1d0+4d0*zc)+1d0)**2/(4d0*zc))
      ps=dlog((dsqrt(1d0+4d0*zs)+1d0)**2/(4d0*zs))
      pt=dlog((dsqrt(1d0+4d0*zt)+1d0)**2/(4d0*zt))
      pb=dlog((dsqrt(1d0+4d0*zb)+1d0)**2/(4d0*zb))
      pw=dlog((dsqrt(1d0+4d0*zw)+1d0)**2/(4d0*zw))

c...MSbar from On-shell mixing angle
      cw2os=mw2/mz2
      sw2os=1d0-cw2os
      sw2cm=sw2os/(1d0
     .     -cw2os/sw2os*(dREAL(SIGFZS(MZ2)/MZ2-SIGFWS(MW2)/MW2)))      

      cw2cm=1d0-sw2cm
      swcm=dsqrt(sw2cm)
      cwcm=dsqrt(cw2cm)
      
      kappaf=-alp2pi/swcm/cwcm*
     .       1d0/3d0*(1d0/2d0-2d0*sw2cm)*(
     .       dlog(me2/mz2)-5d0/3d0+4d0*ze
     .       +(1d0-2d0*ze)*dsqrt(1d0+4d0*ze)*pe)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       1d0/3d0*(1d0/2d0-2d0*sw2cm)*(
     .       dlog(mmy2/mz2)-5d0/3d0+4d0*zmu
     .       +(1d0-2d0*zmu)*dsqrt(1d0+4d0*zmu)*pmu)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       1d0/3d0*(1d0/2d0-2d0*sw2cm)*(
     .       dlog(mtau2/mz2)-5d0/3d0+4d0*ztau
     .       +(1d0-2d0*ztau)*dsqrt(1d0+4d0*ztau)*ptau)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/3d0-8d0/9d0*sw2cm)*(
     .       dlog(mu2/mz2)-5d0/3d0+4d0*zu
     .       +(1d0-2d0*zu)*dsqrt(1d0+4d0*zu)*pu)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/3d0-8d0/9d0*sw2cm)*(
     .       dlog(mc2/mz2)-5d0/3d0+4d0*zc
     .       +(1d0-2d0*zc)*dsqrt(1d0+4d0*zc)*pc)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/3d0-8d0/9d0*sw2cm)*(
     .       dlog(mt2/mz2)-5d0/3d0+4d0*zt
     .       +(1d0-2d0*zt)*dsqrt(1d0+4d0*zt)*pt)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/6d0-2d0/9d0*sw2cm)*(
     .       dlog(md2/mz2)-5d0/3d0+4d0*zd
     .       +(1d0-2d0*zd)*dsqrt(1d0+4d0*zd)*pd)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/6d0-2d0/9d0*sw2cm)*(
     .       dlog(ms2/mz2)-5d0/3d0+4d0*zs
     .       +(1d0-2d0*zs)*dsqrt(1d0+4d0*zs)*ps)
      kappaf=kappaf
     .       -alp2pi/swcm/cwcm*
     .       (1d0/6d0-2d0/9d0*sw2)*(
     .       dlog(mb2/mz2)-5d0/3d0+4d0*zb
     .       +(1d0-2d0*zb)*dsqrt(1d0+4d0*zb)*pb)
      
      if (q2.gt.0.0006d0) then 
        kappab=-alp2pi/swcm/cwcm*(
     .         -(42d0*cw2cm+1d0)/12d0*dlog(cw2cm)+1d0/18d0
     .         -(dsqrt(1d0+4d0*zw)*pw/2d0-1d0)
     .          *((7d0-4d0*zw)*cw2cm+(1d0+4d0*zw)/6d0)
     .         -zw*(3d0/4d0-zw+(zw-3d0/2d0)*dsqrt(1d0+4d0*zw)*pw
     .              +zw*(2d0-zw)*pw*pw
     .             )
     .        )
      else
        kappab=-alp2pi/swcm/cwcm*(
     .         -(42d0*cw2cm+1d0)/12d0*dlog(cw2cm)
     .         +(1d0/3d0-37d0/60d0/zw)*cw2cm
     .         +4d0/9d0-43d0/720d0/zw
     .        )
      endif
      
      kappa=1d0+cwcm/swcm*(kappaf+kappab)
      eprc_sw2_cm=sw2cm*kappa
      
      return
      end

c---------------------------------------------------------------------------
