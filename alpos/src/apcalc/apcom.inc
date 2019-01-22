*     Common for APCALC applications
*     __________________________________________________________________
      IMPLICIT NONE
*     apcalc arrays
      INTEGER NXMAX,NFMAX,NMAX,NSYMAX, NX,NF,INI,ICON,IST,ISP
*     ___________________________
      PARAMETER (NXMAX=5,NFMAX=5)  !           <===< adapt these numbers
      PARAMETER (NMAX=NXMAX+NFMAX,NSYMAX=(NMAX*NMAX+NMAX)/2)
      DOUBLE PRECISION X(NXMAX),VX((NXMAX*NXMAX+NXMAX)/2)
      DOUBLE PRECISION F(NFMAX)
      DOUBLE PRECISION WORK(NSYMAX+5*NMAX+5*NFMAX)          ! work array
      DOUBLE PRECISION STATUS(6)
*     __________________________________________________________________
