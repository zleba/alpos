!> \file
!! apcalc program, subroutines.
!!
!! \author Volker Blobel, University Hamburg, 2014-2015
!!
!! \copyright
!! Copyright (c) 2014 - 2015 Volker Blobel, University Hamburg
!! This library is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Library General Public License as
!! published by the Free Software Foundation; either version 2 of the
!! License, or (at your option) any later version. \n\n
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Library General Public License for more details. \n\n
!! You should have received a copy of the GNU Library General Public
!! License along with this program (see the file COPYING.LIB for more
!! details); if not, write to the Free Software Foundation, Inc.,
!! 675 Mass Ave, Cambridge, MA 02139, USA.
!!

*  see e.g.     http://www.gnu.org/licenses/gpl-3.0.en.html


*_______________________________________________________________________

*                        Constrained least squares "apc"
*                        ===============================
*_______________________________________________________________________

      SUBROUTINE APCALC(INI,NX,X,VX, NF,F,A, DX,STATUS,IP,    W)
      IMPLICIT NONE
*     apply constraints
      INTEGER INI,NX,NF,NDF,IP,MAXIT,MAXITE
      DOUBLE PRECISION X(*),VX(*),DX(NX), F(NF),A(NX,NF),STATUS(6)
      DOUBLE PRECISION W(*)                ! work array
      DOUBLE PRECISION PULL(NX),EPSITE,EPSIT,RATIO,REDUCE
*     Printing routine (APCXVX):
      DOUBLE PRECISION ERRF,ERRU,PULLI !,DEI,DEJ,RHO 
*     internal
      INTEGER I,J,IJ,II,N,NRANK,NXSYM
      INTEGER NSYM,IB,IBF,IDA,INE,IDX
      DATA MAXIT/32/
      DATA EPSIT/1.0D-8/ 
*_______________________________________________________________________
      N    =NX+NF                      ! matrix dimension 
      NSYM =(N*N+N)/2                  ! number of matrix elements
      NXSYM=(NX*NX+NX)/2               ! number of X matrix elements
      NDF  =NF                         ! number of degrees of freedom

*     0                                ! pointer to saved X(.) initial NX   
*     NX                               ! pointer to matrix W       NSYM
      IB =NSYM+NX                      ! pointer to B(.)    length N
      IBF=NSYM+NX+N                    ! pointer to BF(.)   length NF
      IDX=NSYM+N+N                     ! pointer to DX(.)   length NX       
*     NSYM+N+N+NX                      ! words used so far        

      IDA=NSYM+N+N+NX                  ! pointer to DIAG(.)
      INE=NSYM+N+N+N+NX                ! pointer to NEXT(.)
*     NDI=NSYM+N+N+N+N+NX              ! work array dimension
*     NAL=NSYM+N+N+N+N+NX  +5*NF       ! numdif array 
*        =NSYM+5*N+4*NF 

*     copy covariance matrix to work array W(.)

      DO I=1,NXSYM 
       W(NX+I)=-VX(I)                   ! copy -VX(.) into W_11
      END DO
      DO I=1,NX                        ! analyse W_11
       II=(I*I+I)/2
       IF(VX(II).EQ.0.0D0) THEN
          NDF=NDF-1 ! number of degrees of freedom
       ELSE
c          IF(INI.LT.0) VX(II)=2.0D0*VX(II)
       END IF
       W(IB+I)=0.0D0                   ! reset
      END DO 
      INI=ABS(INI)

      IF(INI.EQ.0) THEN
         STATUS(1)=0.0D0 ! initial chisquare=0 
         STATUS(3)=0.0D0 ! initial F-average=0
         STATUS(6)=0.0D0 ! counter
         DO I=1,NX
          W(I) =X(I)     ! initial x-values 
          DX(I)=0.0D0    ! initial correction=0 !!!
         END DO
      END IF
      INI=INI+1

*     build matrix W(.) and right-hand side B(.)=W(IB+.)

      DO J=1,NF                           ! remaining constraints, rhs 
       W(IBF+J)=-F(J)                     ! constraint equations 
       DO I=1,NX
        W(IBF+J)=W(IBF+J)+A(I,J)*DX(I)    ! subtract actual step
       END DO
       W(IB+NX+J)=W(IBF+J)                ! right hand side 
      END DO

*     copy matrix A(NX,NF) to part of symmetric matrix
*     (backward, because of potential storage overlap of A(.,.) and W) 
      IJ=NSYM 
      DO J=NF,1,-1
       DO I=1,J 
        W(NX+IJ)=0.0D0
        IJ=IJ-1
       END DO 
       DO I=NX,1,-1
        W(NX+IJ)=A(I,J)
        IJ=IJ-1 
       END DO
      END DO
      CALL DXMINV(W(NX+1),W(IB+1),NX,NF, NRANK,  W(IDA+1),W(INE+1))
      NDF=NDF-(NX+NF-NRANK)

      DO I=1,NX
       DX(I)=W(IB+I) ! solution
      END DO
      STATUS(2)=STATUS(1) ! previous chi square
      STATUS(4)=STATUS(3) ! previous  F average       
      STATUS(5)=NDF       ! number of degrees of freedom        

      STATUS(1)=0.0D0
      STATUS(3)=0.0D0
      DO J=1,NF
       STATUS(1)=STATUS(1)-W(IB+NX+J)*W(IBF+J)  ! chi^2
       STATUS(3)=STATUS(3)+ABS(F(J))            ! sum |F|
      END DO        
      IF(STATUS(1).LT.0.0D0) STATUS(1)=0.0D0
      STATUS(3)=STATUS(3)/DFLOAT(NF)            ! sum |F|/NF
      IF(IP.EQ.0) RETURN

      IF(INI.EQ.1) THEN
         WRITE(*,*) ' '
         WRITE(*,100) NX,NF
         WRITE(*,101)
         WRITE(*,102) INI,NDF,STATUS(1),STATUS(3)
      ELSE  
         RATIO=LOG10((STATUS(3)+1.0D-8)/(STATUS(4)+1.0D-8))
         WRITE(*,103) INI,NDF,STATUS(1),STATUS(1)-STATUS(2),
     +                        STATUS(3),STATUS(3)-STATUS(4),RATIO 
         IF(INI.GE.MAXIT) THEN
            WRITE(*,*) 'Maximum number of iterations reached!'
         END IF
      END IF
 100  FORMAT('__________________________________________________'/
     +       'Constrained fit by apcalc with NX, NF=',2I6)
 101  FORMAT('_____________________________________',
     +       '_____________________________________'/  
     + ' iter  ndf       chi^2_____delta          ',
     +                       '|F|/NF_________delta       ratio')
 102  FORMAT(I4,I6,2X,F10.3,     16X,G14.4)
 103  FORMAT(I4,I6,2X,F10.3,G14.4,2X,G14.4,G14.4,F8.2)
      RETURN
*     __________________________________________________________________

      ENTRY APCLIN(INI,NX,X,DX,STATUS, W)  ! linear iteration
      
      IF(INI.GT.1) THEN
         IF(ABS(STATUS(1)-STATUS(2)).LT.EPSIT) THEN    ! d chi square
            IF(STATUS(3).LT.EPSIT) THEN                ! |F|
               INI=0                                   !=> convergence
            END IF
         END IF   
      END IF   
      IF(INI.GE.MAXIT) INI=0    ! maximum number of iterations reached
      IF(INI.EQ.0)     RETURN
      DO I=1,NX
       X(I)=W(I)+DX(I)     ! add correction dx to saved x
      END DO
      RETURN

*     __________________________________________________________________
      ENTRY APCNON(INI,NX,X,DX,STATUS, W)  ! non-linear iteration
      REDUCE=1.0D0 
      IF(INI.EQ.1) THEN
         REDUCE=0.25D0
      ELSE IF(INI.GT.1) THEN
         RATIO=(STATUS(3)+EPSIT)/(STATUS(4)+EPSIT)         
         IF(.NOT.(RATIO.GE.0.0D0.OR.RATIO.LE.0.0D0)) THEN
            RATIO=1.0D10
         END IF  
         IF(RATIO.GT.1.01D0) THEN
            REDUCE=0.5D0
            IF(RATIO.GT.1.0D1) REDUCE=0.1D0
            IF(RATIO.GT.1.0D3) REDUCE=0.01D0 
            IF(RATIO.GT.1.0D9) REDUCE=0.001D0 
            DO I=1,NX  
             DX(I)=X(I)-W(I)   ! (previous) step 
            END DO
         END IF
      END IF
      IF(INI.GE.3) THEN ! test for convergence
         IF(ABS(STATUS(1)-STATUS(2)).LT.EPSIT) THEN    ! chi square
            IF(ABS(STATUS(3)-STATUS(4)).LT.EPSIT) THEN ! |F|
               IF(STATUS(3).LT.EPSIT) THEN
                  INI=0
               END IF
            END IF
         END IF   
      END IF
      IF(INI.GE.MAXIT) INI=0 ! maximum number of iterations reached
      IF(INI.EQ.0)     RETURN  
      DO I=1,NX
       DX(I)=REDUCE*DX(I)  ! reduce step
       X(I)=W(I)+DX(I)     ! add correction dx to saved x
      END DO
c      IF(RATIO.EQ.1.0D9) INI=-INI  ! enlarge diagonal for solution
      RETURN

*     __________________________________________________________________

      ENTRY APCPAR(MAXITE,EPSITE) ! define MAXIT and EPSIT
      MAXIT=MAXITE
      IF(MAXIT.LE.0)     MAXIT=10
      EPSIT=EPSITE
      IF(EPSIT.LE.0.0D0) EPSIT=1.0D-8 
      RETURN 

*     __________________________________________________________________
      ENTRY APCPUL(NX,X,VX,PULL, W) ! get pull
      DO I=1,NX
       II=(I*I+I)/2
       PULL(I)=0.0D0
       IF(VX(II).GT.W(NX+II)) PULL(I)=(X(I)-W(I))/SQRT(VX(II)-W(NX+II))
      END DO
      RETURN

*     __________________________________________________________________
      ENTRY APCOVA(NX,VX, W)          ! get covariance matrix
      DO I=1,(NX*NX+NX)/2
       VX(I)=W(NX+I)
      END DO
      RETURN

*     __________________________________________________________________
      ENTRY APCRES(NX,X,VX, W) ! print start X and fitted X  with errors 
*     ...
      WRITE(*,*) ' '
      WRITE(*,201)
      DO I=1,NX
       II=(I*I+I)/2
       ERRF=SQRT(ABS(W(NX+II)))
       ERRU=SQRT(ABS(VX(II)))
       PULLI=0.0D0 
       IF(VX(II).GT.W(NX+II)) PULLI=(X(I)-W(I))/SQRT(VX(II)-W(NX+II))
       IF(PULLI.NE.0.0D0) THEN
          WRITE(*,202) I,X(I),ERRF,W(I),ERRU,PULLI
       ELSE IF(ERRU.NE.0.0D0) THEN
          WRITE(*,202) I,X(I),ERRF,W(I),ERRU
       ELSE
          WRITE(*,202) I,X(I),ERRF,W(I)
       END IF
      END DO     
 201  FORMAT(' par_i  fitted X_i_____+-std.dev.    ',
     +           '  initial X_i_____+-std.dev      pull')
 202  FORMAT(I6,2X,2G14.6,4X,2G14.6,F6.2)          
      RETURN

*     __________________________________________________________________
      ENTRY APCERO(NX,X,VX,NF,F)  ! zero 
      N=NX+NF
      DO I=1,NX
       X(I)=0.0D0
      END DO
      DO I=1,(NX*NX+NX)/2
       VX(I)=0.0D0
      END DO 
      DO J=1,NF
       F(J)=0.0D0
      END DO
      END


*_______________________________________________________________________
*     internal matrix routine for solution and inverse calculation     

      SUBROUTINE DCMINV(W,B,NARG, NRANK, AUX,NEXT) ! matrix inversion
*
*     Obtain solution of a system of
*                   linear equations  V *  X  =  B 
*     with symmetric matrix V, including matrix inversion
*
*                    - - -
*        CALL DCMINV(W,B,N,NRANK,AUX,NEXT)
*                    - -   -----
*
*           W = symmetric N-by-N matrix in symmetric storage mode
*               W(1) = W11, W(2) = W12, W(3) = W22, W(4) = W13, . . .
*               replaced by inverse matrix
*           B = N-vector, replaced by solution vector
*
*           NRANK is detected rank of matrix
*           AUX() is an auxiliary array of length N
*           NEXT() is an auxiliary array of length N. It is not integer, 
*           but double precision to make use of storage.
*
*     Method of solution is by elimination selecting the  pivot  on  the
*     diagonal each stage. The rank of the matrix is returned in  NRANK.
*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
*     matrix V and the corresponding elements of  B  are  set  to  zero.
*     SMINV can be used for a dimension up to 100 (see INVCDR).
*
      IMPLICIT NONE
      INTEGER NARG,NRANK,NX,NF
      DOUBLE PRECISION W(*),B(*),AUX(*),NEXT(*)
      DOUBLE PRECISION VKK,SUM,EPS,VJK
      PARAMETER (EPS=1.0D-6)
      
      INTEGER I,J,K,L,N,JK,JL,LK,KK, JFIRST,JLAST, LAST,MEAS,ME 

      INTEGER M,MM        ! statement function for indices 
      M(I)=(I*I+I)/2      ! diagonal index in symmetric storage mode
      MM(I,J)=(I*I-I)/2+J ! index of off-diagonal element       

*     ...
*     _________________________________________________________________
*     normal inversion entry
 
      N=NARG                    ! dimension parameter 
      JFIRST=1                  ! start flag
      DO I=1,N
       NEXT(I)=I+1              ! set all "next" pointer
       AUX(I)=ABS(W(M(I)))      ! save abs of diagonal elements
      END DO
      NEXT(N)=-1                ! end flag
      NRANK=0 
*     __________________________________________________________________      
*     common code for inversion and solution of matrix equation  

 10   DO I=1,N                  ! loop on all remaining elements
       VKK=0.0D0                ! search for pivot element
       K=0                      ! pivot index
       J=JFIRST                 ! first candidate index
       LAST=0

 20    IF(J.GT.0) THEN          ! test for linearity and zero matrix 
          IF(ABS(W(M(J))).GT.MAX(ABS(VKK),EPS*AUX(J))) THEN
             VKK=W(M(J))        ! largest pivot candidate so far
             K=J                ! index of largest
             L=LAST
          END IF
          LAST=J
          J=NEXT(J)            ! index of next candidate
          GOTO 20
       END IF

       IF(K.NE.0) THEN          ! pivot element found - proceed
          NRANK=NRANK+1         ! increase rank counter
          KK=M(K)               ! diagonal index
          IF(L.EQ.0) THEN
             JFIRST=NEXT(K)     ! new first index
          ELSE
             NEXT(L)=NEXT(K)    ! bridge used index
          END IF
          NEXT(K)=0.0D0         ! reset used index
c          NRANK=NRANK+1         ! increase rank ? again

          VKK  =1.0D0/VKK       ! invert pivot 
          W(KK)=-VKK          
          B(K) =B(K)*VKK
          JK   =KK-K
          JL   =0
          DO J=1,N              ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =W(JK)
              W(JK)=VKK*VJK
              B(J) =B(J)-B(K)*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  W(JL)=W(JL)-W(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE                     ! no pivot candidate found - reset
          DO K=1,N
           IF(NEXT(K).NE.0.0D0) THEN   ! undefined variable
              B(K)=0.0D0 ! clear undefined vector element
              DO J=1,K
               IF(NEXT(J).NE.0.0D0) W(MM(K,J))=0.0D0  ! clear
              END DO
           END IF
          END DO
          GOTO 30               ! end of inversion loop
       END IF
      END DO                    ! end of inversion loop

 30   DO I=1,M(N)
       W(I)=-W(I)               ! finally reverse sign
      END DO
      RETURN 

*     _________________________________________________________________
*     special entry for partial inversion of block matrix equation

      ENTRY DXMINV(W,B,NX,NF, NRANK, AUX,NEXT)
c      WRITE(*,*) 'DXMINV ',NX,NF
      N=NX+NF                   ! dimension parameter 
      DO I=1,N
       NEXT(I)=0.0D0
      END DO
      NRANK=0

*     distinguish between measured and unmeasured variables_____________

      JFIRST=0                  ! first index of measured variable 
      MEAS=0                    ! number of measured variables
      DO I=1,NX
       IF(W(M(I)).LT.0.0D0) THEN! measured variable
          NRANK=NRANK+1 
          IF(JFIRST.EQ.0) THEN
             JFIRST=I           ! first index of measured variable
          ELSE
             NEXT(JLAST)=I     ! insert index at previous index
          END IF
          JLAST=I               ! save index
          MEAS=MEAS+1
       END IF
      END DO ! I=1,NX
c      WRITE(*,*) 'MEAS NRANK ',MEAS,NRANK
c      WRITE(*,*) 'DXMINV JLAST',JLAST 
      IF(JLAST.EQ.0) GOTO 10    ! nothing to do
      NEXT(JLAST)=-1            ! stop index for last measured variable 

*     apply exchange algorithm to sub-matrices__________________________

c      WRITE(*,*) 'Starting loop'
      DO I=NX+1,N               ! loop I over constraint equations
       J=JFIRST                 ! first index of unmeasured variable
       DO ME=1,MEAS             ! already inverted element index J
        SUM=0.0D0
        JK=M(J-1)               ! index of diagonal element before 
        DO K=1,NX
         IF(K.LE.J) JK=JK+1     ! index in j column
         IF(NEXT(K).NE.0.0D0) SUM=SUM+W(JK)*W(MM(I,K)) 
         IF(K.GE.J) JK=JK+K     ! index in j row
        END DO ! K=1,NX
        AUX(J)=SUM              ! = A-row * VX-row/col
        J=NEXT(J)               ! next index of unmeasured variable
       END DO ! ME=1,MEAS

       DO K=I,N 
        SUM=0.0D0
        J=JFIRST                ! first index of unmeasured variable 
        DO ME=1,MEAS            ! already inverted element index J
         SUM=SUM+W(MM(K,J))*AUX(J) ! = A-row * H
         J=NEXT(J)             ! next index of unmeasured variable
        END DO ! ME=1,MEAS
        W(MM(K,I))=W(MM(K,I))+SUM  ! add to diagonal W_22
       END DO ! K=I,N

       J=JFIRST                 ! first index of unmeasured variable
       DO ME=1,MEAS
        W(MM(I,J))=-AUX(J)      ! add to off-diagonal W_22
        J=NEXT(J)              ! next index of unmeasured variable
       END DO ! ME=1,MEAS
      END DO ! I=NX+1,N
c      WRITE(*,*) 'Ending loop'

c      WRITE(*,*) 'inside DXMINV'
c      WRITE(*,777) (W(I),I=1,6)
c      WRITE(*,777) (B(I),I=1,3)
c 777  FORMAT(5F15.5)       


*     set pointer for unmeasured variables______________________________ 

      JFIRST=0
      JLAST=0
c      WRITE(*,*) 'before pointers N=',N
c      WRITE(*,*) 'JFIRST JLAST',JFIRST,JLAST
c      WRITE(*,777) (NEXT(I),I=1,8) 
      DO I=1,N    
c       WRITE(*,*) 'NEXT I=',I,NEXT(I)               
       IF(NEXT(I).EQ.0.0D0) THEN ! unmeasured variable
c          WRITE(*,*) 'is zero'
          IF(JFIRST.EQ.0) THEN  
             JFIRST=I           ! first index of unmeasured variable
          ELSE
             NEXT(JLAST)=I     ! next index of unmeasured variable
          END IF
          JLAST=I
       ELSE
          NEXT(I)=0.0D0        ! reset index for measured variable
       END IF 
      END DO ! I=1,N
      IF(JLAST.NE.0) NEXT(JLAST)=-1 ! end-flag
c      WRITE(*,*) 'JFIRST JLAST',JFIRST,JLAST
c      WRITE(*,777) (NEXT(I),I=1,8) 

      GOTO 10
      END 

      SUBROUTINE APCDER(NX,X,NF,F,AJAC,ICON,W)
*     __________________________________________________________________
*     Numerical derivative calculation, in a loop
*                       X(NX) = vector X 
*                       F(NF) = function (constraint) value F_j(X)
*               AUX(-2:+2,NF) = work array         
*             => AJAC(NX,NF)) = Jacobian matrix
*       
*           Code:
*          [ICON=2 for 5-point method] 
*     10    F( 1)=function of X(.)
*           ...
*           F(NF)=function of X(.)
*           CALL APC(NX,X,NF,F,AJAC,ICON,AUX)   
*           IF(ICON.NE.0) GOTO 10
*     __________________________________________________________________
      IMPLICIT NONE
      DOUBLE PRECISION X(NX),F(NF),AJAC(NX,NF),W(*) !,AUX(-2:+2,NF)
      INTEGER NX,NF,ICON,IREP,MODE,MODEST,J,N,NAL

      DOUBLE PRECISION EPSOPT       ! optimal step factor
      DOUBLE PRECISION STEP         ! step 
      DOUBLE PRECISION H            ! derivative step             
      DOUBLE PRECISION FAC(2)       ! factors for step
      DOUBLE PRECISION ARG          ! argument
      DOUBLE PRECISION XUN(-2:+2)   ! X values for function values
      DOUBLE PRECISION DER(4)       ! central derivatives
      DOUBLE PRECISION REL          ! estimate of rel. error
      PARAMETER (REL=7.0D-17)       ! realistic estimate of rel. error
      DOUBLE PRECISION EPSLON       ! copy (corrected) estimate
      DATA EPSLON/REL/              ! default epslon
      DATA IREP/0/                  ! neutral position
*     __________________________________________________________________
      SAVE
*     ...
      N=NX+NF
      NAL=(N*N+N)/2+4*N+NX-2        ! pointer to -2:+2 array
      IF(IREP.EQ.0) THEN            ! initial entry
         IF(ICON.NE.2) THEN         ! 3-point mode
            MODEST=1                ! [1.0000]
            EPSOPT=(3.0D0*EPSLON)**(1.0D0/3.0D0)
         ELSE IF(ICON.EQ.2) THEN    ! 5-point mode
            MODEST=2                ! [0.7071,1.0000]
            EPSOPT=(13.5D0*EPSLON)**(1.0D0/5.0D0)
            FAC(1)=SQRT(0.5D0)
         END IF
         FAC(MODEST)=1.0D0
         MODE=0                     ! reset counter
         DO J=1,NF   
 !         AUX(0,J)=F(J)            ! save central function values
           W(NAL+5*J)=F(J)          ! save central function values
         END DO
      END IF

      IF(IREP.GT.0) THEN            ! ... after -h step
         DO J=1,NF   
 !         AUX(-MODE,J)=F(J)         ! save function values
           W(NAL+5*J-MODE)=F(J)      ! save function values
         END DO
         X(IREP)=XUN(0)             ! restore variable value
         IF(MODE.EQ.MODEST) THEN    ! all function values are calculated:
            IF(MODEST.EQ.1) THEN    ! 3-point mode
               H=0.5D0*(XUN(1)-XUN(-1))
               DO J=1,NF
 !               AJAC(IREP,J)=(AUX(1,J)-AUX(-1,J))/(XUN(1)-XUN(-1))
                 AJAC(IREP,J)=
     +           (W(NAL+5*J+1)-W(NAL+5*J-1))/(XUN(1)-XUN(-1))
               END DO
            ELSE                    ! 5-point mode: Chebyshev
               DO J=1,NF
 !               CALL PRODER(AUX(-2,J),DER,0.5D0*(XUN(2)-XUN(-2)))
                 CALL PRODER(W(NAL+5*J-2),DER,0.5D0*(XUN(2)-XUN(-2)))
                AJAC(IREP,J)=DER(1) ! first derivative
               END DO
            END IF   
            IF(IREP.EQ.NX) THEN     ! all elements done
               IREP=0               ! reset flag
               ICON=0               ! reset argument
               DO J=1,NF
!                F(J)=AUX(0,J)      ! restore function values
                F(J)=W(NAL+5*J)     ! restore function values
               END DO
               RETURN               ! leave calder
            END IF   
            MODE=0                  ! reset counter
         END IF
      END IF
      IF(MODE.EQ.0) THEN            ! next variable
         IREP=IREP+1                ! next index
         XUN(0)=X(IREP)             ! save X
         STEP=EPSOPT*(1.0D0+ABS(XUN(0))) ! nominal step
      END IF

      IF(IREP.GE.0) THEN            ! initial or +- h done
         MODE=MODE+1 
         XUN(MODE)=XUN(0)+STEP*FAC(MODE)
         X(IREP)=XUN(MODE)          ! x_0 + h
         IREP=-IREP                 ! reverse switch -> negative
      ELSE IF(IREP.LT.0) THEN       ! ... after +h step
         DO J=1,NF   
!          AUX(MODE,J)=F(J)         ! save function values
           W(NAL+5*J+MODE)=F(J)     ! save function values
         END DO
         IREP=-IREP                 ! reverse switch -> positive
         XUN(-MODE)=XUN(0)-STEP*FAC(MODE)
         X(IREP)=XUN(-MODE)         ! x_0 - h
      END IF
      ICON=IREP                     ! continue: flag not zero
      RETURN

      ENTRY APCEPS(ARG)             ! multiply REL by factor
      EPSLON=ARG*REL                ! epsilon
      END

      SUBROUTINE PRODER(F,D,H)
*     __________________________________________________________________
*     calculate 1.st to 4.th derivative from 5 Chebyshev nodes
*        F(.) = function values at -1 -sqrt(1/2) 0 +sqrt(1/2) +1 
*        H    = step size
*     => D(4) = derivatives
*     __________________________________________________________________
      DOUBLE PRECISION F(-2:+2),D(4),H
*     ... 
      D(1)=0.5D0*(F(-2)-SQRT(8.0D0)*(F(-1)-F(+1))-F(+2))/H
      D(2)=(-F(-2)+4.0D0*F(-1)-6.0D0*F(0)+4.0D0*F(+1)-F(+2))/H**2
      D(3)=6.0D0*(-F(-2)+SQRT(2.0D0)*(F(-1)-F(+1))+F(+2))/H**3
      D(4)=24.0D0*(F(-2)-2.0D0*(F(-1)-F(0)+F(+1))+F(+2))/H**4
      END 

      SUBROUTINE APC(NX,X,VX,NF,F,STATUS, ISP,IST,W)
*     __________________________________________________________________
*     derivative calculation, step computation and convergence strategy
*     ISP =  0  quiet,    3-point derivatives     "almost" linear
*           +1  printout
*           +2  non-linear
*           +4  5-point derivatives  
*     CALL APC(NX,X,NF,F, ISP,IST,WORK)
*     IF(IST.GT.0) GOTO 10          ! continue iteration
*     IF(IST.LT.0) GOTO 20          ! continue derivative calculation
*     ... ! IST=0                   ! done
*
*     or
*
*     CALL APC(NX,X,NF,F, ISP,IST,WORK) !  
*     IF(IST.NE.0) GOTO 10          ! continue
*     ...                           ! done  
*
*     IST=1 printout   "almost" linear  3-point derivatives
*     IST=3 printout   "non-linear"     3-point derivative
*     IST=5 printout   "linear"         5-point derivative 
*     IST=7 printout   "non-linear"     5-point derivative 
*     __________________________________________________________________ 
      INTEGER NX,NF,ISP,IST, INI,ICON,N,IDX,IAA  
      DOUBLE PRECISION X(NX),VX(*),F(NF),STATUS(6),W(*)
      DATA INI/0/
      DATA NCOUNT/0/
      SAVE 
*     ... 
      NCOUNT=NCOUNT+1
      IP=MOD(ISP,2)    ! print flag
*     IL=MOD(ISP/2,2)  ! IL=0  non-lin; IL=1 lin
*     ID=ISP/4         ! =1 5-point derivative
      N=NX+NF
      IDX=(N*N+N)/2+N+N 
      IAA=NX+(NX*NX+NX)/2
      IF(ISP/4.EQ.1) ICON=2 
      CALL APCDER(NX,X,NF,F,W(IAA+1),ICON,W) ! derivatives
      IF(ICON.NE.0) THEN                 ! to be continued 
         IST=-ABS(ICON)
      ELSE                               ! derivatives done
         CALL APCALC(INI,NX,X,VX, NF,F,W(IAA+1),W(IDX+1),
     +               STATUS,IP, W)
         IF(MOD(ISP/2,2).EQ.0) THEN      ! linear 
            CALL APCLIN(INI,NX,X,W(IDX+1),STATUS,  W)
         ELSE                            ! non-linear
            CALL APCNON(INI,NX,X,W(IDX+1),STATUS,  W)
         END IF
         IF(INI.EQ.0) THEN          ! end-of-fit
            IF(IP.NE.0) THEN
               NUMDIF=3
               IF(ISP/4.NE.0) NUMDIF=5
               WRITE(*,101) NCOUNT,ISP,NUMDIF
            END IF
            STATUS(6)=NCOUNT
            NCOUNT=0
         END IF 
         IST=INI 
      END IF 
 101  FORMAT(I10,' F-evals (ISP=',I1,', ',
     +I1,'-point num.diff. calculation)')
      END 

*     uilities - not specific for apc


      SUBROUTINE APCOPY(X,VX,IX,Y,VY,JY,NN)  ! copy symmetric matrix
      IMPLICIT NONE
*     __________________________________________________________________
      INTEGER IX,JY,NN,NX,NF
      DOUBLE PRECISION X(*),VX(*),Y(*),VY(*),A(NX,NF)
*     __________________________________________________________________
      INTEGER  I,J,K,L,M,IM,JM,II
      DOUBLE PRECISION ERRF,DEI,DEJ,RHO
      REAL CM(10)
*     ...
      IM=(IX*IX+IX)/2-1 ! start 'from' index
      JM=(JY*JY+JY)/2-1 ! start 'to'   index
      DO K=1,NN 
       Y(JY-1+K)=X(IX-1+K)  
       DO L=1,K
        IM=IM+1
        JM=JM+1
        VY(JM)=VX(IM)   ! copy
       END DO
       IM=IM+IX-1
       JM=JM+JY-1                             
      END DO
      RETURN

*     __________________________________________________________________
      ENTRY APRXVX(NX,X,VX) ! print X, VX with errors
*     ...
      WRITE(*,*) ' '
      WRITE(*,101)
      DO I=1,NX       
       IF(VX(1).LT.0.0D0) THEN  ! missing covariance matrix VX(.)
          WRITE(*,102) I,X(I)
       ELSE
          II=(I*I+I)/2
          ERRF=SQRT(ABS(VX(II)))
          IF(ERRF.NE.0.0D0) THEN
             WRITE(*,102) I,X(I),ERRF
          ELSE
             WRITE(*,102) I,X(I)
          END IF
       END IF 
      END DO
 101  FORMAT(' par_i   __X_i_______+-__std.dev____')
 102  FORMAT(I6,2X,2G14.6)
      RETURN

*     __________________________________________________________________
      ENTRY APCORR(NX,VX)   ! print correlations
*     ...
      WRITE(*,*)
      WRITE(*,*) 'Correlation coefficients between parameters i and j:'
      WRITE(*,203) (I,I=1,MIN(10,NX))
      DO I=1,NX
       DEI=SQRT(VX((I*I+I)/2))
       IF(DEI.GT.0.0D0) THEN 
          M=0
          DO J=1,I
           RHO=0.0D0
           DEJ=SQRT(VX((J*J+J)/2))
           IF(DEI*DEJ.GT.0.0D0) THEN
              RHO=VX(J+(I*I-I)/2)/(DEI*DEJ)
           END IF
           M=M+1
           IF(M.GT.10) WRITE(*,*) 'Heavy warning'
           CM(M)=RHO
           IF(M.LT.0.OR.M.GT.10) WRITE(*,*) 'ALarm ',M
           IF(M.EQ.10) THEN
              IF(J.LE.10) THEN
                 WRITE(*,204) I,CM
              ELSE
                 WRITE(*,205)   CM
              END IF
              M=0
           END IF
          END DO
          IF(M.NE.0) THEN
             IF(J.LE.10) THEN
                WRITE(*,204) I,(CM(L),L=1,M)
             ELSE
                WRITE(*,205)   (CM(L),L=1,M)
             END IF
             M=0
          END IF
C       ELSE   
C          WRITE(*,204) I
       END IF 
      END DO  
c 201  FORMAT(' par_i  fitted X_i_____+-std.dev.    ',
c     +           '  initial X_i_____+-std.dev      pull')
c 202  FORMAT(I6,2X,2G14.6,4X,2G14.6,F6.2B)
 203  FORMAT(' par_i  \ j=',10I6)
 204  FORMAT(I6,6X,10F6.3)
 205  FORMAT(12X,   10F6.3)
      RETURN  

*     __________________________________________________________________
      ENTRY APDER(NX,NF,A)  ! print derivative matrix
      WRITE(*,*) ' '
      WRITE(*,301) 'Derivative matrix of dimension',NX,' x',NF
      DO J=1,NF
       WRITE(*,302) J,(A(I,J),I=1,NX)
      END DO
 301  FORMAT(A,I4,A,I4)
 302  FORMAT(I4,2X,5G14.6,1X/(6X,5G14.6))
      END

