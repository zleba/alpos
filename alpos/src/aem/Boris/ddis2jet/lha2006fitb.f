      subroutine H1evolveb(xin,qin,pdf)
      implicit real*8 (a-h,o-z)
c  done on  28/06/06 at  13.30.12
C PARAMETERs moved along eigenvector  7 by    -0.65019 dchi2=   2.73            
c evolution has been made starting at q2_input =      2.500
c  available for :       1.960 <= q2 <=   30000.000
c            and :  0.001000 <= x <=  0.932603
c*                                                 
c   for x outside limits, the closest limit        
c        is assumed : f2(x>xmax,q2)=f2(xmax,q2)    
c                     f2(x<xmin,q2)=f2(xmin,q2)    
c  for q2 outside limits, the closest limit        
c        is assumed : f2(x,q2>q2max)=f2(x,q2max)   
c                     f2(x,q2<q2min)=f2(x,q2min)   
c*                                                 
c  comments, etc... to C. Pascaud or F. Zomer      
***************************************************
       PARAMETER(n_bin_q2= 76)
       PARAMETER(n_bin_x=100)
      REAL xl_bin(n_bin_x),q2l_bin(n_bin_q2) 
      PARAMETER(ngrid=40)
      REAL*4 f(0:ngrid,8,n_bin_x,n_bin_q2),val(8) 
      real*8 pdf(-6:6)
      real*4 q2in,x,y 
      save
c
      q2in = qin*qin
      x=log(xin)                                   
      y=log(q2in)                                  
      DO i=2,n_bin_x                               
        IF(x.LT.xl_bin(i))  goto 1                 
        IF(xl_bin(i).ge.0.)  goto 1                
      ENDDO                                        
      i=n_bin_x                                    
    1 i=i-1                                        
      DO j=2,n_bin_q2                              
        IF(y.LT.q2l_bin(j))  GOTO 2                
      ENDDO                                        
      j=n_bin_q2                                   
    2 j=j-1                                        
      dx=xl_bin(i+1)-xl_bin(i)                     
      xd=(x-xl_bin(i))/dx                          
      dy=q2l_bin(j+1)-q2l_bin(j)                   
      yd=(y-q2l_bin(j))/dy                         
c                                                  
      do k=1,8                                     
      val(k)=f(imem,k,i,j)+xd*(f(imem,k,i+1,j)-f(imem,k,i,j))
     &+yd*(f(imem,k,i,j+1)-f(imem,k,i,j))          
     &+xd*yd*(f(imem,k,i+1,j+1)+f(imem,k,i,j)      
     &-f(imem,k,i+1,j)-f(imem,k,i,j+1))            
      enddo                                        
      pdf(-6) = 0.0d0                              
       pdf(6) = 0.0d0                              
      pdf(-5) = val(7)                             
       pdf(5) = val(7)                             
      pdf(-4) = val(6)                             
       pdf(4) = val(6)                             
      pdf(-3) = val(5)                             
       pdf(3) = val(5)                             
      pdf(-2) = val(4)                             
       pdf(2) = val(3)+val(4)                      
      pdf(-1) = val(2)                             
       pdf(1) = val(1)+val(2)                      
       pdf(0) = val(8)                             
      return                                       
c                                                  
      entry H1readb                                 
c                                                  
      read(1,*)nmem,ndef                           
      read(1,1000) (xl_bin(nx),nx=1,n_bin_x)       
      read(1,1000) (q2l_bin(nq2),nq2=1,n_bin_q2)   
        do i=1,n_bin_q2                            
          q2l_bin(i)=log(q2l_bin(i))               
        enddo                                      
      do nm = 0,nmem                               
        do jval = 1,8                              
      read(1,1000)((f(nm,jval,nx,nq2),nx=1,n_bin_x),nq2=1,n_bin_q2)   
        enddo                                      
      enddo                                        
      return                                       
c                                                  
c      entry H1alfa(alfas,qalfa)                    
c      CALL alphah1(alfas,Qalfa)                    
c      return                                       
c                                                  
      entry H1initb(Eorder,Q2fit)                   
      return                                       
c                                                  
      entry H1pdfb(mem)                             
      imem = mem                                   
      return                                       
c                                                  
 1000 format(5e13.5)                               
      END                                          
