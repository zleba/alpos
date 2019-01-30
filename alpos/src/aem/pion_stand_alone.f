c======================================================================
c
c subroutine STROWP1(X,SCALE,UPV,DNV,SEA,STR,CHM,GL) 
c taken from PDFLIB provides pion structure function acc.
c to Owens Set 1 ... (type=2, group=1, set=1)
c
c=====================================================================


      SUBROUTINE STROWP1(X,SCALE,UPV,DNV,SEA,STR,CHM,GL)
C :::::::::  OWENS SET 1 PION STRUCTURE FUNCTION  :::::::::

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      integer k
      DOUBLE PRECISION DGAMMA,X
      DOUBLE PRECISION COW(3,5,4),TS(6),XQ(9)
      double precision gg1,gg2,gg3,z(0:270),gg(0:270)
      
C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1),IS=1,5),IP=1,3)/
     1  4.0000D-01,  7.0000D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -6.2120D-02,  6.4780D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     3 -7.1090D-03,  1.3350D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
     1  8.8800D-01,  0.0000D+00,  3.1100D+00,  6.0000D+00,  0.0000D+00,
     2 -1.8020D+00, -1.5760D+00, -1.3170D-01,  2.8010D+00, -1.7280D+01,
     3  1.8120D+00,  1.2000D+00,  5.0680D-01, -1.2160D+01,  2.0490D+01/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
     1  9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -2.4280D-01, -2.1200D-01,  8.6730D-01,  1.2660D+00,  2.3820D+00,
     3  1.3860D-01,  3.6710D-03,  4.7470D-02, -2.2150D+00,  3.4820D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
     1  0.0000D+00, -2.2120D-02,  2.8940D+00,  0.0000D+00,  0.0000D+00,
     2  7.9280D-02, -3.7850D-01,  9.4330D+00,  5.2480D+00,  8.3880D+00,
     3 -6.1340D-02, -1.0880D-01, -1.0852D+01, -7.1870D+00, -1.1610D+01/

       DATA ZEROD/0.D0/, ONED/1.D0/, SIXD/6.D0/
       DATA ALAM/0.2D0/, Q02/4.D0/, QMAX2/2.D3/
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        Q2 = SCALE*SCALE
        Q2IN = MIN( QMAX2,MAX( Q02,Q2))
        SD = LOG( LOG( Q2IN/ALAM**2)/ LOG( Q02/ALAM**2))

C...Calculate structure functions.
        DO 240 KFL=1,4
        DO 230 IS=1,5
  230   TS(IS)=COW(1,IS,KFL)+COW(2,IS,KFL)*SD+
     &  COW(3,IS,KFL)*SD*SD
        IF(KFL.EQ.1) THEN

	do k=0,270,1
	z(k)=(3d0-0.3d0)*dble(k)/dble(270)+0.3d0
	enddo
	
	gg(   0)=   2.992
	gg(   1)=   2.890
	gg(   2)=   2.796
	gg(   3)=   2.707
	gg(   4)=   2.624
	gg(   5)=   2.546
	gg(   6)=   2.473
	gg(   7)=   2.404
	gg(   8)=   2.338
	gg(   9)=   2.277
	gg(  10)=   2.218
	gg(  11)=   2.163
	gg(  12)=   2.110
	gg(  13)=   2.061
	gg(  14)=   2.013
	gg(  15)=   1.968
	gg(  16)=   1.925
	gg(  17)=   1.884
	gg(  18)=   1.845
	gg(  19)=   1.808
	gg(  20)=   1.772
	gg(  21)=   1.738
	gg(  22)=   1.706
	gg(  23)=   1.675
	gg(  24)=   1.645
	gg(  25)=   1.616
	gg(  26)=   1.589
	gg(  27)=   1.562
	gg(  28)=   1.537
	gg(  29)=   1.513
	gg(  30)=   1.489
	gg(  31)=   1.467
	gg(  32)=   1.445
	gg(  33)=   1.424
	gg(  34)=   1.404
	gg(  35)=   1.385
	gg(  36)=   1.366
	gg(  37)=   1.348
	gg(  38)=   1.331
	gg(  39)=   1.314
	gg(  40)=   1.298
	gg(  41)=   1.282
	gg(  42)=   1.267
	gg(  43)=   1.253
	gg(  44)=   1.239
	gg(  45)=   1.225
	gg(  46)=   1.212
	gg(  47)=   1.200
	gg(  48)=   1.187
	gg(  49)=   1.176
	gg(  50)=   1.164
	gg(  51)=   1.153
	gg(  52)=   1.142
	gg(  53)=   1.132
	gg(  54)=   1.122
	gg(  55)=   1.112
	gg(  56)=   1.103
	gg(  57)=   1.094
	gg(  58)=   1.085
	gg(  59)=   1.077
	gg(  60)=   1.069
	gg(  61)=   1.061
	gg(  62)=   1.053
	gg(  63)=   1.046
	gg(  64)=   1.038
	gg(  65)=   1.031
	gg(  66)=   1.025
	gg(  67)=   1.018
	gg(  68)=   1.012
	gg(  69)=   1.006
	gg(  70)=   1.000
	gg(  71)=   0.994
	gg(  72)=   0.989
	gg(  73)=   0.984
	gg(  74)=   0.978
	gg(  75)=   0.974
	gg(  76)=   0.969
	gg(  77)=   0.964
	gg(  78)=   0.960
	gg(  79)=   0.955
	gg(  80)=   0.951
	gg(  81)=   0.947
	gg(  82)=   0.944
	gg(  83)=   0.940
	gg(  84)=   0.936
	gg(  85)=   0.933
	gg(  86)=   0.930
	gg(  87)=   0.927
	gg(  88)=   0.924
	gg(  89)=   0.921
	gg(  90)=   0.918
	gg(  91)=   0.916
	gg(  92)=   0.913
	gg(  93)=   0.911
	gg(  94)=   0.909
	gg(  95)=   0.906
	gg(  96)=   0.904
	gg(  97)=   0.903
	gg(  98)=   0.901
	gg(  99)=   0.899
	gg( 100)=   0.897
	gg( 101)=   0.896
	gg( 102)=   0.895
	gg( 103)=   0.893
	gg( 104)=   0.892
	gg( 105)=   0.891
	gg( 106)=   0.890
	gg( 107)=   0.889
	gg( 108)=   0.889
	gg( 109)=   0.888
	gg( 110)=   0.887
	gg( 111)=   0.887
	gg( 112)=   0.886
	gg( 113)=   0.886
	gg( 114)=   0.886
	gg( 115)=   0.886
	gg( 116)=   0.886
	gg( 117)=   0.886
	gg( 118)=   0.886
	gg( 119)=   0.886
	gg( 120)=   0.886
	gg( 121)=   0.887
	gg( 122)=   0.887
	gg( 123)=   0.888
	gg( 124)=   0.888
	gg( 125)=   0.889
	gg( 126)=   0.890
	gg( 127)=   0.890
	gg( 128)=   0.891
	gg( 129)=   0.892
	gg( 130)=   0.894
	gg( 131)=   0.895
	gg( 132)=   0.896
	gg( 133)=   0.897
	gg( 134)=   0.899
	gg( 135)=   0.900
	gg( 136)=   0.902
	gg( 137)=   0.903
	gg( 138)=   0.905
	gg( 139)=   0.907
	gg( 140)=   0.909
	gg( 141)=   0.911
	gg( 142)=   0.913
	gg( 143)=   0.915
	gg( 144)=   0.917
	gg( 145)=   0.919
	gg( 146)=   0.921
	gg( 147)=   0.924
	gg( 148)=   0.926
	gg( 149)=   0.929
	gg( 150)=   0.931
	gg( 151)=   0.934
	gg( 152)=   0.937
	gg( 153)=   0.940
	gg( 154)=   0.943
	gg( 155)=   0.946
	gg( 156)=   0.949
	gg( 157)=   0.952
	gg( 158)=   0.955
	gg( 159)=   0.958
	gg( 160)=   0.962
	gg( 161)=   0.965
	gg( 162)=   0.969
	gg( 163)=   0.972
	gg( 164)=   0.976
	gg( 165)=   0.980
	gg( 166)=   0.984
	gg( 167)=   0.988
	gg( 168)=   0.992
	gg( 169)=   0.996
	gg( 170)=   1.000
	gg( 171)=   1.004
	gg( 172)=   1.009
	gg( 173)=   1.013
	gg( 174)=   1.018
	gg( 175)=   1.022
	gg( 176)=   1.027
	gg( 177)=   1.032
	gg( 178)=   1.037
	gg( 179)=   1.041
	gg( 180)=   1.046
	gg( 181)=   1.052
	gg( 182)=   1.057
	gg( 183)=   1.062
	gg( 184)=   1.068
	gg( 185)=   1.073
	gg( 186)=   1.079
	gg( 187)=   1.084
	gg( 188)=   1.090
	gg( 189)=   1.096
	gg( 190)=   1.102
	gg( 191)=   1.108
	gg( 192)=   1.114
	gg( 193)=   1.120
	gg( 194)=   1.127
	gg( 195)=   1.133
	gg( 196)=   1.140
	gg( 197)=   1.146
	gg( 198)=   1.153
	gg( 199)=   1.160
	gg( 200)=   1.167
	gg( 201)=   1.174
	gg( 202)=   1.181
	gg( 203)=   1.188
	gg( 204)=   1.196
	gg( 205)=   1.203
	gg( 206)=   1.211
	gg( 207)=   1.218
	gg( 208)=   1.226
	gg( 209)=   1.234
	gg( 210)=   1.242
	gg( 211)=   1.250
	gg( 212)=   1.259
	gg( 213)=   1.267
	gg( 214)=   1.276
	gg( 215)=   1.284
	gg( 216)=   1.293
	gg( 217)=   1.302
	gg( 218)=   1.311
	gg( 219)=   1.320
	gg( 220)=   1.329
	gg( 221)=   1.339
	gg( 222)=   1.348
	gg( 223)=   1.358
	gg( 224)=   1.368
	gg( 225)=   1.378
	gg( 226)=   1.388
	gg( 227)=   1.398
	gg( 228)=   1.408
	gg( 229)=   1.419
	gg( 230)=   1.430
	gg( 231)=   1.440
	gg( 232)=   1.451
	gg( 233)=   1.463
	gg( 234)=   1.474
	gg( 235)=   1.485
	gg( 236)=   1.497
	gg( 237)=   1.509
	gg( 238)=   1.520
	gg( 239)=   1.532
	gg( 240)=   1.545
	gg( 241)=   1.557
	gg( 242)=   1.570
	gg( 243)=   1.582
	gg( 244)=   1.595
	gg( 245)=   1.608
	gg( 246)=   1.622
	gg( 247)=   1.635
	gg( 248)=   1.649
	gg( 249)=   1.662
	gg( 250)=   1.676
	gg( 251)=   1.691
	gg( 252)=   1.705
	gg( 253)=   1.720
	gg( 254)=   1.734
	gg( 255)=   1.749
	gg( 256)=   1.765
	gg( 257)=   1.780
	gg( 258)=   1.796
	gg( 259)=   1.811
	gg( 260)=   1.827
	gg( 261)=   1.844
	gg( 262)=   1.860
	gg( 263)=   1.877
	gg( 264)=   1.894
	gg( 265)=   1.911
	gg( 266)=   1.928
	gg( 267)=   1.946
	gg( 268)=   1.964
	gg( 269)=   1.982
	gg( 270)=   2.000



	do k=0,270,1
	if((TS(1).ge.z(k)).and.(TS(1).lt.z(k+1)))then
	gg1=(gg(k)+gg(k+1))/2d0
	endif
	enddo	

	do k=0,270,1
        if(((TS(2)+1d0).ge.z(k)).and.((TS(2)+1d0).lt.z(k+1)))then
	gg2=(gg(k)+gg(k+1))/2d0
	endif	
	enddo	

	do k=0,270,1
c        if(((TS(1)+TS(2)+1d0).ge.z(k)).and.((TS(1)+TS(2)+1d0).lt.z(k+1)))then
        if(((TS(1)+TS(2)+1d0).ge.z(k)))then
           if((TS(1)+TS(2)+1d0).lt.z(k+1))then
              gg3=(gg(k)+gg(k+1))/2d0
           endif
	endif
	enddo
	
	DENOM=gg1*gg2/gg3

c          DENOM = DGAMMA(TS(1))*DGAMMA(TS(2)+ONED)/
c     +                                          DGAMMA(TS(1)+TS(2)+ONED)
     
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/DENOM
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
  240   CONTINUE

C...Put into output arrays.
        UPV = XQ(1)
        DNV = XQ(1)
        SEA = XQ(3)/SIXD
        STR = XQ(3)/SIXD
        CHM = XQ(4)
        BOT = ZEROD
        TOP = ZEROD
        GL  = XQ(2)
C
        RETURN
        END