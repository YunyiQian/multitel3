      SUBROUTINE BODYW13(iFlt,X,N,DT,IB,IC,H,P,FC0,zqp,ZQS,
     -                  ZI,dely,ND)
*===================================================*
*  Calculate synthetic waveforms                    *
*    Type of body wave  IB= 1/2/3/4: P/SV/SH/PP     *
*    Component          IC= 1/2/3: UD/NS/EW (IB=1/4)*
*                       IC= 1/2  : UD/HR    (IB=2)  *
*                       IC= any  : SH       (IB=3)  *
*
* X   - Output Vector of Time series (OUT)
* N   - Number of Points (IN)
* DT  - Delta or Time between points (IN)
* IB  - Type of Body Wave (IN)
*                        1 = P 
*                        2 = SV
*                        3 = SH
*                        4 = PP
* IC  - Component         (IN)
*                        1 = Z   IF IB = 1 or 4
*                        2 = NS  IF IB = 1 or 4
*                        3 = EW  IF IB = 1 or 4
*                        1 = UD  IF IB = 2
*                        2 = HR  IF IB = 2
*                      any = SH  IF IB = 3      
* F1  - strike       (IN)
* D1  - dip          (IN)
* A1  - slip         (IN)
* H   - source depth (IN)
* AZ  - azimuth      (IN)  
* P   - (IN)
* FC0 - (IN)
* ZQP - P attenuation (IN)
* ZQS - S attenuation (IN)
* ZI  - Instrument response (IN)
* dely - Delay (OUT)
*     
*============================Ver.900715 ===========*
* Modification:                                     *
* 1) An isotropic component of M.T. is added        *
*         for d1 > 360 (degree)        -900531      *
* 2) Source layer # is determined in this subroutine*
*      independently from a reference point -900715 *
*===================================================*
      PARAMETER (NL0=20,RADIUS=6371,PI=3.141593)
      IMPLICIT COMPLEX*8 (Z)
      DIMENSION X(N),Z(ND),ZI(ND),ZQP(ND),ZQS(ND)
     -,  ZR0(ND),ZRPU(ND),ZRPD(ND),ZRSU(ND),ZRSD(ND),ZRPP(ND),ZDM(ND)
      COMMON /STR0/NL ,VP (NL0),VS (NL0),DEN (NL0),DEP (NL0)
      COMMON /STR1/NL1,VP1(NL0),VS1(NL0),DEN1(NL0),DEP1(NL0)
      COMMON /STR2/NL2,VP2(NL0),VS2(NL0),DEN2(NL0),DEP2(NL0)
	 common /sourceRegion/vsrc

         DF=1/(DT*N)
         DW=DF*2*PI
         TL=DT*N
* < Source layer # >
       hl=0.
       do 15 l=1,nl-1
	 hl=hl+dep(l)
	 dh=h-hl
15     if(dh.lt.0.) goto 16
16      ll=l
	hl=hl-dep(ll)
	vsrc=vp(ll)
* < Radiation pattern >
      CALL RADP3(iFlt,P,PD,PU,SVD,SVU,SHD,SHU,VP(LL),VS(LL))
c	write(*,*) 'rad',iFlt,pd,pu,svd,svu,shd,shu
* < Structure effects: Near-source & near-reciever >
         DH=H-HL
       IF(IB.EQ.3) GOTO 1
       CALL REFL(ZRPU,ZRPD,ZRSU,ZRSD,DW,N,P,VP,VS,DEN,DEP,NL,LL,IB,TR1)
       IF(IC.EQ.1)
     - CALL CNVR(ZDM,ZR0,DW,N,P,VP1,VS1,DEN1,DEP1,NL1,IB,TR2)
       IF(IC.NE.1)
     - CALL CNVR(ZR0,ZDM,DW,N,P,VP1,VS1,DEN1,DEP1,NL1,IB,TR2)
*   PP-reflector
       CALL REFL(ZRPP,ZDM,ZDM,ZDM
     -              ,DW,N,P,VP2,VS2,DEN2,DEP2,NL2,NL2,IB,TR0)
              GOTO 2
1     CALL REFLSH(ZRSU,ZRSD,DW,N,P,VS,DEN,DEP,NL,LL,TR1)
      CALL CNVRSH(ZR0,DW,N,P,VS1,DEN1,DEP1,NL1,TR2)
2      CONTINUE
* < Delay for depth phase >
          YS=SQRT(1/VS(LL)**2-P**2)
          TS=YS*DH
       IF(IB.EQ.1.OR.IB.EQ.4) THEN
          YP=SQRT(1/VP(LL)**2-P**2)
          TP=YP*DH
       IF(IB.EQ.1)              DELY=TR1+TR2-TP-5.*dt
* an additional delay of 10 sec is put for PP wave:
       IF(IB.EQ.4)              DELY=TR1+TR2-TP-10.
       ELSE
                                DELY=TR1+TR2-TS-5.*dt
       ENDIF
c               write(*,*) tr1,tr2,ts,tp,dely
      DO 3 I=1,N/2
           W=DW*(I-1)
*--------------------------------------------------------------
* P or PP wave
      IF(IB.EQ.1.OR.IB.EQ.4) THEN
        FC=FC0/(4*PI*DEN(NL)*VP(NL)**3)
       Z(I) =      ZRPD(I)* PD*EXP(CMPLX(0.,+W*TP))
*   Exclude pP & sP phases outside the time window
      IF(LL.EQ.NL.AND.TP*2.0.GE.TL) GOTO 11
           Z(I) = Z(I)+ZRPU(I)* PU*EXP(CMPLX(0.,-W*TP))
     -                -ZRSD(I)*SVD*EXP(CMPLX(0.,+W*TS))
     -                -ZRSU(I)*SVU*EXP(CMPLX(0.,-W*TS))
11      Z(I) = Z(I)*ZQP(I)
       IF(IB.EQ.1) GOTO 50
*  PP-reflector & additional Q & Hilbert-transform
         Z(I) = Z(I)*ZRPP(I)*ZQP(I)*CMPLX(0.,1.)
* SV-wave
      ELSEIF(IB.EQ.2) THEN
        FC=FC0/(4*PI*DEN(NL)*VS(NL)**3)
        Z(I) =      ZRSD(I)*SVD*EXP(CMPLX(0.,+W*TS))
*   Exclude pS & sS phases outside the time window
      IF(LL.EQ.NL.AND.TS*2.0.GE.TL) GOTO 21
           Z(I) = Z(I)-ZRPU(I)* PU*EXP(CMPLX(0.,-W*TP))
     -                -ZRPD(I)* PD*EXP(CMPLX(0.,+W*TP))
     -                +ZRSU(I)*SVU*EXP(CMPLX(0.,-W*TS))
21      Z(I) = Z(I)*ZQS(I)
* SH-wave
      ELSEIF(IB.EQ.3) THEN
        FC=FC0/(4*PI*DEN(NL)*VS(NL)**2*VS(LL))
        Z(I) =      ZRSD(I)*SHD*EXP(CMPLX(0.,+W*TS))
*   Exclude sS phase outside the time window
      IF(LL.EQ.NL.AND.TS*2.0.GE.TL) GOTO 31
           Z(I) = Z(I)+ZRSU(I)*SHU*EXP(CMPLX(0.,-W*TS))
31      Z(I) = Z(I)*ZQS(I)
      ENDIF
*--------------------------------------------------------------
50      Z(I)=FC*Z(I)*ZR0(I)*ZI(I)*EXP(CMPLX(0.,W*DELY))
      IF(I.EQ.1) GOTO 3
       Z(N+2-I)=CONJG(Z(I))
3     CONTINUE
	 z(1)=0
         Z(N/2+1)=0
      CALL CFFT(Z,N,1)
      DO 5 I=1,N
5     X(I)=Z(I)*DF
      END
