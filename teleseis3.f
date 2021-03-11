**===================================================**
**    Calculate  Green's functions W(i;n,k)          **
**
**===================================================**
**  Calculate direct and core-reflected phase        **
**  t* changes with epicentral distances             **
**  Yunyi Qian, 09-12-2015                           **
**===================================================**
      PARAMETER (ND0=8192,NM0=6,LK0=10,NL0=20,PI=3.141593,RAD=.0174533)
      IMPLICIT COMPLEX*8 (Z)
      CHARACTER NAME*40,NAM*4, prefix1*80, prefix2*80,
     -          outfile1*80, outfile2*80
      character*1 cm(9), model*80
      character*1 greenV(9)
      character   sta*10
      real*4 depth
      real*8 gs(500),ps(500),vp00,gscs(500),pscs(500),vs00
      DIMENSION w(ND0),ZI(ND0),ZQScS(ND0),ZQS(ND0),ZQP(ND0)
     -,     ZZ(50),ZP0(50),h(lk0),ZQPcP(ND0)
     -,     ib(3),ic(3),p(3),g(3),time1(3),time2(3)
      COMMON /STR0/NL ,VP (NL0),VS (NL0),DEN (NL0),DEP (NL0)
      COMMON /STR1/NL1,VP1(NL0),VS1(NL0),DEN1(NL0),DEP1(NL0)
      COMMON /STR2/NL2,VP2(NL0),VS2(NL0),DEN2(NL0),DEP2(NL0)
      common /sourceRegion/vsrc
      dimension gcarc(10000,6), para(10000,6)
      dimension geom(10000,6), time(10000,6)
      dimension tstar(10000,6), ref_co(10000,6)
      real*4 parapcp,timepcp,geompcp,tstarpcp,ref_copcp
      real*4 parascs,timescs,geomscs,tstarscs,ref_coscs
      real*4 paras,times,geoms,tstars,ref_cos
      real*4 parap,timep,geomp,tstarp,ref_cop


* < Components and phases >
      data cm/'Z', 'R','T', 'Z', 'R', 'T','Z', 'R', 'T'/
      data greenV/'0','1','2','3','4','5','6','7','8'/
      data ib/1, 1, 3/
      data ic/1, 2, 2/


      nl= 4
      data vp/  5.2, 6.24, 6.58, 7.8, 16*0.0/
      data vs/  3.0, 3.60, 3.80, 4.4, 16*0.0/
      data den/ 2.4, 2.67, 2.80, 3.3, 16*0.0/
      data dep/ 5.5, 9.50, 19.0, 0.0, 16*0.0/

      nl1=3
      data vp1/  5.2, 6.40, 8.0, 17*0.0/
      data vs1/  3.0, 3.75, 4.5, 17*0.0/
      data den1/ 2.4, 2.70, 3.3, 17*0.0/
      data dep1/ 5.0, 30.0, 0.0, 17*0.0/


* < Default Green's function parameters >
      xcorrection = 10**(-5.0)
      tstart=-30.0


c ------------ End initialization ----------------------------
      read(*,'(a80)') model
      if(model(1:4) .ne. 'none') then
         open(2,file=model)
         READ(2,'(a40)') name
         READ(2,*) NL ,(VP (L),VS (L),DEN (L),DEP (L),L=1,NL )
         READ(2,*) NL1,(VP1(L),VS1(L),DEN1(L),DEP1(L),L=1,NL1)
         r0=6371
         close(2)
      endif
        read(*,*) nt, dt
        nstart=-tstart/dt
        
        DF=1/(DT*NT)
        DW=DF*2*PI

* < Set ZI >
        CALL INSTG(ZI,NT,DW,0,ZP0,ZZ,0,0,1.,0)
        read(*,*) depth

* < Travel time >
        call readray(gcarc,para,geom,time,tstar,ref_co)


* < Direct S Synthetics >
        do 100, j=1,1000
           READ(*,*,END=999) az, dist, sta
           read(*,'(a80)') prefix1
           read(*,'(a80)') prefix2
           idx1=index(prefix1,' ')-1
           idx2=index(prefix2,' ')-1
           call interpolation(dist,gcarc,para,geom,time,tstar,ref_co,4
     -                  , paras,geoms,times,tstars,ref_cos)
           call interpolation(dist,gcarc,para,geom,time,tstar,ref_co,1
     -                  , parap,geomp,timep,tstarp,ref_cop)
           time1(1)=timep+tstart
           time1(2)=timep+tstart
           time1(3)=times+tstart
           call qf(zqs,nt,tstars,df)
           call qf(zqp,nt,tstarp,df)
           p(1)=parap/180*pi/111.2*r0/(r0-depth)
           g(1)=geomp
           p(2)=parap/180*pi/111.2*r0/(r0-depth)
           g(2)=geomp
           p(3)=paras/180*pi/111.2*r0/(r0-depth)
           g(3)=geoms
           DO 99 K=1,3  
                do 98 iFlt=1,3
                   CALL BODYW13(iFlt,w,NT,DT,IB(k),IC(k),
     -             depth,p(k),g(k),zqp,zqs,ZI,dely,nt)
                   outfile1=prefix1(1:idx1)//greenV((iFlt-1)*3+K)
     -             //char(0)
                   do 96, m=nt,nstart+1,-1
                      w(m)=w(m-nstart)-w(1)
                      w(m)= w(m)*xcorrection
96                 continue
                   do 95, m=1,nstart
                       w(m)=0.0
95                 continue
                   do 94, m=nt,2,-1
                      w(m) = (w(m)-w(m-1))/dt
94                 continue
                   call wrtsac1(outfile1,dt,nt,time1(k),
     -                 time1(1),time1(3),dist,az,xstart,w)
98              continue   
99         continue
          
           call interpolation(dist,gcarc,para,geom,time,tstar,ref_co,6
     -                  ,parascs,geomscs,timescs,tstarscs,ref_coscs)
           call interpolation(dist,gcarc,para,geom,time,tstar,ref_co,3
     -                  ,parapcp,geompcp,timepcp,tstarpcp,ref_copcp)
           time2(1)=timepcp+tstart
           time2(2)=timepcp+tstart
           time2(3)=timescs+tstart          
           call qf(zqscs,nt,tstarscs,df)
           call qf(zqpcp,nt,tstarpcp,df)
           p(1)=parapcp/180*pi/111.2*r0/(r0-depth)
           g(1)=geompcp*ref_copcp
           p(2)=parapcp/180*pi/111.2*r0/(r0-depth)
           g(2)=geompcp*ref_copcp
           p(3)=parascs/180*pi/111.2*r0/(r0-depth)
           g(3)=geomscs
           DO 199 K=1,3
                do 198 iFlt=1,3
                   CALL BODYW13(iFlt,w,NT,DT,IB(k),IC(k),
     -             depth,p(k),g(k),zqpcp,zqscs,ZI,dely,nt)
                   outfile2=prefix2(1:idx2)//greenV((iFlt-1)*3+K)
     -                          //char(0)
                   do 196, m=nt,nstart+1,-1
                      w(m)=w(m-nstart)-w(1)
                      w(m)= w(m)*xcorrection
196                continue
                   do 195, m=1,nstart
                       w(m)=0.0
195                continue
                   do 194, m=nt,2,-1
                      w(m) = (w(m)-w(m-1))/dt
194                continue
                   call wrtsac1(outfile2,dt,nt,time2(k),
     -                   time1(1),time1(3),dist,az,xstart,w)
198             continue
199        continue
100     continue
999     CONTINUE
      END
