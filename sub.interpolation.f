c        program interpolation
        subroutine interpolation(dist,gcarc,para,geom,time,tstar,ref_co,
     &              raytype,outpara,outgeom,outtime,outtstar,outref_co)
c       subroutine for interpolatation
c       IN:
c         dist - great circle distance
c         gcarc - 2-Dimension array which stores gcarc read by readray()
c         para  - 2-Dimension array which stores parameters
c         geom  - 2-Dimension array which stores geometrical spreading
c                     factor (1/R in 1/km)
c         time  - 2-Dimension array which stores traveltimes(sec)
c         tstar - 2-Dimension array which stores t* in attenuation
c         ref_co - 2-Dimension array which stores reflection coefficient
c         raytype - 1 for P wave, 2 for PP, 3 for PcP, 4 for S, 5 for SS
c                   6 for ScS
c       OUT:
c         outpara - parameter of given gcarc (dist)
c         outgeom - geometrical spreading factor
c         outtime - traveltime
c         outtstar - t* in attenuation
c         outref_co - reflection coefficient
c
        implicit none
        integer max_ray,raykind
        parameter (max_ray = 10000, raykind = 6)

        real dist
        integer raytype
        integer i
        real gcarc(max_ray,raykind), para(max_ray,raykind)
        real time(max_ray,raykind), geom(max_ray,raykind)
        real tstar(max_ray,raykind), ref_co(max_ray,raykind)

        real factor
        real outpara, outtime, outgeom
        real outtstar, outref_co

        do i = 1, max_ray
          if (abs(dist-gcarc(i,raytype))>0.8) cycle
 
          if ((dist >= gcarc(i-1,raytype) 
     &    .and. dist < =gcarc(i,raytype))
     &   .or. (dist < =gcarc(i-1,raytype) 
     &   .and. dist >= gcarc(i,raytype)))
     &    then
                
            factor = (dist - gcarc(i-1,raytype))
            factor = factor/(gcarc(i,raytype)-gcarc(i-1,raytype))

            outpara   = para(i-1,raytype)
     &                + factor*(para(i,raytype)-para(i-1,raytype))

            outtime   = time(i-1,raytype) 
     &                + factor*( time(i,raytype)-time(i-1,raytype) )

            outgeom   = geom(i-1,raytype) 
     &                + factor*( geom(i,raytype)-geom(i-1,raytype) )

            outtstar  = tstar(i-1,raytype) 
     &                + factor*(tstar(i,raytype)-tstar(i-1,raytype))

            outref_co = ref_co(i-1,raytype) 
     &                + factor*(ref_co(i,raytype)-ref_co(i-1,raytype))

            exit
          endif
        enddo

        end 
