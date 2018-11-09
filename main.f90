!                                 Notice
!  This program was prepared by the University of California (University)
!  under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
!  All rights in the program are reserved by DOE on behalf of the Government
!  and the University pursuant to the contract. You are authorized to use
!  this program for Government purposes but it is not to be released or
!  distributed to the public.
!  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
!  NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
!  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
!  OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
!  ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
!  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
      program main
!     subroutine main
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! QWIC-URB is a multi-building flow solver using 3 dimensional    
! successive overrelaxation method solver and is based on the  
! work of Rockle (1990) and Kaplan and Dinar (1996).     
!                       
! p1 and p2 are cell centered quantites   and called Lagrange  
! Multipliers.                   
! u,v,w and uo,vo,wo are cell face quantites       
! uo,vo, and wo are the inital non-mass conserved velocities   
! u,v,w is the final mass conserved velocity field 
!                       
! bcsetup.f does much of the work in the code, here the ground 
! and buildings are defined, the empircal flow parameterizations
! input and the boundary conditions for the SOR solver set up. 
!                       
!  - the final velocity field is written out in euler.f  
!  - the Lagrangde multiplier field is writtein out in sor3d.f 
!                       
! UPDATES                     
! Eric Pardyjak                     
! QWIC-URBv1.00Beta September 2002        
!  September 2002 Update:           
!  This version has an updated input file. The height of zref  
!  for the power inlet velocity profile has been added.     
! QWIC-URBv1.00 October 2002           
! This version has multiple updates that include      
!  1. new output format of output files & new output files  
!  2. input.dat has a new line for a rooftop recirculation  
!  3. new coordinate system               
!  4. fixes in array out of bounds in sor3d.f         
!
! QWIC-URBv1.00f90 January 2003
!  1. Winds out of the North and South can be handle (ERP 12/17)
!  2. A bug in the street canyon routine was fixed (MDW 1/8)
!  3. Allocatable arrays are now deallocated allowing the qwicurb
!        to be run multiple time in the GUI. (TWH 1/8)
!
! Note if this version of the code is being used with a Matlab GUI
! and Visual Compaq Fortran Compiler 6.6, it will be necessary to do the
! following:
!  1. Download the 6.6B update (VFRUN66BI.exe). It can be found at 
!      the following location: 
!        http://h18009.www1.hp.com/fortran/visual/redist.html
!  2. Once patch has been installed move the file dformd.dll from
!     the C:\Windows\system32 on XP machines or 
!        C:\WINNT\system32 on Windows2000 machines folder to 
!        C:\MATLAB6p5\bin\win32
!  3. Be sure to run mex -setup so that the new df66opts.bat is
!      updated.
!
! erp 1/29/04
! This verson of QUICURB has Petra Kastner-Klein's street canyon
! models available. See subroutine init.f90 and bcsetup.f90
!
! ERP 10/04/04
! This version of QUICURB allows for multiple runs to simulate
! "quasi-steady" time dependent flows
! ERP 3/8/05
! - This version of the code contains the new building sorting logic to sort
! building from small to tall so that the algorithms are applied in that order
! - THis version also uses a local velocity to determine the orientation of a
! street canyon vortex
! ERP 6/17/05
! This version of QUICURB has the basic data assimilation algorithms based
! on Barnes objective mapping as implemented by Tom Booth
! new subroutines include: sensorinit
!
!The following file unit numbers are currently in use in QUIC4.0 
! 8/8/2005 ERP
! File numbering usage:
! unit   name           open location  close location
! 28     uofield.dat       outfile.f90       main.f90
! 33     QU_celltype.dat      bcsetup.f90       main.f90
! 34     QU_velocity.dat      outfile.f90       main.f90
! 35     QU_simparams.inp  init.f90       init.f90
! 36     QU_metparams.inp  init.f90       main.f90
! 37     QU_buildings.inp  init.f90       init.f90
! 38     QU_velocity.bin      outfile.f90       main.f90
! 39     QU_celltype.bin      bcsetup.f90       main.f90
! 46     QP_buildout.inp      init.f90       main.f90
! 47     QU_fileoptions.inp   init.f90       init.f90
! 52     f_name            sensorinit.f90    sensorinit.f90 f_name is character string variable
! 55     QU_veg.inp        sensorinit.f90    sensorinit.f90
! 60     uoutmatu.dat      outfile.f90
! 61     
! 62
! 74        QU_intersect.dat  street_intersect.f90
! 75        QU_icellflagtest.datstreet_intersect.f90
! 
! test cvs
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         use datamodule ! make data from module "datamodule" visible
         implicit none
         integer  idate(6)
!         real t1

!  start timing sequence
!  t1=secnds(0.0)

         alpha1=1.  !correlation of wind components term
         alpha2=1.  !correlation of wind components term
         eta=(alpha1/alpha2)**2 ! ratio of gaussian precision moduli

         ! since omega is also used for an angle in bcsetup

         omegarelax=1.78      !acceleration term in SOR

!   iteration parameters
!         itermax=3000     !max number of iterations
!         eps=1.e-8     !max error
! read input file
         call init
! sort buildings by height small to tall  
         call sort
         if(damageflag .eq. 1)then
            call building_damage
         endif
         if(diffusion_flag .gt. 0)then
            itermax=itermax/diffstep
         endif
!         time=start_time
!     begin multirun pseudo time loop
         do i_time = 1,num_time_steps
!     read met input file for each time step 
            call unix2c(time(i_time)+utc_offset*3600,idate)
            write(*,14001)idate(2),idate(3),idate(1),idate(4),idate(5),idate(6)
14001 format('Calculating winds for ',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2)
            call sensorinit
            max_velmag=0.
            do i=1,nx
               do j=1,ny
                  max_velmag=max(max_velmag,sqrt((uo(i,j,nz)**2.)+(vo(i,j,nz)**2.)))
               enddo
            enddo
            max_velmag=1.2*max_velmag
            ! print*,max_velmag
!     call boundary condition set up routine
            call building_parameterizations
! MAN 06/04/2007 check to see if there are buildings, vegetation, or multiple profiles
            if(inumbuild .eq. 0 .and. num_sites .eq. 1 .and. inumcanopy .eq. 0 &
                  .and. (landuse_flag .eq. 0 .or. (landuse_veg_flag .eq. 0 .and. landuse_urb_flag .eq. 0)))then
               diffusion_flag=0
               itermax=0
            endif
            print*,'Solving for Mass Consistency'
! end MAN 06/04/2007
            if(diffusion_flag .gt. 1)then
               do diffiter= 1,diffstep
                  call divergence
!  call sor3d routine
                  call sor3d
!  call Euler equations to get new updated velocities
                  call euler
!   call Diffusion operator
                  call diffusion
               enddo
            else
!  call divergence routine to calculate divergence of uo field
               call divergence
!  call sor3d routine
               call sor3d
!  call Euler equations to get new updated velocities
!  note that Euler call outfile.f
               call euler
            endif
            call outfile
!            time=time + time_incr*i_time
         enddo

! deallocate allocatable arrays - f90 specific
         
         deallocate(z,zm,dz_array) !MAN 7/21/2008 stretched grid
         deallocate(icellflag,ibldflag, bldnum, bldtype)         ! twh - added this line
         deallocate(bldgeometry,bldstartidx,bldstopidx,numpolygons) !MAN 08/10/2010 polygon buildings
         deallocate(LrNode,LrFace)
         deallocate(uo, vo, wo)                         ! twh - added this line
         deallocate(ufarwake,vfarwake)
         deallocate(p1, p2, r)                          ! twh - added this line
         deallocate(e, f, g, h, m, n, denom, u, v, w) ! twh - added this line
         deallocate(Ht, Wti, Lti, aa, bb)               ! erp 1/31/2003
         deallocate(xfo, yfo, zfo, gamma)                      ! twh - added this line
         deallocate(bld_damage)
         if(inumcanopy .gt. 0)then
            deallocate(canopy_ktop,canopy_top,canopy_atten)
            deallocate(canopy_zo,canopy_ustar,canopy_d)
         endif
         !MAN 8/30/2005 stacked building fix
         deallocate(zfo_actual)
         if(diffusion_flag .gt. 0)deallocate(Fxd,Fyd,Fzd,visc)
         deallocate(u_prof,v_prof)                 !TMB 2/25/05
         deallocate(uoint,voint)                   !TMB 2/25/05
         deallocate(wm,wms)                           !TMB 2/25/05
!ERP 8/17/05
         deallocate(group_id)
!erp 6/08/04 deallocate new variables
         deallocate(uo_bl,vo_bl)
         deallocate(time)
         if(landuse_flag .eq. 1)then
            deallocate(landuse)
            if(landuse_veg_flag .eq. 1 .or. landuse_urb_flag .eq. 1)then
               deallocate(landuse_height,landuse_atten)
               if(inumcanopy .eq. 0)then
                  deallocate(canopy_ktop,canopy_top,canopy_atten)
                  deallocate(canopy_zo,canopy_ustar,canopy_d)
               endif
            endif
         endif
         if(uofield_flag.eq.1) close(28)
         if(format_flag.eq.1 .or. format_flag.eq.3)close(33) !erp 3/02/05
         if(format_flag.eq.1 .or. format_flag.eq.3)close(34) !erp 3/02/05
         if(format_flag.eq.2 .or. format_flag.eq.3)close(38) !erp 3/02/05
         if(format_flag.eq.2 .or. format_flag.eq.3)close(39) !erp 3/02/05
!  close(28)      !close uofield.dat

         close(36)      !close QU_metparams.inp
         close(46)      !close QP_buildout.inp
         if(staggered_flag .eq. 1)close(99)
         if(canopy_flag .gt. 0)close(100)
         !close(48)     !close QP_buildorder.inp

         !TMB 3/11/05 I open/close these files in sensorinit, but I want to have a complete list of all
         !the files opened/closed listed here so there would be no confusion later on

         !close(51)     !close uosensorfield.dat   this is now done in sensorinit
         !close(52)     !close the individual profile files with this is also done in sensorinit
   
         !TMB end

!  end timing sequence
!  print*,secnds(t1)

!ifndef PGI_DEBUG
!else
    stop
!endif


!  return
      end
      
      subroutine dateToUnix(idate, utime)
         implicit none
         integer(8),intent(in)  :: idate(6)
         integer(8),intent(out) :: utime
         integer(8) i,days(12)
      
         days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         !initially, set the timestamp to be zero
         utime = 0
         
         !determine if the current year is a leap year, which will potentially modify the 
         !number of days in February
         if(mod(idate(1),400) == 0) then
            days(2) = 29
         elseif(mod(idate(1),100) == 0) then
            days(2) = 28
         elseif(mod(idate(1),4) == 0) then 
            days(2) = 29
         endif
         
         !then check if the input values are valid
         if(idate(1) .lt. 1901 .or. idate(1) .gt. 2038)then
            print*, "dateToUnix: invalid year specified", idate(1)
            return
         elseif(idate(2) .lt. 1 .or. idate(2) .gt. 12)then
            print*, "dateToUnix: invalid month specified", idate(2)
            return
         elseif(idate(3) .lt. 1 .or. idate(3) .gt. days(idate(2)))then
            print*, "dateToUnix: invalid day specified", idate(3)
            return
         elseif(idate(4) .lt. 0 .or. idate(4) .gt. 23)then
            print*, "dateToUnix: invalid hour specified", idate(4)
            return
         elseif(idate(5) .lt. 0 .or. idate(5) .gt. 59)then
            print*, "dateToUnix: invalid minute specified", idate(5)
            return
         elseif(idate(6) .lt. 0 .or. idate(6) .gt. 59)then
            print*, "dateToUnix: invalid second specified", idate(6)
            return
         endif
         
         !determine if the specified date is before or after the beginning of the epoch
         if(idate(1) .lt. 1970)then
            !then calculate the number of seconds between jan 1, 1970, 12:00am and the beginning
            !of the specified year, moving backwards (before beginning of epoch)
            do i=1970, idate(1) + 1, -1
               !leap years are years which are multiples of 400 or years that are multiples of 
               ! 4 and not 100
               if(mod(i,400) == 0) then
                  utime = utime - 31622400 !366 * 24 * 60 * 60 = 31622400s
               elseif(mod(i,100) == 0) then
                  utime = utime - 31536000 !365 * 24 * 60 * 60 = 31536000s
               elseif(mod(i,4) == 0) then 
                  utime = utime - 31622400 !366 * 24 * 60 * 60 = 31622400s
               else
                  utime = utime - 31536000 !365 * 24 * 60 * 60 = 31536000s
               endif
            enddo
            if(days(2) == 29)then
               utime = utime - 86400 ! on a leap year, the beginning of the year is a day further away from 0
            endif
         else
            !then calculate the number of seconds between jan 1, 1970, 12:00am and the beginning
            !of the specified year (after beginning of epoch)
            do i=1970, idate(1) - 1
               !leap years are years which are multiples of 400 or years that are multiples of 
               ! 4 and not 100
               if(mod(i,400) == 0) then
                  utime = utime + 31622400 !366 * 24 * 60 * 60 = 31622400s
               elseif(mod(i,100) == 0) then
                  utime = utime + 31536000 !365 * 24 * 60 * 60 = 31536000s
               elseif(mod(i,4) == 0) then 
                  utime = utime + 31622400 !366 * 24 * 60 * 60 = 31622400s
               else
                  utime = utime + 31536000 !365 * 24 * 60 * 60 = 31536000s
               endif
            enddo
         endif
         
         !then calculate the number of seconds between the beginning of the month specified
         !and the beginning of the year specified
         do i=1,idate(2) - 1
            utime = utime + days(i) * 86400 ! 24 * 60 * 60 = 86400s
         enddo
         
         !From here on out, everything should be standard.  no leap year stuff to deal with, 
         !and leap seconds are not accounted for in UNIX time
         utime = utime + (idate(3) - 1) * 86400 ! 24 * 60 * 60 = 86400s
         utime = utime + idate(4) * 3600 ! 60 * 60 = 3600s
         utime = utime + idate(5) * 60
         utime = utime + idate(6) 
      end subroutine
      
      ! unix2c Converts Unix system time to date/time integer array.
      subroutine unix2c(utime, idate)
         implicit none
         integer utime, idate(6), days(12)
         ! utime  input  Unix system time, seconds since 1970.0
         ! idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
         ! -Author  Clive Page, Leicester University, UK.   1995-MAY-2
         integer mjday, nsecs
         real day
         ! Note the MJD algorithm only works from years 1901 to 2099.
         mjday    = int(utime/86400 + 40587)
         idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
         day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
         idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
         idate(3) = 1 + int(mod(day,30.6))
         nsecs    = mod(utime, 86400)
         idate(6) = mod(nsecs, 60)
         nsecs    = nsecs / 60
         idate(5) = mod(nsecs, 60)
         idate(4) = nsecs / 60
         
         days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         !initially, set the timestamp to be zero
         utime = 0
         
         !determine if the current year is a leap year, which will potentially modify the 
         !number of days in February
         if(mod(idate(1),400) == 0) then
            days(2) = 29
         elseif(mod(idate(1),100) == 0) then
            days(2) = 28
         elseif(mod(idate(1),4) == 0) then 
            days(2) = 29
         endif
         
         !normalize seconds
         do while(idate(6) .lt. 0)
            idate(5) = idate(5) - 1
            idate(6) = idate(6) + 60
         enddo
         do while(idate(6) .gt. 59)
            idate(5) = idate(5) + 1
            idate(6) = idate(6) - 60
         enddo
         !normalize minutes
         do while(idate(5) .lt. 0)
            idate(4) = idate(4) - 1
            idate(5) = idate(5) + 60
         enddo
         do while(idate(5) .gt. 59)
            idate(4) = idate(4) + 1
            idate(5) = idate(5) - 60
         enddo
         !normalize hours
         do while(idate(4) .lt. 0)
            idate(3) = idate(3) - 1
            idate(4) = idate(4) + 24
         enddo
         do while(idate(4) .gt. 23)
            idate(3) = idate(3) + 1
            idate(4) = idate(4) - 24
         enddo
         !normalize days
         !note that day is a special case.  It can increment or decrement month, and in order
         !to get the number of days to add/subract, we need to know the normalized month.  As
         !such, the month has to be normalized in the middle of normalizing the day
         do while(idate(3) .lt. 1)
            idate(2) = idate(2) - 1
            !normalize month
            do while(idate(2) .lt. 1)
               idate(1) = idate(1) - 1
               idate(2) = idate(2) + 12
            enddo
            do while(idate(2) .gt. 12)
               idate(1) = idate(1) + 1
               idate(2) = idate(2) - 12
            enddo
            idate(3) = idate(3) + days(idate(2))
         enddo
         do while(idate(3) .gt. days(idate(2)))
            idate(2) = idate(2) + 1
            !normalize month
            do while(idate(2) .lt. 1)
               idate(1) = idate(1) - 1
               idate(2) = idate(2) + 12
            enddo
            do while(idate(2) .gt. 12)
               idate(1) = idate(1) + 1
               idate(2) = idate(2) - 12
            enddo
            idate(3) = idate(3) - days(idate(2))
         enddo
         !normalize month 
         do while(idate(2) .lt. 1)
            idate(1) = idate(1) - 1
            idate(2) = idate(2) + 12
         enddo
         do while(idate(2) .gt. 12)
            idate(1) = idate(1) + 1
            idate(2) = idate(2) - 12
         enddo
      end
