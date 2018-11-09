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
      subroutine sensorinit
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subroutine to initialize the initial velocity field (uo,vo,wo) for each time step
! 
! This program interpolates irregularly spaced data onto a 
! uniform grid using the Barnes Objective Map Analysis  
! scheme as implemented in Koch et al., 1983
! 
! This program has been modified from Eric Pardyjak's original 2D code to work with 
! quicurb 3D      
!
! this subroutine uses a quasi-3D barnes objective mapping technique 
! quasi-3D is just using sensor height (zref) and sensor wind speed (uref)
! to extrapolate a vertical velocity profile at each sensor location 
! to get a velocity at every height at each location
! from these hieght varying velocities, a regular 2D barnes objective map analysis
! is done at each planar height throughout the whole domain.
! 
! Called by met_init.f90
! 
! 
! 
! Tom Booth 2/17/04
! 
! Variable information:
!  site_xcoord,site_ycoord,site_zcoord - the coordinates of each site (meters)
!  t_site - is the time step for each site
!  dir_site - is the direction of the wind for each site at each time step
!  vel_site - is the magnitude of the wind for each site at each time step
!  total_time_step - is the total number of time steps in a 24 hr period
!  num_sites - is the total number of measurement sites
!  sgamma - numerical convergence parameter
!  lamda - weight parameter (ko)
!  deln - average Radius (i.e. computed data spacing)
! TMB/ERP 9/20/05 
!  Added a logarithmic interpolation below the the lowest input data point for
!  input files that read in wind profile data

!  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule
         implicit none

         integer kk,jj !,site_num
         integer iwork, jwork
         integer met_input_flag


         ! real site_time_incr
         real site_zref, theta_site, mag_site, umult_site, vmult_site

         real yc, xc, rc, rcsum, rcval
         real sumw, sumwv, sumwu
         real sgamma, lamda, deln
         real dxx, dyy, u12, u34, v12, v34
! MAN 04/05/2007 Time varying profile parameters
         real bisect
         real xtemp,psi_m !erp 9/18/2006 stability variables
         
         real, allocatable :: x(:,:), y(:,:)   ! Locations of grid cell centers in meters
         ! MAN 10/23/2013 advanced data point interpolation
         integer iter,logflag
         real zonew,zolow,zohigh,unew,unewl,unewh
         real dwinddir,zval,dwsdzu,dwsdzl
         real sumprofile1,sumprofile2,profileratio,a1,a2,a3
         real utot,utotu,utotl,z1,z2,z3,denoml,denomm,denomu
         integer ku,kl,km
         ! real, allocatable :: profile1(:),profile2(:),profiletheta(:)
         ! real, allocatable :: site_dwsdz(:)
         
         allocate(x(nx,ny),y(nx,ny))
         ! allocate(profile1(nz),profile2(nz),profiletheta(nz))
         ! allocate(site_dwsdz(num_vert_points))
         
         ! read in whole file on first quasi-time step
         if (i_time.eq.1) then
            vk=0.4 !Von Karman's constant
            read(36,*) ! QUIC version header line
            read(36,*)met_input_flag
! MAN 4/17/2007 moved met data read statements to separate subroutines and added architecture for new 3D vectorfield input
            if(met_input_flag .le. 1)then
               call read_quic_met
            elseif(met_input_flag .eq. 2)then
               call read_ittmm5_met
            else
               call read_hotmac_met
            endif
         endif ! if i_time = 1
         
         ! ---------------------------------------------------------------------------------
         ! This do loop defines each vertical velocity profile at each sensor site so the
         ! Barnes Mapping interpolation scheme has a velocity at each site for each height
lp003:   do kk=1,num_sites
            if(site_blayer_flag(kk,i_time) .lt. 4)then
               theta_site=(270.-site_wd_data(kk,i_time,1))*pi/180.
               mag_site = site_ws_data(kk,i_time,1)
               site_zref = site_z_data(kk,i_time,1)
               umult_site=cos(theta_site)
               vmult_site=sin(theta_site)
! power law profile
               if(site_blayer_flag(kk,i_time) .eq. 2)then
                  ! MAN 07/25/2008 stretched vertical grid
                  do k=2,nz      ! loop through vertical cell indices
                     u_prof(kk,k)=umult_site*mag_site*((zm(k)/site_zref)**site_pp(kk,i_time))
                     v_prof(kk,k)=vmult_site*mag_site*((zm(k)/site_zref)**site_pp(kk,i_time))
                  enddo
               endif !erp 2/6/2003 end power law profile
! logrithmic velocity profile
               if(site_blayer_flag(kk,i_time) .eq. 1)then
! MAN 05/15/2007 adjust for stability
                  do k=2,nz      ! loop through vertical cell indices
                     if(k .eq. 2)then
                        if(site_zref*site_rL(kk,i_time) .ge. 0)then
                           psi_m=4.7*site_zref*site_rL(kk,i_time)
                        else
                           xtemp=(1.-15*site_zref*site_rL(kk,i_time))**(0.25)
                           psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                        endif
                        ustar=mag_site*vk/(log((site_zref+site_pp(kk,i_time))/site_pp(kk,i_time))+psi_m)
                     endif
! end MAN 05/15/2007        
                     ! MAN 07/25/2008 stretched vertical grid             
                     if(zm(k)*site_rL(kk,i_time) .ge. 0)then
                        psi_m=4.7*zm(k)*site_rL(kk,i_time)
                     else
                        xtemp=(1.-15*zm(k)*site_rL(kk,i_time))**(0.25)
                        psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                     endif
                     u_prof(kk,k)=(umult_site*ustar/vk)*(log((zm(k)+site_pp(kk,i_time)) &
                           /site_pp(kk,i_time))+psi_m)
                     v_prof(kk,k)=(vmult_site*ustar/vk)*(log((zm(k)+site_pp(kk,i_time)) &
                           /site_pp(kk,i_time))+psi_m)
                  enddo
               endif ! erp 2/6/2003 end log law profile
! Canopy profile
               if(site_blayer_flag(kk,i_time) .eq. 3)then
                  do k=2,nz      ! loop through vertical cell indices
                     if(k .eq. 2)then ! only calculate d once
! MAN 05/15/2007 adjust for stability                        
                        if(site_zref*site_rL(kk,i_time) .ge. 0)then
                           psi_m=4.7*site_zref*site_rL(kk,i_time)
                        else
                           xtemp=(1.-15*site_zref*site_rL(kk,i_time))**(0.25)
                           psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                        endif
                        ustar = mag_site*vk/(log(site_zref/site_pp(kk,i_time))+psi_m)                        
                        d = bisect(ustar,site_pp(kk,i_time),site_H(kk,i_time),site_ac(kk,i_time),vk,psi_m)
! end MAN 05/15/2007
                        if(site_H(kk,i_time)*site_rL(kk,i_time) .ge. 0)then
                           psi_m=4.7*(site_H(kk,i_time)-d)*site_rL(kk,i_time)
                        else
                           xtemp=(1.-15*(site_H(kk,i_time)-d)*site_rL(kk,i_time))**(0.25)
                           psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                        endif
                        uH = (ustar/vk)*(log((site_H(kk,i_time)-d)/site_pp(kk,i_time))+psi_m);
                        if(site_zref .le. site_H(kk,i_time))then
                           mag_site=mag_site/(uH*exp(a*(site_zref/site_H(kk,i_time) -1.)))
                        else
                           if(site_zref*site_rL(kk,i_time) .ge. 0)then
                              psi_m=4.7*(site_zref-d)*site_rL(kk,i_time)
                           else
                              xtemp=(1.-15*(site_zref-d)*site_rL(kk,i_time))**(0.25)
                              psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                           endif
                           mag_site=mag_site/((ustar/vk)*(log((site_zref-d)/zo)+psi_m))
                        endif
                        ustar=mag_site*ustar
                        uH=mag_site*uH
                     endif
                     ! MAN 07/25/2008 stretched vertical grid
                     if(zm(k) .le. site_H(kk,i_time))then ! lower canopy profile
                        u_prof(kk,k) = umult_site * uH*exp(site_ac(kk,i_time)&
                                       *(zm(k)/site_H(kk,i_time) -1))
                        v_prof(kk,k) = vmult_site * uH*exp(site_ac(kk,i_time)&
                                       *(zm(k)/site_H(kk,i_time) -1))
                     endif
                     if(zm(k) .gt. site_H(kk,i_time))then ! upper canopy profile
                        if(zm(k)*site_rL(kk,i_time) .ge. 0)then
                           psi_m=4.7*(zm(k)-d)*site_rL(kk,i_time)
                        else
                           xtemp=(1.-15*(zm(k)-d)*site_rL(kk,i_time))**(0.25)
                           psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
                        endif
                        u_prof(kk,k)=(umult_site*ustar/vk)*&
                                  (log((zm(k)-d)/site_pp(kk,i_time))+psi_m)
                                  
                        v_prof(kk,k)=(vmult_site*ustar/vk)*&
                                  (log((zm(k)-d)/site_pp(kk,i_time))+psi_m)
                     endif !end urban canopy TMB 6/16/03
                  enddo
               endif
            endif
! new 2/7/2005 velocity profile entry
            if(site_blayer_flag(kk,i_time) .eq. 4)then  ! data entry profile
               ii=0
               theta_site=(270.-site_wd_data(kk,i_time,1))*pi/180.
               do k=2,nz ! loop through vertical cell indices
                  if(zm(k) .lt. site_z_data(kk,i_time,1) .or. site_nz_data(kk,i_time) .eq. 1)then
                     u_prof(kk,k)= (site_ws_data(kk,i_time,1)*cos(theta_site)/(log((site_z_data(kk,i_time,1)+site_pp(kk,i_time)) &
                        /site_pp(kk,i_time))))*log((zm(k)+site_pp(kk,i_time))/site_pp(kk,i_time))
                     v_prof(kk,k)= (site_ws_data(kk,i_time,1)*sin(theta_site)/(log((site_z_data(kk,i_time,1)+site_pp(kk,i_time)) &
                        /site_pp(kk,i_time))))*log((zm(k)+site_pp(kk,i_time))/site_pp(kk,i_time))
                  else
                     if(ii .lt. site_nz_data(kk,i_time)-1 .and. zm(k) .ge. site_z_data(kk,i_time,ii+1))then
                        ii=ii+1
                        if(abs(site_wd_data(kk,i_time,ii+1)-site_wd_data(kk,i_time,ii)) .gt. 180.)then
                           if(site_wd_data(kk,i_time,ii+1) .gt. site_wd_data(kk,i_time,ii))then
                              dwinddir=(site_wd_data(kk,i_time,ii+1)-360.-site_wd_data(kk,i_time,ii)) &
                                 /(site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii))
                           else
                              dwinddir=(site_wd_data(kk,i_time,ii+1)-site_wd_data(kk,i_time,ii)+360.) &
                                 /(site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii))
                           endif
                        else
                           dwinddir=(site_wd_data(kk,i_time,ii+1)-site_wd_data(kk,i_time,ii)) &
                              /(site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii))
                        endif
                        zohigh=20.
                        ustar=vk*site_ws_data(kk,i_time,ii)/log((site_z_data(kk,i_time,ii)+zohigh)/zohigh)
                        unewh=(ustar/vk)*log((site_z_data(kk,i_time,ii)+zohigh)/zohigh)
                        zolow=1.e-9
                        ustar=vk*site_ws_data(kk,i_time,ii)/log((site_z_data(kk,i_time,ii)+zolow)/zolow)
                        unewl=(ustar/vk)*log((site_z_data(kk,i_time,ii+1)+zolow)/zolow)
                        if(site_ws_data(kk,i_time,ii+1) .gt. unewl .and. site_ws_data(kk,i_time,ii+1) .lt. unewh)then
                           logflag=1
                           iter=0
                           zonew=site_pp(kk,i_time)
                           ustar=vk*site_ws_data(kk,i_time,ii)/log((site_z_data(kk,i_time,ii)+zonew)/zonew)
                           unew=(ustar/vk)*log((site_z_data(kk,i_time,ii+1)+zonew)/zonew)
                           do while(iter .lt. 200 .and. abs(unew-site_ws_data(kk,i_time,ii)) &
                                 .gt. 0.0001*site_ws_data(kk,i_time,ii))
                              iter=iter+1
                              zonew=0.5*(zolow+zohigh)
                              ustar=vk*site_ws_data(kk,i_time,ii)/log((site_z_data(kk,i_time,ii)+zonew)/zonew)
                              unew=(ustar/vk)*log((site_z_data(kk,i_time,ii+1)+zonew)/zonew)
                              if(unew .gt. site_ws_data(kk,i_time,ii+1))then
                                  zohigh=zonew
                              else
                                  zolow=zonew
                              endif
                           enddo
                        else
                           logflag=0
                           if(ii .lt. site_nz_data(kk,i_time)-1)then
                              a1=((site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii)) &
                                 *(site_ws_data(kk,i_time,ii+2)-site_ws_data(kk,i_time,ii)) &
                                 +(site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,ii+2)) &
                                 *(site_ws_data(kk,i_time,ii+1)-site_ws_data(kk,i_time,ii))) &
                                 /((site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii)) &
                                 *((site_z_data(kk,i_time,ii+2)**2)-(site_z_data(kk,i_time,ii)**2)) &
                                 +((site_z_data(kk,i_time,ii+1)**2)-(site_z_data(kk,i_time,ii)**2)) &
                                 *(site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,ii+2)))
                           else
                              a1=0
                           endif
                           a2=((site_ws_data(kk,i_time,ii+1)-site_ws_data(kk,i_time,ii)) &
                              -a1*((site_z_data(kk,i_time,ii+1)**2)-(site_z_data(kk,i_time,ii)**2))) &
                              /(site_z_data(kk,i_time,ii+1)-site_z_data(kk,i_time,ii))
                           a3=site_ws_data(kk,i_time,ii)-a1*(site_z_data(kk,i_time,ii)**2) &
                              -a2*site_z_data(kk,i_time,ii);
                        endif
                     endif
                     if(logflag .eq. 1)then
                        mag_site=(ustar/vk)*log((zm(k)+zonew)/zonew);
                     else
                        mag_site=a1*(zm(k)**2)+a2*zm(k)+a3
                     endif
                     theta_site=(270.-(site_wd_data(kk,i_time,ii)+ &
                        dwinddir*(zm(k)-site_z_data(kk,i_time,ii))))*pi/180.
                     u_prof(kk,k)= mag_site*cos(theta_site)
                     v_prof(kk,k)= mag_site*sin(theta_site)
                  endif
               enddo
               ! site_dwsdz(:)=0.
               ! profile1(:)=0.
               ! profile2(:)=0.
               ! profiletheta(:)=0.
               ! k=2
               ! do ii=2,site_nz_data(kk,i_time)
               !    zohigh=20.
               !    ustar=vk*site_ws_data(kk,i_time,ii-1)/log((site_z_data(kk,i_time,ii-1)+zohigh)/zohigh)
               !    unewh=(ustar/vk)*log((site_z_data(kk,i_time,ii)+zohigh)/zohigh)
               !    zolow=1.e-9
               !    ustar=vk*site_ws_data(kk,i_time,ii-1)/log((site_z_data(kk,i_time,ii-1)+zolow)/zolow)
               !    unewl=(ustar/vk)*log((site_z_data(kk,i_time,ii)+zolow)/zolow)
               !    if(site_ws_data(kk,i_time,ii) .gt. unewl .and. site_ws_data(kk,i_time,ii) .lt. unewh)then
               !       iter=0
               !       zonew=site_pp(kk,i_time)
               !       ustar=vk*site_ws_data(kk,i_time,ii-1)/log((site_z_data(kk,i_time,ii-1)+zonew)/zonew)
               !       unew=(ustar/vk)*log((site_z_data(kk,i_time,ii)+zonew)/zonew)
               !       do while(iter .lt. 200 .and. abs(unew-site_ws_data(kk,i_time,ii)) &
               !             .gt. 0.0001*site_ws_data(kk,i_time,ii))
               !          iter=iter+1
               !          zonew=0.5*(zolow+zohigh)
               !          ustar=vk*site_ws_data(kk,i_time,ii-1)/log((site_z_data(kk,i_time,ii-1)+zonew)/zonew)
               !          unew=(ustar/vk)*log((site_z_data(kk,i_time,ii)+zonew)/zonew)
               !          if(unew .gt. site_ws_data(kk,i_time,ii))then
               !              zohigh=zonew
               !          else
               !              zolow=zonew
               !          endif
               !       enddo
               !       site_dwsdz(ii)=(ustar/vk)/(site_z_data(kk,i_time,ii)+zonew)
               !       if(ii .eq. 2)then
               !          site_dwsdz(ii-1)=(ustar/vk)/(site_z_data(kk,i_time,ii-1)+zonew)
               !       endif
               !    else
               !       if(ii .gt. k .and. k .lt. site_nz_data(kk,i_time)-1)k=ii
               !       if(site_nz_data(kk,i_time) .gt. 2)then
               !          site_dwsdz(ii)=site_ws_data(kk,i_time,k-1) &
               !                *(2*site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k+1)) &
               !                /((site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k)) &
               !                *(site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k+1))) &
               !                +site_ws_data(kk,i_time,k) &
               !                *(2*site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k+1)) &
               !                /((site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k-1)) &
               !                *(site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k+1))) &
               !                +site_ws_data(kk,i_time,k+1) &
               !                *(2*site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k)) &
               !                /((site_z_data(kk,i_time,k+1)-site_z_data(kk,i_time,k-1)) &
               !                *(site_z_data(kk,i_time,k+1)-site_z_data(kk,i_time,k)))
               !          if(ii .eq. 2)then
               !             site_dwsdz(ii-1)=site_ws_data(kk,i_time,k-1) &
               !                *(2*site_z_data(kk,i_time,ii-1)-site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k+1)) &
               !                /((site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k)) &
               !                *(site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k+1))) &
               !                +site_ws_data(kk,i_time,k) &
               !                *(2*site_z_data(kk,i_time,ii-1)-site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k+1)) &
               !                /((site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k-1)) &
               !                *(site_z_data(kk,i_time,k)-site_z_data(kk,i_time,k+1))) &
               !                +site_ws_data(kk,i_time,k+1) &
               !                *(2*site_z_data(kk,i_time,ii-1)-site_z_data(kk,i_time,k-1)-site_z_data(kk,i_time,k)) &
               !                /((site_z_data(kk,i_time,k+1)-site_z_data(kk,i_time,k-1)) &
               !                *(site_z_data(kk,i_time,k+1)-site_z_data(kk,i_time,k)))
               !          endif
               !       else
               !          site_dwsdz(ii)=(site_ws_data(kk,i_time,ii)-site_ws_data(kk,i_time,ii-1)) &
               !                /(site_z_data(kk,i_time,ii)-site_z_data(kk,i_time,ii-1))
               !          site_dwsdz(ii-1)=site_dwsdz(ii)
               !       endif
               !    endif
               !    if(site_z_data(kk,i_time,ii) .gt. z(nz))exit
               ! enddo
               ! theta_site=(270.-site_wd_data(kk,i_time,1))*pi/180.
               ! mag_site=site_ws_data(kk,i_time,1)
               ! if(site_nz_data(kk,i_time) .eq. 2)then
               !    jj=1
               !    if(abs(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)) .gt. 180.)then
               !       if(site_wd_data(kk,i_time,jj+1) .gt. site_wd_data(kk,i_time,jj))then
               !          dwinddir=(site_wd_data(kk,i_time,jj+1)-360.-site_wd_data(kk,i_time,jj)) &
               !                /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !       else
               !          dwinddir=(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)+360.) &
               !                /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !       endif
               !    else
               !       dwinddir=(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)) &
               !                /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !    endif
               ! else
               !    jj=0
               ! endif
               ! do k=2,nz      ! loop through vertical cell indices
               !    if(zm(k) .lt. site_z_data(kk,i_time,1))then
               !       profile1(k)=(site_ws_data(kk,i_time,1)/(log((site_z_data(kk,i_time,1)+site_pp(kk,i_time)) &
               !             /site_pp(kk,i_time))))*log((zm(k)+site_pp(kk,i_time))/site_pp(kk,i_time))
               !       profile2(k)=profile1(k)
               !       profiletheta(k)=theta_site
               !       dwsdzl=(site_ws_data(kk,i_time,1)/log((site_z_data(kk,i_time,1)+site_pp(kk,i_time)) &
               !             /site_pp(kk,i_time)))/(zm(k)+site_pp(kk,i_time))
               !    else
               !       mag_site=profile2(k-1)
               !       if(site_nz_data(kk,i_time) .gt. 2)then
               !          if(jj .lt. site_nz_data(kk,i_time)-1 .and. zm(k) .ge. site_z_data(kk,i_time,jj+1))then
               !             jj=jj+1
               !             if(abs(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)) .gt. 180.)then
               !                if(site_wd_data(kk,i_time,jj+1) .gt. site_wd_data(kk,i_time,jj))then
               !                   dwinddir=(site_wd_data(kk,i_time,jj+1)-360.-site_wd_data(kk,i_time,jj)) &
               !                         /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !                else
               !                   dwinddir=(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)+360.) &
               !                         /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !                endif
               !             else
               !                dwinddir=(site_wd_data(kk,i_time,jj+1)-site_wd_data(kk,i_time,jj)) &
               !                         /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))
               !             endif
               !          endif
               !          theta_site=(270.-(site_wd_data(kk,i_time,jj)+ &
               !             dwinddir*(zm(k)-site_z_data(kk,i_time,jj))))*pi/180.
               !          dwsdzu=((site_dwsdz(jj+1)-site_dwsdz(jj)) &
               !             /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))) &
               !             *(zm(k)-site_z_data(kk,i_time,jj))+site_dwsdz(jj)
               !          mag_site=mag_site+0.5*(dwsdzu+dwsdzl)*(zm(k)-zm(k-1))
               !          profile2(k)=mag_site
               !          profiletheta(k)=theta_site
               !          dwsdzl=dwsdzu
               !          profile1(k)=((site_ws_data(kk,i_time,jj+1)-site_ws_data(kk,i_time,jj)) &
               !                   /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))) & 
               !                   *(zm(k)-site_z_data(kk,i_time,jj))+site_ws_data(kk,i_time,jj)
               !       elseif(site_nz_data(kk,i_time) .eq. 2)then
               !          theta_site=(270.-(site_wd_data(kk,i_time,jj)+ &
               !             dwinddir*(zm(k)-site_z_data(kk,i_time,jj))))*pi/180.
               !          dwsdzu=((site_dwsdz(jj+1)-site_dwsdz(jj)) &
               !             /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))) &
               !             *(zm(k)-site_z_data(kk,i_time,jj))+site_dwsdz(jj)
               !          mag_site=mag_site+0.5*(dwsdzu+dwsdzl)*(zm(k)-zm(k-1))
               !          profile2(k)=mag_site
               !          profiletheta(k)=theta_site
               !          dwsdzl=dwsdzu
               !          profile1(k)=((site_ws_data(kk,i_time,jj+1)-site_ws_data(kk,i_time,jj)) &
               !                   /(site_z_data(kk,i_time,jj+1)-site_z_data(kk,i_time,jj))) & 
               !                   *(zm(k)-site_z_data(kk,i_time,jj))+site_ws_data(kk,i_time,jj)
               !       else
               !          profile1(k)=(site_ws_data(kk,i_time,1)/(log((site_z_data(kk,i_time,1)+site_pp(kk,i_time)) &
               !             /site_pp(kk,i_time))))*log((zm(k)+site_pp(kk,i_time))/site_pp(kk,i_time))
               !          profile2(k)=profile1(k)
               !          profiletheta(k)=theta_site
               !       endif
               !    endif
               ! enddo
               ! sumprofile1=0.
               ! sumprofile2=0.
               ! do k=2,nz
               !    sumprofile1=sumprofile1+profile1(k) ! *dz_array(k)
               !    sumprofile2=sumprofile2+profile2(k) ! *dz_array(k)
               ! enddo
               ! profileratio=sumprofile1/sumprofile2
               ! do k=2,nz
               !    u_prof(kk,k)= profileratio*profile2(k)*cos(profiletheta(k))
               !    v_prof(kk,k)= profileratio*profile2(k)*sin(profiletheta(k))
               ! enddo
            endif ! erp 2/6/2003 end data entry
         enddo   lp003       ! num_sites kk
! end MAN 04/05/2007
         ! ---------------------------------------------------------------------------------
         ! find average distance of each measuring station relative to all other
         ! measuring stations
         
         if(num_sites .eq. 1)then
            do k=2,nz
               uo(:,:,k)=u_prof(1,k)
               vo(:,:,k)=v_prof(1,k)
            enddo
         else ! Barnes Mapping Scheme
            rcsum=0.
            do kk=1,num_sites
               rcval=1000000. ! ignore anything over 1000 kilometers
               do k=1,num_sites
                  xc=site_xcoord(k)-site_xcoord(kk)
                  yc=site_ycoord(k)-site_ycoord(kk)
                  rc=sqrt(xc**2+yc**2)
                  if(rc .lt. rcval .and. k .ne. kk) rcval=rc ! shortest distance
               enddo
               rcsum=rcval+rcsum ! sum of shortest distances
            enddo
            deln=rcsum/real(num_sites)  ! average Radius (i.e. computed data spacing)
            lamda=5.052*(2.*deln/pi)**2 ! weight parameter
            ! numerical convergence parameter
            sgamma = 0.2         ! gamma=.2 max detail, gamma=1 min detail
! MAN 7/7/2005 var dz conversion
! calculate grid cell center locations in meters
            !$omp parallel do private(i)
            do j=1,ny
               do i=1,nx
                  x(i,j)=real(i-1)*dx-.5*dx ! calculating in meters  ddx=(meters/grid)
                  y(i,j)=real(j-1)*dy-.5*dy
               enddo
            enddo
            !$omp end parallel do
! end MAN 7/7/2005

            ! ------------------------------------------------------------------------------
            ! first and second barnes pass done for each cell level in z direction
            
            ! compute weight of each site on point (i,j)
            !$omp parallel do default(shared) private(i,j,kk)
            do j=1,ny				
               do i=1,nx
                  do kk=1,num_sites
                     wm(kk,i,j)=exp(-1/lamda*(site_xcoord(kk)-x(i,j))**2   &
                               -1/lamda*(site_ycoord(kk)-y(i,j))**2)        
                     wms(kk,i,j)=exp(-1/(sgamma*lamda)*(site_xcoord(kk)-x(i,j))**2   &
                                    -1/(sgamma*lamda)*(site_ycoord(kk)-y(i,j))**2)
                  enddo
                  
                  if(sum(wm(:,i,j)) .eq. 0.)then
                     wm(:,i,j)=1e-20
                  endif
               enddo
            enddo
            !$omp end parallel do
lp004:      do k=2,nz
               ! interpolate onto the grid
               ! do first Barnes pass
               !$omp parallel do private(i,sumwu,sumwv,sumw)
               do j=1,ny
                  do i=1,nx
                     sumwu=sum(wm(:,i,j)*u_prof(:,k))
                     sumwv=sum(wm(:,i,j)*v_prof(:,k))
                     sumw=sum(wm(:,i,j))
                     uo(i,j,k)=sumwu/sumw
                     vo(i,j,k)=sumwv/sumw
                  enddo ! i=nx
               enddo ! j=ny
               !$omp end parallel do
! before doing the 2nd pass for the Barnes Method     
! use a 4-point bilinear interpolation       
! scheme to get estimated values at measured point (+)  
! using the 1st pass calculated data at grid points (*) 
!                       
!     *       *            1        2     
!                    
!          +                    +        !definition of points
!                       
!     *       *            3        4
!     |    |  |      
!     | dxx|  |  
!     |       |  !definition of measurements
!     |  ddx  |
!        
               do kk=1,num_sites
                  if(site_xcoord(kk) .gt. 0. .and. site_xcoord(kk) .lt. real(nx-1)*dx .and. &
                        site_ycoord(kk) .gt. 0. .and. site_ycoord(kk) .lt. real(ny-1)*dy)then
                     do j=1,ny
                        !find closest grid location on lower side of site
                        
!!!!!!!!!               if(y(1,j).lt.site_ycoord(kk)) jwork=j
                        if(y(1,j).lt.site_ycoord(kk))then 
							jwork=j
						endif	
!!!!!!!!!                 
                     enddo
                     do i=1,nx
                        !find closest grid location on left side of site
!!!!!!!!!               if(x(i,1).lt.site_xcoord(kk)) iwork=i
                        if(x(i,1).lt.site_xcoord(kk))then 
							iwork=i
                        endif
!!!!!!!!!                                             
                     enddo !i=nx
                     ! distance to site point from lower and left sides
                     dxx=site_xcoord(kk)-x(iwork,jwork)
                     dyy=site_ycoord(kk)-y(iwork,jwork)
! MAN 7/7/2005 fixed interpolation of velocities and var dz conversion  
                     ! upper u interpolated velocity
                     u12 = (1-dxx/dx)*uo(iwork,jwork+1,k)+(dxx/dx)*uo(iwork+1,jwork+1,k)
                     ! lower u interplotaed velocity
                     u34 = (1-dxx/dx)*uo(iwork,jwork,k)+(dxx/dx)*uo(iwork+1,jwork,k)
                     ! total interpolated u velocity
                     uoint(kk)=(dyy/dy)*u12+(1-dyy/dy)*u34
                     
                     ! upper v interpolated velocity
                     v12 = (1-dxx/dx)*vo(iwork,jwork+1,k)+(dxx/dx)*vo(iwork+1,jwork+1,k)
                     ! lower v interplotaed velocity
                     v34 = (1-dxx/dx)*vo(iwork,jwork,k)+(dxx/dx)*vo(iwork+1,jwork,k)
                     ! total interpolated u velocity
                     voint(kk)=(dyy/dy)*v12+(1-dyy/dy)*v34         
                  else
                     uoint(kk)=u_prof(kk,k)
                     voint(kk)=v_prof(kk,k)
                  endif
! end MAN 7/7/2005 
               enddo ! kk=num_sites
! end bilinear interpolation section
! Begin 2nd Barnes pass
               !$omp parallel do private(i,sumwu,sumwv,sumw)
               do j=1,ny
                  do i=1,nx
                     sumwu=sum(wms(:,i,j)*(u_prof(:,k)-uoint(:)))
                     sumwv=sum(wms(:,i,j)*(v_prof(:,k)-voint(:)))
                     sumw=sum(wms(:,i,j))
                     if(sumw .ne. 0.)then
                        uo(i,j,k)=uo(i,j,k) + sumwu/sumw
                        vo(i,j,k)=vo(i,j,k) + sumwv/sumw
                     endif
                  enddo
               enddo
               !$omp end parallel do
            enddo   lp004       ! k=1:nz
         endif
         ! if this is the last time loop then deallocate this subroutines variables
         ! deallocate arrays
         ! deallocate(u_prof,v_prof)
         deallocate(x,y)
         ! deallocate(profile1,profile2,profiletheta)
         ! deallocate(site_dwsdz)
         ! initialize vector fields TMB 7/11/03
         if (i_time .eq. 1) then
            open(unit=200,file="QU_velocity_o.bin",status="unknown",form='unformatted')
         endif
         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
         do k=1,nz-1
            do j=1,ny-1
               do i=1,nx-1
                  u(i,j,k)=(0.5*(uo(i,j,k)+uo(i+1,j,k)))
                  v(i,j,k)=(0.5*(vo(i,j,k)+vo(i,j+1,k)))
                  if(uo(i,j,k) .ne. uo(i,j,k) .or. vo(i,j,k) .ne. vo(i,j,k))then
                     print*,'Interpolation scheme NaN at ',i,j,k
                  endif
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         write(200)(((u(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1), &
            (((v(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
         if(i_time .eq. num_time_steps)then
            close(200)
         endif
         u(1:nx,1:ny,1:nz)=0.
         v(1:nx,1:ny,1:nz)=0.
         w(1:nx,1:ny,1:nz)=0.
! erp initialize upwind array erp 6/18/04
         uo_bl(1:nz)=uo(1,1,1:nz)
         vo_bl(1:nz)=vo(1,1,1:nz)
         print*,'finished interpolating winds'
         return
      end
