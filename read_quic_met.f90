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
      subroutine read_quic_met

         use datamodule ! make data from module "datamodule" visible
! operations  done a priori to speed up the code AAG & IS  07/03/06
         implicit none
         integer kk,tt,coordFlag, utm_zone
         character*128 f_name    ! filename for a data profile at one of the sites.
         real*8 convergence, lon, lat, utm_x, utm_y, theta
         
         ! blayer_flag has already been read in
         read(36,*)num_sites			! total number of sites
         read(36,*)num_vert_points
         ! allocate arrays
         allocate(site_blayer_flag(num_sites,num_time_steps))
         allocate(site_pp(num_sites,num_time_steps))
         allocate(site_H(num_sites,num_time_steps))
         allocate(site_rL(num_sites,num_time_steps))
         allocate(site_ac(num_sites,num_time_steps))
         allocate(site_xcoord(num_sites))
         allocate(site_ycoord(num_sites))
!         allocate(site_zcoord(num_sites))
!         allocate(t_site(num_sites,num_time_steps))
!         allocate(dir_site(num_sites,num_time_steps))
!         allocate(vel_site(num_sites,num_time_steps))
         allocate(u_prof(num_sites,nz),v_prof(num_sites,nz))
         allocate(uoint(num_sites),voint(num_sites))
         allocate(wm(num_sites,nx,ny),wms(num_sites,nx,ny))
         allocate(site_nz_data(num_sites,num_time_steps))
         ! MAN 4/5/2007 number of data points to allocate in data points profiles
         allocate(site_z_data(num_sites,num_time_steps,num_vert_points))
         allocate(site_ws_data(num_sites,num_time_steps,num_vert_points))
         allocate(site_wd_data(num_sites,num_time_steps,num_vert_points))
         
!         allocate(site_ustar(num_sites),site_ac(num_sites))
!         allocate(site_lc(num_sites),site_d(num_sites))


         ! initialize velocity profiles at each site
         u_prof(1:num_sites,1:nz) = 0
         v_prof(1:num_sites,1:nz) = 0
         
         theta=real(domain_rotation*pi/180.,8)

         ! initialize the vectors
         site_xcoord(1:num_sites) = 0
         site_ycoord(1:num_sites) = 0
         site_nz_data(:,:)=1
!         site_zcoord(1:num_sites) = 0
!         t_site(1:num_sites,1:num_time_steps) = -99
!         dir_site(1:num_sites,1:num_time_steps) = 0
!         vel_site(1:num_sites,1:num_time_steps) = 0

lp001:   do kk=1,num_sites
            read(36,*)						! The site name
            read(36,*)						! Description line
            read(36,*)f_name				! File name of the individual file
            ! Reading each profile/sensor file
            open(unit=52,file=f_name,status='old')
            read(52,*)						! The site name
            read(52,*)coordFlag
            !print*,"coordFlag is",coordFlag
            read(52,*)site_xcoord(kk)		! x coordinate of site location (meters)
            !print*,"site_xcoord is",site_xcoord(kk)
            read(52,*)site_ycoord(kk)		! y coordinate of site location (meters)
            !print*,"site_ycoord is",site_ycoord(kk)
            select case(coordFlag)
               case(2) ! UTM coordinates
                  read(52,*) utm_x      ! UTMX position of site location (meters)
                    !print*,"utm_x is",utm_x
                  read(52,*) utm_y      ! UTMY position of site location (meters)
                    !print*,"utm_y is",utm_y
                  read(52,*) utm_zone   ! UTM zone number of site location
                    !print*,"utm_zone is",utm_zone
                  read(52,*) ! UTM zone letter of site location (1=A, 2=B, etc.)
               case(3)
                  read(52,*) lat ! latitude of site location (dd)
                    !print*,"lat is",lat
                  read(52,*) lon ! longitude of site location (dd)
                    !print*,"lon is",lon
            end select
            convergence = 0.0
            if( utmx .ne. 0. .and. utmy .ne. 0.) then
            	select case(coordFlag)
            		case(1) ! QUIC coordinates
            			utm_x = real(site_xcoord(kk),8)*dcos(theta) + real(site_ycoord(kk),8)*dsin(theta) + real(utmx,8)
            			utm_y = -real(site_xcoord(kk),8)*dsin(theta) + real(site_ycoord(kk),8)*dcos(theta) + real(utmy,8)
            			call utm_geo(lon, lat, utm_x, utm_y, utmzone, 1)
                  case(2) ! UTM coordinates
                     call utm_geo(lon, lat, utm_x, utm_y, utm_zone, 1)
               end select
               call getConvergence(lon, lat, utmzone, convergence)
            endif
            
            
            do tt=1,num_time_steps
               read(52,*) !time stamp
               read(52,*)site_blayer_flag(kk,tt) ! boundary layer flag for each site (1=log,2=exp,3=canopy,4=data)e
                !print*,"site_blayer_flag is",site_blayer_flag(kk,tt)
	            read(52,*)site_pp(kk,tt)			! if blayer = 2 site_pp = exp else site_pp = zo
                !print*,"site_pp is",site_pp(kk,tt)
               select case(site_blayer_flag(kk,tt))
                  case(1)! logarithmic profile
                     read(52,*)site_rL(kk,tt)			! reciprocal Monin-Obukhov length
                    !print*,"log profile has been selected, site_rL is",site_rL(kk,tt)
                  case(3)! urban canopy
                     read(52,*)site_rL(kk,tt)			! reciprocal Monin-Obukhov length
                     read(52,*)site_H(kk,tt)			   ! canopy height
                     read(52,*)site_ac(kk,tt)			! atenuation coefficient
                        !print*,"urban canopy has been selected, site_rL (recip M-O length) is",site_rL(kk,tt)
                        !print*,"site_H (canopy height) is",site_H(kk,tt)
                        !print*,"site_ac (att. coeff.) is",site_ac(kk,tt)
                  case(4)! data points
                     read(52,*)site_nz_data(kk,tt)		! number of data points in vertical wind profile
                        !print*,"profile from data has been selected"
               end select
               read(52,*)! skip line			!"height  direction   magnitude"  Label
               do ii=1,site_nz_data(kk,tt)
                  read(52,*)site_z_data(kk,tt,ii),site_ws_data(kk,tt,ii),site_wd_data(kk,tt,ii)
                    !print*,"size_z_data is",site_z_data(kk,tt,ii)
                    !print*,"size_ws_data is",site_ws_data(kk,tt,ii)
                    !print*,"size_wd_data is",site_wd_data(kk,tt,ii)
! MAN 02/05/2007 Domain Rotation
                  ! Scot - Modify here... 
                  site_wd_data(kk,tt,ii)=site_wd_data(kk,tt,ii) - domain_rotation - real(convergence)
                  if(site_wd_data(kk,tt,ii).lt. 0.)then
                       site_wd_data(kk,tt,ii)=site_wd_data(kk,tt,ii)+360.
                  elseif(site_wd_data(kk,tt,ii).ge. 360.)then
                       site_wd_data(kk,tt,ii)=site_wd_data(kk,tt,ii)-360.
                  endif
               enddo
            enddo
            close(52)	! closing the profile/sensor data file
         enddo   lp001       !kk = num_sites
! MAN 10/10/2007 if there is only one measurement make sure it is in the domain
         if(num_sites .eq. 1)then
            site_xcoord(1)=0.5*real(nx-1)*dx
            site_ycoord(1)=0.5*real(ny-1)*dy
         endif
! end MAN 10/10/2007
         return
      end


subroutine getConvergence(lon, lat, zone, convergence)
    implicit none
    integer zone
    double precision, parameter :: PI = 3.141592653589793d0
    double precision, parameter :: degrad=PI/180.d0, raddeg=180.d0/PI
    double precision lat, lon, tempLon, convergence
    
    tempLon = (zone * 6.0) - 183.0 - lon
    
    convergence = datan(dtan(tempLon * degrad) * dsin(lat * degrad)) * raddeg
end subroutine
