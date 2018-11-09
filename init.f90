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
      subroutine init
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subroutine to initialize main array variables          
! reads in data from the input.dat file
! init.f90 subroutine is called by main.f90
! init.f90 calls zwake.f90 function       
! ERP 2001     
! Variable information:
! uo,vo,wo - are the initial velocity that are prescribed prior to mass 
!            conservation.
! u,v,w - final velocity field
! xfo,yfo - denote the x and y locations that represent the 
!           left side center of the building
! nx,ny, and nz are the number of cells in the x,y and z directions                 
! theta - wind angle in standard meteorological format. 0 degrees is    
!    out of the north, 90 degrees is a wind out of the east, etc     
! inumbuild - number of buildings in building array      
! 
! Building type designations
!  bldtype = 1 Regular building 
!  bldtype = 2 Vegetation, requires an attenuation coefficient (see Cionco, 1965)
!  bldtype = 3 Bridge
!  bldtype = 5 Relief
! 
! 
! * note that the velocity field is initialized at the end of the    
!   subroutine.   
! erp 6/2/03 modifications to allow for data input velocity profiles. For
!     example wind direction can be varied as a function of height
! erp 6/2/03 modifications to allow for variable grid resolutions to
!     be entered via the ddx variable. ddx is specified in meters
! erp 6/5/03 added empirical parameterization subdomain. A subdomain box
!     may be defined in which the empirical parameterizations are 
!     applied. Outside of this domain conservation of mass is applied only
!     to an incoming boundary layer profile.
! erp 7/25/03   This version of qwicurb has the added array zfo(ibuild) 
!     which allows buildings of different sizes to be stacked on one another. 
! erp 8/14/03 This version of the code incorporates Tom Booth's upwind urban
!     boundary layer work. As such it calls the function zwake.
! erp 2/11/04 This version has been modified for Nilesh's upwind vortex
! erp 6/08/04 Modifications to canyon parameterization 
! erp 10/05/04  This version removes the meteorological input information
!     and puts it in the subroutine met_init.f90 to allow for multi-run
!     capability. 
! erp 6/30/05 This version adds variable dx,dy and dz capability based on the
!     work of Tau.
!           
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         use datamodule
         implicit none
! erp 2/17/04 NLB modifications

!  real Az,Bz
!  real logvel,Az,Bz
! erp 4/20/04 end changes logvel moved to datamodule

!  real zw, ustar, lc, d, vk, ac !TMB 7/10/03 canopy variables
!         real zwake
         integer count,numnodes ! MAN 05-24-2007       
         integer landuse_params(3)
         integer, allocatable:: bin_int_read(:)
         real, allocatable:: bin_real_read(:)
         pi=4.*atan(1.0)
   
! open and read input file input.dat
         open(unit=35,file='QU_simparams.inp',status='old')
         open(unit=36,file='QU_metparams.inp',status='old')
         open(unit=37,file='QU_buildings.inp',status='old')
         open(unit=47,file='QU_fileoptions.inp',status='old')
         open(unit=46,file='QP_buildout.inp',status='unknown')
         
         !open(unit=48,file='QP_buildorder.inp',status='unknown')

! read file read writing option data from QU_options.dat
         read(47,*) ! QUIC version header line
         read(47,*)format_flag   !output format flag 1=ascii,2=binary,3=both
         read(47,*)uofield_flag  !write out uofield.dat 1=yes, 0=no
         read(47,*)uosensor_flag !write out flag for sensor velocities 1=yes, 0=no
         read(47,*)staggered_flag !write out flag for staggered velocities 1=yes, 0=no

! read input data from QU_domain.inp
         read(35,*) ! QUIC version header line
         read(35,*)nx      !nx defined in input file
         read(35,*)ny      !ny defined in input file
         read(35,*)nz      !nz defined in input file
         read(35,*)dx      !dx defined in input file
         read(35,*)dy      !dy defined in input file
         
         dxy=min(dx,dy)
! man 1/14/05 account for difference in grid cell definitions
         nx=nx+1
         ny=ny+1
         nz=nz+2
! end man 1/14/05 account for difference in grid cell definitions
! MAN 07/25/2008 stretched vertical grid
         allocate(z(nz),zm(nz),dz_array(nz))
         read(35,*)stretchgridflag   !Stretched grid flag (0= dz constant with z)
         z(:)=0.
         zm(:)=0.
         dz_array(:)=0.
         select case(stretchgridflag)
            case(0) !uniform
               read(35,*)dz      !dz defined in input file
               dz_array(:)=dz
            case(1) !custom
               read(35,*)
               do k=2,nz-1
                  read(35,*)dz_array(k)
               enddo
            case(2,3,4) !parabolic or exponential
               read(35,*)
               read(35,*)
               read(35,*)
               do k=2,nz-1
                  read(35,*)dz_array(k)
               enddo
         endselect
         dz_array(1)=dz_array(2)
         dz_array(nz)=dz_array(nz-1)
         zm(1)=-0.5*dz_array(1)
         z(1)=0.0
         do k=2,nz
            z(k)=z(k-1)+dz_array(k)
            zm(k)=z(k)-0.5*dz_array(k)
         enddo
         
         dz=minval(dz_array)
! erp 6/30.05 coefficients for sor solver
         A=dx**2/dy**2
         B=eta*(dx**2/dz**2)
! MAN 09/02/2008 time inforation
!         read(35,*)start_time      !decimal start time
!         read(35,*)time_incr         !time increment
         read(35,*)num_time_steps    !total time increments
         allocate(time(num_time_steps))
         ! read(35,*) ! day of year (not necessary with unix time stamp)
         read(35,*)utc_offset ! UTC conversion
         read(35,*) ! header line
         do i_time=1,num_time_steps
            read(35,*)time(i_time)
         enddo
! building parameterization flags
         read(35,*)roofflag   ! rooftop recirc flag
         read(35,*)upwindflag ! upwind cavity flag
         read(35,*)streetcanyonflag    ! street canyon initialization method, added PKK 05/12/03
         read(35,*)intersectionflag    ! MAN 7/11/2006
         read(35,*)wakeflag ! MAN 06/29/2007 added wake flag to QU_simparams.inp
         read(35,*)sidewall_flag ! MAN 04/24/2013 added sidewall flag to QU_simparams.inp
! MAN 7/10/2006 convergence criteria
         read(35,*)itermax ! max number of iterations
         read(35,*)residual_reduction     ! MAN 09/26/2006 added residual reduction to input
! AAG 08/25/2006 turbulent diffusion parameters
         read(35,*)diffusion_flag   ! turns on diffusion
         read(35,*)diffstep      ! diffusion iterations
! MAN 02/05/2007 Geo-referencing parameters
         !!! Scot - look here
         read(35,*)domain_rotation
         read(35,*)utmx
         read(35,*)utmy
         read(35,*)utmzone
         read(35,*)utmzoneletter
! MAN 05-24-2007 QUIC-CFD flag turns off the SOR solver        
         read(35,*)qcfd_flag
         read(35,*)damageflag
         if(qcfd_flag .gt. 0)then
            itermax=0
            diffusion_flag=0
         endif
! input from QU_buildings.inp
         read(37,*) ! QUIC version header line
         read(37,*)zo   ! MAN 8-19-2005 Updated input output file structures
         read(37,*)inumbuild  ! number of buildings
         read(37,*)inumpolygon ! number of polygon nodes

         allocate(uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz))
         allocate(ufarwake(nx,ny),vfarwake(nx,ny))
! need to move to met_init.f90
         allocate(uo_bl(nz),vo_bl(nz)) ! erp 6/08/04

         allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz))
         allocate(p1(nx-1,ny-1,nz-1),p2(nx-1,ny-1,nz-1),r(nx-1,ny-1,nz-1))
         ! MAN 8/30/2005 stacked building fix
         allocate(zfo_actual(inumbuild))
         allocate(bldgeometry(inumbuild),numpolygons(inumbuild)) ! MAN 08/10/2010
         allocate(bldstartidx(inumbuild),bldstopidx(inumbuild),bldroof(inumbuild)) ! MAN 08/10/2010
         allocate(bldx(inumpolygon),bldy(inumpolygon),bldwall(inumbuild)) ! MAN 08/10/2010
         allocate(LrNode(inumpolygon),LrFace(inumpolygon),FaceRelWindDir(inumpolygon))
         allocate(bldcx(inumbuild),bldcy(inumbuild)) ! MAN 08/10/2010
         allocate(xfo(inumbuild),yfo(inumbuild),zfo(inumbuild))   ! erp 7/25/03
         allocate(gamma(inumbuild)) ! erp 7/26/03
         allocate(aa(inumbuild),bb(inumbuild))              ! erp 1/31/2003
         allocate(icellflag(nx-1,ny-1,nz-1),ibldflag(nx-1,ny-1,nz-1))
         allocate(e(nx-1,ny-1,nz-1),f(nx-1,ny-1,nz-1),g(nx-1,ny-1,nz-1))
         allocate(h(nx-1,ny-1,nz-1),m(nx-1,ny-1,nz-1),n(nx-1,ny-1,nz-1))
         allocate(denom(nx-1,ny-1,nz-1))
         allocate(Ht(inumbuild),Wti(inumbuild),Lti(inumbuild))
         allocate(bldnum(inumbuild),bldtype(inumbuild),group_id(inumbuild))
         allocate(atten(inumbuild)) ! erp 1/3/2006
         allocate(rooftop_flag(inumbuild)) ! AAG 09/13/06 
         rooftop_flag(:)=0 ! AAG 09/13/06  intialized rooftop flag to 1
         allocate(bld_damage(inumbuild)) ! MAN 08/17/2009
         bld_damage(:)=0 ! MAN 08/17/2009
         bldroof(:)=0
         Wti(:)=0.
         Lti(:)=0.
         if(diffusion_flag .gt. 0)allocate(Fxd(nx,ny,nz),Fyd(nx,ny,nz),Fzd(nx,ny,nz),visc(nx,ny,nz))
! Read in the building number, type, height, width, length ,xfo,yfo,zfo,gamma and atten
! atten - attenuation coefficient for vegetation
! erp 1/31/2003
! note that for now if the building is cylindrical, enter Lti = 0.
         count=0
         do i=1,inumbuild
            read(37,*)
            bldnum(i)=i
            read(37,*)group_id(i)
            read(37,*)bldgeometry(i)
            read(37,*)bldtype(i)
            if(bldtype(i) .eq. 2)read(37,*)atten(i)
            if(bldgeometry(i) .eq. 4 .or. bldgeometry(i) .eq. 5)then
               read(37,*)bldwall(i)
               read(37,*)bldroof(i)
            endif
            read(37,*)Ht(i)
            read(37,*)zfo(i)
            read(37,*)bldcx(i)
            read(37,*)bldcy(i)
            if(bldgeometry(i) .eq. 6)then
               read(37,*)numpolygons(i)
               do j=1,numpolygons(i)
                  read(37,*)
                  read(37,*)numnodes
                  read(37,*)
                  if(j .eq. 1)bldstartidx(i)=count+1
                  do k=1,numnodes
                     count=count+1
                     read(37,*)bldx(count),bldy(count)
                  enddo
                  read(37,*)
               enddo
               bldstopidx(i)=count
            else
               read(37,*)xfo(i)
               read(37,*)yfo(i)
               read(37,*)Lti(i)
               read(37,*)Wti(i)
               read(37,*)gamma(i)
               if(bldgeometry(i) .eq. 2 .or. bldgeometry(i) .eq. 5)then    !if the building is a cylinder/ellipse
                  bb(i)=Wti(i)/2.         !set minor axis to input Width
                  aa(i)=Lti(i)/2.         !set major axis to input Lenth
               endif
               if(bldgeometry(i) .eq. 3)then      ! if the building is a Pentagon
                  bb(i)=Wti(i)/2.                 ! Radius Pentagon is inscribed in
                  xfo(i)=xfo(i)-bb(i)
               endif
            endif
            read(37,*)
            Ht(i)=Ht(i)+zfo(i)
            gamma(i)=gamma(i)*pi/180.  !erp 7/25/03
         enddo

   
! erp 1/3/2006 check to see if building is actually vegetation
         inumcanopy = 0
         inumbuildneg = 0
         inumgarage = 0
         do i=1,inumbuild
            if(bldtype(i) .eq. 2)then
               inumcanopy=inumcanopy+1       ! total number of vegative canopies
            endif
            if(bldtype(i) .eq. 3)then
               inumgarage=inumgarage+1       ! total number of garage buildings
            endif
            if(bldtype(i) .eq. 0)then
               inumbuildneg=inumbuildneg+1         ! total number of negative buildings
            endif
         enddo
         
! if vegetation exist rename variables to be consistent with plantinit.f90 rountine  
         if(inumcanopy.gt.0)then
            allocate(canopy_ktop(nx-1,ny-1),canopy_top(nx-1,ny-1),canopy_atten(nx-1,ny-1,nz-1))
            allocate(canopy_zo(nx-1,ny-1),canopy_ustar(nx-1,ny-1),canopy_d(nx-1,ny-1))
         endif
! end 1/3/2006
     
! MAN 8/25/2009 read in land use data
         open(unit=66,file='QU_landuse.inp',form='unformatted',status='old')
         read(66)landuse_params
         landuse_flag=landuse_params(1)
         landuse_veg_flag=landuse_params(2)
         landuse_urb_flag=landuse_params(3)
         if(landuse_flag .eq. 1)then
            allocate(landuse(nx-1,ny-1),bin_int_read((nx-1)*(ny-1)))
            read(66)bin_int_read
            count = 1
            do i=1,nx-1
               do j=1,ny-1
                  landuse(i,j)=bin_int_read(count)
                  count=count+1
               enddo
            enddo
            deallocate(bin_int_read)
            if(landuse_veg_flag .eq. 1 .or. landuse_urb_flag .eq. 1)then
               lu_canopy_flag=1
               allocate(landuse_height(nx-1,ny-1),landuse_atten(nx-1,ny-1),bin_real_read((nx-1)*(ny-1)))
               read(66)bin_real_read
               count = 1
               do i=1,nx-1
                  do j=1,ny-1
                     landuse_height(i,j)=bin_real_read(count)
                     count=count+1
                  enddo
               enddo
               read(66)bin_real_read
               count = 1
               do i=1,nx-1
                  do j=1,ny-1
                     landuse_atten(i,j)=bin_real_read(count)
                     count=count+1
                  enddo
               enddo
               deallocate(bin_real_read)
               if(inumcanopy .eq. 0)then
                  allocate(canopy_ktop(nx-1,ny-1),canopy_top(nx-1,ny-1),canopy_atten(nx-1,ny-1,nz-1))
                  allocate(canopy_zo(nx-1,ny-1),canopy_ustar(nx-1,ny-1),canopy_d(nx-1,ny-1))
               endif
            endif
         else
            landuse_veg_flag=0
            landuse_urb_flag=0
            lu_canopy_flag=0
         endif
         close(66)
         
         if(inumcanopy .gt. 0 .or. lu_canopy_flag .gt. 0)then
            canopy_flag=1
         else
            canopy_flag=0
         endif
! end 8/25/2009

! erp 10/2004 move to met_init.f
! erp 1/31/2003
! convert from real world units to grid units
! 
!  zref=zref/ddx
!  uin=uin/ddx
!  !modify zo in grid terms
!  if(blayer_flag.eq.2.or.blayer_flag.eq.4)pp=pp/ddx  
! end erp 10/2004 move to met_init.f

! calculate domain Length width and height erp 1/30/2003
         Lx=(nx-1)*dx
         Ly=(ny-1)*dy
         ! MAN 07/25/2008 stretched vertical grid
         Lz=z(nz-1)
! calculate domain Length width and height erp 1/30/2003


!       initialize p1 and p2 - they are Lagrange multipliers
! TMB using vectors instead of DO loops 
 
         p1(1:nx-1,1:ny-1,1:nz-1)=0.1
         p2(1:nx-1,1:ny-1,1:nz-1)=0.
         r(1:nx-1,1:ny-1,1:nz-1)=0.


! erp 1/17/03 initilize all arrays at zero to start   
         uo(1:nx,1:ny,1:nz)=0.
         vo(1:nx,1:ny,1:nz)=0.
         wo(1:nx,1:ny,1:nz)=0.
   
         close(35)
         close(37)
         close(47)
         
         return
      end
