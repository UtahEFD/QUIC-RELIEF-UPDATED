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
      subroutine building_parameterizations
!************************************************************************
! bcsetup - boundary condition setup program for qwic urb      
!    - called by main.f90                 
!    - calls defbuild.f90, upwind.f90, street_intersec.f90  
! The empirical parameterizations are applied in the following order:
!  1. uninterupted boundary layer
!  2. upwind vortex cavity (calls upwind.f90)
!  3. rooftop recirculation
!  4. wakes
! 
!                 
! THIS VERSION OF QWIC-URB CONTAINS MULTIPLE BUILDING CAPABILITY  
! all velocities are cell face quantities          
! icellflag,celltype,Lagrange multipliers (p1,p2) and their    
! coeficients  (e,f,g,m,n,o,p,q,h,p) are a cell centered quantities. 
!                          
! Lfx is the length of the front eddy on a building in xdirection 
! Lfy is the length of the front eddy on a building in ydirection 
! Lr is the length of the rear vortex cavity behind  a building      
! theta is the angle of mean wind speed (meters/sec)        
! inumbuild is the number of buildings in the city array    
! xfo is the x coord. left front center of bld
! yfo is the y coord. left front center of bld
! zfo is the z coord of building base 
! xlo is the x coord. lower front center of bld
! xuo is the x coord. upper front center of bld
!                          
! REFERENCE FRAME DEFINITIONS:                  
! below  is defined as  lesser  z value (ie z(k) is below z(k+1)  
! above  is defined as  greater z value (ie z(k+1) is above z(k)  
! right  is defined as  greater x value (ie x(i+1) is right of x(i)  
! left   is defined as  lesser  x value (ie x(i) is right of x(i+1)  
! front  is defined as  greater y value (ie y(j+1) is in front of y(j)  
! behind is defined as  lesser  y value (ie y(j) is behind y(j+1) 
! ERP dec/2000                      
! moving  over various wind angle routine from earlier version of code  
! ERP August/2002                   
! Most recent Modification:
! Changing of coordintate system to reflect true stagard grid   
! ERP Sept/Oct 2002
! erp 12/17 changes fixing north/south flows
! mdw 1/8/2003 changes include fix to jstartb
! erp 7/23/03 modifications to allow for building blocks to be stacked
!  this involves modifications to k loops,zfo and zbo
! erp 11/13/03 fixing Lfx bug
! erp 11/13/03 fixing rooftop vortex bug
! erp 1/29/04 added Petra Kastner-Klein's finalized street canyon technique
!  that she idependantly tested and verified
!  an option is included to use the CPB potential flow formulas for the velocity 
!  field initialization in the street canyon
!  the input file was modified, in line 9 the initialization method is chosen
!    streetcanyonflag=1 ==> original Roeckle approach
!    streetcanyonflag=2 ==> CPB approach
!  potential flow formulas are not applied near the lateral edges of the canyon
!  depth of lateral vortex zone equal to canyon heigth, u is set to zero, v=u0(z)
! erp 2/11/04 added Nilesh's upwind vortex 
! erp 03/09/04 added canyon flag check to remove wake for skimming flow?
! NLB 02/10/04 Added Upwind Vortex Parameterizations for Perpendicular and Varying Incident Wind Angles
! NLB 10/11/04 Added Rooftop Parameterizations for Perpendicular and Varying Incident Wind Angles 
! erp 03/02/05 This subroutine now can writeout the building grid locations in
!     both a binary (QU_celltype.bin) and ASCII format (celltype2.dat). If
!     format_flag = 1, ASCII only. If format_flag=2, binary only and if
!     format_flag = 3, both the ASCII and binaries are written out
!
! ERP 6/8/2006 this version includes rooftop fixes for both off angle and normal
!  angle calculations (Suhas Pols implementation of off angle fixes)
! MATR 04/02/18 Adding subroutines call for relief type (skipping wake and rooftop param over relief,
! while computing variables needed for the rest of the code to work)
!
! Cellflag designations 
!  icellflag = 0  building
!  icellflag = 1  fluid cell BL parameterization
!  icellflag = 2  upwind cavity parameterization
!  icellflag = 3  rooftop parameterization
!  icellflag = 4  near wake cavity parameterization
!  icellflag = 5  far wake cavity parameterization
!  icellflag = 6  street canyon parameterization
!  icellflag = 8  vegetation parameterization
!  icellflag = 9  street intersection parameterization
!  icellflag = 10  parking garage parameterization
!
!
!************************************************************************
         use datamodule ! make data from module "datamodule" visible

         implicit none
         real vegvelfrac,avg_atten,num_atten
         integer, allocatable:: bld_num_array(:,:,:)
		
		
		 
!      integer icelltemp

!erp 3/1/2006 Rooftop

         

!erp 3/1/2006 Rooftop
         
         if (i_time .eq.1)then
! allocate arrays which are needed in bcsetup, PKK 10/01/02
! they are not passed to other subroutines and will be deallocated at the end of bcsetup 
   
            allocate(istart(inumbuild),iend(inumbuild))
            allocate(jstart(inumbuild),jend(inumbuild))
            allocate(kstart(inumbuild),kend(inumbuild))
            !allocate(istart_canyon_N(inumbuild),iend_canyon_N(inumbuild)) !MAN7/6/2006
            !allocate(jstart_canyon_E(inumbuild),jend_canyon_E(inumbuild)) !MAN7/6/2006
            !allocate(istart_canyon_S(inumbuild),iend_canyon_S(inumbuild)) !MAN7/6/2006
            !allocate(jstart_canyon_W(inumbuild),jend_canyon_W(inumbuild)) !MAN7/6/2006
            !allocate(kend_canyon_W(inumbuild),kend_canyon_E(inumbuild)) !MAN7/6/2006
            !allocate(kend_canyon_N(inumbuild),kend_canyon_S(inumbuild)) !MAN7/6/2006
!            allocate(f_flag(inumbuild),w_flag(inumbuild),f_flagchk(inumbuild))
            !allocate(c_flag_E(inumbuild),c_flag_W(inumbuild),c_flag_N(inumbuild),c_flag_S(inumbuild))
            allocate(Lf(inumbuild),Lr(inumbuild))
            allocate(Weff(inumbuild),Leff(inumbuild))
            allocate(Wt(inumbuild),Lt(inumbuild))
!            allocate(Lfx(inumbuild),Lfy(inumbuild))
!            allocate(Lfx1(inumbuild),Lfy1(inumbuild))  !NLB 02/10/04
!            allocate(Roofcx(nx,ny,nz))   !NLB 10/11/04 For Rooftop
            allocate(Rscale(inumbuild),Rcx(inumbuild))  !NLB 10/11/04 For Rooftop
            allocate(vo_roof(nx,ny,nz),uo_roof(nx,ny,nz)) !NLB 10/10/05
lp001:      do ibuild=1,inumbuild
               if(bld_damage(ibuild) .eq. 2)cycle
               do k=2,nz-1
                  kstart(ibuild)=k
                  if(zfo_actual(ibuild) .le. zm(k))exit
               enddo
               do k=kstart(ibuild),nz-1
                  kend(ibuild)=k
                  if(Ht(ibuild) .lt. zm(k+1))exit
               enddo
            enddo   lp001
         endif    !end 1st time through if
! added AAG 09/20/06  for multiple time steps, uo_roof was not atllocated  
        
!end veg 

lp002:   do i=1,inumbuild
            Lf(i)=-999.
            Leff(i)=0.
            Weff(i)=0.
            Wt(i)=0.5*Wti(i)
            Lt(i)=0.5*Lti(i)
         enddo    lp002

!erp 1/30/2003
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate buildings
! this now includes a building generation loop that allows
! for multiple buildings
! calculate building spacing s and set wake flags
! erp 3/9/05 Note that defbuild calls pentagon
         if(i_time .eq. 1)then
            print*,'Importing Building Data'
            call defbuild  !erp 1/05/05
            call wallbc
         else
            do k=2,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     if(icellflag(i,j,k) .gt. 1)icellflag(i,j,k)=1
                  enddo
               enddo
            enddo
            do ibuild=1,inumbuild
               if(bldgeometry(ibuild) .eq. 3)call pentagon
            enddo
         endif
         
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate Canopies         
			if(inumcanopy.gt.0. .or. (landuse_flag .eq. 1 .and. &
					(landuse_veg_flag .eq. 1 .or. landuse_urb_flag .eq. 1)))then
				print*,'Applying Vegetation Parameterizations'
				call plantinit
			endif
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Call upwind - routine to generate the upwind cavity on a building
! MATR 04/02/18 Cycling when bd type = 5 (relief) 
        if(upwindflag .gt. 0)then
            print*,'Applying Upwind Rotor Parameterizations'
			do ibuild=1,inumbuild
				if(bld_damage(ibuild) .eq. 2)cycle
				if(bldtype(ibuild) .eq. 2 .or. bldtype(ibuild) .eq. 5)cycle
				select case(bldgeometry(ibuild))
					case(1,4,6)
						if(zfo(ibuild) .eq. 0.)call upwind
				end select
			enddo
        endif
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! wake section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! MATR 04/02/18 Calling ReliefWake when bdtype = 5 (relief)
			if(wakeflag .ne. 0)then
				print*,'Applying Building Wake Parameterizations'
				do ibuild=1,inumbuild
				   ! print*,bldnum(ibuild)
					if(bld_damage(ibuild) .eq. 2)cycle
					if(bldtype(ibuild) .eq. 2)cycle
					if(bldtype(ibuild) .eq. 4)then
						call bridgewake
					elseif(bldtype(ibuild) .eq. 5)then
						call reliefwake
					else
						select case(bldgeometry(ibuild))
							case(1,4)
								call rectanglewake
							case(2,5)
								call cylinderwake
							case(6)
								call polywake
						endselect						
					endif
				enddo
			endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! street canyon section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
         if(streetcanyonflag .ne. 0)then
			print*,'Applying Street Canyon Parameterizations'
			call streetcanyon
         endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! sidewall section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(sidewall_flag .ne. 0)then
            print*,'Applying Sidewall Rotor Parameterizations'
            do ibuild=1,inumbuild
				if(bld_damage(ibuild) .eq. 2)cycle
				if(bldtype(ibuild) .eq. 1 .or. bldtype(ibuild) .eq. 3)then
					select case(bldgeometry(ibuild))
						case(1,4,6)
							call sidewall
					end select
				endif
            enddo
         endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! rooftop or courtyard section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! MATR 04/02/18 Calling reliefrooftop when bdtype = 5 (relief)
         if(roofflag .gt. 0)then
            uo_roof=uo
            vo_roof=vo
            print*,'Applying Rooftop and/or Courtyard Parameterizations'
         endif
			do ibuild=1,inumbuild
				if(bld_damage(ibuild) .eq. 2)cycle
 				if(bldtype(ibuild) .eq. 5)then
 					call reliefrooftop
 				endif
				if(bldtype(ibuild) .eq. 1 .or. bldtype(ibuild) .eq. 3 )then
					select case(bldgeometry(ibuild))
						case(1,2)
							call rooftop
						case(4,5)
							call courtyard
						case(6)
							call rooftop
							if(numpolygons(ibuild) .gt. 1)then
								call courtyard
							endif
					endselect
				endif
			enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! parking garage section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc 
			if(inumgarage .gt. 0)then
			   print*,'Applying Garage Parameterizations' 
				do ibuild=1,inumbuild
					if(bld_damage(ibuild) .eq. 2)cycle
					if(bldtype(ibuild) .eq. 3)call parking_garage
				enddo
			endif
         
! End Building parameterization section

!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! street intersection section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
			if(streetcanyonflag .ne. 0 .and. intersectionflag .eq. 1)then
			    print*,'Applying Blended Region Parameterizations'
				if(inumbuild.gt.0)then
					call street_intersect
					call poisson ! AG 04/06/2007 blends street intersection winds
				endif
			endif
! MAN 4/12/2007 Moved Poisson boundary conditions to new subroutine "wallbc"
!ANU 01/04/2006vegetation parameterization
         if(inumcanopy.gt.0. .or. (landuse_flag .eq. 1 .and. &
               (landuse_veg_flag .eq. 1 .or. landuse_urb_flag .eq. 1)))then
            !$omp parallel do private(i,k,avg_atten,num_atten,vegvelfrac)
            do j=1,ny-1
               do i=1,nx-1
                  if(canopy_top(i,j) .gt. 0.)then
                     do k=2,canopy_ktop(i,j)
                        if(canopy_atten(i,j,k) .gt. 0. .and. icellflag(i,j,k) .ne. 8)then
                           avg_atten = canopy_atten(i,j,k)
                           if(canopy_atten(i,j,k+1) .ne. canopy_atten(i,j,k) &
                                 .or. canopy_atten(i,j,k-1) .ne. canopy_atten(i,j,k))then
                              num_atten=1.
                              if(canopy_atten(i,j,k+1) .gt. 0.)then
                                 avg_atten = avg_atten + canopy_atten(i,j,k+1)
                                 num_atten=num_atten+1.
                              endif
                              if(canopy_atten(i,j,k-1) .gt. 0.)then
                                 avg_atten = avg_atten + canopy_atten(i,j,k-1)
                                 num_atten=num_atten+1.
                              endif
                              avg_atten=avg_atten/num_atten
                           endif
                           vegvelfrac=log(canopy_top(i,j)/canopy_zo(i,j))*&
                                 exp(avg_atten*((zm(k)/canopy_top(i,j))-1.))/&
                                 log(zm(k)/canopy_zo(i,j))
                           if(vegvelfrac .lt. 1. .and. vegvelfrac .ge. 0.)then
                              uo(i,j,k)=uo(i,j,k)*vegvelfrac
                              vo(i,j,k)=vo(i,j,k)*vegvelfrac
                              if(j .lt. ny-1)then
                                 if(canopy_atten(i,j+1,k) .eq. 0.)then
                                    vo(i,j+1,k)=vo(i,j+1,k)*vegvelfrac
                                 endif
                              endif
                              if(i .lt. nx-1)then
                                 if(canopy_atten(i+1,j,k) .eq. 0.)then
                                    uo(i+1,j,k)=uo(i+1,j,k)*vegvelfrac
                                 endif
                              endif
                           endif
                           if(icellflag(i,j,k) .gt. 0)then
                              icellflag(i,j,k)=8
                           endif
                        endif
                     enddo
                  endif
               enddo
            enddo
            !$omp end parallel do
         endif
            
            

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!Celltype Coeficient Section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin defining celltypes using the cellflags based on boundary cells
         !$omp parallel do private(i,j)
         do k=1,nz-1
            do j=1,ny-1
               do i=1,nx-1
                  if(icellflag(i,j,k) .eq. 0)then ! MAN 7/8/2005 Celltype definition change
! all cells that are buildings have a zero velocity within them
                     uo(i,j,k)  =0.
                     uo(i+1,j,k)=0.
                     vo(i,j,k)  =0.
                     vo(i,j+1,k)=0.
                     wo(i,j,k)  =0.
                     wo(i,j,k+1)=0.
                  endif
               enddo
            enddo
         enddo
         !omp end parallel do
! MAN 7/8/2005 Celltype definition change



! Tecplot format output
! commented out for qwicurb GUI
!        open(unit=29,file="celltype.dat",status="unknown")
! man 7/12/2005 changes to celltype
         if(i_time .eq. 1)then
            allocate(bld_num_array(nx-1,ny-1,nz-1))
            bld_num_array(:,:,:)=0
            !$omp parallel do private(i,j)
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     if(ibldflag(i,j,k) .gt. 0)then
                        bld_num_array(i,j,k)=bldnum(ibldflag(i,j,k))
                     endif
                  enddo
               enddo
            enddo
            !$omp end parallel do
         endif
         if(format_flag.eq.1 .or. format_flag.eq.3)then
            if(i_time.eq.1)then
               open(unit=33,file="QU_celltype.dat",status="unknown")
               open(unit=67,file="QU_buildflag.dat",status="unknown")
               do k=1,nz-1
                  do j=1,ny-1
                     do i=1,nx-1
                        write(67,71)(real(i)-.5)*dx,(real(j)-.5)*dy,zm(k),bld_num_array(i,j,k)
                     enddo
                  enddo
               enddo
            endif
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     write(33,71)(real(i)-.5)*dx,(real(j)-.5)*dy,zm(k),icellflag(i,j,k)
                  enddo
               enddo
            enddo
         endif
!erp 3/02/2005 lines added to write out unformatted for binary read into Matlab
         if(format_flag.eq.2 .or. format_flag.eq.3)then  !erp 3/2/05
            if(i_time.eq.1)then   
                  open(unit=39,file="QU_celltype.bin",form='unformatted',status="unknown")
                  open(unit=68,file="QU_buildflag.bin",form='unformatted',status="unknown")
                  write(68)(((bld_num_array(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
            endif
            write(39)(((icellflag(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
         endif !end format_flag if erp 3/03/05
         if(i_time .eq. 1)then
            deallocate(bld_num_array)
            if(format_flag.eq.1 .or. format_flag.eq.3)close(67)
            if(format_flag.eq.2 .or. format_flag.eq.3)close(68)
         endif


!  close(29)
!  close(33)
!  close(46)

 71      format(3(1x,f8.3),i5)
!  erp multirun statements

! new print out removed from upwind and now located here because of the
! removal of the upwind cavity under certain conidtions
! erp 12/06/04
         if(i_time .eq. 1)then 
            write(46,*)inumbuild,' ! total number of buildings'
            write(46,*)inumcanopy+inumbuildneg, ' ! total number of vegitative canopies'
            write(46,*)inumpolygon,' ! total number of polygon nodes'
         endif

         do ibuild=1,inumbuild
            if(bldtype(ibuild) .eq. 1 .or. bldtype(ibuild) .eq. 3 &
                  .or. bldtype(ibuild) .eq. 4)then
! write to buildout.dat erp 1/15/2004
               write(46,*)bldnum(ibuild),' !number'
               write(46,*)bldgeometry(ibuild),' !geometry'
               write(46,*)bldtype(ibuild),' !type'
               write(46,*)Ht(ibuild)-zfo_actual(ibuild),' !height'
               write(46,*)zfo_actual(ibuild),' !zfo'
               write(46,*)bld_damage(ibuild),' !damage'
               if(bldgeometry(ibuild) .eq. 3)then
                  write(46,*)Wti(ibuild),' !Weff'
                  write(46,*)Lti(ibuild),' !Leff'
               else
                  write(46,*)Weff(ibuild),' !Weff'
                  write(46,*)Leff(ibuild),' !Leff'
               endif
               write(46,*)Lf(ibuild),' !Lf'
               write(46,*)Lr(ibuild),' !Lr'
               if(bldgeometry(ibuild) .eq. 6)then
                  write(46,*)bldcx(ibuild),' !xc'
                  write(46,*)bldcy(ibuild),' !yc'
                  write(46,*)bldstartidx(ibuild),' !start'
                  write(46,*)bldstopidx(ibuild),' !stop'
                  write(46,*)numpolygons(ibuild),' !num polygons'
                  do i=bldstartidx(ibuild),bldstopidx(ibuild)
                     write(46,*)bldx(i),bldy(i)
                  enddo
               else
                  write(46,*)xfo(ibuild),' !xfo'
                  write(46,*)yfo(ibuild),' !yfo'
                  write(46,*)Lti(ibuild),' !length'
                  write(46,*)Wti(ibuild),' !width'
                  write(46,*)gamma(ibuild)*180/pi,' !gamma'
                  if(bldgeometry(ibuild) .gt. 3)then
                     write(46,*)bldwall(ibuild),' !wall'
                     write(46,*)bldroof(ibuild),' !roof'
                  endif
               endif
            endif
         enddo
         !$omp parallel do private(i,j)
         do k=2,nz
            do j=1,ny
               do i=1,nx
                  if((uo(i,j,k) .ne. uo(i,j,k)) .or. (vo(i,j,k) .ne. vo(i,j,k)) .or. (wo(i,j,k) .ne. wo(i,j,k)))then
                     print*,'NaN found at ',i,j,k,' Cellflag = ',icellflag(i,j,k)                     
                     if(uo(i,j,k) .ne. uo(i,j,k))then
                        uo(i,j,k)=0.
                     endif
                     if(vo(i,j,k) .ne. vo(i,j,k))then
                        vo(i,j,k)=0.
                     endif
                     if(wo(i,j,k) .ne. wo(i,j,k))then
                        wo(i,j,k)=0.
                     endif
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do
!end change
         if(i_time.eq.num_time_steps)then

            deallocate(istart,iend)             !twh - added this line 01/08/03
            deallocate(jstart,jend)             !twh - added this line 01/08/03
            deallocate(kstart,kend)             !twh - added this line 01/08/03
            deallocate(Lf,Lr)
            deallocate(Weff,Leff)               !twh - added this line 01/08/03
            deallocate(Wt,Lt)            !twh - added this line 01/08/03
            deallocate(Rscale,Rcx)              !NLB - added this line 10/11/04
         endif

         return
      end
