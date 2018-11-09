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
      subroutine defbuild
!************************************************************************
! defbuild- define building geometries    
!    - called by bcsetup.f       
!    - calls pentagon.f                   
! this is a new subroutine as of 7/26/03 that allows buildings that are
! none orthagonal to the main coordinate system to be used.
! gamma - building angle, is valid +/- 45 degrees
!
!ERP 2003         
!************************************************************************
         use datamodule ! make data from module "datamodule" visible

         implicit none
         real x1,x2,x3,x4,y1,y2,y3,y4,xL1,xL2,yL3,yL4,x_c,y_c,z_c
         real x1in,x2in,x3in,x4in,y1in,y2in,y3in,y4in,xL1in,xL2in,yL3in,yL4in
         real xfoin,yfoin,court_frac
         real slope,xmin,ymin,xmax,ymax
         !real testcany,testcanx,sstaryN,sstaryS,sstarxW,sstarxE,LoverH ! MAN 03/09/2007
         real yco,xco
!         real chk,chk2,chkH,Hlow_east,Hlow_west,Hlow_north,Hlow_south !MAN 7/6/2006
!         real Lup_north,Lup_south,Wup_north,Wup_south,Lup_east,Lup_west,Wup_east,Wup_west ! MAN 03/09/2007
         real radius,thetacell,radius_out,radius_in,r_c,wall_thickness,roof_ratio,roof_zfo
         real z_c_x,z_c_y,rayintersect
         !integer, allocatable :: bldq(:)
         integer kroof,ilevel,ivert,startpoly,numcrossing,omega
         !real det
         
         pi=4.*atan(1.0) !NLB,pi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! this is just a building generation part
! define cellflags & initialize all cells as fluid cells
! make solid gound, ie the floor
         icellflag(:,:,:)=1
         icellflag(:,:,1)=0 ! MAN 7/8/2005 Celltype definition change
         ibldflag(:,:,:)=0 ! MAN 8/29/2007 building flags
         if(lu_canopy_flag .gt. 0 .or. inumcanopy.gt.0)then
            canopy_top(:,:)=0.
            canopy_atten(:,:,:)=0.
         endif
         !if(inumpolygon .gt. 0)allocate(bldq(inumpolygon))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate buildings
! this now includes a building generation loop that allows
! for multiple buildings
! calculate building spacing s and set wake flags
         !$omp parallel do
         do ibuild=1,inumbuild
            if(bldgeometry(ibuild) .lt. 6)then
               istart(ibuild)=nint(xfo(ibuild)/dx)+1     ! front of the building  
               iend(ibuild)=nint((xfo(ibuild)+Lti(ibuild))/dx) ! back of the bld  
               jend(ibuild)=nint((yfo(ibuild)+Wt(ibuild))/dy)  ! far side of bld  
               jstart(ibuild)=nint((yfo(ibuild)-Wt(ibuild))/dy)+1! close side of bld
            endif
         enddo
         !$omp end parallel do
         if(lu_canopy_flag .gt. 0)then
            !$omp parallel do private(i,k)
            do j=1,ny-1
               do i=1,nx-1
                  canopy_top(i,j)=landuse_height(i,j)
                  if(canopy_top(i,j) .gt. 0.)then
                     do k=2,nz-1
                        if(k .eq. 2 .and. canopy_top(i,j) .lt. zm(k))then
                           canopy_top(i,j)=0.
                           exit
                        endif
                        canopy_atten(i,j,k)=landuse_atten(i,j)
                        if(canopy_top(i,j) .le. zm(k))exit
                     enddo
                  endif
               enddo
            enddo
            !$omp end parallel do
         endif
!building Loop2Loop2Loop2!Loop2Loop2Loop2!Loop2Loop2Loop2!Loop2Loop2 begin
!ccccc need an even number of cells in building
         !$omp parallel do private(i,j,k,ilevel,x1,x2,x3,x4,y1,y2,y3,y4, &
         !$omp xmin,xmax,ymin,ymax,x_c,y_c,slope,xL1,xL2,yL3,yL4,xco,yco, &
         !$omp thetacell,roof_ratio,roof_zfo,kroof,z_c,court_frac, &
         !$omp xfoin,yfoin,x1in,x2in,x3in,x4in,y1in,y2in,y3in,y4in, &
         !$omp xL1in,xL2in,yL3in,yL4in,z_c_x,z_c_y,wall_thickness, &
         !$omp radius_out,radius_in,r_c,ivert,numcrossing,startpoly,rayintersect)
lp002:   do ibuild=1,inumbuild         !begin building gen loop 2
! set non-fluid cell type flags
! icellflag = 0 is a solid cell
            if(bld_damage(ibuild) .eq. 2)then
               cycle
            endif
            ilevel=0
            if(lu_canopy_flag .gt. 0 .and. bldtype(ibuild) .eq. 2 &
                  .and. kend(ibuild) .lt. 2)then
               kend(ibuild)=2
            endif
! if the building is orthogonal to the main coordinate system
! just set the cellflags as follows
            select case(bldgeometry(ibuild))
               case(1) !Rectangular buildings
                  if(gamma(ibuild) .eq. 0)then
!erp  do k=int(zfo(ibuild)),kend(ibuild)
! int changed to nint on next line 8-14-06   
                     do j=jstart(ibuild),jend(ibuild)
                        do i=istart(ibuild),iend(ibuild)
                           select case(bldtype(ibuild))
                              case(0)
                                 icellflag(i,j,kstart(ibuild):kend(ibuild))=1
                                 ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                              case(2)
                                 do k=kstart(ibuild),kend(ibuild)
                                    if(icellflag(i,j,k) .ne. 0)then
                                       if(lu_canopy_flag .gt. 0)then
                                          if(canopy_top(i,j) .eq. landuse_height(i,j))then
                                             if(k .eq. 2)canopy_atten(i,j,:)=0.
                                             if(Ht(ibuild) .lt. 0.5*dz_array(1))then
                                                canopy_top(i,j)=0.
                                             else
                                                canopy_top(i,j)=Ht(ibuild)
                                                canopy_atten(i,j,k)=atten(ibuild)
                                             endif
                                          elseif(Ht(ibuild) .gt. canopy_top(i,j))then
                                             canopy_top(i,j)=Ht(ibuild)
                                          endif
                                       else
                                          if(Ht(ibuild) .gt. canopy_top(i,j))then
                                             canopy_top(i,j)=Ht(ibuild)
                                          endif
                                          canopy_atten(i,j,k)=atten(ibuild)
                                       endif
                                    endif
                                 enddo
                              case(3)
                                 ilevel=0
                                 do k=kstart(ibuild),kend(ibuild)
                                    ilevel=ilevel+1
                                    if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                    icellflag(i,j,k)=0
                                    ibldflag(i,j,k)=ibuild
                                 enddo                             
                              case default
                                 icellflag(i,j,kstart(ibuild):kend(ibuild))=0
                                 ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                           endselect
                        enddo
                     enddo
! if the building is NON-orthoganol to the main coordinate system
! use the following algorithm
                  else
! calculate corner coordinates of the building
                     x1=xfo(ibuild)+Wt(ibuild)*sin(gamma(ibuild))
                     y1=yfo(ibuild)-Wt(ibuild)*cos(gamma(ibuild))
                     x2=x1+Lti(ibuild)*cos(gamma(ibuild))
                     y2=y1+Lti(ibuild)*sin(gamma(ibuild))
                     x4=xfo(ibuild)-Wt(ibuild)*sin(gamma(ibuild))
                     y4=yfo(ibuild)+Wt(ibuild)*cos(gamma(ibuild))
                     x3=x4+Lti(ibuild)*cos(gamma(ibuild))
                     y3=y4+Lti(ibuild)*sin(gamma(ibuild))
 271                 format(8f8.3)
                     if(gamma(ibuild).gt.0)then
                        xmin=x4
                        xmax=x2
                        ymin=y1
                        ymax=y3
                     endif
                     if(gamma(ibuild).lt.0)then
                        xmin=x1
                        xmax=x3
                        ymin=y2
                        ymax=y4
                     endif
                     istart(ibuild)=nint(xmin/dx)
                     iend(ibuild)=nint(xmax/dx)
                     jstart(ibuild)=nint(ymin/dy)
                     jend(ibuild)=nint(ymax/dy)
!erp  do k=int(zfo(ibuild)),kend(ibuild)  
!erp        do j=int(ymin),int(ymax)
!erp     do i=int(xmin),int(xmax)
!erp     x_c=real(i) !x coordinate to be checked
!erp     y_c=real(j) !y coordinate to be checked
! changed int to nint in next three lines 8-14-06
                     do j=nint(ymin/dy)+1,nint(ymax/dy)+1   !convert back to real world unit, TZ 10/29/04
                        do i=nint(xmin/dx)+1,nint(xmax/dx)+1   !convert back to real world unit, TZ 10/29/04
                           x_c=(real(i)-0.5)*dx !x coordinate to be checked   !convert back to real world unit, TZ 10/29/04
                           y_c=(real(j)-0.5)*dy !y coordinate to be checked   !convert back to real world unit, TZ 10/29/04
!calculate the equations of the lines making up the 4 walls of the
!building
						   if( x4 .eq. x1)x4=x4+.0001
                           slope = (y4-y1)/(x4-x1) !slope of L1
                           xL1 = x4 + (y_c-y4)/slope
                           if( x3 .eq. x2)x3=x3+.0001
                           slope = (y3-y2)/(x3-x2) !slope of L2
                           xL2 = x3 + (y_c-y3)/slope
                           if( x2 .eq. x1)x2=x2+.0001
                           slope = (y2-y1)/(x2-x1) !slope of L3
                           yL3 = y1 + slope*(x_c-x1)
                           if( x3 .eq. x4)x3=x3+.0001
                           slope = (y3-y4)/(x3-x4) !slope of L4
                           yL4 = y4 + slope*(x_c-x4)
                           if(x_c.gt.xL1.and.x_c.lt.xL2.and.y_c.gt.yL3.and.y_c.lt.yL4)then
                              select case(bldtype(ibuild))
                                 case(0)
                                    icellflag(i,j,kstart(ibuild):kend(ibuild))=1
                                    ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                                 case(2)
                                    do k=kstart(ibuild),kend(ibuild)
                                       if(icellflag(i,j,k) .ne. 0)then
                                          if(lu_canopy_flag .gt. 0)then
                                             if(canopy_top(i,j) .eq. landuse_height(i,j))then
                                                if(k .eq. 2)canopy_atten(i,j,:)=0.
                                                if(Ht(ibuild) .lt. 0.5*dz_array(1))then
                                                   canopy_top(i,j)=0.
                                                else
                                                   canopy_top(i,j)=Ht(ibuild)
                                                   canopy_atten(i,j,k)=atten(ibuild)
                                                endif
                                             elseif(Ht(ibuild) .gt. canopy_top(i,j))then
                                                canopy_top(i,j)=Ht(ibuild)
                                             endif
                                          else
                                             if(Ht(ibuild) .gt. canopy_top(i,j))then
                                                canopy_top(i,j)=Ht(ibuild)
                                             endif
                                             canopy_atten(i,j,k)=atten(ibuild)
                                          endif
                                       endif
                                    enddo
                                 case(3)
                                    ilevel=0
                                    do k=kstart(ibuild),kend(ibuild)
                                       ilevel=ilevel+1
                                       if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                       icellflag(i,j,k)=0
                                       ibldflag(i,j,k)=ibuild
                                    enddo                                 
                                 case default
                                    icellflag(i,j,kstart(ibuild):kend(ibuild))=0
                                    ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                              endselect
                           endif
                        enddo
                     enddo
                  endif
! generate cylindrical buildings
! need to specify a and b as the major and minor axis of
! the ellipse
! xco and yco are the coordinates of the center of the ellipse
               case(2)
                  if(aa(ibuild) .gt. 0. .and. bb(ibuild) .gt. 0.)then
                     if(gamma(ibuild) .ne. 0.)then
                        xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))
                        yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
                        istart(ibuild)=nint((xco-max(Lt(ibuild),Wt(ibuild)))/dx)
                        iend(ibuild)=nint((xco+max(Lt(ibuild),Wt(ibuild)))/dx)
                        jstart(ibuild)=nint((yco-max(Lt(ibuild),Wt(ibuild)))/dy)
                        jend(ibuild)=nint((yco+max(Lt(ibuild),Wt(ibuild)))/dy)
                     else
                        xco = xfo(ibuild) + Lt(ibuild)
                        yco = yfo(ibuild)
                     endif
                     if(istart(ibuild) .le. 0)istart(ibuild)=1
                     if(iend(ibuild) .le. 0)iend(ibuild)=nx-1
                     if(jstart(ibuild) .le. 0)jstart(ibuild)=1
                     if(jend(ibuild) .le. 0)jend(ibuild)=ny-1
!erp 7/23/03 do k=1,kend(ibuild)
!erp  do k=int(zfo(ibuild)),kend(ibuild)  !erp 7/23/03
! int changed to nint in next line 8-14-06
                     do j=jstart(ibuild),jend(ibuild)
                        do i=istart(ibuild),iend(ibuild)
                           x_c=(real(i)-0.5)*dx-xco
                           y_c=(real(j)-0.5)*dy-yco
                           thetacell=atan2(y_c,x_c)
                           if(sqrt(x_c**2.+y_c**2.) .le. radius(aa(ibuild),bb(ibuild),&
                                 thetacell,gamma(ibuild)))then
                              select case(bldtype(ibuild))
                                 case(0)
                                    icellflag(i,j,kstart(ibuild):kend(ibuild))=1
                                    ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                                 case(2)
                                    do k=kstart(ibuild),kend(ibuild)
                                       if(icellflag(i,j,k) .ne. 0)then
                                          if(lu_canopy_flag .gt. 0)then
                                             if(canopy_top(i,j) .eq. landuse_height(i,j))then
                                                if(k .eq. 2)canopy_atten(i,j,:)=0.
                                                if(Ht(ibuild) .lt. 0.5*dz_array(1))then
                                                   canopy_top(i,j)=0.
                                                else
                                                   canopy_top(i,j)=Ht(ibuild)
                                                   canopy_atten(i,j,k)=atten(ibuild)
                                                endif
                                             elseif(Ht(ibuild) .gt. canopy_top(i,j))then
                                                canopy_top(i,j)=Ht(ibuild)
                                             endif
                                          else
                                             if(Ht(ibuild) .gt. canopy_top(i,j))then
                                                canopy_top(i,j)=Ht(ibuild)
                                             endif
                                             canopy_atten(i,j,k)=atten(ibuild)
                                          endif
                                       endif
                                    enddo
                                 case(3)
                                    ilevel=0
                                    do k=kstart(ibuild),kend(ibuild)
                                       ilevel=ilevel+1
                                       if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                       icellflag(i,j,k)=0
                                       ibldflag(i,j,k)=ibuild
                                    enddo                                 
                                 case default
                                    icellflag(i,j,kstart(ibuild):kend(ibuild))=0
                                    ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                              endselect
                           endif
                        enddo
                     enddo
                  endif
! build a Pentagon shaped building
               case(3)
                  ! pentagon will be done after all other buildings are finished
! Rectangular Stadium Building
               case(4)
!calculate corner coordinates of the building
                  x1=xfo(ibuild)+Wt(ibuild)*sin(gamma(ibuild))
                  y1=yfo(ibuild)-Wt(ibuild)*cos(gamma(ibuild))
                  x2=x1+Lti(ibuild)*cos(gamma(ibuild))
                  y2=y1+Lti(ibuild)*sin(gamma(ibuild))
                  x4=xfo(ibuild)-Wt(ibuild)*sin(gamma(ibuild))
                  y4=yfo(ibuild)+Wt(ibuild)*cos(gamma(ibuild))
                  x3=x4+Lti(ibuild)*cos(gamma(ibuild))
                  y3=y4+Lti(ibuild)*sin(gamma(ibuild))
                  xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))
                  yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
                  if(gamma(ibuild) .gt. 0.)then
                     xmin=x4
                     xmax=x2
                     ymin=y1
                     ymax=y3
                  elseif(gamma(ibuild) .lt. 0.)then
                     xmin=x1
                     xmax=x3
                     ymin=y2
                     ymax=y4
                  else
                     xmin=x1
                     xmax=x3
                     ymin=y2
                     ymax=y4
                  endif
                  istart(ibuild)=nint(xmin/dx)
                  iend(ibuild)=nint(xmax/dx)
                  jstart(ibuild)=nint(ymin/dy)
                  jend(ibuild)=nint(ymax/dy)
                  if(bldroof(ibuild) .eq. 0)then
                     roof_ratio=1
                     roof_zfo=Ht(ibuild)
                  else
                     roof_ratio=0.8
                     roof_zfo=(Ht(ibuild)-zfo_actual(ibuild))*roof_ratio+zfo_actual(ibuild)
                  endif
                  ! MAN 07/25/2008 stretched vertical grid
                  do k=kstart(ibuild),nz-1
                     kroof=k
                     if(roof_zfo .le. z(k-1))exit
                  enddo
! changed int to nint in next three lines 8-14-06
                  do k=kstart(ibuild),kroof
                     if(bldtype(ibuild) .eq. 3)then
                        ilevel=ilevel+1
                        if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                     endif
                     z_c=zm(k)-zfo_actual(ibuild)  !z coordinate to be checked
                     court_frac=bldwall(ibuild)*(1.-z_c/(roof_zfo-zfo_actual(ibuild)))
                     xfoin=xfo(ibuild)+court_frac*cos(gamma(ibuild))
                     yfoin=yfo(ibuild)+court_frac*sin(gamma(ibuild))
                     x1in=xfoin+(Wt(ibuild)-court_frac)*sin(gamma(ibuild))
                     y1in=yfoin-(Wt(ibuild)-court_frac)*cos(gamma(ibuild))
                     x2in=x1in+(Lti(ibuild)-2.*court_frac)*cos(gamma(ibuild))
                     y2in=y1in+(Lti(ibuild)-2.*court_frac)*sin(gamma(ibuild))
                     x4in=xfoin-(Wt(ibuild)-court_frac)*sin(gamma(ibuild))
                     y4in=yfoin+(Wt(ibuild)-court_frac)*cos(gamma(ibuild))
                     x3in=x4in+(Lti(ibuild)-2.*court_frac)*cos(gamma(ibuild))
                     y3in=y4in+(Lti(ibuild)-2.*court_frac)*sin(gamma(ibuild))
                     do j=nint(ymin/dy)+1,nint(ymax/dy)+1   
                        y_c=(real(j)-0.5)*dy !y coordinate to be checked
                        do i=nint(xmin/dx)+1,nint(xmax/dx)+1   
                           x_c=(real(i)-0.5)*dx !x coordinate to be checked   
!calculate the equations of the lines making up the 4 walls of the
!building
                           if(gamma(ibuild) .eq. 0.)then
                              xL1 = x1
                              xL1in = x1in
                              xL2 = x2
                              xL2in = x2in
                              yL3 = y1
                              yL3in = y1in
                              yL4 = y3
                              yL4in = y3in
                           else
                              slope = (y4-y1)/(x4-x1) !slope of L1
                              xL1 = x4 + (y_c-y4)/slope
                              xL1in = x4in + (y_c-y4in)/slope
                              slope = (y3-y2)/(x3-x2) !slope of L2
                              xL2 = x3 + (y_c-y3)/slope
                              xL2in = x3in + (y_c-y3in)/slope
                              slope = (y2-y1)/(x2-x1) !slope of L3
                              yL3 = y1 + slope*(x_c-x1)
                              yL3in = y1in + slope*(x_c-x1in)
                              slope = (y3-y4)/(x3-x4) !slope of L4
                              yL4 = y4 + slope*(x_c-x4)
                              yL4in = y4in + slope*(x_c-x4in)
                           endif
                           if(abs(xL1-xL1in) .lt. dx)xL1in=xL1+dx
                           if(abs(xL2-xL2in) .lt. dx)xL2in=xL2-dx
                           if(abs(yL3-yL3in) .lt. dy)yL3in=yL3+dy
                           if(abs(yL4-yL4in) .lt. dx)yL4in=yL4-dy
                           if(x_c .gt. xL1 .and. x_c .lt. xL2 .and. &
                                 y_c .gt. yL3 .and. y_c .lt. yL4)then
                              icellflag(i,j,k)=0
                              ibldflag(i,j,k)=ibuild ! MAN 8/29/2007 building flags
                              if(x_c .gt. xL1in .and. x_c .lt. xL2in .and. &
                                     y_c .gt. yL3in .and. y_c .lt. yL4in)then
                                 icellflag(i,j,k)=1 ! MAN 7/8/2005 Celltype definition change
                                 ibldflag(i,j,k)=0 ! MAN 8/29/2007 building flags
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  if(bldroof(ibuild) .gt. 0)then
                     court_frac=bldwall(ibuild)
                     do j=jstart(ibuild),jend(ibuild)
                        do i=istart(ibuild),iend(ibuild)
                           x_c=((real(i)-0.5)*dx-xco)*cos(gamma(ibuild)) + ((real(j)-0.5)*dy-yco)*sin(gamma(ibuild))
                           y_c=-((real(i)-0.5)*dx-xco)*sin(gamma(ibuild)) + ((real(j)-0.5)*dy-yco)*cos(gamma(ibuild))
                           if(abs(x_c) .lt. Lt(ibuild) .and. abs(y_c) .lt. Wt(ibuild))then
                              if(abs(x_c) .gt. Lt(ibuild)-court_frac .or. abs(y_c) .gt. Wt(ibuild)-court_frac)then
                                 z_c_x=roof_zfo+(Ht(ibuild)-roof_zfo)*sqrt((Lt(ibuild)-abs(x_c))/court_frac)
                                 z_c_y=roof_zfo+(Ht(ibuild)-roof_zfo)*sqrt((Wt(ibuild)-abs(y_c))/court_frac)
                                 z_c=min(z_c_x,z_c_y)
                                 ! MAN 07/25/2008 stretched vertical grid
                                 do k=kstart(ibuild),nz-1
                                    kroof=k
                                    if(z_c .le. z(k-1))exit
                                 enddo
                                 ibldflag(i,j,kroof)=ibuild ! MAN 8/29/2007 building flags
                                 icellflag(i,j,kroof)=0
                              endif
                           endif
                        enddo
                     enddo
                  endif
! Elliptical stadium building
               case(5)
!calculate corner coordinates of the building
                  if(aa(ibuild) .gt. 0. .and. bb(ibuild) .gt. 0.)then
                     if(gamma(ibuild) .ne. 0.)then
                        xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))
                        yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
                        istart(ibuild)=nint((xco-max(Lt(ibuild),Wt(ibuild)))/dx)
                        iend(ibuild)=nint((xco+max(Lt(ibuild),Wt(ibuild)))/dx)
                        jstart(ibuild)=nint((yco-max(Lt(ibuild),Wt(ibuild)))/dy)
                        jend(ibuild)=nint((yco+max(Lt(ibuild),Wt(ibuild)))/dy)
                     else
                        xco = xfo(ibuild) + Lt(ibuild)
                        yco = yfo(ibuild)
                     endif
                     if(bldroof(ibuild) .eq. 0)then
                        roof_ratio=1
                        roof_zfo=Ht(ibuild)
                     else
                        roof_ratio=0.8
                        roof_zfo=(Ht(ibuild)-zfo_actual(ibuild))*roof_ratio+zfo_actual(ibuild)
                     endif
                     if(istart(ibuild) .le. 0)istart(ibuild)=1
                     if(iend(ibuild) .le. 0)iend(ibuild)=nx-1
                     if(jstart(ibuild) .le. 0)jstart(ibuild)=1
                     if(jend(ibuild) .le. 0)jend(ibuild)=ny-1
                     wall_thickness=max(dx,dy)
                     ! MAN 07/25/2008 stretched vertical grid
                     do k=kstart(ibuild),nz-1
                        kroof=k
                        if(roof_zfo .le. z(k-1))exit
                     enddo
!erp 7/23/03 do k=1,kend(ibuild)
!erp  do k=int(zfo(ibuild)),kend(ibuild)  !erp 7/23/03
! int changed to nint in next line 8-14-06
                     do k=kstart(ibuild),kroof
                        if(bldtype(ibuild) .eq. 3)then
                           ilevel=ilevel+1
                           if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                        endif
                        z_c=zm(k)-zfo_actual(ibuild)
                        court_frac=bldwall(ibuild)*(1.-z_c/(roof_zfo-zfo_actual(ibuild)))
                        do j=jstart(ibuild),jend(ibuild)
                           y_c=(real(j)-0.5)*dy-yco
                           do i=istart(ibuild),iend(ibuild)
                              x_c=(real(i)-0.5)*dx-xco
                              thetacell=atan2(y_c,x_c)
                              radius_out=radius(aa(ibuild),bb(ibuild),thetacell,gamma(ibuild))
                              radius_in=radius(aa(ibuild)-court_frac,bb(ibuild)-court_frac,thetacell,gamma(ibuild))
                              if(radius_out-radius_in .le. 1.25*wall_thickness)radius_in=radius_out-1.25*wall_thickness
                              r_c=sqrt(x_c**2.+y_c**2.)
                              if(r_c .le. radius_out .and. r_c .gt. radius_in)then
                                 ibldflag(i,j,k)=ibuild ! MAN 8/29/2007 building flags
                                 icellflag(i,j,k)=0
                              endif
                           enddo
                        enddo
                     enddo
                     if(bldroof(ibuild) .gt. 0)then
                        court_frac=bldwall(ibuild)
                        do j=jstart(ibuild),jend(ibuild)
                           y_c=(real(j)-0.5)*dy-yco
                           do i=istart(ibuild),iend(ibuild)
                              x_c=(real(i)-0.5)*dx-xco
                              thetacell=atan2(y_c,x_c)
                              radius_out=radius(aa(ibuild),bb(ibuild),thetacell,gamma(ibuild))
                              radius_in=radius(aa(ibuild)-court_frac,bb(ibuild)-court_frac,thetacell,gamma(ibuild))
                              r_c=sqrt(x_c**2.+y_c**2.)
                              if(r_c .gt. radius_in .and. r_c .lt. radius_out)then
                                 z_c=roof_zfo+(Ht(ibuild)-roof_zfo)*sqrt((radius_out-r_c)/(radius_out-radius_in))
                                 ! MAN 07/25/2008 stretched vertical grid
                                 do k=kstart(ibuild),nz-1
                                    kroof=k
                                    if(z_c .le. z(k-1))exit
                                 enddo
                                 ibldflag(i,j,kroof)=ibuild ! MAN 8/29/2007 building flags
                                 icellflag(i,j,kroof)=0
                              endif
                           enddo
                        enddo
                     endif
                  endif
               case(6)
                  x1=bldx(bldstartidx(ibuild))
                  x2=x1
                  y1=bldy(bldstartidx(ibuild))
                  y2=y1
                  do ivert=bldstartidx(ibuild)+1,bldstopidx(ibuild)
                     if(bldx(ivert) .eq. bldx(bldstartidx(ibuild)) .and. &
                           bldy(ivert) .eq. bldy(bldstartidx(ibuild)))exit
                     if(bldx(ivert) .lt. x1)x1=bldx(ivert)
                     if(bldx(ivert) .gt. x2)x2=bldx(ivert)
                     if(bldy(ivert) .lt. y1)y1=bldy(ivert)
                     if(bldy(ivert) .gt. y2)y2=bldy(ivert)
                  enddo
                  istart(ibuild)=int(x1/dx)+1      !front of the building  
                  iend(ibuild)=int(x2/dx)+1  !back of the bld  
                  jend(ibuild)=int(y2/dy)+1  !far side of bld  
                  jstart(ibuild)=int(y1/dy)+1!close side of bld
                  do j=jstart(ibuild),jend(ibuild)
                     y_c=(real(j)-0.5)*dy
                     do i=istart(ibuild),iend(ibuild)
                        x_c=(real(i)-0.5)*dx
                        ivert=bldstartidx(ibuild)
                        startpoly=ivert
                        numcrossing=0
                        !Based on Wm. Randolph Franklin, "PNPOLY - Point Inclusion in Polygon Test" Web Page (2000)
                        do while(ivert .lt. bldstopidx(ibuild))
                           if(((bldy(ivert) .le. y_c) .and. (bldy(ivert+1) .gt. y_c)) &
                                 .or. ((bldy(ivert) .gt. y_c) .and. (bldy(ivert+1) .le. y_c)))then
                              rayintersect=(y_c-bldy(ivert))/(bldy(ivert+1)-bldy(ivert))
                              if(x_c .lt. bldx(ivert)+rayintersect*(bldx(ivert+1)-bldx(ivert)))then
                                 numcrossing=numcrossing+1
                              endif
                           endif
                           ivert=ivert+1
                           if(bldx(ivert) .eq. bldx(startpoly) .and. &
                                 bldy(ivert) .eq. bldy(startpoly))then
                              ivert=ivert+1
                              startpoly=ivert
                           endif
                        enddo
                        if(numcrossing/2 .ne. ceiling(0.5*real(numcrossing)))then
                        
                        !!Based on Hormann and Agathos (2001)
                        !do while(ivert .lt. bldstopidx(ibuild))
                        !   if(bldx(ivert) .gt. x_c .and. bldy(ivert) .ge. y_c)bldq(ivert)=0
                        !   if(bldx(ivert) .le. x_c .and. bldy(ivert) .gt. y_c)bldq(ivert)=1
                        !   if(bldx(ivert) .lt. x_c .and. bldy(ivert) .le. y_c)bldq(ivert)=2
                        !   if(bldx(ivert) .ge. x_c .and. bldy(ivert) .lt. y_c)bldq(ivert)=3
                        !   ivert=ivert+1
                        !   if(bldx(ivert) .eq. bldx(startpoly) .and. &
                        !         bldy(ivert) .eq. bldy(startpoly))then
                        !      bldq(ivert)=bldq(startpoly)
                        !      ivert=ivert+1
                        !      startpoly=ivert
                        !   endif
                        !enddo
                        !ivert=bldstartidx(ibuild)
                        !startpoly=ivert
                        !omega=0
                        !do while(ivert .lt. bldstopidx(ibuild))
                        !   select case(bldq(ivert+1)-bldq(ivert))
                        !      case(-3)
                        !         omega=omega+1
                        !      case(3)
                        !         omega=omega-1
                        !      case(-2)
                        !         if(det(bldx(ivert),bldy(ivert),bldx(ivert+1),bldy(ivert+1),x_c,y_c) .gt. 0)then
                        !            omega=omega+1
                        !         endif
                        !      case(2)
                        !         if(det(bldx(ivert),bldy(ivert),bldx(ivert+1),bldy(ivert+1),x_c,y_c) .lt. 0)then
                        !            omega=omega-1
                        !         endif
                        !   endselect
                        !   ivert=ivert+1
                        !   if(bldx(ivert) .eq. bldx(startpoly) .and. &
                        !         bldy(ivert) .eq. bldy(startpoly))then
                        !      ivert=ivert+1
                        !      startpoly=ivert
                        !   endif
                        !enddo
                        !if(omega .ne. 0)then
                           select case(bldtype(ibuild))
                              case(0)
                                 icellflag(i,j,kstart(ibuild):kend(ibuild))=1
                                 ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                              case(2)
                                 do k=kstart(ibuild),kend(ibuild)
                                    if(icellflag(i,j,k) .ne. 0)then
                                       if(lu_canopy_flag .gt. 0)then
                                          if(canopy_top(i,j) .eq. landuse_height(i,j))then
                                             if(k .eq. 2)canopy_atten(i,j,:)=0.
                                             if(Ht(ibuild) .lt. 0.5*dz_array(1))then
                                                canopy_top(i,j)=0.
                                             else
                                                canopy_top(i,j)=Ht(ibuild)
                                                canopy_atten(i,j,k)=atten(ibuild)
                                             endif
                                          elseif(Ht(ibuild) .gt. canopy_top(i,j))then
                                             canopy_top(i,j)=Ht(ibuild)
                                          endif
                                       else
                                          if(Ht(ibuild) .gt. canopy_top(i,j))then
                                             canopy_top(i,j)=Ht(ibuild)
                                          endif
                                          canopy_atten(i,j,k)=atten(ibuild)
                                       endif
                                    endif
                                 enddo
                              case(3)
                                 ilevel=0
                                 do k=kstart(ibuild),kend(ibuild)
                                    ilevel=ilevel+1
                                    if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                    icellflag(i,j,k)=0
                                    ibldflag(i,j,k)=ibuild
                                 enddo
                              case(5)
                                 icellflag(i,j,kstart(ibuild):kend(ibuild))=1
                                 ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild                                 
                              case default
                                 icellflag(i,j,kstart(ibuild):kend(ibuild))=0
                                 ibldflag(i,j,kstart(ibuild):kend(ibuild))=ibuild
                           endselect
                        endif
                     enddo
                  enddo
            endselect
! erp 1/31/2003
         enddo   lp002
         !$omp end parallel do
         do ibuild=1,inumbuild
            if(bldgeometry(ibuild) .eq. 3)call pentagon
         enddo
         !if(inumpolygon .gt. 0)deallocate(bldq)
! end building generation 2a
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         return
      end


!      real function det(Pxi,Pyi,Pxip1,Pyip1,Rx,Ry)
!         implicit none
!         real Pxi,Pxip1,Pyi,Pyip1,Rx,Ry
!         det=(Pxi-Rx)*(Pyip1-Ry)-(Pxip1-Rx)*(Pyi-Ry)
!         return
!      end
