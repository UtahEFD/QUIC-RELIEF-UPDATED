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
      subroutine streetcanyon
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! STREETCANYON applies the canyon vortices for building separations
! that are less than the critical value for skimming flow.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule ! make data from module "datamodule" visible
         implicit none
         
         integer perpendicular_flag,x_idx,y_idx,x_idx_max,canyon_flag
         integer k_ref,top_flag,dbuild,ic,jc,circle_flag,iu,ju,iv,jv
         integer x_idx_min,reverse_flag
         real uo_h,vo_h,upwind_dir,upwind_rel,downwind_rel,xco,yco
         real beta,xd,yd,thetad,xcd,ycd,rd,along_dir,cross_dir
         real velmag,canyon_dir,along_mag,cross_mag,usign,xc,yc
         real x1,y1,x2,y2,x3,y3,x4,y4,xwall,xpos,S
         real xu,yu,xv,yv,xw,yw,xp,yp,xwallu,xwallv,xwallw
         real xw1,yw1,xw2,yw2,xw3,yw3,tol
         real thetamax,thetamin,thetai,xnorm_bisect,radius
         real upwind_norm,upwind_norm1,upwind_norm2,angle_tol
         real ucomponent,vcomponent,numu,numv
         integer ivert,jvert
         real segmentLength
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Begin Street Canyon vortex subsection for skimming 
! flow regime 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
         angle_tol=3*pi/4
         do ibuild=1,inumbuild
            if(bld_damage(ibuild) .eq. 2)cycle
            if(bldgeometry(ibuild) .eq. 6)then
               xco = bldcx(ibuild)
               yco = bldcy(ibuild)
            else
               xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))!CENTER of building in QUIC domain coordinates
               yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
            endif
            ! find upwind direction and determine the type of flow regime
            uo_h=uo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
            vo_h=vo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
            upwind_norm=0.
            upwind_norm1=0.
            upwind_norm2=0.
            upwind_dir=atan2(vo_h,uo_h)
            upwind_rel=upwind_dir-gamma(ibuild)
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
            select case(bldgeometry(ibuild))
               case(1,4)
                  tol=0.01*pi/180.
                  x1=xfo(ibuild)+Wt(ibuild)*sin(gamma(ibuild))-xco
                  y1=yfo(ibuild)-Wt(ibuild)*cos(gamma(ibuild))-yco
                  x2=x1+Lti(ibuild)*cos(gamma(ibuild))
                  y2=y1+Lti(ibuild)*sin(gamma(ibuild))
                  x4=xfo(ibuild)-Wt(ibuild)*sin(gamma(ibuild))-xco
                  y4=yfo(ibuild)+Wt(ibuild)*cos(gamma(ibuild))-yco
                  x3=x4+Lti(ibuild)*cos(gamma(ibuild))
                  y3=y4+Lti(ibuild)*sin(gamma(ibuild))
                  if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
                     xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     xw2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     upwind_norm1=pi+gamma(ibuild)
                     upwind_norm2=0.5*pi+gamma(ibuild)
                     perpendicular_flag=0
                     usign=1
                  elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
                     xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     upwind_norm=0.5*pi+gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
                     xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     upwind_norm1=0.5*pi+gamma(ibuild)
                     upwind_norm2=gamma(ibuild)
                     perpendicular_flag=0
                     usign=-1
                  elseif(abs(upwind_rel) .le. tol)then
                     xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     upwind_norm=gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
                     xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     upwind_norm1=gamma(ibuild)
                     upwind_norm2=-0.5*pi+gamma(ibuild)
                     perpendicular_flag=0
                     usign=1
                  elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
                     xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     upwind_norm=-0.5*pi+gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
                     xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     upwind_norm1=-0.5*pi+gamma(ibuild)
                     upwind_norm2=-pi+gamma(ibuild)
                     perpendicular_flag=0
                     usign=-1
                  else
                     xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     upwind_norm=-pi+gamma(ibuild)
                     perpendicular_flag=1
                  endif
                  if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                  if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                  if(upwind_norm1 .gt. pi)upwind_norm1=upwind_norm1-2*pi
                  if(upwind_norm1 .le. -pi)upwind_norm1=upwind_norm1+2*pi
                  if(upwind_norm2 .gt. pi)upwind_norm2=upwind_norm2-2*pi
                  if(upwind_norm2 .le. -pi)upwind_norm2=upwind_norm2+2*pi
                  do y_idx=0,2*ceiling((yw1-yw3)/dxy)
                     yc=0.5*real(y_idx)*dxy+yw3
                     ! print*,yc
                     top_flag=0
                     if(perpendicular_flag .gt. 0)then
                        xwall=xw1
                     elseif(yc.ge.yw2)then
                        xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
                        upwind_norm=upwind_norm1
                     else
                        xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
                        upwind_norm=upwind_norm2
                     endif
                     do k=kend(ibuild),2,-1
                        canyon_flag=0
                        S=0.
                        x_idx_min=-1
                        reverse_flag=0
                        do x_idx=1,2*ceiling(Lr(ibuild)/dxy)
                           xc=0.5*real(x_idx)*dxy
                           i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                           j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                           if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                              exit
                           endif
                           if(icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                              x_idx_min=x_idx
                           endif
                           if(icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                              canyon_flag=1
                              x_idx_max=x_idx-1
                              S=0.5*real(x_idx_max-x_idx_min)*dxy
                              if(top_flag .eq. 0 .and. S .gt. 0.)then
                                 k_ref=k+1
                                 ic=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                 jc=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                 if(icellflag(ic,jc,k_ref) .ne. 0)then
                                    numu=0.
                                    numv=0.
                                    ucomponent=0.
                                    vcomponent=0.
                                    if(icellflag(ic-1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+uo(ic,jc,k_ref)
                                    endif
                                    if(icellflag(ic+1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+uo(ic+1,jc,k_ref)
                                    endif
                                    if(icellflag(ic,jc-1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+vo(ic,jc,k_ref)
                                    endif
                                    if(icellflag(ic,jc+1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+vo(ic,jc+1,k_ref)
                                    endif
                                    if( ucomponent .ne. 0. .and. numu > 0.) then
                                       ucomponent=ucomponent/numu
                                    else
                                       ucomponent=0.
                                    endif
                                    if( vcomponent .ne. 0. .and. numv > 0.) then
                                       vcomponent=vcomponent/numv
                                    else
                                       vcomponent=0.
                                    endif
                                    if(numu .eq. 0. .and. numv .eq. 0.)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    elseif(numu .gt. 0 .and. numv .gt. 0.)then
                                       velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                       canyon_dir=atan2(vcomponent,ucomponent)
                                    elseif(numu .gt. 0)then
                                       velmag=abs(ucomponent)
                                       if(ucomponent .gt. 0.)then
                                          canyon_dir=0.
                                       else
                                          canyon_dir=pi
                                       endif
                                    else
                                       velmag=abs(vcomponent)
                                       if(vcomponent .gt. 0.)then
                                          canyon_dir=0.5*pi
                                       else
                                          canyon_dir=-0.5*pi
                                       endif
                                    endif
                                    top_flag=1
                                    if(Ht(ibldflag(i,j,k)) .lt. Ht(ibuild) .and. z(k)/S .lt. 0.65)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    endif
                                 else
                                    canyon_flag=0
                                    top_flag=0
                                    S=0.
                                    exit
                                 endif
                                 if(velmag .gt. max_velmag)then
                                    print*,'Parameterized velocity exceeds max in street canyon',&
                                       velmag,max_velmag,i,j,k,yc
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                              endif
                              dbuild=ibldflag(i,j,k)
! Find the along canyon and cross canyon directions
                              i=ceiling(((xc-0.5*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                              j=ceiling(((xc-0.5*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                              if(bldtype(dbuild) .eq. 4)then
                                 canyon_flag=0
                                 S=0.
                                 top_flag=0
                                 exit
                              endif
                              select case(bldgeometry(dbuild))
                                 case(1,4)
                                    beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                    downwind_rel=canyon_dir-gamma(dbuild)
                                    if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                    if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                    xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                    ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                    xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                    yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                       if(downwind_rel .le. 0.)then
                                          if(downwind_rel .le. -0.5*pi)then
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                       if(abs(downwind_rel) .ge. 0.5*pi)then
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                       if(downwind_rel .ge. 0.)then
                                          if(downwind_rel .le. 0.5*pi)then
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    else
                                       if(abs(downwind_rel) .lt. 0.5*pi)then
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    endif
                                 case(2,5)
                                    downwind_rel=canyon_dir-gamma(dbuild)
                                    if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                    if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                    xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                    ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                    xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                    yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                       ((Wt(dbuild)*cos(thetad))**2.))
                                    along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                       (Lt(dbuild)**2.)*rd*sin(thetad))
                                    if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       else
                                          along_dir=along_dir-pi+gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       endif
                                    else
                                       reverse_flag=1
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       else
                                          along_dir=along_dir-pi+gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       endif
                                    endif
                                 case(6)
                                    do jvert=bldstartidx(dbuild),bldstopidx(dbuild)
                                       cross_dir=atan2(bldy(jvert+1)-bldy(jvert),bldx(jvert+1)-bldx(jvert))+0.5*pi
                                       if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                       xcd=0.5*(bldx(jvert+1)+bldx(jvert))
                                       ycd=0.5*(bldy(jvert+1)+bldy(jvert))
                                       xd=((real(i)-0.5)*dx-xcd)*cos(cross_dir)+&
                                          ((real(j)-0.5)*dy-ycd)*sin(cross_dir)
                                       yd=-((real(i)-0.5)*dx-xcd)*sin(cross_dir)+&
                                          ((real(j)-0.5)*dy-ycd)*cos(cross_dir)
                                       if(abs(xd) .lt. dxy)then
                                          segmentLength=sqrt(((bldy(jvert+1)-bldy(jvert))**2.) &
                                             +((bldx(jvert+1)-bldx(jvert))**2.))
                                          if(abs(yd) .le. 0.5*segmentLength)then
                                             downwind_rel=canyon_dir-cross_dir
                                             if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                             if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir-0.5*pi
                                                else
                                                   along_dir=cross_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=0
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir+0.5*pi
                                                else
                                                   along_dir=cross_dir-0.5*pi
                                                endif
                                             endif
                                             if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                             if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                             exit
                                          endif
                                       endif
                                       if(bldx(jvert+1) .eq. bldx(bldstartidx(dbuild)) &
                                             .and. bldy(jvert+1) .eq. bldy(bldstartidx(dbuild)))exit
                                    enddo
                              end select
                              if(along_dir .gt. pi)along_dir=along_dir-2*pi
                              if(along_dir .le. -pi)along_dir=along_dir+2*pi
                              if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                              if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                              if(reverse_flag .eq. 1)then
                                 if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              else
                                 if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              endif
                              exit
                           endif
                        enddo
                        ! if(top_flag .eq. 1)print*,yc,along_dir,cross_dir
                        if(canyon_flag .eq. 1 .and. S .gt. 0.9*dxy)then
                           along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                              log(zm(k)/zo)/log(zm(k_ref)/zo)
                           cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                           if(abs(along_mag) .gt. max_velmag)then
                              print*,'Along canyon exceeds max in street canyon',&
                                 along_mag,max_velmag,i,j,k
                           endif
                           if(abs(cross_mag) .gt. max_velmag)then
                              print*,'Cross canyon exceeds max in street canyon',&
                                 cross_mag,max_velmag,i,j,k
                           endif
                           do x_idx=x_idx_min,x_idx_max
                              xc=0.5*real(x_idx)*dxy
                              i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                              j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                              
                              if(icellflag(i,j,k) .ne. 0)then ! icellflag(i,j,k) .ne. 6 .and. 
! u component
                                 iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                                 ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                 xp=real(iu-1)*dx-xco
                                 yp=(real(ju)-0.5)*dy-yco
                                 xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(perpendicular_flag .gt. 0)then
                                    xwallu=xw1
                                 elseif(yu.ge.yw2)then
                                    xwallu=((xw2-xw1)/(yw2-yw1))*(yu-yw1)+xw1
                                 else
                                    xwallu=((xw3-xw2)/(yw3-yw2))*(yu-yw2)+xw2
                                 endif
                                 xpos=xu-xwallu
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    uo(iu,ju,k)=along_mag*cos(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                    if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                                       print*,'Parameterized U exceeds max in street canyon',&
                                          uo(iu,ju,k),velmag,max_velmag,iu,ju,k
                                    endif
                                 endif
! v component
                                 iv=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                                 jv=nint(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                                 xp=(real(iv)-0.5)*dx-xco
                                 yp=real(jv-1)*dy-yco
                                 xv=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yv=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(perpendicular_flag .gt. 0)then
                                    xwallv=xw1
                                 elseif(yv.ge.yw2)then
                                    xwallv=((xw2-xw1)/(yw2-yw1))*(yv-yw1)+xw1
                                 else
                                    xwallv=((xw3-xw2)/(yw3-yw2))*(yv-yw2)+xw2
                                 endif
                                 xpos=xv-xwallv
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    vo(iv,jv,k)=along_mag*sin(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                    if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                                       print*,'Parameterized V exceeds max in street canyon',&
                                          vo(iv,jv,k),velmag,max_velmag,iv,jv,k
                                    endif
                                 endif
! w component
                                 xp=(real(i)-0.5)*dx-xco
                                 yp=(real(j)-0.5)*dy-yco
                                 xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(perpendicular_flag .gt. 0)then
                                    xwallw=xw1
                                 elseif(yw.ge.yw2)then
                                    xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                                 else
                                    xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                                 endif
                                 xpos=xw-xwallw
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    if(icellflag(i,j,k-1) .ne. 0)then
                                       if(reverse_flag .eq. 0)then
                                          wo(i,j,k)=-abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       else
                                          wo(i,j,k)=abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       endif
                                       if(abs(wo(i,j,k)) .gt. max_velmag)then
                                          print*,'Parameterized W exceeds max in street canyon',&
                                             wo(i,j,k),velmag,max_velmag,i,j,k
                                       endif
                                    endif
                                    icellflag(i,j,k)=6
                                 endif
                              endif
                           enddo
                        endif
                     enddo
                  enddo
               case(2,5)
                  tol=0.01*min(dx,dy)
                  if(abs(aa(ibuild)-bb(ibuild)) .lt. tol)then
                     circle_flag=1
                  else
                     circle_flag=0
                  endif
                  if(circle_flag .eq. 1)then
                     thetamin=-0.5*pi
                     thetamax=0.5*pi
                     yw1=Lt(ibuild)
                     yw3=-Lt(ibuild)
                  else
                     y1=0.
                     y2=0.
                     do i=1,180
                        thetai=real(180-i)*pi/180.
                        y1=radius(aa(ibuild),bb(ibuild),thetai,&
                            gamma(ibuild)-upwind_dir)*sin(thetai)
                        if(y1 .lt. y2)then
                           exit
                        endif
                        y2=y1
                     enddo
                     thetamax=thetai+pi/180.
                     thetamin=thetamax-pi
                     yw1=y2
                     yw3=-y2
                  endif
                  do y_idx=1,2*int((yw1-yw3)/dxy)
                     yc=0.5*real(y_idx)*dxy+yw3
                     top_flag=0
                     if(circle_flag .eq. 1)then
                        xwall=sqrt((Lt(ibuild)**2.)-(yc**2.))
                        upwind_norm=atan2(xwall*sin(upwind_dir)+yc*cos(upwind_dir),&
                              xwall*cos(upwind_dir)-yc*sin(upwind_dir))
                     else
                        xwall=xnorm_bisect(aa(ibuild),bb(ibuild),&
                              gamma(ibuild)-upwind_dir,yc,thetamin,thetamax,dxy)
                        xd=xwall*cos(gamma(ibuild))+yc*sin(gamma(ibuild))
                        yd=-xwall*sin(gamma(ibuild))+yc*cos(gamma(ibuild))
                        thetad=atan2(yd,xd)
                        rd=Lt(ibuild)*Wt(ibuild)/sqrt(((Lt(ibuild)*sin(thetad))**2.)+&
                           ((Wt(ibuild)*cos(thetad))**2.))
                        upwind_norm=atan2(-(Wt(ibuild)**2.)*rd*cos(thetad),&
                           (Lt(ibuild)**2.)*rd*sin(thetad))+gamma(ibuild)+0.5*pi
                        if(cos(upwind_norm-upwind_dir) .le. 0.)upwind_norm=upwind_norm-pi
                        if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                        if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                     endif
                     do k=kend(ibuild),kstart(ibuild),-1
                        canyon_flag=0
                        S=0.
                        x_idx_min=-1
                        reverse_flag=0
                        do x_idx=0,2*int(Lr(ibuild)/dxy)+1
                           xc=0.5*real(x_idx)*dxy
                           i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                           j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                           if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                              exit
                           endif
                           if(icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                              x_idx_min=x_idx
                           endif
                           if(icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                              canyon_flag=1
                              x_idx_max=x_idx-1
                              S=0.5*real(x_idx_max-x_idx_min)*dxy
                              if(top_flag .eq. 0)then
                                 k_ref=k+1
                                 ic=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                 jc=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                 if(icellflag(ic,jc,k_ref) .ne. 0)then
                                    numu=0.
                                    numv=0.
                                    ucomponent=0.
                                    vcomponent=0.
                                    if(icellflag(ic-1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+uo(ic,jc,k_ref)
                                    endif
                                    if(icellflag(ic+1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+uo(ic+1,jc,k_ref)
                                    endif
                                    if(icellflag(ic,jc-1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+vo(ic,jc,k_ref)
                                    endif
                                    if(icellflag(ic,jc+1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+vo(ic,jc+1,k_ref)
                                    endif
                                    if(numu > 0.)ucomponent=ucomponent/numu
                                    if(numv > 0.)vcomponent=vcomponent/numv
                                    if(numu == 0. .and. numv == 0.)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    elseif(numu > 0. .and. numv > 0.)then
                                       velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                       canyon_dir=atan2(vcomponent,ucomponent)
                                    elseif(numu > 0.)then
                                       velmag=abs(ucomponent)
                                       if(ucomponent .gt. 0.)then
                                          canyon_dir=0.
                                       else
                                          canyon_dir=pi
                                       endif
                                    else
                                       velmag=abs(vcomponent)
                                       if(vcomponent .gt. 0.)then
                                          canyon_dir=0.5*pi
                                       else
                                          canyon_dir=-0.5*pi
                                       endif
                                    endif
                                    top_flag=1
                                    if(S <= 0.)then
                                    	canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    elseif(Ht(ibldflag(i,j,k)) .lt. Ht(ibuild) .and. z(k)/S .lt. 0.65)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    endif
                                 else
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                                 if(velmag .gt. max_velmag)then
                                    print*,'Parameterized velocity exceeds max in street canyon',&
                                       velmag,max_velmag,i,j,k
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                              endif
                              dbuild=ibldflag(i,j,k)
! Find the along canyon and cross canyon directions
                              i=ceiling(((xc-0.5*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                              j=ceiling(((xc-0.5*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                              if(bldtype(dbuild) .eq. 4)then
                                 canyon_flag=0
                                 S=0.
                                 top_flag=0
                                 exit
                              endif
                              select case(bldgeometry(dbuild))
                                 case(1,4)
                                    beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                    downwind_rel=canyon_dir-gamma(dbuild)
                                    xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                    ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                    xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                    yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                       if(downwind_rel .le. 0.)then
                                          if(downwind_rel .le. -0.5*pi)then
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                       if(abs(downwind_rel) .ge. 0.5*pi)then
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                       if(downwind_rel .ge. 0.)then
                                          if(downwind_rel .le. 0.5*pi)then
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    else
                                       if(abs(downwind_rel) .lt. 0.5*pi)then
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-0.5*pi+gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    endif
                                 case(2,5)
                                    downwind_rel=canyon_dir-gamma(dbuild)
                                    xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                    ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                    xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                    yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                       ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                       ((Wt(dbuild)*cos(thetad))**2.))
                                    along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                       (Lt(dbuild)**2.)*rd*sin(thetad))
                                    if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       else
                                          along_dir=along_dir-pi+gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       endif
                                    else
                                       reverse_flag=1
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       else
                                          along_dir=along_dir-pi+gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       endif
                                    endif
                                 case(6)
                                    do jvert=bldstartidx(dbuild),bldstopidx(dbuild)
                                       cross_dir=atan2(bldy(jvert+1)-bldy(jvert),bldx(jvert+1)-bldx(jvert))+0.5*pi
                                       if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                       xcd=0.5*(bldx(jvert+1)+bldx(jvert))
                                       ycd=0.5*(bldy(jvert+1)+bldy(jvert))
                                       xd=((real(i)-0.5)*dx-xcd)*cos(cross_dir)+&
                                          ((real(j)-0.5)*dy-ycd)*sin(cross_dir)
                                       yd=-((real(i)-0.5)*dx-xcd)*sin(cross_dir)+&
                                          ((real(j)-0.5)*dy-ycd)*cos(cross_dir)
                                       if(abs(xd) .lt. dxy)then
                                          segmentLength=sqrt(((bldy(jvert+1)-bldy(jvert))**2.) &
                                             +((bldx(jvert+1)-bldx(jvert))**2.))
                                          if(abs(yd) .le. 0.5*segmentLength)then
                                             downwind_rel=canyon_dir-cross_dir
                                             if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                             if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir-0.5*pi
                                                else
                                                   along_dir=cross_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=0
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir+0.5*pi
                                                else
                                                   along_dir=cross_dir-0.5*pi
                                                endif
                                             endif
                                             if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                             if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                             exit
                                          endif
                                       endif
                                       if(bldx(jvert+1) .eq. bldx(bldstartidx(dbuild)) &
                                             .and. bldy(jvert+1) .eq. bldy(bldstartidx(dbuild)))exit
                                    enddo
                              end select
                              if(along_dir .gt. pi)along_dir=along_dir-2*pi
                              if(along_dir .le. -pi)along_dir=along_dir+2*pi
                              if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                              if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                              if(reverse_flag .eq. 1)then
                                 if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              else
                                 if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              endif
                              exit
                           endif
                        enddo
                        if(canyon_flag .eq. 1 .and. S .gt. 0.9*dxy)then
                           along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                              log(zm(k)/zo)/log(zm(k_ref)/zo)
                           cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                           if(abs(along_mag) .gt. max_velmag)then
                              print*,'Along canyon exceeds max in street canyon',&
                                 along_mag,max_velmag,i,j,k
                           endif
                           if(abs(cross_mag) .gt. max_velmag)then
                              print*,'Cross canyon exceeds max in street canyon',&
                                 cross_mag,max_velmag,i,j,k
                           endif
                           do x_idx=x_idx_min,x_idx_max
                              xc=0.5*real(x_idx)*dxy
                              i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                              j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                              if(icellflag(i,j,k) .ne. 0)then ! icellflag(i,j,k) .ne. 6 .and. 
! u component
                                 iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                                 ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                 xp=real(iu-1)*dx-xco
                                 yp=(real(ju)-0.5)*dy-yco
                                 xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(circle_flag .eq. 1)then
                                    if(abs(yu) .gt. Lt(ibuild))then
                                       cycle
                                    else
                                       xwallu=sqrt((Lt(ibuild)**2.)-(yu**2.))
                                    endif
                                 else
                                    xwallu=xnorm_bisect(aa(ibuild),bb(ibuild),&
                                          gamma(ibuild)-upwind_dir,yu,thetamin,thetamax,dxy)
                                 endif
                                 xpos=xu-xwallu
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    uo(iu,ju,k)=along_mag*cos(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                    if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                                       print*,'Parameterized U exceeds max in street canyon',&
                                          uo(iu,ju,k),max_velmag,iu,ju,k
                                    endif
                                 endif
! v component
                                 iv=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                 jv=nint(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                                 xp=(real(iv)-0.5)*dx-xco
                                 yp=real(jv-1)*dy-yco
                                 xv=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yv=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(circle_flag .eq. 1)then
                                    if(abs(yv) .gt. Lt(ibuild))then
                                       cycle
                                    else
                                       xwallv=sqrt((Lt(ibuild)**2.)-(yv**2.))
                                    endif
                                 else
                                    xwallv=xnorm_bisect(aa(ibuild),bb(ibuild),&
                                          gamma(ibuild)-upwind_dir,yv,thetamin,thetamax,dxy)
                                 endif
                                 xpos=xv-xwallv
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    vo(iv,jv,k)=along_mag*sin(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                    if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                                       print*,'Parameterized V exceeds max in street canyon',&
                                          vo(iv,jv,k),max_velmag,iv,jv,k
                                    endif
                                 endif
! w component
                                 xp=(real(i)-0.5)*dx-xco
                                 yp=(real(j)-0.5)*dy-yco
                                 xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(circle_flag .eq. 1)then
                                    if(abs(yw) .gt. Lt(ibuild))then
                                       cycle
                                    else
                                       xwallw=sqrt((Lt(ibuild)**2.)-(yw**2.))
                                    endif
                                 else
                                    xwallw=xnorm_bisect(aa(ibuild),bb(ibuild),&
                                          gamma(ibuild)-upwind_dir,yw,thetamin,thetamax,dxy)
                                 endif
                                 xpos=xw-xwallw
                                 if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                    if(icellflag(i,j,k-1) .ne. 0)then
                                       if(reverse_flag .eq. 0)then
                                          wo(i,j,k)=-abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       else
                                          wo(i,j,k)=abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       endif
                                       if(abs(wo(i,j,k)) .gt. max_velmag)then
                                          print*,'Parameterized W exceeds max in street canyon',&
                                             wo(i,j,k),max_velmag,i,j,k
                                       endif
                                    endif
                                    icellflag(i,j,k)=6
                                 endif
                              endif
                           enddo
                        endif
                     enddo
                  enddo
               case(6)
                  tol=0.01*pi/180.
                  do ivert=bldstartidx(ibuild),bldstopidx(ibuild)
                     xw1=(bldx(ivert)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*sin(upwind_dir)
                     yw1=-(bldx(ivert)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*cos(upwind_dir)
                     xw3=(bldx(ivert+1)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*sin(upwind_dir)
                     yw3=-(bldx(ivert+1)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*cos(upwind_dir)
                     upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
                     if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
                     if(abs(upwind_rel) .lt. 0.5*pi)then
                        if(abs(upwind_rel) .gt. pi-tol .or. abs(upwind_rel) .lt. tol)then
                           perpendicular_flag=1
                        else
                           perpendicular_flag=0
                        endif
                        yw2=min(yw1,yw3)
                        upwind_norm=atan2(bldy(ivert+1)-bldy(ivert),bldx(ivert+1)-bldx(ivert))+0.5*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                        do y_idx=0,2*ceiling(abs(yw1-yw3)/dxy)
                           yc=0.5*real(y_idx)*dxy+yw2
                           top_flag=0
                           if(perpendicular_flag .gt. 0)then
                              xwall=xw1
                           else
                              xwall=((xw3-xw1)/(yw3-yw1))*(yc-yw1)+xw1
                           endif
                           do k=kend(ibuild),2,-1
                              canyon_flag=0
                              S=0.
                              x_idx_min=-1
                              reverse_flag=0
                              do x_idx=1,2*ceiling(Lr(ibuild)/dxy)
                                 xc=0.5*real(x_idx)*dxy
                                 i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                 j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                 if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                                    exit
                                 endif
                                 if(icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                                    x_idx_min=x_idx
                                 endif
                                 if(icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                                    canyon_flag=1
                                    x_idx_max=x_idx-1
                                    S=0.5*real(x_idx_max-x_idx_min)*dxy
                                    if(top_flag .eq. 0)then
                                       k_ref=k+1
                                       ic=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                       jc=ceiling(((0.5*real(x_idx_max)*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                       if(icellflag(ic,jc,k_ref) .ne. 0)then
                                          numu=0.
                                          numv=0.
                                          ucomponent=0.
                                          vcomponent=0.
                                          if(icellflag(ic-1,jc,k_ref) .ne. 0)then
                                             numu=numu+1.
                                             ucomponent=ucomponent+uo(ic,jc,k_ref)
                                          endif
                                          if(icellflag(ic+1,jc,k_ref) .ne. 0)then
                                             numu=numu+1.
                                             ucomponent=ucomponent+uo(ic+1,jc,k_ref)
                                          endif
                                          if(icellflag(ic,jc-1,k_ref) .ne. 0)then
                                             numv=numv+1.
                                             vcomponent=vcomponent+vo(ic,jc,k_ref)
                                          endif
                                          if(icellflag(ic,jc+1,k_ref) .ne. 0)then
                                             numv=numv+1.
                                             vcomponent=vcomponent+vo(ic,jc+1,k_ref)
                                          endif
                                          if( ucomponent .ne. 0. .and. numu > 0.) then
                                             ucomponent=ucomponent/numu
                                          else
                                             ucomponent=0.
                                          endif
                                          if( vcomponent .ne. 0. .and. numv > 0.) then
                                             vcomponent=vcomponent/numv
                                          else
                                             vcomponent=0.
                                          endif
                                          if(numu == 0. .and. numv == 0.)then
                                             canyon_flag=0
                                             top_flag=0
                                             S=0.
                                             exit
                                          elseif(numu > 0. .and. numv > 0.)then
                                             velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                             canyon_dir=atan2(vcomponent,ucomponent)
                                          elseif(numu > 0.)then
                                             velmag=abs(ucomponent)
                                             if(ucomponent .gt. 0.)then
                                                canyon_dir=0.
                                             else
                                                canyon_dir=pi
                                             endif
                                          else
                                             velmag=abs(vcomponent)
                                             if(vcomponent .gt. 0.)then
                                                canyon_dir=0.5*pi
                                             else
                                                canyon_dir=-0.5*pi
                                             endif
                                          endif
                                          top_flag=1
                                          if(abs(S) .gt. 0.) then
                                             if(Ht(ibldflag(i,j,k)) .lt. Ht(ibuild) .and. z(k)/S .lt. 0.65)then
                                               canyon_flag=0
                                               top_flag=0
                                               S=0.
                                               exit
                                             endif
                                          endif
                                       else
                                          canyon_flag=0
                                          top_flag=0
                                          S=0.
                                          exit
                                       endif
                                       if(velmag .gt. max_velmag)then
                                          print*,'Parameterized velocity exceeds max in street canyon',&
                                             velmag,max_velmag,i,j,k
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                          exit
                                       endif
                                    endif
                                    dbuild=ibldflag(i,j,k)
! Find the along canyon and cross canyon directions
                                    i=ceiling(((xc-0.5*dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                    j=ceiling(((xc-0.5*dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                    if(bldtype(dbuild) .eq. 4)then
                                       canyon_flag=0
                                       S=0.
                                       top_flag=0
                                       exit
                                    endif
                                    select case(bldgeometry(dbuild))
                                       case(1,4)
                                          beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                          downwind_rel=canyon_dir-gamma(dbuild)
                                          if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                          if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                          xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                          ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                          xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                             ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                          yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                             ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                          thetad=atan2(yd,xd)
                                          if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                             if(downwind_rel .le. 0.)then
                                                if(downwind_rel .le. -0.5*pi)then
                                                   along_dir=-pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.5*pi)then
                                                   along_dir=-pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             endif
                                          elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                             if(abs(downwind_rel) .ge. 0.5*pi)then
                                                if(downwind_rel .lt. 0)then
                                                   along_dir=-0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .lt. 0)then
                                                   along_dir=-0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             endif
                                          elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                             if(downwind_rel .ge. 0.)then
                                                if(downwind_rel .le. 0.5*pi)then
                                                   along_dir=gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=-pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.5*pi)then
                                                   along_dir=gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=-pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             endif
                                          else
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=-0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=-0.5*pi+gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             endif
                                          endif
                                       case(2,5)
                                          downwind_rel=canyon_dir-gamma(dbuild)
                                          if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                          if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                          xcd=xfo(dbuild)+Lt(dbuild)*cos(gamma(dbuild))
                                          ycd=yfo(dbuild)+Lt(dbuild)*sin(gamma(dbuild))
                                          xd=((real(i)-0.5)*dx-xcd)*cos(gamma(dbuild))+&
                                             ((real(j)-0.5)*dy-ycd)*sin(gamma(dbuild))
                                          yd=-((real(i)-0.5)*dx-xcd)*sin(gamma(dbuild))+&
                                             ((real(j)-0.5)*dy-ycd)*cos(gamma(dbuild))
                                          thetad=atan2(yd,xd)
                                          rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                             ((Wt(dbuild)*cos(thetad))**2.))
                                          along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                             (Lt(dbuild)**2.)*rd*sin(thetad))
                                          if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                             if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                                along_dir=along_dir+gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=along_dir-pi+gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                                along_dir=along_dir+gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=along_dir-pi+gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          endif
                                       case(6)
                                          do jvert=bldstartidx(dbuild),bldstopidx(dbuild)
                                             cross_dir=atan2(bldy(jvert+1)-bldy(jvert),bldx(jvert+1)-bldx(jvert))+0.5*pi
                                             if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                             xcd=0.5*(bldx(jvert+1)+bldx(jvert))
                                             ycd=0.5*(bldy(jvert+1)+bldy(jvert))
                                             xd=((real(i)-0.5)*dx-xcd)*cos(cross_dir)+&
                                                ((real(j)-0.5)*dy-ycd)*sin(cross_dir)
                                             yd=-((real(i)-0.5)*dx-xcd)*sin(cross_dir)+&
                                                ((real(j)-0.5)*dy-ycd)*cos(cross_dir)
                                             if(abs(xd) .lt. 0.75*dxy)then
                                                segmentLength=sqrt(((bldy(jvert+1)-bldy(jvert))**2.) &
                                                   +((bldx(jvert+1)-bldx(jvert))**2.))
                                                if(abs(yd) .le. 0.5*segmentLength)then
                                                   downwind_rel=canyon_dir-cross_dir
                                                   if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                                   if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                                   if(abs(downwind_rel) .lt. 0.5*pi)then
                                                      reverse_flag=1
                                                      if(downwind_rel .ge. 0.)then
                                                         along_dir=cross_dir-0.5*pi
                                                      else
                                                         along_dir=cross_dir+0.5*pi
                                                      endif
                                                   else
                                                      reverse_flag=0
                                                      if(downwind_rel .ge. 0.)then
                                                         along_dir=cross_dir+0.5*pi
                                                      else
                                                         along_dir=cross_dir-0.5*pi
                                                      endif
                                                   endif
                                                   if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                                   if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                                   exit
                                                endif
                                             endif
                                             if(bldx(jvert+1) .eq. bldx(bldstartidx(dbuild)) &
                                                   .and. bldy(jvert+1) .eq. bldy(bldstartidx(dbuild)))exit
                                          enddo
                                    end select
                                    if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                    if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                    if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                                    if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                                    if(reverse_flag .eq. 1)then
                                       if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                       endif
                                    else
                                       if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                       endif
                                    endif
                                    exit
                                 endif
                              enddo
                              if(canyon_flag .eq. 1 .and. S .gt. 0.9*dxy)then
                                 along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                                    log(zm(k)/zo)/log(zm(k_ref)/zo)
                                 cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                                 if(abs(along_mag) .gt. max_velmag)then
                                    print*,'Along canyon exceeds max in street canyon',&
                                       along_mag,max_velmag,i,j,k
                                 endif
                                 if(abs(cross_mag) .gt. max_velmag)then
                                    print*,'Cross canyon exceeds max in street canyon',&
                                       cross_mag,max_velmag,i,j,k
                                 endif
                                 do x_idx=x_idx_min,x_idx_max
                                    xc=0.5*real(x_idx)*dxy
                                    i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                    j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                    
                                    if(icellflag(i,j,k) .ne. 0)then ! icellflag(i,j,k) .ne. 6 .and. 
! u component
                                       iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                                       ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                                       xp=real(iu-1)*dx-xco
                                       yp=(real(ju)-0.5)*dy-yco
                                       xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                       yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                       if(perpendicular_flag .gt. 0)then
                                          xwallu=xw1
                                       else
                                          xwallu=((xw3-xw1)/(yw3-yw1))*(yu-yw1)+xw1
                                       endif
                                       xpos=xu-xwallu
                                       if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                          uo(iu,ju,k)=along_mag*cos(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                          if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                                             print*,'Parameterized U exceeds max in street canyon',&
                                                uo(iu,ju,k),velmag,max_velmag,iu,ju,k
                                          endif
                                       endif
! v component
                                       iv=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                                       jv=nint(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                                       xp=(real(iv)-0.5)*dx-xco
                                       yp=real(jv-1)*dy-yco
                                       xv=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                       yv=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                       if(perpendicular_flag .gt. 0)then
                                          xwallv=xw1
                                       else
                                          xwallv=((xw3-xw1)/(yw3-yw1))*(yv-yw1)+xw1
                                       endif
                                       xpos=xv-xwallv
                                       if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                          vo(iv,jv,k)=along_mag*sin(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                          if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                                             print*,'Parameterized V exceeds max in street canyon',&
                                                vo(iv,jv,k),velmag,max_velmag,iv,jv,k
                                          endif
                                       endif
! w component
                                       xp=(real(i)-0.5)*dx-xco
                                       yp=(real(j)-0.5)*dy-yco
                                       xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                       yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                       if(perpendicular_flag .gt. 0)then
                                          xwallw=xw1
                                       else
                                          xwallw=((xw3-xw1)/(yw3-yw1))*(yw-yw1)+xw1
                                       endif
                                       xpos=xw-xwallw
                                       if(xpos .le. S .and. xpos .gt. -0.5*dxy)then
                                          if(icellflag(i,j,k-1) .ne. 0)then
                                             if(reverse_flag .eq. 0)then
                                                wo(i,j,k)=-abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             else
                                                wo(i,j,k)=abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             endif
                                             if(abs(wo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized W exceeds max in street canyon',&
                                                   wo(i,j,k),velmag,max_velmag,i,j,k
                                             endif
                                          endif
                                          icellflag(i,j,k)=6
                                       endif
                                    endif
                                 enddo
                              endif
                           enddo
                        enddo
                     endif
                     if(bldx(ivert+1) .eq. bldx(bldstartidx(ibuild)) &
                        .and. bldy(ivert+1) .eq. bldy(bldstartidx(ibuild)))exit
                  enddo
            end select
         enddo
         return
      end
