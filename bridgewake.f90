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
      subroutine bridgewake
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule ! make data from module "datamodule" visible
         implicit none
         integer perpendicular_flag,uwakeflag,vwakeflag,wwakeflag
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real xw1,yw1,xw2,yw2,xw3,yw3,xf2,yf2,tol,zb,ynorm
         real farwake_exponent,farwake_factor,farwake_velocity
         real cav_fac,wake_fac,beta,LoverH,upwind_rel_norm
         real bridge_thickness,xc,yc,dNu,dNv,dNw,xwall
         real xu,yu,xv,yv,xp,yp,xwallu,xwallv,xwallw,xw,yw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,iw,jw
         integer ktop,kbottom
         
         xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))!CENTER of building in QUIC domain coordinates
         yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
         ! find upwind direction and determine the type of flow regime
         uo_h=uo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         vo_h=vo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         upwind_rel=upwind_dir-gamma(ibuild)
         if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
         if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
         upwind_rel_norm=upwind_rel+0.5*pi
         if(upwind_rel_norm .gt. pi)upwind_rel_norm=upwind_rel_norm-2*pi

         tol=0.01*pi/180.
         farwake_exponent=1.5
         farwake_factor=3
         !Location of corners relative to the center of the building
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
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yf2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            perpendicular_flag=0
         elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
            xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
            xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yf2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            perpendicular_flag=0
         elseif(abs(upwind_rel) .le. tol)then
            xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
            xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xf2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yf2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            perpendicular_flag=0
         elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
            xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
            xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yf2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            perpendicular_flag=0
         else
            xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
         endif
         if(wakeflag .gt. 1)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         select case(wakeflag)
            case(1)
               Weff(ibuild)=abs(yw3-yw1)
               if(perpendicular_flag .eq. 1)then
                  Leff(ibuild)=abs(xf2-xw1)
               else
                  Leff(ibuild)=abs(xf2-xw2)
               endif
            case(2)
               beta=abs(atan2(Lti(ibuild),Wti(ibuild)))
               if(abs(upwind_rel) .gt. 0.5*pi-beta .and. &
                     abs(upwind_rel) .lt. 0.5*pi+beta)then
                  Leff(ibuild)=abs(Wti(ibuild)/sin(upwind_rel))
               else
                  Leff(ibuild)=abs(Lti(ibuild)/cos(upwind_rel))
               endif
               if(abs(upwind_rel_norm) .gt. 0.5*pi-beta .and. &
                     abs(upwind_rel_norm) .lt. 0.5*pi+beta)then
                  Weff(ibuild)=abs(Wti(ibuild)/sin(upwind_rel_norm))
               else
                  Weff(ibuild)=abs(Lti(ibuild)/cos(upwind_rel_norm))
               endif
            case(3)
               Leff(ibuild)=Wti(ibuild)*Lti(ibuild)/abs(yw3-yw1)
               if(perpendicular_flag .eq. 1)then
                  Weff(ibuild)=Wti(ibuild)*Lti(ibuild)/abs(xf2-xw1)
               else
                  Weff(ibuild)=Wti(ibuild)*Lti(ibuild)/abs(xf2-xw2)
               endif
         endselect
         
         bridge_thickness=0.5*(Ht(ibuild)-zfo(ibuild))
         LoverH=Leff(ibuild)/bridge_thickness
         if(LoverH.gt.3.)LoverH=3.
         if(LoverH.lt.0.3)LoverH=0.3
         Lr(ibuild)=0.9*Weff(ibuild)/((LoverH**(0.3))*   &
                        (1+0.24*Weff(ibuild)/bridge_thickness))
         
         do k=2,kstart(ibuild)
            kbottom=k
            if(zfo(ibuild) .le. zm(k))exit
         enddo
         do k=kstart(ibuild),nz-1
            ktop=k
            if(Ht(ibuild) .lt. zm(k+1))exit
         enddo
lp003:   do k=ktop,kbottom,-1
            zb=zm(k)-(bridge_thickness+zfo(ibuild))
            ufarwake(:,:)=0
            vfarwake(:,:)=0
            !$omp parallel do private(yc,xwall,ynorm,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,xwallv,dNv,iw,jw,xw,yw,xwallw,dNw,wwakeflag)
lp002:      do y_idx=1,2*int((yw1-yw3)/dxy)-1
               yc=0.5*real(y_idx)*dxy+yw3
               if(perpendicular_flag .gt. 0)then
                  xwall=xw1
               elseif(yc.ge.yw2)then
                  xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
               else
                  xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
               endif
               if(yc .ge. 0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif    
               x_idx_min=-1
lp001:         do x_idx=0,2*ceiling(farwake_factor*Lr(ibuild)/dxy)
                  uwakeflag=1
                  vwakeflag=1
                  wwakeflag=1
                  xc=0.5*real(x_idx)*dxy
                  i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                  j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                  if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(icellflag(i,j,k) .eq. 0)then
                     if(x_idx_min .ge. 0)then
                        if(ibldflag(i,j,k) .eq. ibuild)then
                           x_idx_min=-1
                        else
                           exit
                        endif
                     endif
                  endif
! u values
! Far wake
                  if(icellflag(i,j,k) .ne. 0 )then
                     iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                     ju=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
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
                     xu=xu-xwallu
                     dNu=sqrt((1.-(yu/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(Lr(ibuild))**2)
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .eq. dNu .and. uwakeflag .eq. 1 .and. icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then
                           farwake_velocity=ufarwake(iu,ju)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                           ! if(icellflag(iu,ju,k) .ne. 4)then
                              uo(iu,ju,k)=farwake_velocity
                              wo(i,j,k)=0.
                           ! endif
! Cavity                   
                        else
                           uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                              print*,'Parameterized U exceeds max in rectangle wake',&
                                 uo(iu,ju,k),max_velmag,iu,ju,k
                           endif
                           wo(i,j,k)=0.
                        endif
                     endif
! v values
! Far wake
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
                     xv=xv-xwallv
                     dNv=sqrt((1.-(yv/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(Lr(ibuild))**2)
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .eq. dNv .and. vwakeflag .eq. 1 .and. icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           farwake_velocity=vfarwake(iv,jv)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                           ! if(icellflag(iv,jv,k) .ne. 4)then
                              vo(iv,jv,k)=farwake_velocity
                              wo(i,j,k)=0.
                           ! endif
! Cavity                   
                        else
                           vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                              print*,'Parameterized V exceeds max in rectangle wake',&
                                 vo(iv,jv,k),max_velmag,iv,jv,k
                           endif
                           wo(iv,jv,k)=0.
                        endif  
                     endif
! check cell centers to mark cell flags
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                     xp=(real(iw)-0.5)*dx-xco
                     yp=(real(jw)-0.5)*dy-yco
                     xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                     yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                     if(perpendicular_flag .gt. 0)then
                        xwallw=xw1
                     elseif(yw.ge.yw2)then
                        xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                     else
                        xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                     endif
                     xw=xw-xwallw
                     dNw=sqrt((1.-(yw/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(Lr(ibuild))**2)
                     if(xw .gt. farwake_factor*dNw)wwakeflag=0
                     if(dNw .eq. dNw .and. wwakeflag .eq. 1 .and. icellflag(iw,jw,k) .ne. 0)then
                        if(xw .gt. dNw)then
                           if(icellflag(iw,jw,k) .eq. 4)then
                              icellflag(iw,jw,k)=12
                           else  
                              icellflag(iw,jw,k)=5
                           endif
! Cavity                   
                        else
                            icellflag(iw,jw,k)=4
                        endif  
                     endif
                  endif
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit
               enddo   lp001      
            enddo   lp002 
            !$omp end parallel do
         enddo   lp003
         return
      end
