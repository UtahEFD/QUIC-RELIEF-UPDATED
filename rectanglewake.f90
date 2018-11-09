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
      subroutine rectanglewake
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
         real cav_fac,wake_fac,beta,LoverH,WoverH,upwind_rel_norm,eff_height
         real canyon_factor,xc,yc,dNu,dNv,xwall,xu,yu,xv,yv,xp,yp,xwallu,xwallv,xwallw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,iw,jw
         real vd,hd,Bs,BL,shell_height,xw,yw,dNw
         integer roof_perpendicular_flag,ns_flag
         integer ktop,kbottom,nupwind
         real LrRect(3),LrLocal,LrLocalu,LrLocalv,LrLocalw
         real epsilon
         
         epsilon = 10e-10
         
         if(bldgeometry(ibuild) .eq. 4 .and. bldroof(ibuild) .gt. 0)then
            eff_height=0.8*(Ht(ibuild)-zfo_actual(ibuild))+zfo_actual(ibuild)
         else
            eff_height=Ht(ibuild)
         endif
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
         LoverH=Leff(ibuild)/eff_height
         WoverH=Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         if((bldtype(ibuild) .eq. 1 .or. bldtype(ibuild) .eq. 10) .and. roofflag .eq. 2)then
            tol=30*pi/180.
            roof_perpendicular_flag=0
            if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
               roof_perpendicular_flag=0
            elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
               roof_perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
               roof_perpendicular_flag=0
            elseif(abs(upwind_rel) .le. tol)then
               roof_perpendicular_flag=1
               ns_flag=0
            elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
               roof_perpendicular_flag=0
            elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
               roof_perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
               roof_perpendicular_flag=0
            else
               roof_perpendicular_flag=1
               ns_flag=0
            endif
            rooftop_flag(ibuild)=1
            k=kend(ibuild)
            nupwind=0
            do y_idx=1,2*int((yw1-yw3)/dxy)-1
               yc=0.5*real(y_idx)*dxy+yw3
               if(perpendicular_flag .gt. 0)then
                  xwall=xf2
               elseif(yc.ge.yf2)then
                  xwall=((xf2-xw1)/(yf2-yw1))*(yc-yw1)+xw1
               else
                  xwall=((xw3-xf2)/(yw3-yf2))*(yc-yf2)+xf2
               endif
               if(yc.ge.0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif
               do x_idx=int(Lr(ibuild)/dxy)+1,1,-1
                  xc=-real(x_idx)*dxy
                  i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                  j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                  if(i .ge. 1 .and. i .le. nx-1 .and. j .ge. 1 .and. j .le. ny-1)then
                     !print*,i,j,k,icellflag(i,j,k)
                     if(icellflag(i,j,k) .eq. 0)then
                        nupwind=nupwind+1
                        exit
                     endif
                  endif
               enddo
            enddo
            if(nupwind .ge. int((yw1-yw3)/dxy))rooftop_flag(ibuild)=0
            !print*,ibuild,xfo(ibuild),rooftop_flag(ibuild),nupwind,int((yw1-yw3)/dxy),Lr(ibuild),xwall
            if(roof_perpendicular_flag .eq. 1 .and. rooftop_flag(ibuild) .eq. 1)then
               if(ns_flag .eq. 1)then
                  hd=Wti(ibuild)
               else
                  hd=Lti(ibuild)
               endif
               Bs=min(Weff(ibuild),Ht(ibuild))
               BL=max(Weff(ibuild),Ht(ibuild))
               Rscale(ibuild) = ((Bs**(2./3.))*(BL**(1./3.)))
               Rcx(ibuild)=(0.9*Rscale(ibuild))
               vd= 0.5*0.22*Rscale(ibuild)
               if(hd .lt. Rcx(ibuild))then
                  shell_height=vd*sqrt(1-((0.5*Rcx(ibuild)-hd)/(0.5*Rcx(ibuild)))**2.)
                  if(shell_height .gt. 0)then
                     eff_height=eff_height+shell_height
                  endif
               endif
            endif
         endif
         LoverH=Leff(ibuild)/eff_height
         WoverH=Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         tol=0.01*pi/180.
         if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
            LrRect(1)=Lr(ibuild)*abs(cos(upwind_rel))
            LrRect(3)=Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
         elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
            LrRect(1)=Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
            LrRect(3)=Lr(ibuild)*abs(cos(upwind_rel))
         elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
            LrRect(1)=Lr(ibuild)*abs(cos(upwind_rel))
            LrRect(3)=Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
         elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
            LrRect(1)=Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
            LrRect(3)=Lr(ibuild)*abs(cos(upwind_rel))
         endif
         if(perpendicular_flag .eq. 0)then
            LrRect(2)=((yw1-yw2)*LrRect(1)+(yw2-yw3)*LrRect(3))/(yw1-yw3)
            Lr(ibuild)=((yw1-yw2)*0.5*(LrRect(1)+LrRect(2))+(yw2-yw3)*0.5*(LrRect(2)+LrRect(3)))/(yw1-yw3)
         endif
         if(zfo_actual(ibuild) .gt. 0)then
            call building_connect
         endif
         ! MAN 07/25/2008 stretched vertical grid
         do k=2,kstart(ibuild)
            kbottom=k
            if(zfo(ibuild) .le. zm(k))exit
         enddo
         do k=kstart(ibuild),nz-1
            ktop=k
            if(eff_height .lt. zm(k+1))exit
         enddo
         do k=kstart(ibuild),kend(ibuild)
            kk=k
            if(0.75*(Ht(ibuild)-zfo_actual(ibuild))+zfo_actual(ibuild) .le. zm(k))exit
         enddo
lp003:   do k=ktop,kbottom,-1
            zb=zm(k)-zfo(ibuild)
            ufarwake(:,:)=uo(:,:,k)
            vfarwake(:,:)=vo(:,:,k)
            !$omp parallel do private(yc,LrLocal,xwall,ynorm,canyon_factor,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,LrLocalu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,LrLocalv,xwallv,dNv,iw,jw,xw,yw,LrLocalw,xwallw,dNw,wwakeflag)
lp002:      do y_idx=1,2*int((yw1-yw3)/dxy)-1
               yc=0.5*real(y_idx)*dxy+yw3
               if(perpendicular_flag .gt. 0)then
                  xwall=xw1
                  LrLocal=Lr(ibuild)
               elseif(yc .ge. yw2)then
                  xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
                  LrLocal=LrRect(1)+(yc-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
               else
                  xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
                  LrLocal=LrRect(2)+(yc-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
               endif
               if(yc .ge. 0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif
               !check for building that will disrupt the wake
               canyon_factor=1.
               x_idx_min=-1
               do x_idx=1,int(Lr(ibuild)/dxy)+1
                  xc=real(x_idx)*dxy
                  i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                  j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                  if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(icellflag(i,j,kk) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(icellflag(i,j,kk) .eq. 0 .and. ibldflag(i,j,kk) .ne. ibuild .and. x_idx_min .gt. 0)then
                     canyon_factor=xc/Lr(ibuild)
                     exit
                  endif
               enddo
               x_idx_min=-1
lp001:         do x_idx=0,2*ceiling(farwake_factor*Lr(ibuild)/dxy)
                  uwakeflag=1
                  vwakeflag=1
                  wwakeflag=1
                  xc=0.5*real(x_idx)*dxy !
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
                        elseif(canyon_factor .lt. 1.)then
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
                        LrLocalu=Lr(ibuild)
                     elseif(yu.ge.yw2)then
                        xwallu=((xw2-xw1)/(yw2-yw1))*(yu-yw1)+xw1
                        LrLocalu=LrRect(1)+(yu-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallu=((xw3-xw2)/(yw3-yw2))*(yu-yw2)+xw2
                        LrLocalu=LrRect(2)+(yu-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xu=xu-xwallu
                     if(abs(yu) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNu=sqrt((1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalu)**2)
                     else
                     	dNu = 0.
                     endif
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then
                           farwake_velocity=ufarwake(iu,ju)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then ! .and. icellflag(iu,ju,k) .ne. 4
                              uo(iu,ju,k)=farwake_velocity
                              wo(i,j,k)=0.
                              ! icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
						         if(abs(uo(iu,ju,k)) .gt. max_velmag)then
							        print*,'Parameterized U exceeds max in rectangle wake',&
							      	 uo(iu,ju,k),max_velmag,iu,ju,k
						         endif
						         wo(i,j,k)=0.
						         ! icellflag(i,j,k)=4
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
                        LrLocalv=Lr(ibuild)
                     elseif(yv.ge.yw2)then
                        xwallv=((xw2-xw1)/(yw2-yw1))*(yv-yw1)+xw1
                        LrLocalv=LrRect(1)+(yv-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallv=((xw3-xw2)/(yw3-yw2))*(yv-yw2)+xw2
                        LrLocalv=LrRect(2)+(yv-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xv=xv-xwallv
                     if(abs(yv) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNv=sqrt((1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalv)**2)
                     else 
                     	dNv = 0.
                     endif
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           farwake_velocity=vfarwake(iv,jv)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then ! .and. icellflag(iv,jv,k) .ne. 4
                              vo(iv,jv,k)=farwake_velocity
                              wo(i,j,k)=0.
                              ! icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                              print*,'Parameterized V exceeds max in rectangle wake',&
                                 vo(iv,jv,k),max_velmag,iv,jv,k
                           endif
                           wo(iv,jv,k)=0.
                           ! icellflag(i,j,k)=4
                        endif  
                     endif
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                     xp=(real(iw)-0.5)*dx-xco
                     yp=(real(jw)-0.5)*dy-yco
                     xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                     yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                     if(perpendicular_flag .gt. 0)then
                        xwallw=xw1
                        LrLocalw=Lr(ibuild)
                     elseif(yw .ge. yw2)then
                        xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                        LrLocalw=LrRect(1)+(yv-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                        LrLocalw=LrRect(2)+(yv-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xw=xw-xwallw
                     if(abs(yw) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNw=sqrt((1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalw)**2)
                     else 
                     	dNw = 0.
                     endif
                     if(xw .gt. farwake_factor*dNw)wwakeflag=0
                     if(dNw .gt. 0. .and. wwakeflag .eq. 1 .and. icellflag(iw,jw,k) .ne. 0)then
                        if(xw .gt. dNw)then
                           if(canyon_factor .eq. 1.)then
                              if(icellflag(iw,jw,k) .eq. 4)then
                                 icellflag(iw,jw,k)=12
                              else
                                 icellflag(iw,jw,k)=5
                              endif
                           endif
! Cavity                   
                        else
                           icellflag(iw,jw,k)=4
                        endif  
                     endif
                  endif
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit !  .and. wwakeflag .eq. 0
               enddo   lp001      
            enddo   lp002
            !$omp end parallel do      
         enddo   lp003
         return
      end
