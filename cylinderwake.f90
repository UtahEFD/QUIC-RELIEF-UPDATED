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
      subroutine cylinderwake
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule ! make data from module "datamodule" visible
         implicit none
         integer circle_flag,uwakeflag,vwakeflag,wwakeflag
         real uo_h,vo_h,upwind_dir,xco,yco,tol,zb
         real farwake_exponent,farwake_factor,farwake_velocity
         real thetamax,thetamin,thetai,LoverH,WoverH
         real ynorm,radius,xnorm_bisect,cav_fac,wake_fac,eff_height
         real canyon_factor,xc,yc,yw1,yw3,y1,y2,dNu,dNv,xwall
         real xu,yu,xv,yv,xp,yp,xwallu,xwallv
         real xw,yw,dNw,xwallw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,iw,jw
         integer ktop,kbottom
         real epsilon
         
         epsilon = 10e-10
         
         if(bldgeometry(ibuild) .eq. 5 .and. bldroof(ibuild) .gt. 0)then
            eff_height=0.8*(Ht(ibuild)-zfo_actual(ibuild))+zfo_actual(ibuild)
         else
            eff_height=Ht(ibuild)
         endif
         farwake_exponent=1.5
         farwake_factor=3
         xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))!CENTER of building in QUIC domain coordinates
         yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
         uo_h=uo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         vo_h=vo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
! find upwind direction and determine the type of flow regime
         upwind_dir=atan2(vo_h,uo_h)
         tol=0.01*dxy
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
            Weff(ibuild)=Lti(ibuild)
            Leff(ibuild)=Lti(ibuild)
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
            Weff(ibuild)=2*yw1
            Leff(ibuild)=2*radius(aa(ibuild),bb(ibuild),upwind_dir,gamma(ibuild))
         endif
         if(wakeflag .eq. 2)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         LoverH=Leff(ibuild)/eff_height
         WoverH=Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         Lr(ibuild)=0.9*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         ynorm=yw1
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
            !$omp parallel do private(yc,xwall,canyon_factor,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,xwallv,dNv,iw,jw,xw,yw,xwallw,dNw,wwakeflag)
lp002:      do y_idx=1,2*int((yw1-yw3)/dxy)-1
               yc=0.5*real(y_idx)*dxy+yw3
               if(circle_flag .eq. 1)then
                  if(abs(yc) .gt. Lt(ibuild))then
                     cycle
                  else
                     xwall=sqrt((Lt(ibuild)**2.)-(yc**2.))
                  endif
               else
                  xwall=xnorm_bisect(aa(ibuild),bb(ibuild),&
                        gamma(ibuild)-upwind_dir,yc,thetamin,thetamax,dxy)
               endif
               ! check for building that will disrupt the wake
               canyon_factor=1.
               x_idx_min=-1
               do x_idx=1,ceiling(Lr(ibuild)/dxy)
                  xc=real(x_idx)*dxy
                  i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                  j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
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
                     xu=xu-xwallu

                     
                     if(ynorm > epsilon .and. abs(yu) < ynorm &
                           .and. eff_height > epsilon  .and. zb < eff_height)then
                        dNu = (1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*Lr(ibuild))**2
                        dNu=sqrt(dNu)
                     else
                        dNu = 0.0
                     endif
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then
                           farwake_velocity=ufarwake(iu,ju)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then ! .and. icellflag(iu,ju,k) .ne. 4
                              uo(iu,ju,k)=farwake_velocity
                              wo(i,j,k)=0.
                              ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                              print*,'Parameterized U exceeds max in cylinder wake',&
                                 uo(iu,ju,k),max_velmag,iu,ju,k
                           endif
                           wo(i,j,k)=0.
                           ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=4
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
                     xv=xv-xwallv
                     
                     if(ynorm > epsilon .and. abs(yv) < ynorm &
                           .and. eff_height > epsilon .and. zb < eff_height)then
                        dNv = (1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*Lr(ibuild))**2
                        dNv=sqrt(dNv)
                     else
                        dNv = 0.0
                     endif
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           farwake_velocity=vfarwake(iv,jv)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then !  .and. icellflag(iv,jv,k) .ne. 4
                              vo(iv,jv,k)=farwake_velocity
                              wo(i,j,k)=0.
                              ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                              print*,'Parameterized V exceeds max in cylinder wake',&
                                 vo(iv,jv,k),max_velmag,iv,jv,k
                           endif
                           wo(iv,jv,k)=0.
                           ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=4
                        endif  
                     endif
! celltype values
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)
                     xp=(real(iw)-0.5)*dx-xco
                     yp=(real(jw)-0.5)*dy-yco
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
                     xw=xw-xwallw
                     if(ynorm > epsilon .and. abs(yw) < ynorm &
                           .and. eff_height > epsilon .and. zb < eff_height)then
                        dNw = (1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*Lr(ibuild))**2
                        dNw=sqrt(dNw)
                     else
                        dNw = 0.0
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
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit
               enddo   lp001      
            enddo   lp002
            !$omp end parallel do
         enddo   lp003
         return
      end
      
      
      real function radius(a,b,gamma,theta)
         implicit none
         real a,b,gamma,theta
         radius=a*b/sqrt( (a*sin(theta-gamma))**2. + (b*cos(theta-gamma))**2. )
         return
      end
      
      
      real function xnorm_bisect(a,b,gamma,y,thetamin,thetamax,dxy)
         implicit none
         integer i
         real, intent(in) :: a,b,gamma,y,dxy,thetamin,thetamax
         real yguess,yguess_low,eps,prod,theta,rad,radius,thetalow,thetahigh
         i=0
         eps=a+b
         thetalow=thetamin
         thetahigh=thetamax
         do while (i .lt. 100 .and. eps .gt. 0.1*dxy)
            i=i+1
            theta=0.5*(thetalow+thetahigh)
            rad=radius(a,b,gamma,theta)
            yguess=rad*sin(theta)-y
            yguess_low=radius(a,b,gamma,thetalow)*sin(thetalow)-y
            eps=abs(yguess)
            prod=yguess*yguess_low
            if(prod .lt. 0)then
               thetahigh=theta
            elseif(prod .gt. 0)then
               thetalow=theta
            else
               eps=0.
            endif
         enddo
         xnorm_bisect=rad*cos(theta)
         return
      end