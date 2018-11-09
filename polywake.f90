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
      subroutine polywake
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule ! make data from module "datamodule" visible
         implicit none
         integer perpendicular_flag,uwakeflag,vwakeflag,wwakeflag
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,x2,polyarea
         real y1,y2,y3,y4,xw1,yw1,xw3,yw3,tol,zb,ynorm,ynormp,ynormm
         real farwake_exponent,farwake_factor,farwake_velocity,LrAve,seglength,totseglength
         real cav_fac,wake_fac,LoverH,WoverH,eff_height,LrLocal,LrLocalu,LrLocalv
         real canyon_factor,xc,yc,dNu,dNv,xwall,xu,yu,xv,yv,xp,yp,xwallu,xwallv
         real xw,yw,dNw,xwallw,LrLocalw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,stop_idx,iw,jw
         integer ktop,kbottom,ivert,endFace,lastFace,nextFace, test1, test2
         real epsilon
         
         epsilon = 10e-10
         
         
         eff_height=Ht(ibuild)
         xco = bldcx(ibuild)!CENTER of building in QUIC domain coordinates
         yco = bldcy(ibuild)
         ! find upwind direction and determine the type of flow regime
         uo_h=uo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         vo_h=vo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         tol=0.01*pi/180.
         farwake_exponent=1.5
         farwake_factor=3
         x1=0.
         x2=0.
         y1=0.
         y2=0.
         polyarea=0.
         stop_idx = 0
         do j=bldstartidx(ibuild),bldstopidx(ibuild)
            polyarea=polyarea+bldx(j)*bldy(j+1)-bldx(j+1)*bldy(j)
            xp=(bldx(j)-bldcx(ibuild))*cos(upwind_dir)+(bldy(j)-bldcy(ibuild))*sin(upwind_dir)
            yp=-(bldx(j)-bldcx(ibuild))*sin(upwind_dir)+(bldy(j)-bldcy(ibuild))*cos(upwind_dir)
            if(xp .lt. x1)x1=xp
            if(xp .gt. x2)x2=xp
            if(yp .lt. y1)y1=yp
            if(yp .gt. y2)y2=yp
            if(bldx(j+1) .eq. bldx(bldstartidx(ibuild)) .and. bldy(j+1) .eq. bldy(bldstartidx(ibuild)))exit
         enddo
         polyarea=0.5*abs(polyarea)
         Weff(ibuild)=polyarea/(x2-x1)
         Leff(ibuild)=polyarea/(y2-y1)
         if(wakeflag .gt. 1)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         LrNode(:)=0.
         LrFace(:)=-1.
         ynormp=0.
         ynormm=0.
         do ivert=bldstartidx(ibuild),bldstopidx(ibuild)
            y1=-(bldx(ivert)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*cos(upwind_dir)
            if(y1 .lt. ynormm)ynormm=y1
            if(y1 .gt. ynormp)ynormp=y1
            if(bldx(ivert+1) .eq. bldx(bldstartidx(ibuild)) &
                  .and. bldy(ivert+1) .eq. bldy(bldstartidx(ibuild)))exit
         enddo
         LoverH=Leff(ibuild)/eff_height
         WoverH=Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         do ivert=bldstartidx(ibuild),bldstopidx(ibuild)
            xw1=(bldx(ivert)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*sin(upwind_dir)
            yw1=-(bldx(ivert)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*cos(upwind_dir)
            xw3=(bldx(ivert+1)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*sin(upwind_dir)
            yw3=-(bldx(ivert+1)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*cos(upwind_dir)
            upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(abs(upwind_rel) .lt. 0.5*pi)then
               LrFace(ivert)=Lr(ibuild)*cos(upwind_rel)
            endif
            if(bldx(ivert+1) .eq. bldx(bldstartidx(ibuild)) &
                  .and. bldy(ivert+1) .eq. bldy(bldstartidx(ibuild)))then
               endFace=ivert
               exit
            endif
         enddo
         LrAve=0.
         totseglength=0.
         do ivert=bldstartidx(ibuild),bldstopidx(ibuild)
            if(LrFace(ivert) .gt. 0.)then
               if(ivert .eq. 1)then
                  lastFace=endFace
                  nextFace=ivert+1
               elseif(ivert .eq. endFace)then
                  lastFace=ivert-1
                  nextFace=bldstartidx(ibuild)
               else
                  lastFace=ivert-1
                  nextFace=ivert+1
               endif
               y1=-(bldx(lastFace)-bldcx(ibuild))*sin(upwind_dir)+(bldy(lastFace)-bldcy(ibuild))*cos(upwind_dir)
               y2=-(bldx(ivert)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*cos(upwind_dir)
               y3=-(bldx(nextFace)-bldcx(ibuild))*sin(upwind_dir)+(bldy(nextFace)-bldcy(ibuild))*cos(upwind_dir)
               y4=-(bldx(nextFace+1)-bldcx(ibuild))*sin(upwind_dir)+(bldy(nextFace+1)-bldcy(ibuild))*cos(upwind_dir)
               if(LrFace(lastFace) .lt. 0. .and. LrFace(nextFace) .lt. 0.)then
                  LrNode(ivert)=LrFace(ivert)
                  LrNode(ivert+1)=LrFace(ivert)
               elseif(LrFace(lastFace) .lt. 0.)then
                  LrNode(ivert)=LrFace(ivert)
                  LrNode(ivert+1)=((y3-y4)*LrFace(nextFace)+(y2-y3)*LrFace(ivert))/(y2-y4)
               elseif(LrFace(nextFace) .lt. 0.)then
                  LrNode(ivert)=((y2-y3)*LrFace(ivert)+(y1-y2)*LrFace(lastFace))/(y1-y3)
                  LrNode(ivert+1)=LrFace(ivert)
               else
                  LrNode(ivert)=((y2-y3)*LrFace(ivert)+(y1-y2)*LrFace(lastFace))/(y1-y3)
                  LrNode(ivert+1)=((y3-y4)*LrFace(nextFace)+(y2-y3)*LrFace(ivert))/(y2-y4)
               endif
               seglength=y2-y3
               LrAve=LrAve+LrFace(ivert)*seglength
               totseglength=totseglength+seglength
            endif
            
            if(bldx(ivert+1) .gt. bldx(bldstartidx(ibuild)) - .1 .and. bldx(ivert+1) .lt. bldx(bldstartidx(ibuild)) + .1 &
            .and. bldy(ivert+1) .gt. bldy(bldstartidx(ibuild)) - .1 .and. bldy(ivert+1) .lt. bldy(bldstartidx(ibuild)) + .1)then
               stop_idx=ivert

               exit
            endif
         enddo
         Lr(ibuild)=LrAve/totseglength
         if(zfo_actual(ibuild) .gt. 0)then
            call building_connect
         endif
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
         xco=bldcx(ibuild)
         yco=bldcy(ibuild)
lp004:   do k=ktop,kbottom,-1
            zb=zm(k)-zfo(ibuild)
            ufarwake(:,:)=uo(:,:,k)
            vfarwake(:,:)=vo(:,:,k)
            test1 = ibuild
            test2 = bldstartidx(ibuild)
            !!$omp parallel do private(xw1,yw1,xw3,yw3,upwind_rel,perpendicular_flag, &
            !!$omp y_idx,yc,LrLocal,xwall,ynorm,canyon_factor,x_idx_min,x_idx, &
            !!$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,LrLocalu,xwallu,dNu,farwake_velocity, &
            !!$omp iv,jv,xv,yv,LrLocalv,xwallv,dNv) &
            !!$omp firstprivate(ibuild)
            
            
            
lp003:      do ivert=bldstartidx(ibuild),stop_idx
               xw1=(bldx(ivert)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*sin(upwind_dir)
               yw1=-(bldx(ivert)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert)-bldcy(ibuild))*cos(upwind_dir)
               xw3=(bldx(ivert+1)-bldcx(ibuild))*cos(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*sin(upwind_dir)
               yw3=-(bldx(ivert+1)-bldcx(ibuild))*sin(upwind_dir)+(bldy(ivert+1)-bldcy(ibuild))*cos(upwind_dir)
               upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
               if(upwind_rel .gt. pi)upwind_rel=upwind_rel-2*pi
               if(abs(upwind_rel) .lt. 0.5*pi)then
                  if(abs(upwind_rel) .lt. tol)then
                     perpendicular_flag=1
                  else
                     perpendicular_flag=0
                  endif
lp002:            do y_idx=0,2*ceiling(abs(yw1-yw3)/dxy)
                     yc=yw1-0.5*real(y_idx)*dxy
                     LrLocal=LrNode(ivert)+(yc-yw1)*(LrNode(ivert+1)-LrNode(ivert))/(yw3-yw1)
                     if(perpendicular_flag .gt. 0)then
                        xwall=xw1
                     else
                        xwall=((xw3-xw1)/(yw3-yw1))*(yc-yw1)+xw1
                     endif
                     if(yc .ge. 0.)then
                        ynorm=ynormp
                     else
                        ynorm=ynormm
                     endif
                     canyon_factor=1.
                     x_idx_min=-1
                     do x_idx=1,ceiling(LrLocal/dxy)
                        xc=real(x_idx)*dxy
                        i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                        j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                        if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                           exit
                        endif
                        if(icellflag(i,j,kk) .ne. 0 .and. x_idx_min .lt. 0)then
                           x_idx_min=x_idx
                        endif
                        if(icellflag(i,j,kk) .eq. 0 .and. x_idx_min .gt. 0)then
                           canyon_factor=xc/Lr(ibuild)
                           exit
                        endif
                     enddo
                     x_idx_min=-1
lp001:               do x_idx=1,2*ceiling(farwake_factor*LrLocal/dxy)
                        uwakeflag=1
                        vwakeflag=1
                        wwakeflag=1
                        xc=0.5*real(x_idx)*dxy
                        i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                        j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
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
                              elseif(icellflag(i,j,kk) .eq. 0)then
                                 exit
                              endif
                           endif
                        endif
! u values
! Far wake
                        if(icellflag(i,j,k) .ne. 0)then
                           iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/dx)+1
                           ju=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/dy)+1
                           xp=real(iu-1)*dx-xco
                           yp=(real(ju)-0.5)*dy-yco
                           xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                           yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                           LrLocalu=LrNode(ivert)+(yu-yw1)*(LrNode(ivert+1)-LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallu=xw1
                           else
                              xwallu=((xw3-xw1)/(yw3-yw1))*(yu-yw1)+xw1
                           endif
                           xu=xu-xwallu
                           
                           if(abs(yu) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                          	   dNu = (1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalu)**2
                          	   dNu=sqrt(dNu)
                           else
                           	 dNu=0.
                           endif
                           if(xu .gt. farwake_factor*dNu)uwakeflag=0
                           if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. yu .le. yw1 .and. yu .ge. yw3 .and. &
                                 icellflag(iu,ju,k) .ne. 0)then
                              if(xu .gt. dNu)then
                                 farwake_velocity=ufarwake(iu,ju)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                                 if(canyon_factor .eq. 1.)then ! icellflag(iu,ju,k) .ne. 4
                                    uo(iu,ju,k)=farwake_velocity
                                    wo(i,j,k)=0.
                                    ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=5
                                 endif
! Cavity                   
                              else
                                 uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 if(abs(uo(iu,ju,k)) .gt. max_velmag)then
                                    print*,'Parameterized U exceeds max in polygon wake',&
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
                           LrLocalv=LrNode(ivert)+(yv-yw1)*(LrNode(ivert+1)-LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallv=xw1
                           else
                              xwallv=((xw3-xw1)/(yw3-yw1))*(yv-yw1)+xw1
                           endif
                           xv=xv-xwallv
                           
                           if(abs(yv) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                           	dNv = (1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalv)**2
                           	dNv=sqrt(dNv)
                           else
                           	dNv=0.
                           endif
                           if(xv .gt. farwake_factor*dNv)vwakeflag=0
                           if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. yv .le. yw1 .and. yv .ge. yw3 .and. &
                                 icellflag(iv,jv,k) .ne. 0)then
                              if(xv .gt. dNv)then
                                 farwake_velocity=vfarwake(iv,jv)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                                 if(canyon_factor .eq. 1.)then ! .and. icellflag(iv,jv,k) .ne. 4
                                    vo(iv,jv,k)=farwake_velocity
                                    wo(i,j,k)=0.
                                    ! if(icellflag(i,j,k) .ne. 0)icellflag(i,j,k)=5
                                 endif
! Cavity                   
                              else
                                 vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                                 if(abs(vo(iv,jv,k)) .gt. max_velmag)then
                                    print*,'Parameterized V exceeds max in polygon wake',&
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
                           LrLocalw=LrNode(ivert)+(yw-yw1)*(LrNode(ivert+1)-LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallw=xw1
                           else
                              xwallw=((xw3-xw1)/(yw3-yw1))*(yw-yw1)+xw1
                           endif
                           xw=xw-xwallw
                           
                           if(abs(yw) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                           	dNw = (1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalw)**2
                           	dNw=sqrt(dNw)
                           else
                           	dNw=0.
                           endif
                           if(xw .gt. farwake_factor*dNw)wwakeflag=0
                           if(dNw .gt. 0. .and. wwakeflag .eq. 1 .and. yw .le. yw1 .and. yw .ge. yw3 .and. &
                                 icellflag(iw,jw,k) .ne. 0)then
                              if(xw .gt. dNw)then
                                 if(canyon_factor .eq. 1.)then
                                    if(icellflag(iv,jv,k) .eq. 4)then
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
               endif
            enddo   lp003
            !!$omp end parallel do
         enddo   lp004
         return
      end
