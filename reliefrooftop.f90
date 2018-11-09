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
      subroutine reliefrooftop
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MATR 04/02/18 - Based on Rooftop subroutine used for regular building

         use datamodule ! make data from module "datamodule" visible
         implicit none
         real perpendicular_flag,ns_flag
         integer roofflag_temp,uflag,vflag,wflag,k_ref,kendv
         integer ivert,k_shellu,k_shellv,k_shellw
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real tol,xfront,yfront
         real zr,x_u,y_u,x_v,y_v,x_w,y_w
         real vd,Bs,BL,roofangle,hx,hy,hdu,hdv,hdwx,hdwy
         real xnorm,ynorm,vel_mag,zref2,zref2u,zref2v
         real shell_heightu,shell_heightv,rotationAngle
         real sinRotationAngle,cosRotationAngle,vel_magu,vel_magv
         real cosUpwindDir,sinUpwindDir,denomu,denomv,invzo,invvd
         real hxu,hyu,hxv,hyv,hxw,hyw,tanRoofangle
         real cosXnorm,sinXnorm,cosXnormpPi,sinXnormpPi
         real cosYnorm,sinYnorm,cosYnormpPi,sinYnormpPi
         real shell_heightu_part, shell_heightv_part
         
         if(roofflag .eq. 0)then
            rooftop_flag(ibuild)=0
            return
         endif
         if(bldgeometry(ibuild) .eq. 6)then
            xco=bldcx(ibuild)
            yco=bldcy(ibuild)
         else
            xco = xfo(ibuild) + Lt(ibuild)*cos(gamma(ibuild))! CENTER of building in QUIC domain coordinates
            yco = yfo(ibuild) + Lt(ibuild)*sin(gamma(ibuild))
         endif
         ! find upwind direction and determine the type of flow regime
         uo_h=uo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         vo_h=vo(nint(xco/dx),nint(yco/dy),kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         cosUpwindDir=cos(upwind_dir)
         sinUpwindDir=sin(upwind_dir)
         tol=30*pi/180.
         invzo=1./zo
         if(bldgeometry(ibuild) .eq. 6)then
            perpendicular_flag=0
            xfront=0.
            yfront=0.
            do ivert=bldstartidx(ibuild),bldstopidx(ibuild)
               x1=(bldx(ivert)-bldcx(ibuild))*cosUpwindDir+(bldy(ivert)-bldcy(ibuild))*sinUpwindDir
               y1=-(bldx(ivert)-bldcx(ibuild))*sinUpwindDir+(bldy(ivert)-bldcy(ibuild))*cosUpwindDir
               if(x1 .lt. xfront)then
                  xfront=x1
                  yfront=y1
               endif
               if(bldx(ivert+1) .eq. bldx(bldstartidx(ibuild)) &
                     .and. bldy(ivert+1) .eq. bldy(bldstartidx(ibuild)))exit
            enddo
            ns_flag=0
            rotationAngle=upwind_dir
         else
            rotationAngle=gamma(ibuild)
            upwind_rel=upwind_dir-gamma(ibuild)
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
            ! Location of corners relative to the center of the building
            x1=xfo(ibuild)+Wt(ibuild)*sin(gamma(ibuild))-xco
            y1=yfo(ibuild)-Wt(ibuild)*cos(gamma(ibuild))-yco
            x2=x1+Lti(ibuild)*cos(gamma(ibuild))
            y2=y1+Lti(ibuild)*sin(gamma(ibuild))
            x4=xfo(ibuild)-Wt(ibuild)*sin(gamma(ibuild))-xco
            y4=yfo(ibuild)+Wt(ibuild)*cos(gamma(ibuild))-yco
            x3=x4+Lti(ibuild)*cos(gamma(ibuild))
            y3=y4+Lti(ibuild)*sin(gamma(ibuild))
            perpendicular_flag=0
            if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
               xfront=Lt(ibuild)
               yfront=-Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(0.5*pi-upwind_rel)-2*abs(0.75*pi-upwind_rel)))
               xnorm=gamma(ibuild)!+roofangle
               ynorm=gamma(ibuild)-0.5*pi!-roofangle
            elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
               xfront=Lt(ibuild)
               yfront=-Wt(ibuild)
               perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
               xfront=-Lt(ibuild)
               yfront=-Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(upwind_rel)-2*abs(0.25*pi-upwind_rel)))
               xnorm=gamma(ibuild)+pi!-roofangle
               ynorm=gamma(ibuild)-0.5*pi!+roofangle
            elseif(abs(upwind_rel) .le. tol)then
               xfront=-Lt(ibuild)
               yfront=-Wt(ibuild)
               perpendicular_flag=1
               ns_flag=0
            elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
               xfront=-Lt(ibuild)
               yfront=Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(upwind_rel)-2.0*abs(-0.25*pi-upwind_rel)))
               xnorm=gamma(ibuild)-pi!+roofangle
               ynorm=gamma(ibuild)+0.5*pi!-roofangle
            elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
               xfront=-Lt(ibuild)
               yfront=Wt(ibuild)
               perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
               xfront=Lt(ibuild)
               yfront=Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(-0.5*pi-upwind_rel)-2.0*abs(-0.75*pi-upwind_rel)))
               xnorm=gamma(ibuild)!-roofangle
               ynorm=gamma(ibuild)+0.5*pi!+roofangle
            else
               xfront=Lt(ibuild)
               yfront=Wt(ibuild)
               perpendicular_flag=1
               ns_flag=0
            endif
            if(perpendicular_flag .lt. 1)then
               tanRoofangle=tan(roofangle)
               cosXnorm=cos(xnorm)
               sinXnorm=sin(xnorm)
               cosXnormpPi=cos(xnorm+pi)
               sinXnormpPi=sin(xnorm+pi)
               cosYnorm=cos(ynorm)
               sinYnorm=sin(ynorm)
               cosYnormpPi=cos(ynorm+pi)
               sinYnormpPi=sin(ynorm+pi)
            endif
         endif
         sinRotationAngle=sin(rotationAngle)
         cosRotationAngle=cos(rotationAngle)
         ! MAN 07/25/2008 stretched vertical grid
         do k=kend(ibuild)+1,nz-1
            k_ref=k
            if(1.5*Ht(ibuild) .lt. z(k))exit
         enddo
         if(k_ref .lt. nz)then
            Bs=min(Weff(ibuild),Ht(ibuild))
            BL=max(Weff(ibuild),Ht(ibuild))
            Rscale(ibuild) = ((Bs**(2./3.))*(BL**(1./3.)))
            Rcx(ibuild)=(0.9*Rscale(ibuild))
            vd= 0.5*0.22*Rscale(ibuild)
            zref2=(vd/sqrt(0.5*Rcx(ibuild)))
            invvd=1./vd
            ! MAN 07/25/2008 stretched vertical grid
            do k=kend(ibuild)+1,k_ref
               kendv=k
               if(Ht(ibuild)+vd .lt. zm(k))exit
            enddo
                  if(perpendicular_flag .gt. 0)then
                     !$omp parallel do private(i,k,uflag,vflag,wflag,x_u,y_u,x_v,y_v,hx,hy,hdu,hdv, &
                     !$omp zref2u,zref2v,k_shellu,k_shellv,shell_heightu,shell_heightv,denomu,denomv, &
                     !$omp vel_magu,vel_magv,zr, shell_heightu_part, shell_heightv_part)
lp005:               do j=jstart(ibuild),jend(ibuild)
lp004:                  do i=istart(ibuild),iend(ibuild)
                           uflag=0
                           vflag=0
                           wflag=0
                           !check to see if velocity vector is above the building or in a street canyon cell
                           if(icellflag(i,j,kend(ibuild)) .eq. 0)then
                              uflag=1
                              vflag=1
                              wflag=1
                           else
                              if(icellflag(i-1,j,kend(ibuild)) .eq. 0)uflag=1
                              if(icellflag(i,j-1,kend(ibuild)) .eq. 0)vflag=1
                           endif
                           if(icellflag(i,j,kend(ibuild)+1) .eq. 6 .or. icellflag(i,j,kend(ibuild)+1) .eq. 4)then
                              uflag=0
                              vflag=0
                              wflag=0
                           endif
                           if(uflag+vflag+wflag .gt. 0 .and. icellflag(i,j,kend(ibuild)+1) .gt. 0)then
                              x_u=((real(i)-1)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_u=-((real(i)-1)*dx-xco)*sinRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              x_v=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-1)*dy-yco)*sinRotationAngle
                              y_v=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-1)*dy-yco)*cosRotationAngle
                              hx=abs(x_u-xfront)
                              hy=abs(y_u-yfront)
                              hdu=(ns_flag*hy+(1-ns_flag)*hx)
                              hx=abs(x_v-xfront)
                              hy=abs(y_v-yfront)
                              hdv=(ns_flag*hy+(1-ns_flag)*hx)
                              zref2u=zref2*sqrt(hdu)
                              zref2v=zref2*sqrt(hdv)
                              k_shellu=0
                              k_shellv=0
                              do k=kend(ibuild),nz
                                 if(zref2u+Ht(ibuild) .lt. zm(k+1) .and. k_shellu .lt. 1)then
                                    k_shellu=k
                                 endif
                                 if(zref2v+Ht(ibuild) .lt. zm(k+1) .and. k_shellv .lt. 1)then
                                    k_shellv=k
                                 endif
                                 if(k_shellu .gt. 0 .and. k_shellv .gt. 0)exit
                              enddo
                              shell_heightu_part = 1-((0.5*Rcx(ibuild)-hdu)/(0.5*Rcx(ibuild)))**2.
                              shell_heightv_part = 1-((0.5*Rcx(ibuild)-hdv)/(0.5*Rcx(ibuild)))**2.
                              if(shell_heightu_part .gt. 0)then
                              	shell_heightu=vd*sqrt(1-((0.5*Rcx(ibuild)-hdu)/(0.5*Rcx(ibuild)))**2.)
                              else
                              	shell_heightu=0.
                              endif
                              if(shell_heightv_part .gt. 0.) then
                              	shell_heightv=vd*sqrt(1-((0.5*Rcx(ibuild)-hdv)/(0.5*Rcx(ibuild)))**2.)
                              else
                              	shell_heightv=0.
                              endif
                              if(k_shellu .le. kend(ibuild))then
                                 uflag=0
                                 denomu=1.
                              else
                                 denomu=1./log((zm(k_shellu)-Ht(ibuild))*invzo)
                              endif
                              if(k_shellv .le. kend(ibuild))then
                                 vflag=0
                                 denomv=1.
                              else
                                 denomv=1./log((zm(k_shellv)-Ht(ibuild))*invzo)
                              endif
                              vel_magu=uo_roof(i,j,k_shellu)
                              vel_magv=vo_roof(i,j,k_shellv)
lp003:                        do k=kend(ibuild)+1,kendv
                                 k_shellu=0
                                 k_shellv=0
                                 zr=zm(k)-Ht(ibuild)
                                 if(icellflag(i,j,k) .lt. 1)then
                                    exit
                                 else
                                    if(uflag .eq. 1)then
                                       if(zr .le. zref2u)then
                                          uo(i,j,k)=vel_magu*log(zr*invzo)*denomu
                                          if(abs(uo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized U exceeds max in rooftop',&
                                                uo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellu=1
                                          if(wflag .eq. 1)then
                                             icellflag(i,j,k)=3
                                          endif
                                       endif
                                       if(hdu .lt. Rcx(ibuild) .and. zr .le. shell_heightu)then
                                          uo(i,j,k)=-uo_roof(i,j,k)*abs((shell_heightu-zr)*invvd)
                                          if(abs(uo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized U exceeds max in rooftop',&
                                                uo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellu=1
                                          if(wflag .eq. 1)then
                                             icellflag(i,j,k)=3
                                          endif
                                       endif
                                    endif
                                    if(vflag .eq. 1)then
                                       if(zr .le. zref2v)then
                                          vo(i,j,k)=vel_magv*log(zr*invzo)*denomv
                                          if(abs(vo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized V exceeds max in rooftop',&
                                                vo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellv=1
                                          if(wflag .eq. 1)then
                                             icellflag(i,j,k)=3
                                          endif
                                       endif
                                       if(hdv .lt. Rcx(ibuild) .and. zr .le. shell_heightv)then
                                          vo(i,j,k)=-vo_roof(i,j,k)*abs((shell_heightv-zr)*invvd)
                                          if(abs(vo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized V exceeds max in rooftop',&
                                                vo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellv=1
                                          if(wflag .eq. 1)then
                                             icellflag(i,j,k)=3
                                          endif
                                       endif
                                    endif
                                 endif
                                 if(k_shellu+k_shellv .lt. 1)exit
                              enddo  lp003
                           endif
                        enddo   lp004      
                     enddo   lp005
                     !$omp end parallel do
                  else !delta wing vortex
                     !$omp parallel do private(i,k,uflag,vflag,wflag,x_u,y_u,x_v,y_v,x_w,y_w, &
                     !$omp hxu,hyu,hdu,hxv,hyv,hdv,hxw,hyw,hdwx,hdwy,k_shellu,k_shellv,k_shellw, &
                     !$omp vel_mag,zr)
lp008:               do j=jstart(ibuild),jend(ibuild)
lp007:                  do i=istart(ibuild),iend(ibuild)
                           uflag=0
                           vflag=0
                           wflag=0
                           !check to see if velocity vector is above the building or in a street canyon cell
                           if(icellflag(i,j,kend(ibuild)) .eq. 0)then
                              uflag=1
                              vflag=1
                              wflag=1
                           else
                              if(icellflag(i-1,j,kend(ibuild)) .eq. 0)uflag=1
                              if(icellflag(i,j-1,kend(ibuild)) .eq. 0)vflag=1
                           endif
                           if(icellflag(i,j,kend(ibuild)+1) .eq. 6 .or. icellflag(i,j,kend(ibuild)+1) .eq. 4)then
                              uflag=0
                              vflag=0
                              wflag=0
                           endif
                           if(uflag+vflag+wflag .gt. 0 .and. icellflag(i,j,kend(ibuild)+1) .gt. 0)then
                              x_u=((real(i)-1)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_u=-((real(i)-1)*dx-xco)*sinRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              x_v=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-1)*dy-yco)*sinRotationAngle
                              y_v=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-1)*dy-yco)*cosRotationAngle
                              x_w=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_w=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              hxu=abs(x_u-xfront)
                              hyu=abs(y_u-yfront)
                              hdu=min(hxu,hyu)
                              hxv=abs(x_v-xfront)
                              hyv=abs(y_v-yfront)
                              hdv=min(hxv,hyv)
                              hxw=abs(x_w-xfront)
                              hyw=abs(y_w-yfront)
                              hdwx=hyw*tanRoofangle
                              hdwy=hxw*tanRoofangle
lp006:                        do k=kend(ibuild)+1,k_ref
                                 k_shellu=0
                                 k_shellv=0
                                 k_shellw=0
                                 if(icellflag(i,j,k) .lt. 1)then
                                    exit
                                 else
                                    zr=zm(k)-Ht(ibuild)
                                    vel_mag=sqrt((uo_roof(nint(xco/dx),nint(yco/dy),k)**2.)+&
                                       (vo_roof(nint(xco/dx),nint(yco/dy),k)**2.))
                                    if(uflag .eq. 1)then
                                       if(hxu .le. min(Rcx(ibuild),2*hyu*tanRoofangle))then
                                          if(zr .le. min(Rcx(ibuild),hyu*tanRoofangle))then
                                             uo(i,j,k)=vel_mag*cosXnorm
                                             if(abs(uo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized U exceeds max in rooftop',&
                                                   uo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                          if(zr .le. min(Rcx(ibuild),2*hyu*tanRoofangle))then
                                             uo(i,j,k)=vel_mag*cosXnormpPi
                                             if(abs(uo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized U exceeds max in rooftop',&
                                                   uo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                       endif
                                       if(hyu .le. min(Rcx(ibuild),2*hxu*tanRoofangle))then
                                          if(zr .le. min(Rcx(ibuild),hxu*tanRoofangle))then
                                             uo(i,j,k)=vel_mag*cosYnorm
                                             if(abs(uo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized U exceeds max in rooftop',&
                                                   uo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                          if(zr .le. min(Rcx(ibuild),2*hxu*tanRoofangle))then
                                             uo(i,j,k)=vel_mag*cosYnormpPi
                                             if(abs(uo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized U exceeds max in rooftop',&
                                                   uo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                       endif
                                    endif
                                    if(vflag .eq. 1)then
                                       if(hxv .le. min(Rcx(ibuild),2*hyv*tanRoofangle))then
                                          if(zr .le. min(Rcx(ibuild),hyv*tanRoofangle))then
                                             vo(i,j,k)=vel_mag*sinXnorm
                                             if(abs(vo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized V exceeds max in rooftop',&
                                                   vo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                          if(zr .le. min(Rcx(ibuild),2*hyv*tanRoofangle))then
                                             vo(i,j,k)=vel_mag*sinXnormpPi
                                             if(abs(vo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized V exceeds max in rooftop',&
                                                   vo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                       endif
                                       if(hyv .le. min(Rcx(ibuild),2*hxv*tanRoofangle))then
                                          if(zr .le. min(Rcx(ibuild),hxv*tanRoofangle))then
                                             vo(i,j,k)=vel_mag*sinYnorm
                                             if(abs(vo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized V exceeds max in rooftop',&
                                                   vo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                          if(zr .le. min(Rcx(ibuild),2*hxv*tanRoofangle))then
                                             vo(i,j,k)=vel_mag*sinYnormpPi
                                             if(abs(vo(i,j,k)) .gt. max_velmag)then
                                                print*,'Parameterized V exceeds max in rooftop',&
                                                   vo(i,j,k),max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                       endif
                                    endif
                                    if(wflag .eq. 1)then
                                       if(hxw .le. min(Rcx(ibuild),2*hdwx) .and. zr .le. min(Rcx(ibuild),2*hdwx))then
                                          wo(i,j,k)=0.1*vel_mag*((hdwx-hxw)/hdwx)*(1-abs((zr-hdwx)/hdwx))
                                          if(abs(wo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized W exceeds max in rooftop',&
                                                wo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellw=1
                                          icellflag(i,j,k)=3
                                       endif
                                       if(hyw .le. min(Rcx(ibuild),2*hdwy) .and. zr .le. min(Rcx(ibuild),2*hdwy))then
                                          wo(i,j,k)=0.1*vel_mag*((hdwy-hyw)/hdwy)*(1-abs((zr-hdwy)/hdwy))
                                          if(abs(wo(i,j,k)) .gt. max_velmag)then
                                             print*,'Parameterized W exceeds max in rooftop',&
                                                wo(i,j,k),max_velmag,i,j,k
                                          endif
                                          k_shellw=1
                                          icellflag(i,j,k)=3
                                       endif
                                    endif
                                 endif
                                 if(k_shellu+k_shellv+k_shellw .lt. 1)exit
                              enddo  lp006
                           endif
                        enddo   lp007      
                     enddo   lp008
                  endif
            
         endif
         return
      end
