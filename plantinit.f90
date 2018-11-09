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
      subroutine plantinit

         !this function initializes the velocity profile in the specified
         !vegetative area area using the MacDonald (2000) approach.
         !ca - attenuation coefficient
         !vk - von Karmen constant

         !TMB 7/18/03 internal plant canopy variables
         !ANU 6/2005 implemented

         use datamodule
         implicit none
   
         real bisect  !bisection function
         real vegvelfrac,avg_atten,num_atten,canopy_slope_match
         canopy_ktop(:,:)=0
         canopy_ustar(:,:)=0.
         canopy_zo(:,:)=0.
         canopy_d(:,:)=0.
         vk=.4 !VonKarmen Constant
!erp  add subroutine to calculate ustar and zo  using a least squares
! regression of the current initialization
         
         call regress
         do j=1,ny-1
            do i=1,nx-1
               if(canopy_top(i,j) .gt. 0.)then
                  if(damageflag .eq. 1)then
                     Rbuild=sqrt((((real(i)-0.5)-explosionx)**2.)+(((real(j)-0.5)-explosiony)**2.))
                     if(Rbuild .le. Rdestroyed)then
                        canopy_top(i,j)=0.
                        do k=2,nz
                           canopy_atten(i,j,k)=0.
                        enddo
                        cycle
                     endif
                  endif
                  canopy_d(i,j) = bisect(canopy_ustar(i,j),canopy_zo(i,j), &
                        canopy_top(i,j),canopy_atten(i,j,canopy_ktop(i,j)),vk,0.)
                 ! if(canopy_d(i,j) .gt. 0.99*canopy_top(i,j))then
                  if(canopy_d(i,j) .eq. 10000)then
                    print*,"bisect.f90 failed to converge"
                    print*,"Arguments to canopy_slope_match are:"
                    print*,"zo =",canopy_zo(i,j)
                    print*,"H = ",canopy_top(i,j)
                    print*,"a = ",canopy_atten(i,j,canopy_ktop(i,j))
                    canopy_d(i,j)=canopy_slope_match(canopy_zo(i,j),canopy_top(i,j), &
                    canopy_atten(i,j,canopy_ktop(i,j)))
                  !      print*,"d passed from canopy_slope_match to plantinit is",canopy_d(i,j)
                  endif

                  !    canopy_d(i,j)=0.7*canopy_top(i,j)

                  !   canopy_zo(i,j)=min(0.1*canopy_top(i,j),0.5*dz_array(1))
                  !endif ! d = 10000 is the flag that indicates bisect.f90 didn't converge (LDUe)

                  uH = (canopy_ustar(i,j)/vk)*log((canopy_top(i,j)-canopy_d(i,j))/canopy_zo(i,j))
                  do k=2,nz
                     if(zm(k) .le. canopy_top(i,j))then
                        if(canopy_atten(i,j,k) .gt. 0.)then
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
                           vegvelfrac=log((canopy_top(i,j)-canopy_d(i,j))/canopy_zo(i,j))*&
                                 exp(avg_atten*((zm(k)/canopy_top(i,j))-1.))/&
                                 log(zm(k)/canopy_zo(i,j))
                           if(vegvelfrac .gt. 1. .or. vegvelfrac .lt. 0.)vegvelfrac=1
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
                           if(icellflag(i,j,k) .gt. 0)then
                              icellflag(i,j,k)=8
                           endif
                        endif
                     else
                        vegvelfrac=log((zm(k)-canopy_d(i,j))/canopy_zo(i,j))/&
                                 log(zm(k)/canopy_zo(i,j))
                        if(vegvelfrac .gt. 1. .or. vegvelfrac .lt. 0.)vegvelfrac=1
                        uo(i,j,k)=uo(i,j,k)*vegvelfrac
                        vo(i,j,k)=vo(i,j,k)*vegvelfrac
                        if(j .lt. ny-1)then
                           if(canopy_atten(i,j+1,canopy_ktop(i,j)) .eq. 0.)then
                              vo(i,j+1,k)=vo(i,j+1,k)*vegvelfrac
                           endif
                        endif
                        if(i .lt. nx-1)then
                           if(canopy_atten(i+1,j,canopy_ktop(i,j)) .eq. 0.)then
                              uo(i+1,j,k)=uo(i+1,j,k)*vegvelfrac
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
         return
      end
