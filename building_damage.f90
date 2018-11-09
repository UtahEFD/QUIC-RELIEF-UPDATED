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
      subroutine building_damage
         use datamodule
         implicit none
         
         real xco,yco
         
         open(unit=52,file='QP_source.inp',status='old')
         read(52,*)
         read(52,*)
         read(52,*)
         
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         read(52,*)
         
         read(52,*)explosionx
         read(52,*)explosiony
         read(52,*)
         read(52,*)hemass
         
         close(52)
         
         Rdestroyed=255*((hemass/907185)**(1./3.))
         Rdamaged=1425*((hemass/907185)**(1./3.))
         
         do ibuild=1,inumbuild
            if(bldtype(ibuild) .eq. 2)cycle
            if(bldgeometry(ibuild) .eq. 3)then
               xco = xfo(ibuild)
               yco = yfo(ibuild)
            else
               if(gamma(ibuild) .ne. 0.)then
                  xco = xfo(ibuild) + 0.5*Lti(ibuild)*cos(gamma(ibuild))
                  yco = yfo(ibuild) + 0.5*Lti(ibuild)*sin(gamma(ibuild))
               else
                  xco = xfo(ibuild) + 0.5*Lti(ibuild)
                  yco = yfo(ibuild)
               endif
            endif
            Rbuild=sqrt(((xco-explosionx)**2.)+((yco-explosiony)**2.))
            if(Rbuild .le. Rdestroyed)then
               bld_damage(ibuild)=2
            elseif(Rbuild .le. Rdamaged)then
               bld_damage(ibuild)=1
            endif
         enddo
         do ibuild=1,inumbuild
            if(bldtype(ibuild) .eq. 2)cycle
            if(bld_damage(ibuild) .eq. 1 .and. zfo_actual(ibuild) .gt. 0.)then
               do j=1,inumbuild
                  if(group_id(ibuild) .eq. group_id(j) .and. ibuild .ne. j &
                        .and. zfo_actual(j) .lt. zfo_actual(ibuild) .and. bld_damage(j) .eq. 2 &
                        .and. bldtype(j) .ne. 2)then
                     bld_damage(ibuild)=2
                  endif
               enddo
            endif
         enddo
         return
      end