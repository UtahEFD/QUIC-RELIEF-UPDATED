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
      subroutine building_connect
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subroutine building_connect - this subroutine 
!	- called by bcsetup.f90
!	- calls none
! ERP 8/17/05
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         use datamodule ! make data from module "datamodule" visible
         implicit none

         integer jj,grpid_chk
         real Pc
         real debug1, debug2, debug3

!	Coeffeicent
         Pc = 1

         grpid_chk = group_id(ibuild)
	
         do jj=1,inumbuild
            if(jj.ne.ibuild)then
               if(group_id(jj).eq.grpid_chk)then
                  if(Ht(jj).eq.zfo_actual(ibuild)) then
!MAN 8/30/2005 stacked building fix
                     debug1 = Weff(jj)
                     debug2 = Weff(ibuild)
                     debug3 = Ht(jj)
                     if(Weff(jj) > 0.0)then
                        zfo(ibuild) = Ht(jj)*(1 - (Weff(ibuild)/Weff(jj))**Pc)
                     else
                        zfo(ibuild) = 0.0
                     endif
                     zfo(ibuild) = max(zfo(ibuild),0.)
                  endif
               endif
            endif
         enddo



!	stop

         return
      end
