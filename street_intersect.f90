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
!
!
! Subroutine to determine the location of a street intersection
! erp 1/2006
! NOTE: Need to make modification to handle multiple runs 1/24/2006

      subroutine street_intersect
         use datamodule
         implicit none

         integer changeflag,intersect_flag,istart_intflag,jstart_intflag,NS_flag !,EW_flag
         integer, allocatable :: intersect(:,:,:),intersect_1(:,:,:),intersect_2(:,:,:),intersect_1opp(:,:,:),intersect_2opp(:,:,:) 
		 integer, allocatable :: E_W_flag(:,:,:),W_E_flag(:,:,:),N_S_flag(:,:,:),S_N_flag(:,:,:)  !SUP
         
         allocate(intersect(nx-1,ny-1,nz-1),intersect_1(nx-1,ny-1,nz-1),intersect_2(nx-1,ny-1,nz-1))
         allocate(intersect_1opp(nx-1,ny-1,nz-1),intersect_2opp(nx-1,ny-1,nz-1))
	     allocate(E_W_flag(nx-1,ny-1,nz-1),W_E_flag(nx-1,ny-1,nz-1),N_S_flag(nx-1,ny-1,nz-1),S_N_flag(nx-1,ny-1,nz-1)) !SUP
         
         intersect_flag=0
! make sure that all cells that are buildings have a zero velocity within them
         do k=1,nz-1
            do j=1,ny-1
               do i=1,nx-1
                  if(icellflag(i,j,k) .eq. 0)then ! MAN 7/8/2005 celltype definition change
                     uo(i,j,k)  =0.
                     uo(i+1,j,k)=0.
                     vo(i,j,k)  =0.
                     vo(i,j+1,k)=0.
                     wo(i,j,k)  =0.
                     wo(i,j,k+1)=0.
                  endif
               enddo
            enddo
         enddo
         intersect(1:nx-1,1:ny-1,1:nz-1) = 0
         intersect_1(1:nx-1,1:ny-1,1:nz-1) = 0 !x-direction flag
         intersect_2(1:nx-1,1:ny-1,1:nz-1) = 0 !y-direction flag
         intersect_1opp(1:nx-1,1:ny-1,1:nz-1)=0
         intersect_2opp(1:nx-1,1:ny-1,1:nz-1)=0
         E_W_flag=0
         W_E_flag=0
         N_S_flag=0
         S_N_flag=0
         changeflag = 0
! sweep through (x) to find intersections	
         do k=1,nz-1
            do j=1,ny-1
!SUP sweep through +x
               do i=2,nx-1
!determine where the street interesection begins
                  if(icellflag(i-1,j,k) .eq. 6 .and. icellflag(i,j,k) .ne. 6 .and. icellflag(i,j,k) .ne. 0)then
                     changeflag=1
                     istart_intflag = i
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 6 .or.	&	!run into another street canyon   
                        changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 0 .or.&       !run into another building 
						      changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 1)then		!run into free atm.
                     changeflag=0
                  endif
                  intersect_1(i,j,k) = changeflag
               enddo
!if we get to the end of a row and changeflag = 1, then no SI exists reset those 
               if(changeflag .eq. 1)then
                  intersect_1(istart_intflag:nx-1,j,k) = 0
               endif
               changeflag = 0	!reset flag
!SUP sweep through -x
               do i=nx-2,1,-1
!determine where the street interesection begins
                  if(icellflag(i+1,j,k) .eq. 6 .and. icellflag(i,j,k) .ne. 6 .and. icellflag(i,j,k) .ne. 0)then
                     changeflag=1
                     istart_intflag = i
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 6 .or. &	!run into another street canyon   
                        changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 0 .or. &           !run into another building
						      changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 1 )then		!run into free atm.
                     changeflag=0
                  endif
                  intersect_1opp(i,j,k) = changeflag
               enddo
!if we get to the end of a row and changeflag = 1, then no SI exists reset those 
               if(changeflag .eq. 1)then
                  intersect_1opp(nx-1:istart_intflag:-1,j,k) = 0
               endif
               changeflag = 0	!reset flag
            enddo
         enddo
! now sweep in the j direction
         changeflag = 0
         do k=1,nz-1
            do i=1,nx-1
               do j=2,ny-1
                  if(icellflag(i,j-1,k) .eq. 6 .and. icellflag(i,j,k) .ne. 6 .and. icellflag(i,j,k) .ne. 0)then
                     changeflag=1
                     jstart_intflag = j
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 6 .or. &	!run into another street canyon   
                        changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 0 .or. &           !run into another building
						      changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 1)then		!run into free atm.
                     changeflag=0
                  endif
                  intersect_2(i,j,k) = changeflag
               enddo
!if we get to the end of a row and changeflag = 1, then no SI exists reset those 
               if(changeflag.eq.1)then
                  intersect_2(i,jstart_intflag:ny-1,k) = 0 !SUP changed intersect_1 to _2
               endif
               changeflag = 0
!SUP sweep through -y
			   do j=ny-2,1,-1
                  if(icellflag(i,j+1,k) .eq. 6 .and. icellflag(i,j,k) .ne. 6 .and. icellflag(i,j,k) .ne. 0)then
                     changeflag=1
                     jstart_intflag = j
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 6 .or. &	!run into another street canyon   
                        changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 0 .or.  &          !run into another building   
						      changeflag .eq. 1 .and. icellflag(i,j,k) .eq. 1 )then		!run into free atm.
                     changeflag=0
                  endif
                  intersect_2opp(i,j,k) = changeflag
               enddo
!if we get to the end of a row and changeflag = 1, then no SI exists reset those 
               if(changeflag.eq.1)then
                  intersect_2opp(i,ny-1:jstart_intflag:-1,k) = 0
               endif
               changeflag = 0
            enddo
         enddo


         do k=1,nz-1
            do j=1,ny-1
               do i=1,nx-1
                   if((intersect_1(i,j,k) .eq. 1 .or. intersect_1opp(i,j,k) .eq. 1) .and. & 
				             (intersect_2(i,j,k) .eq. 1 .or. intersect_2opp(i,j,k) .eq. 1))intersect(i,j,k)=1
               enddo
            enddo
         enddo
!SUP looking to make sure that there are street canyons on 2 or more adjacent sides
         do k=1,nz-1
            do j=2,ny-1
			      NS_flag=0
               do i=2,nx-1
                  if(intersect(i,j,k) .eq. 1 .and. icellflag(i-1,j,k) .eq. 6) NS_flag=1
				      if(intersect(i,j,k) .ne. 1 .and. NS_flag .eq. 1) NS_flag=0
                  if(NS_flag .eq. 1) E_W_flag(i,j,k)=1
               enddo
               NS_flag=0
			   do i=nx-2,1,-1
                  if(intersect(i,j,k) .eq. 1 .and. icellflag(i+1,j,k) .eq. 6) NS_flag=1
				      if(intersect(i,j,k) .ne. 1 .and. NS_flag .eq. 1) NS_flag=0
                  if(NS_flag .eq. 1) W_E_flag(i,j,k)=1
			   enddo
            enddo
         enddo
         do k=1,nz-1
            do i=2,nx-2
			      NS_flag=0
               do j=2,ny-1
                  if(intersect(i,j,k) .eq. 1 .and. icellflag(i,j-1,k) .eq. 6) NS_flag=1
				      if(intersect(i,j,k) .ne. 1 .and. NS_flag .eq. 1) NS_flag=0
                  if(NS_flag .eq. 1) S_N_flag(i,j,k)=1
               enddo
               NS_flag=0
			      do j=ny-2,1,-1
                  if(intersect(i,j,k) .eq. 1 .and. icellflag(i,j+1,k) .eq. 6) NS_flag=1
				      if(intersect(i,j,k) .ne. 1 .and. NS_flag .eq. 1) NS_flag=0
                  if(NS_flag .eq. 1) N_S_flag(i,j,k)=1
			      enddo
            enddo
         enddo
         do k=1,nz-1
            do j=1,ny-1
               do i=1,nx-1
!                  if(intersect_1(i,j,k).eq.1 .and. intersect_2(i,j,k).eq.1)icellflag(i,j,k)=9 !SUP
                  if((E_W_flag(i,j,k) .eq. 1 .or. W_E_flag(i,j,k) .eq. 1).and.&
				            (S_N_flag(i,j,k) .eq. 1 .or. N_S_flag(i,j,k) .eq. 1))icellflag(i,j,k)=9
               enddo
            enddo
         enddo
         deallocate(intersect_1,intersect_2,intersect)
		   deallocate(E_W_flag,W_E_flag,N_S_flag,S_N_flag)
         return
      end
