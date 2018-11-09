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
      subroutine wallbc

         use datamodule ! make data from module "datamodule" visible
! operations  done a priori to speed up the code AAG & IS  07/03/06
         implicit none
         ! non boundary cells
!erp 6/12/2006 the following includes the redefinition of b.c. based on 
!C. Bathke work
         e(:,:,:)=1.
         e(nx-1,:,:)=0.
         f(:,:,:)=1.
         f(1,:,:)=0.
         g(:,:,:)=1.
         g(:,ny-1,:)=0.
         h(:,:,:)=1.
         h(:,1,:)=0.
         m(:,:,:)=1.
         m(:,:,nz-1)=0.
         n(:,:,:)=1.
         n(:,:,1)=0.
         !$omp parallel do private(i,j)
         do k=2,nz-2
            do j=2,ny-2
               do i=2,nx-2
                  if(icellflag(i,j,k) .ne. 0)then
                     if(icellflag(i,j,k-1) +icellflag(i,j,k+1)+icellflag(i-1,j,k) &
                           +icellflag(i+1,j,k)+icellflag(i,j-1,k)+icellflag(i,j+1,k) .lt. 1)then
                        icellflag(i,j,k)=0
                     else
                        if(icellflag(i,j,k-1) .eq. 0)n(i,j,k)=0.
                        if(icellflag(i,j,k+1) .eq. 0)m(i,j,k)=0.
                        if(icellflag(i-1,j,k) .eq. 0)f(i,j,k)=0.
                        if(icellflag(i+1,j,k) .eq. 0)e(i,j,k)=0.
                        if(icellflag(i,j-1,k) .eq. 0)h(i,j,k)=0.
                        if(icellflag(i,j+1,k) .eq. 0)g(i,j,k)=0.
                     endif
                  endif
               enddo     
            enddo      
         enddo
         !$omp end parallel do
         !$omp parallel do
         do k=2,nz-1
             e(:,:,k)=e(:,:,k)/(dx*dx)
             f(:,:,k)=f(:,:,k)/(dx*dx)
             g(:,:,k)=g(:,:,k)/(dy*dy)
             h(:,:,k)=h(:,:,k)/(dy*dy)
             m(:,:,k)=m(:,:,k)/(dz_array(k)*0.5*(dz_array(k)+dz_array(k+1)))
             n(:,:,k)=n(:,:,k)/(dz_array(k)*0.5*(dz_array(k)+dz_array(k-1)))
             denom(:,:,k)=omegarelax/(e(:,:,k)+f(:,:,k)+g(:,:,k)+h(:,:,k)+m(:,:,k)+n(:,:,k))
          enddo
          !$omp end parallel do
         return
      end
