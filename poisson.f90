      subroutine poisson
         use datamodule ! make data from module "datamodule" visible
         implicit none
         integer iter
         real invomegarelax
         invomegarelax=1./omegarelax
         do iter=1,10
            !$omp parallel do private(i,j)
            do k= 2,nz-1
               do j=2,ny-1
                  do i= 2,nx-1
                     if(icellflag(i,j,k) .eq. 9 .and. icellflag(i-1,j,k) .eq. 9)then
                        uo(i,j,k)=invomegarelax*denom(i,j,k)*((e(i,j,k)*uo(i+1,j,k)+f(i,j,k)*uo(i-1,j,k)) &
                                  +(g(i,j,k)*uo(i,j+1,k)+h(i,j,k)*uo(i,j-1,k)) &
                                  +(m(i,j,k)*uo(i,j,k+1)+n(i,j,k)*uo(i,j,k-1)))
                     endif
                     if(icellflag(i,j,k) .eq. 9 .and. icellflag(i,j-1,k) .eq. 9)then
                        vo(i,j,k)=invomegarelax*denom(i,j,k)*((e(i,j,k)*vo(i+1,j,k)+f(i,j,k)*vo(i-1,j,k)) &
                                    +(g(i,j,k)*vo(i,j+1,k)+h(i,j,k)*vo(i,j-1,k)) &
                                    +(m(i,j,k)*vo(i,j,k+1)+n(i,j,k)*vo(i,j,k-1)))
                     endif
                     if(icellflag(i,j,k) .eq. 9 .and. icellflag(i,j,k-1) .eq. 9)then
                        wo(i,j,k)=invomegarelax*denom(i,j,k)*((e(i,j,k)*wo(i+1,j,k)+f(i,j,k)*wo(i-1,j,k)) &
                                    + (g(i,j,k)*wo(i,j+1,k)+h(i,j,k)*wo(i,j-1,k)) &
                                    + (m(i,j,k)*wo(i,j,k+1)+n(i,j,k)*wo(i,j,k-1)))
                     endif
                  enddo
               enddo
            enddo
            !$omp end parallel do
         enddo
         return
      end
