function canopy_slope_match(zo,H,a)result(d)
! If bisect.f90 fails to converge, this function uses the bisection method to find the displacement height that makes the slopes of the in-canopy and above-canopy velocity profiles match.
implicit none
integer iter
real zo, H, a, tol
real d, d1, d2, f

!this function uses the bisection method to find the root of the specified
!equation

tol=zo/100
f=tol*10

!initial guess (fi)

if(zo .lt. H)then
d1 = zo
elseif(zo .ge. H)then
d1 = 0.1
endif
d2 = H;
d = 0.5*(d1 + d2)


!fi = log((H-d)/zo)-(H/(a*(H-d)))
if(a .gt. 0.)then

    iter = 0;
    do while(iter .lt. 200 .and. abs(f) .gt. tol .and. d .lt. H .and. d .gt. zo)
    iter = iter + 1
    d = 0.5*(d1 + d2) ! Algorithm for bisect method
    f = log((H-d)/zo)-(H/(a*(H-d)))

        if(f.gt.0) then
        d1 = d
        elseif(f.lt.0)then
        d2 = d
        endif
    enddo
    !print*,"d from canopy_slope_match is",d
    if(d .ge. H)then
        print*,"d is outside of range, d =",d
        d=0.7*H
        !d = 10000; ! This functions as a flag indicating the bisection didn't converge. 10000 is an arbitrary value, high enough that no QUIC simulation is likely to yield this normally (LDU)
    endif
else
d=10000
endif

end
