program euler

    implicit none

    character(1000), parameter :: folder = 'data/'
    character(1000), parameter :: file = 'density'
    character(1000) :: file_name
    integer, parameter :: nx = 200, ny = 200
    integer, parameter :: nxL = int(.21*nx), nxR = int(.30*nx)
    integer, parameter :: nyB = int(.50*ny), nyT = int(.51*ny)
    integer, parameter :: nint = 10
    real, parameter :: Lx = 2., Ly = 2.
    real, parameter :: dx = Lx/(nx - 1), dy = Ly/(ny - 1)
    real, parameter :: dt = .1*min(dx, dy)
    real, parameter :: t = 10.
    real, parameter :: g = 1.4

    integer :: n, i, j, k
    real, dimension(nx, ny) :: r, u, v, p
    real, dimension(nx, ny, 4) :: Q, Qn, Qs, Qe, Qw, Qnb, Qsb, Qeb, Qwb
    real, dimension(4) :: Gn, Gs, Fe, Fw
    real :: dn, ds, de, dw, di, dj, ri, rj, xi, xj
    real :: rn, rs, re, rw, rij
    real :: un, us, ue, uw, uij
    real :: vn, vs, ve, vw, vij
    real :: pn, ps, pe, pw, pij
    real :: En, Es, Ee, Ew, Eij

    ! initial conditions
    r = 1.
    u = 1.4
    v = 0.
    p = 1./g

    Q(:, :, 1) = r
    Q(:, :, 2) = r*u
    Q(:, :, 3) = r*v
    Q(:, :, 4) = .5*r*(u*u + v*v) + p/(g - 1)

    ! boundary conditions
    ! outlet (transmissive)
    Q(nx, 2:ny - 1, :) = Q(nx - 1, 2:ny - 1, :)
    ! bottom (reflective)
    Q(2:nx - 1, 1, 1) = Q(2:nx - 1, 2, 1)
    Q(2:nx - 1, 1, 2) = Q(2:nx - 1, 2, 2)
    Q(2:nx - 1, 1, 3) = -Q(2:nx - 1, 2, 3)
    Q(2:nx - 1, 1, 4) = Q(2:nx - 1, 2, 4)
    ! top (reflective)
    Q(2:nx - 1, ny, 1) = Q(2:nx - 1, ny - 1, 1)
    Q(2:nx - 1, ny, 2) = Q(2:nx - 1, ny - 1, 2)
    Q(2:nx - 1, ny, 3) = -Q(2:nx - 1, ny - 1, 3)
    Q(2:nx - 1, ny, 4) = Q(2:nx - 1, ny - 1, 4)
    ! body (reflective)
    ! body left
    Q(nxL, nyB:nyT, 1) = Q(nxL - 1, nyB:nyT, 1)
    Q(nxL, nyB:nyT, 2) = -Q(nxL - 1, nyB:nyT, 2)
    Q(nxL, nyB:nyT, 3) = Q(nxL - 1, nyB:nyT, 3)
    Q(nxL, nyB:nyT, 4) = Q(nxL - 1, nyB:nyT, 4)
    ! body right
    Q(nxR, nyB:nyT, 1) = Q(nxR + 1, nyB:nyT, 1)
    Q(nxR, nyB:nyT, 2) = -Q(nxR + 1, nyB:nyT, 2)
    Q(nxR, nyB:nyT, 3) = Q(nxR + 1, nyB:nyT, 3)
    Q(nxR, nyB:nyT, 4) = Q(nxR + 1, nyB:nyT, 4)
    ! body bottom
    Q(nxL:nxR, nyB, 1) = Q(nxL:nxR, nyB - 1, 1)
    Q(nxL:nxR, nyB, 2) = Q(nxL:nxR, nyB - 1, 2)
    Q(nxL:nxR, nyB, 3) = -Q(nxL:nxR, nyB - 1, 3)
    Q(nxL:nxR, nyB, 4) = Q(nxL:nxR, nyB - 1, 4)
    ! body top
    Q(nxL:nxR, nyT, 1) = Q(nxL:nxR, nyT + 1, 1)
    Q(nxL:nxR, nyT, 2) = Q(nxL:nxR, nyT + 1, 2)
    Q(nxL:nxR, nyT, 3) = -Q(nxL:nxR, nyT + 1, 3)
    Q(nxL:nxR, nyT, 4) = Q(nxL:nxR, nyT + 1, 4)

    do n = 1, int(t/dt)
        print *, n, int(t/dt)
        
        ! Step 1 : Reconstruction

        do concurrent (i = 2:nx - 1)
            do concurrent (j = 2:ny - 1)
                if (.not. in_square(i, j)) then
                    do concurrent (k = 1:4)

                        ! Step 1a : Compute slope vectors

                        dn = Q(i, j + 1, k) - Q(i, j, k)
                        ds = Q(i, j, k) - Q(i, j - 1, k)
                        de = Q(i + 1, j, k) - Q(i, j, k)
                        dw = Q(i, j, k) - Q(i - 1, j, k)
                        di = .5*(de + dw)
                        dj = .5*(dn + ds)
                        
                        ! Step 1b : Compute slope limiters

                        if (dn /= 0.) then
                            rj = ds/dn
                            if (rj <= 0.) then
                                xj = 0.
                            else if (rj > 0. .and. rj <= .5) then
                                xj = 2*rj
                            else if (rj > .5 .and. rj <= 1.) then
                                xj = 1.
                            else
                                xj = min(rj, 2./(1. + rj), 2.)
                            end if
                        else
                            xj = 0.
                        end if

                        if (de /= 0.) then
                            ri = dw/de
                            if (ri <= 0.) then
                                xi = 0.
                            else if (ri > 0. .and. ri <= .5) then
                                xi = 2*ri
                            else if (ri > .5 .and. ri <= 1.) then
                                xi = 1.
                            else
                                xi = min(ri, 2./(1. + ri), 2.)
                            end if
                        else
                            xi = 0.
                        end if

                        ! Step 1c : Compute boundary extrapolated values

                        Qn(i, j, k) = Q(i, j, k) + .5*xj*dj
                        Qs(i, j, k) = Q(i, j, k) - .5*xj*dj
                        Qe(i, j, k) = Q(i, j, k) + .5*xi*di
                        Qw(i, j, k) = Q(i, j, k) - .5*xi*di

                    end do
                end if
            end do
        end do

        ! Step 2 : Evolution

        do concurrent (i = 2:nx - 1)
            do concurrent (j = 2:ny - 1)
                if (.not. in_square(i, j)) then

                    rn = Qn(i, j, 1)
                    un = Qn(i, j, 2)/rn
                    vn = Qn(i, j, 3)/rn
                    En = Qn(i, j, 4)
                    pn = (g - 1)*(En - .5*rn*(un*un + vn*vn))

                    rs = Qs(i, j, 1)
                    us = Qs(i, j, 2)/rs
                    vs = Qs(i, j, 3)/rs
                    Es = Qs(i, j, 4)
                    ps = (g - 1)*(Es - .5*rs*(us*us + vs*vs))

                    re = Qe(i, j, 1)
                    ue = Qe(i, j, 2)/re
                    ve = Qe(i, j, 3)/re
                    Ee = Qe(i, j, 4)
                    pe = (g - 1)*(Ee - .5*re*(ue*ue + ve*ve))

                    rw = Qw(i, j, 1)
                    uw = Qw(i, j, 2)/rw
                    vw = Qw(i, j, 3)/rw
                    Ew = Qw(i, j, 4)
                    pw = (g - 1)*(Ew - .5*rw*(uw*uw + vw*vw))

                    Gn = (/rn*vn, rn*un*vn, rn*vn*vn + pn, vn*(En + pn)/)
                    Gs = (/rs*vs, rs*us*vs, rs*vs*vs + ps, vs*(Es + ps)/)
                    Fe = (/re*ue, re*ue*ue + pe, re*ue*ve, ue*(Ee + pe)/)
                    Fw = (/rw*uw, rw*uw*uw + pw, rw*uw*vw, uw*(Ew + pw)/)

                    Qnb(i, j, :) = Qn(i, j, :) + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
                    Qsb(i, j, :) = Qs(i, j, :) + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
                    Qeb(i, j, :) = Qe(i, j, :) + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
                    Qwb(i, j, :) = Qw(i, j, :) + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
                
                end if
            end do
        end do

        ! boundary conditions
        ! inlet (Dirichlet)
        Qeb(1, 2:ny - 1, :) = Q(1, 2:ny - 1, :)
        ! outlet (transmissive)
        Qwb(nx, 2:ny - 1, :) = Qeb(nx - 1, 2:ny - 1, :)
        ! bottom (reflective)
        Qnb(2:nx - 1, 1, 1) = Qsb(2:nx - 1, 2, 1)
        Qnb(2:nx - 1, 1, 2) = Qsb(2:nx - 1, 2, 2)
        Qnb(2:nx - 1, 1, 3) = -Qsb(2:nx - 1, 2, 3)
        Qnb(2:nx - 1, 1, 4) = Qsb(2:nx - 1, 2, 4)
        ! top (reflective)
        Qsb(2:nx - 1, ny, 1) = Qnb(2:nx - 1, ny - 1, 1)
        Qsb(2:nx - 1, ny, 2) = Qnb(2:nx - 1, ny - 1, 2)
        Qsb(2:nx - 1, ny, 3) = -Qnb(2:nx - 1, ny - 1, 3)
        Qsb(2:nx - 1, ny, 4) = Qnb(2:nx - 1, ny - 1, 4)
        ! body (reflective)
        ! body left
        Qwb(nxL, nyB:nyT, 1) = Qeb(nxL - 1, nyB:nyT, 1)
        Qwb(nxL, nyB:nyT, 2) = -Qeb(nxL - 1, nyB:nyT, 2)
        Qwb(nxL, nyB:nyT, 3) = Qeb(nxL - 1, nyB:nyT, 3)
        Qwb(nxL, nyB:nyT, 4) = Qeb(nxL - 1, nyB:nyT, 4)
        ! body right
        Qeb(nxR, nyB:nyT, 1) = Qwb(nxR + 1, nyB:nyT, 1)
        Qeb(nxR, nyB:nyT, 2) = -Qwb(nxR + 1, nyB:nyT, 2)
        Qeb(nxR, nyB:nyT, 3) = Qwb(nxR + 1, nyB:nyT, 3)
        Qeb(nxR, nyB:nyT, 4) = Qwb(nxR + 1, nyB:nyT, 4)
        ! body bottom
        Qsb(nxL:nxR, nyB, 1) = Qnb(nxL:nxR, nyB - 1, 1)
        Qsb(nxL:nxR, nyB, 2) = Qnb(nxL:nxR, nyB - 1, 2)
        Qsb(nxL:nxR, nyB, 3) = -Qnb(nxL:nxR, nyB - 1, 3)
        Qsb(nxL:nxR, nyB, 4) = Qnb(nxL:nxR, nyB - 1, 4)
        ! body top
        Qnb(nxL:nxR, nyT, 1) = Qsb(nxL:nxR, nyT + 1, 1)
        Qnb(nxL:nxR, nyT, 2) = Qsb(nxL:nxR, nyT + 1, 2)
        Qnb(nxL:nxR, nyT, 3) = -Qsb(nxL:nxR, nyT + 1, 3)
        Qnb(nxL:nxR, nyT, 4) = Qsb(nxL:nxR, nyT + 1, 4)

        ! Step 3 : Solution

        do concurrent (i = 2:nx - 1)
            do concurrent (j = 2:ny - 1)
                if (.not. in_square(i, j)) then

                    rij = r(i, j)
                    uij = u(i, j)
                    vij = v(i, j)
                    pij = p(i, j)
                    Eij = .5*rij*(uij*uij + vij*vij) + pij/(g - 1)

                    Fw = get_xflux(i, j, Qeb, Qwb)
                    Fe = get_xflux(i + 1, j, Qeb, Qwb)
                    Gs = get_yflux(i, j, Qnb, Qsb)
                    Gn = get_yflux(i, j + 1, Qnb, Qsb)

                    Q(i, j, :) = (/rij, rij*uij, rij*vij, Eij/) + (dt/dx)*(Fw - Fe) + (dt/dy)*(Gs - Gn)
                    
                end if
            end do
        end do
        
        ! boundary conditions
        ! outlet (transmissive)
        Q(nx, 2:ny - 1, :) = Q(nx - 1, 2:ny - 1, :)
        ! bottom (reflective)
        Q(2:nx - 1, 1, 1) = Q(2:nx - 1, 2, 1)
        Q(2:nx - 1, 1, 2) = Q(2:nx - 1, 2, 2)
        Q(2:nx - 1, 1, 3) = -Q(2:nx - 1, 2, 3)
        Q(2:nx - 1, 1, 4) = Q(2:nx - 1, 2, 4)
        ! top (reflective)
        Q(2:nx - 1, ny, 1) = Q(2:nx - 1, ny - 1, 1)
        Q(2:nx - 1, ny, 2) = Q(2:nx - 1, ny - 1, 2)
        Q(2:nx - 1, ny, 3) = -Q(2:nx - 1, ny - 1, 3)
        Q(2:nx - 1, ny, 4) = Q(2:nx - 1, ny - 1, 4)
        ! body (reflective)
        ! body left
        Q(nxL, nyB:nyT, 1) = Q(nxL - 1, nyB:nyT, 1)
        Q(nxL, nyB:nyT, 2) = -Q(nxL - 1, nyB:nyT, 2)
        Q(nxL, nyB:nyT, 3) = Q(nxL - 1, nyB:nyT, 3)
        Q(nxL, nyB:nyT, 4) = Q(nxL - 1, nyB:nyT, 4)
        ! body right
        Q(nxR, nyB:nyT, 1) = Q(nxR + 1, nyB:nyT, 1)
        Q(nxR, nyB:nyT, 2) = -Q(nxR + 1, nyB:nyT, 2)
        Q(nxR, nyB:nyT, 3) = Q(nxR + 1, nyB:nyT, 3)
        Q(nxR, nyB:nyT, 4) = Q(nxR + 1, nyB:nyT, 4)
        ! body bottom
        Q(nxL:nxR, nyB, 1) = Q(nxL:nxR, nyB - 1, 1)
        Q(nxL:nxR, nyB, 2) = Q(nxL:nxR, nyB - 1, 2)
        Q(nxL:nxR, nyB, 3) = -Q(nxL:nxR, nyB - 1, 3)
        Q(nxL:nxR, nyB, 4) = Q(nxL:nxR, nyB - 1, 4)
        ! body top
        Q(nxL:nxR, nyT, 1) = Q(nxL:nxR, nyT + 1, 1)
        Q(nxL:nxR, nyT, 2) = Q(nxL:nxR, nyT + 1, 2)
        Q(nxL:nxR, nyT, 3) = -Q(nxL:nxR, nyT + 1, 3)
        Q(nxL:nxR, nyT, 4) = Q(nxL:nxR, nyT + 1, 4)

        r = Q(:, :, 1)
        u = Q(:, :, 2)/r
        v = Q(:, :, 3)/r
        p = (g - 1)*(Q(:, :, 4) - .5*r*(u*u + v*v))

        if (mod(n, nint) == 0) then
            write(file_name, '(a, i10.10, a)') trim(file), n, '.txt'
            open(10, file=trim(folder)//trim(file_name))
            do concurrent (i = 1:nx)
                do concurrent (j = 1:ny)
                    write(10, *) i, j, Q(i, j, 1)
                end do
            end do
        end if

    end do

contains

    pure function in_square(i, j) result(check)
        
        integer, intent(in) :: i, j

        logical :: check

        if (nxL <= i .and. i <= nxR .and. nyB <= j .and. j <= nyT) then
            check = .true.
        else
            check = .false.
        end if
    
    end function in_square

    pure function get_xflux(i, j, Qeb, Qwb) result(FLR)

        integer, intent(in) :: i, j
        real, intent(in) :: Qeb(:, :, :), Qwb(:, :, :)

        real :: rL, uL, vL, pL, EL
        real :: rR, uR, vR, pR, ER
        real :: aL, aR, sL, sR
        real, dimension(4) :: QL, QR, FL, FR
        real, dimension(4) :: FLR

        rL = Qeb(i - 1, j, 1)
        uL = Qeb(i - 1, j, 2)/rL
        vL = Qeb(i - 1, j, 3)/rL
        EL = Qeb(i - 1, j, 4)
        pL = (g - 1)*(EL - .5*rL*(uL*uL + vL*vL))

        rR = Qwb(i, j, 1)
        uR = Qwb(i, j, 2)/rR
        vR = Qwb(i, j, 3)/rR
        ER = Qwb(i, j, 4)
        pR = (g - 1)*(ER - .5*rR*(uR*uR + vR*vR))
    
        ! compute wave speeds
        aL = sqrt(g*pL/rL)
        aR = sqrt(g*pR/rR)
        sL = min(uL - aL, uR - aR)
        sR = max(uL + aL, uR + aR)
        
        ! compute intercell flux
        QL = (/rL, rL*uL, rL*vL, EL/)
        QR = (/rR, rR*uR, rR*vR, ER/)
        FL = (/rL*uL, rL*uL*uL + pL, rL*uL*vL, uL*(EL + pL)/)
        FR = (/rR*uR, rR*uR*uR + pR, rR*uR*vR, uR*(ER + pR)/)
    
        if (sL > 0.) then
            FLR = FL
        else if (sL <= 0. .and. sR >= 0.) then
            FLR = (sR*FL - sL*FR + sL*sR*(QR - QL))/(sR - sL)
        else if (sR < 0.) then
            FLR = FR
        end if

    end function get_xflux

    pure function get_yflux(i, j, Qnb, Qsb) result(GBT)

        integer, intent(in) :: i, j
        real, intent(in) :: Qnb(:, :, :), Qsb(:, :, :)

        real :: rB, uB, vB, pB, EB
        real :: rT, uT, vT, pT, ET
        real :: aB, aT, sB, sT
        real, dimension(4) :: QB, QT, GB, GT
        real, dimension(4) :: GBT

        rB = Qnb(i, j - 1, 1)
        uB = Qnb(i, j - 1, 2)/rB
        vB = Qnb(i, j - 1, 3)/rB
        EB = Qnb(i, j - 1, 4)
        pB = (g - 1)*(EB - .5*rB*(uB*uB + vB*vB))
    
        rT = Qsb(i, j, 1)
        uT = Qsb(i, j, 2)/rT
        vT = Qsb(i, j, 3)/rT
        ET = Qsb(i, j, 4)
        pT = (g - 1)*(ET - .5*rT*(uT*uT + vT*vT))
        
        ! compute wave speeds
        aB = sqrt(g*pB/rB)
        aT = sqrt(g*pT/rT)
        sB = min(vB - aB, vT - aT)
        sT = max(vB + aB, vT + aT)

        ! compute intercell flux
        QB = (/rB, rB*uB, rB*vB, EB/)
        QT = (/rT, rT*uT, rT*vT, ET/)
        GB = (/rB*vB, rB*uB*vB, rB*vB*vB + pB, vB*(EB + pB)/)
        GT = (/rT*vT, rT*uT*vT, rT*vT*vT + pT, vT*(ET + pT)/)

        if (sB > 0.) then
            GBT = GB
        else if (sB <= 0. .and. sT >= 0.) then
            GBT = (sT*GB - sB*GT + sB*sT*(QT - QB))/(sT - sB)
        else if (sT < 0.) then
            GBT = GT
        end if

    end function get_yflux

end program euler
