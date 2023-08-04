program main

    implicit none

    character(100), parameter ::                    mesh        = 'airfoil.su2'                     ! mesh file
    character(100), parameter ::                    data        = 'data/'                           ! data folder
    integer, parameter ::                           precision   = 8                                 ! real number precision
    integer, parameter ::                           interval    = 100                               ! write interval
    real(precision), parameter ::                   g           = 1.4                               ! heat capacity ratio
    real(precision), parameter ::                   M           = .8                                ! freestream Mach
    real(precision), parameter ::                   a           = 1.25                              ! angle of attack (degrees)
    real(precision), parameter ::                   kH          = .05                               ! Harten correction factor
    real(precision), parameter ::                   kM          = 1.                                ! initial freestream Mach fraction
    real(precision), parameter ::                   kt          = .01                               ! final freestream Mach total time fraction
    real(precision), parameter ::                   tt          = 50.                               ! total time
    real(precision), parameter ::                   CFL         = 2.                                ! Courant-Friedrichs-Lewy number
    real(precision), dimension(4), parameter ::     RK          = [.0833, .2069, .4265, 1.]         ! Runge-Kutta scheme
    ! real(precision), dimension(5), parameter ::     RK          = [.0533, .1263, .2375, .4414, 1.]  ! Runge-Kutta scheme

    ! program start

    integer :: n_cells, n_faces, n_nodes, n_freestream_faces, n_slipwall_faces, n_boundary_faces, n_interior_faces
    integer, dimension(:), allocatable :: freestream_faces, slipwall_faces, boundary_faces, interior_faces
    integer, dimension(:, :), allocatable :: cells, faces, cell_faces, face_cells, face_faces, cell_cells
    real(precision), dimension(:, :), allocatable :: nodes

    real(precision), dimension(:), allocatable :: Sf, Vc, wf
    real(precision), dimension(:, :), allocatable :: nf, cf, cc, sc
    real(precision), dimension(:, :, :), allocatable :: rc

    integer :: it
    real(precision) :: t, dt, r0, u0, v0, p0
    real(precision), dimension(4) :: Q0
    real(precision), dimension(:, :), allocatable :: Q, Qf, F, W, W0
    real(precision), dimension(:, :, :), allocatable :: QL, dQ

    ! create mesh

    call set_mesh()

    ! allocate variables

    allocate(Q(n_cells, 4))
    allocate(Qf(n_faces, 4))
    allocate(dQ(n_cells, 4, 2))
    allocate(QL(n_cells, 4, 3))
    allocate(F(n_faces, 4))
    allocate(W(n_cells, 4))
    allocate(W0(n_cells, 4))

    ! initial conditions

    r0 = 1.
    u0 = M*cosd(a)*kM
    v0 = M*sind(a)*kM
    p0 = 1./g
    Q0 = [r0, u0, v0, p0]

    Q(:, 1) = r0
    Q(:, 2) = u0
    Q(:, 3) = v0
    Q(:, 4) = p0
    call set_boundary()
    call set_gradient()
    call set_reconstruction()
    call set_flux()
    call set_W()
    
    ! time iterations

    it = 0
    t = 0.
    call set_timestep()
    
    do while (t < tt)

        call set_iteration()

        print *, it, t, dt, norm2([u0, v0])
        ! print *, 'r', minval(Q(:, 1)), maxval(Q(:, 1))
        ! print *, 'u', minval(Q(:, 2)), maxval(Q(:, 2))
        ! print *, 'v', minval(Q(:, 3)), maxval(Q(:, 3))
        ! print *, 'p', minval(Q(:, 4)), maxval(Q(:, 4))
        print *, '-----'

        it = it + 1
        t = t + dt
        call set_timestep()

        call save_M()  ! save Mach number

    end do

    call save_Cp()  ! save pressure coefficient

contains

    subroutine set_mesh()

        character(100) :: line
        integer :: i, j, n1, n2, n3, f1, f2, fi, fj, c1, c2, ci, cj
        integer, dimension(:), allocatable :: l1, l2, l3
        real(precision) :: nx, ny, n1x, n1y, n2x, n2y, n3x, n3y, ccx, ccy, cfx, cfy, scj, d1, d2
        real(precision), dimension(2) :: nfj
        
        open(10, file=trim(mesh), action='read')

        read(10, *)
        read(10, *) line, n_cells
        allocate(cells(n_cells, 3))
        do i = 1, n_cells
            read(10, *) line, n1, n2, n3
            cells(i, :) = [n1 + 1, n2 + 1, n3 + 1]
        end do

        read(10, *) line, n_nodes
        allocate(nodes(n_nodes, 2))
        do i = 1, n_nodes
            read(10, *) nx, ny
            nodes(i, :) = [nx, ny]
        end do

        read(10, *)
        read(10, *)
        read(10, *) line, n_freestream_faces
        l1 = [integer ::]
        l2 = [integer ::]
        do i = 1, n_freestream_faces
            read(10, *) line, n1, n2
            l1 = [l1, n1 + 1]
            l2 = [l2, n2 + 1]
        end do

        read(10, *)
        read(10, *) line, n_slipwall_faces
        do i = 1, n_slipwall_faces
            read(10, *) line, n1, n2
            l1 = [l1, n2 + 1]  ! swap
            l2 = [l2, n1 + 1]  ! swap
        end do

        do i = 1, n_cells
            n1 = cells(i, 1)
            n2 = cells(i, 2)
            n3 = cells(i, 3)
    
            call set_face(l1, l2, n1, n2)
            call set_face(l1, l2, n2, n3)
            call set_face(l1, l2, n3, n1)
        end do
        n_faces = size(l1)
        allocate(faces(n_faces, 2))
        faces(:, 1) = l1
        faces(:, 2) = l2

        fi = 1
        fj = n_freestream_faces
        freestream_faces = [(i, i = fi, fj)]

        fi = fi + n_freestream_faces
        fj = fj + n_slipwall_faces
        slipwall_faces = [(i, i = fi, fj)]

        n_boundary_faces = n_freestream_faces + n_slipwall_faces
        n_interior_faces = n_faces - n_boundary_faces

        fi = 1
        fj = n_boundary_faces
        boundary_faces = [(i, i = fi, fj)]

        fi = n_boundary_faces + 1
        fj = n_faces
        interior_faces = [(i, i = fi, fj)]

        ! cell_faces

        allocate(cell_faces(n_cells, 3))
        do i = 1, n_cells
            n1 = cells(i, 1)
            n2 = cells(i, 2)
            n3 = cells(i, 3)

            do j = 1, n_faces
                if (check_face(l1(j), l2(j), n1, n2)) then
                    cell_faces(i, 1) = j
                end if
                if (check_face(l1(j), l2(j), n2, n3)) then
                    cell_faces(i, 2) = j
                end if
                if (check_face(l1(j), l2(j), n3, n1)) then
                    cell_faces(i, 3) = j
                end if
            end do
        end do

        ! face_cells

        allocate(face_cells(n_faces, 2))
        do i = 1, n_faces
            l3 = [integer ::]
            do j = 1, n_cells
                if (any(cell_faces(j, :) == i)) then
                    l3 = [l3, j]
                end if
            end do        
            if (size(l3) == 1) then
                face_cells(i, 1) = l3(1)
                face_cells(i, 2) = l3(1)
            else
                face_cells(i, 1) = l3(1)
                face_cells(i, 2) = l3(2)
            end if
        end do
        
        ! Sf, nf, cf

        allocate(Sf(n_faces))
        allocate(nf(n_faces, 2))
        allocate(cf(n_faces, 2))
        allocate(face_faces(n_faces, 2))
        do i = 1, n_faces
            n1 = faces(i, 1)
            n2 = faces(i, 2)
            n1x = nodes(n1, 1)
            n1y = nodes(n1, 2)
            n2x = nodes(n2, 1)
            n2y = nodes(n2, 2)

            Sf(i) = sqrt((n2x - n1x)**2 + (n2y - n1y)**2)
            nf(i, :) = [(n2y - n1y)/Sf(i), -(n2x - n1x)/Sf(i)]
            cf(i, :) = [.5*(n2x + n1x), .5*(n2y + n1y)]

            c1 = face_cells(i, 1)
            c2 = face_cells(i, 2)

            l3 = [integer ::]
            do j = 1, 3  ! separate loop
                if (cell_faces(c1, j) == i) then
                    l3 = [l3, j]
                end if
            end do
            do j = 1, 3  ! separate loop
                if (cell_faces(c2, j) == i) then
                    l3 = [l3, j]
                end if
            end do
            face_faces(i, :) = l3
        end do

        ! Vc, cc, rc, sc

        allocate(Vc(n_cells))
        allocate(cc(n_cells, 2))
        allocate(rc(n_cells, 3, 2))
        allocate(sc(n_cells, 3))
        allocate(cell_cells(n_cells, 3))
        do i = 1, n_cells
            n1 = cells(i, 1)
            n2 = cells(i, 2)
            n3 = cells(i, 3)
            n1x = nodes(n1, 1)
            n1y = nodes(n1, 2)
            n2x = nodes(n2, 1)
            n2y = nodes(n2, 2)
            n3x = nodes(n3, 1)
            n3y = nodes(n3, 2)

            Vc(i) = .5*abs(n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))
            cc(i, :) = [(n1x + n2x + n3x)/3, (n1y + n2y + n3y)/3]

            do j = 1, 3
                fj = cell_faces(i, j)
                cfx = cf(fj, 1)
                cfy = cf(fj, 2)
                ccx = cc(i, 1)
                ccy = cc(i, 2)
                
                rc(i, j, :) = [cfx - ccx, cfy - ccy]

                if (i == face_cells(fj, 1)) then
                    sc(i, j) = 1.
                    cell_cells(i, j) = face_cells(fj, 2)
                else
                    sc(i, j) = -1.
                    cell_cells(i, j) = face_cells(fj, 1)
                end if
            end do
        end do

        ! wf

        allocate(wf(n_faces))
        do i = 1, n_faces
            c1 = face_cells(i, 1)
            c2 = face_cells(i, 2)

            d1 = 1/norm2(cc(c1, :) - cf(i, :))  ! reciprocal
            d2 = 1/norm2(cc(c2, :) - cf(i, :))  ! reciprocal
            wf(i) = d1/(d1 + d2)
        end do

        ! check mesh

        if (.false.) then
            fi = 26
            ci = 84

            c1 = face_cells(fi, 1)
            c2 = face_cells(fi, 2)
            f1 = face_faces(fi, 1)
            f2 = face_faces(fi, 2)

            print *, fi, c1, c2, f1, f2, cell_faces(c1, :), cell_faces(c2, :)

            do j = 1, 3
                fj = cell_faces(ci, j)
                cj = cell_cells(ci, j)
                c1 = face_cells(fj, 1)
                c2 = face_cells(fj, 2)

                scj = sc(ci, j)
                nfj = nf(fj, :)

                print *, ci, fj, cj, c1, c2, scj, nfj
            end do
        end if

    end subroutine set_mesh

    subroutine set_face(la, lb, na, nb)

        integer, dimension(:), allocatable, intent(in out) :: la, lb
        integer, intent(in) :: na, nb
        
        integer, dimension(:), allocatable :: lna, lnb
        logical :: in
        
        allocate(lna, mold=la)
        allocate(lnb, mold=lb)
        lna(:) = na
        lnb(:) = nb

        in = any((la == lna) .and. (lb == lnb)) .or. any((la == lnb) .and. (lb == lna))
        if (.not. in) then
            la = [la, na]
            lb = [lb, nb]
        end if

    end subroutine set_face

    pure function check_face(la, lb, na, nb)

        integer, intent(in) :: la, lb, na, nb
        logical :: check_face

        if ((na == la .and. nb == lb) .or. (na == lb .and. nb == la)) then
            check_face = .true.
        else
            check_face = .false.
        end if
    
    end function check_face

    subroutine set_boundary()

        integer :: i, fi, c1, c2, ci
        real(precision), dimension(2) :: Vi, Vj, nfi
        real(precision), dimension(4) :: Qi
        
        do i = 1, n_freestream_faces
            fi = freestream_faces(i)
            ci = face_cells(fi, 1)
            Qi = Q(ci, :)
            
            Qf(fi, :) = .5*(Qi + Q0)
        end do
        
        do i = 1, n_slipwall_faces
            fi = slipwall_faces(i)
            ci = face_cells(fi, 1)
            Vi = [Q(ci, 2), Q(ci, 3)]
            nfi = nf(fi, :)
            Vj = Vi - dot_product(Vi, nfi)*nfi  ! face value (not dummy value)
            
            Qf(fi, 1) = Q(ci, 1)
            Qf(fi, 2) = Vj(1)
            Qf(fi, 3) = Vj(2)
            Qf(fi, 4) = Q(ci, 4)
        end do

        do i = 1, n_interior_faces
            fi = interior_faces(i)
            c1 = face_cells(fi, 1)
            c2 = face_cells(fi, 2)

            Qf(fi, :) = wf(fi)*Q(c1, :) + (1 - wf(fi))*Q(c2, :)
        end do
        
    end subroutine set_boundary

    subroutine set_gradient()

        integer :: i, j, k, fj
        real(precision) :: Vci, Sfj, scj, Qfj
        real(precision), dimension(2) :: nfj

        do i = 1, n_cells
            Vci = Vc(i)

            do k = 1, 4
                dQ(i, k, :) = [0., 0.]

                do j = 1, 3
                    fj = cell_faces(i, j)
                    Sfj = Sf(fj)
                    nfj = nf(fj, :)
                    scj = sc(i, j)
                    Qfj = Qf(fj, k)

                    dQ(i, k, :) = dQ(i, k, :) + (1./Vci)*Qfj*Sfj*scj*nfj
                end do
            end do
        end do

    end subroutine set_gradient

    subroutine set_reconstruction()

        integer :: i, j, k, fj, cj
        real(precision) :: phi_i, psi_phi_i
        real(precision), dimension(2) :: d_phi_i, rcj, Vi, Vj, nfj
        real(precision), dimension(3) :: phi_j, d_phi_j
        real(precision), dimension(4) :: Qj

        do i = 1, n_cells
            do k = 1, 4
                phi_i = Q(i, k)
                d_phi_i = dQ(i, k, :)

                do j = 1, 3
                    fj = cell_faces(i, j)
                    cj = cell_cells(i, j)
                    rcj = rc(i, j, :)

                    if (any(fj == freestream_faces)) then
                        phi_j(j) = Q0(k)
                        d_phi_j(j) = dot_product(d_phi_i, rcj)
                    else if (any(fj == slipwall_faces)) then
                        Vi = [Q(i, 2), Q(i, 3)]
                        nfj = nf(fj, :)
                        Vj = Vi - 2.*dot_product(Vi, nfj)*nfj  ! dummy value
                        Qj = [Q(i, 1), Vj(1), Vj(2), Q(i, 4)]

                        phi_j(j) = Qj(k)
                        d_phi_j(j) = dot_product(d_phi_i, rcj)
                    else  ! interior faces
                        phi_j(j) = Q(cj, k)
                        d_phi_j(j) = dot_product(d_phi_i, rcj)
                    end if  
                end do

                psi_phi_i = get_limiter(phi_i, phi_j, d_phi_j)

                do j = 1, 3
                    QL(i, k, j) = phi_i + psi_phi_i*d_phi_j(j)
                end do
            end do
        end do

    end subroutine set_reconstruction

    pure function get_limiter(phi_i, phi_j, d_phi_j)
        
        real(precision), intent(in) :: phi_i
        real(precision), dimension(3), intent(in) :: phi_j, d_phi_j
        
        integer :: j
        real(precision) :: phi_max, phi_min, get_limiter
        real(precision), dimension(3) :: psi_phi

        phi_max = max(phi_i, maxval(phi_j))
        phi_min = min(phi_i, minval(phi_j))

        do j = 1, 3
            if (d_phi_j(j) > 0.) then
                psi_phi(j) = min(1., (phi_max - phi_i)/d_phi_j(j))
            else if (d_phi_j(j) < 0.) then
                psi_phi(j) = min(1., (phi_min - phi_i)/d_phi_j(j))
            else
                psi_phi(j) = 1.
            end if
        end do

        get_limiter = minval(psi_phi)

    end function get_limiter

    subroutine set_flux()

        integer :: i, ci, cj, fi, fj
        real(precision) :: rL, uL, vL, pL, EL, HL, rR, uR, vR, pR, ER, HR
        real(precision) :: rt, ut, vt, Ht, ct, Wt, dr, du, dv, dp, dW
        real(precision), dimension(2) :: nfi, Vi, Vj
        real(precision), dimension(4) :: Qi, Qj, FL, FR, F1, F2, F3

        do i = 1, n_faces
            ci = face_cells(i, 1)
            fi = face_faces(i, 1)
            Qi = QL(ci, :, fi)
            nfi = nf(i, :)
            
            ! left state

            rL = Qi(1)
            uL = Qi(2)
            vL = Qi(3)
            pL = Qi(4)
            EL = (uL*uL + vL*vL)/2. + pL/rL/(g - 1.)
            HL = EL + pL/rL
            
            FL = get_F(Qi, nfi)
            
            ! right state

            if (any(i == freestream_faces)) then
                rR = Q0(1)
                uR = Q0(2)
                vR = Q0(3)
                pR = Q0(4)
                ER = (uR*uR + vR*vR)/2. + pR/rR/(g - 1.)
                HR = ER + pR/rR
                
                FR = get_F(Q0, nfi)
            else if (any(i == slipwall_faces)) then
                Vi = [uL, vL]
                Vj = Vi - 2.*dot_product(Vi, nfi)*nfi  ! dummy value

                rR = rL
                uR = Vj(1)
                vR = Vj(2)
                pR = pL
                ER = (uR*uR + vR*vR)/2. + pR/rR/(g - 1.)
                HR = ER + pR/rR

                Qj = [rR, uR, vR, pR]

                FR = get_F(Qj, nfi)
            else ! interior faces
                cj = face_cells(i, 2)
                fj = face_faces(i, 2)  ! local index
                Qj = QL(cj, :, fj)
                
                rR = Qj(1)
                uR = Qj(2)
                vR = Qj(3)
                pR = Qj(4)
                ER = (uR*uR + vR*vR)/2. + pR/rR/(g - 1.)
                HR = ER + pR/rR
                
                FR = get_F(Qj, nfi)
            end if
            
            ! average state

            rt = sqrt(rL*rR)
            ut = (uL*sqrt(rL) + uR*sqrt(rR))/(sqrt(rL) + sqrt(rR))
            vt = (vL*sqrt(rL) + vR*sqrt(rR))/(sqrt(rL) + sqrt(rR))
            Ht = (HL*sqrt(rL) + HR*sqrt(rR))/(sqrt(rL) + sqrt(rR))
            ct = sqrt((g - 1.)*(Ht - (ut*ut + vt*vt)/2.))
            Wt = dot_product([ut, vt], nfi)
            
            ! state difference

            dr = rR - rL
            du = uR - uL
            dv = vR - vL
            dp = pR - pL
            dW = dot_product([uR, vR], nfi) - dot_product([uL, vL], nfi)

            F1 = get_correction(Wt - ct, kH*ct)*((dp - rt*ct*dW)/2./ct/ct) &
            *[real(1., precision), ut - ct*nfi(1), vt - ct*nfi(2), Ht - ct*Wt]
            F2 = get_correction(Wt, kH*ct)*((dr - dp/ct/ct)*[real(1., precision), ut, vt, (ut*ut + vt*vt)/2.] &
            + rt*[real(0., precision), du - dW*nfi(1), dv - dW*nfi(2), ut*du + vt*dv - Wt*dW])
            F3 = get_correction(Wt + ct, kH*ct)*((dp + rt*ct*dW)/2./ct/ct) &
            *[real(1., precision), ut + ct*nfi(1), vt + ct*nfi(2), Ht + ct*Wt]

            F(i, :) = .5*(FL + FR - (F1 + F2 + F3))
        end do

    end subroutine set_flux

    pure function get_F(Qi, nfi)

        real(precision), dimension(4), intent(in) :: Qi
        real(precision), dimension(2), intent(in) :: nfi

        real(precision), dimension(4) :: get_F
        real(precision) :: ri, ui, vi, pi, Ei, Hi, Wi

        ri = Qi(1)
        ui = Qi(2)
        vi = Qi(3)
        pi = Qi(4)
        Ei = (ui*ui + vi*vi)/2. + pi/ri/(g - 1.)
        Hi = Ei + pi/ri
        Wi = dot_product([ui, vi], nfi)

        get_F(1) = ri*Wi
        get_F(2) = ri*ui*Wi + nfi(1)*pi
        get_F(3) = ri*vi*Wi + nfi(2)*pi
        get_F(4) = ri*Hi*Wi

    end function get_F

    pure function get_correction(Li, di)

        real(precision), intent(in) :: Li, di

        real(precision) :: get_correction

        if (abs(Li) > di) then
            get_correction = abs(Li)
        else
            get_correction = (Li*Li + di*di)/2/di
        end if

    end function get_correction
    
    subroutine set_timestep()

        integer :: i, j, fj
        real(precision) :: Vci, lci, Sfj, scj, rfj, ufj, vfj, pfj, afj
        real(precision), dimension(2) :: nfj
        
        dt = 1.
        do i = 1, n_cells
            Vci = Vc(i)

            lci = 0.
            do j = 1, 3
                fj = cell_faces(i, j)
                Sfj = Sf(fj)
                nfj = nf(fj, :)
                scj = sc(i, j)
                rfj = Qf(fj, 1)
                ufj = Qf(fj, 2)
                vfj = Qf(fj, 3)
                pfj = Qf(fj, 4)
                afj = sqrt(g*pfj/rfj)

                lci = lci + (abs(dot_product([ufj, vfj], scj*nfj)) + afj)*Sfj
            end do

            dt = min(dt, CFL*Vci/lci)
        end do

    end subroutine set_timestep

    subroutine set_W()
        
        integer :: i
        real(precision) :: r, u, v, p

        do i = 1, n_cells
            r = Q(i, 1)
            u = Q(i, 2)
            v = Q(i, 3)
            p = Q(i, 4)

            W(i, 1) = r
            W(i, 2) = r*u
            W(i, 3) = r*v
            W(i, 4) = r*(u*u + v*v)/2. + p/(g - 1.)
        end do

    end subroutine set_W

    subroutine set_Q()

        integer :: i
        real(precision) :: r, u, v, p, E

        do i = 1, n_cells
            r = W(i, 1)
            u = W(i, 2)/W(i, 1)
            v = W(i, 3)/W(i, 1)
            E = W(i, 4)/W(i, 1)
            p = (g - 1.)*r*(E - (u*u + v*v)/2.)

            Q(i, 1) = r
            Q(i, 2) = u
            Q(i, 3) = v
            Q(i, 4) = p
        end do

    end subroutine set_Q

    subroutine set_iteration()

        integer :: i, j, k, fj
        real(precision) :: ak, Vci, Sfj, scj
        real(precision), dimension(4) :: Rki
        
        ! mach ramp

        if (t < kt*tt) then
            u0 = M*cosd(a)*(kM + (1 - kM)/(kt*tt)*t)
            v0 = M*sind(a)*(kM + (1 - kM)/(kt*tt)*t)
            Q0 = [r0, u0, v0, p0]
            call set_boundary()
            call set_gradient()
            call set_reconstruction()
            call set_flux()
        end if

        ! runge-kutta

        W0(:, :) = W(:, :)  ! set W0

        do k = 1, size(RK)
            ak = RK(k)

            do i = 1, n_cells
                Vci = Vc(i)

                Rki = [0., 0., 0., 0.]
                do j = 1, 3
                    fj = cell_faces(i, j)
                    Sfj = Sf(fj)
                    scj = sc(i, j)

                    Rki = Rki + scj*F(fj, :)*Sfj
                end do

                W(i, :) = W0(i, :) - ak*dt/Vci*Rki
            end do

            call set_Q()
            call set_boundary()
            call set_gradient()
            call set_reconstruction()
            call set_flux()
        end do

    end subroutine set_iteration

    subroutine save_M()

        character(100) :: file
        integer :: i

        if (mod(it, interval) == 0) then
            write(file, '(a, i10.10, a)') 'M', it, '.txt'
            open(100, file=trim(data)//trim(file))
            do i = 1, n_cells
                write(100, *) cc(i, 1), cc(i, 2), norm2([Q(i, 2), Q(i, 3)])
            end do
        end if

    end subroutine save_M

    subroutine save_Cp()

        integer :: i, fi

        open(100, file=trim(data)//'Cp.txt')
        do i = 1, n_slipwall_faces
            fi = slipwall_faces(i)
            write(100, *) cf(fi, 1), cf(fi, 2), (Qf(fi, 4) - p0)/(.5*r0*(u0*u0 + v0*v0))
        end do

    end subroutine save_Cp

end program main
