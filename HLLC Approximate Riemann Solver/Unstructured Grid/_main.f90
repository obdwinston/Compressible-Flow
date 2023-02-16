program euler
    
    implicit none

    real, parameter :: t = 10.
    real, parameter :: dt = 1e-4
    integer, parameter :: nt = 100 ! save interval

    real, parameter :: g = 1.4
    real, parameter :: r = 1.
    real, parameter :: u = 1.
    real, parameter :: v = 0.
    real, parameter :: p = 1./g

    character(1000), parameter :: mesh = 'cylinder.su2'
    character(1000), parameter :: folder = 'data/'
    character(1000), parameter :: file = 'density'

! program start

    character(1000), dimension(10) :: line
    character(1000) :: file_name
    integer :: file_unit
    integer :: i, j, k, n

    integer :: n_cells, n_nodes
    integer, allocatable, dimension(:, :) :: cells
    real, allocatable, dimension(:, :) :: nodes

    integer :: n1, n2, n3
    integer, allocatable, dimension(:) :: l1, l2, l3
    integer :: n_faces, n_wall_faces, n_body_faces, n_inlet_faces, n_outlet_faces, n_interior_faces
    integer :: start, finish
    integer, allocatable, dimension(:) :: m_wall_faces, m_body_faces, m_inlet_faces, m_outlet_faces, m_interior_faces
    integer, allocatable, dimension(:, :) :: faces, wall_faces, body_faces, inlet_faces, outlet_faces, interior_faces

    integer :: n_boundary_cells, n_interior_cells
    integer, allocatable, dimension(:) :: boundary_cells, interior_cells
    integer, allocatable, dimension(:, :) :: face_cells, cell_cells, cell_faces, cell_adjacent_faces

    integer :: ci, cj, fj, jj
    real :: p1x, p1y, p2x, p2y, p3x, p3y, sign
    real, allocatable, dimension(:) :: Vc
    real, allocatable, dimension(:, :) :: cc
    real, allocatable, dimension(:, :) :: Sf, gf
    real, allocatable, dimension(:, :, :) :: nf, cf, df, ef

    real, dimension(4) :: Ffh, Qfh
    real, dimension(4) :: Ff, QfL, QfR
    real, allocatable, dimension(:, :) :: Q
    real, allocatable, dimension(:, :, :) :: grad_Q
    real, allocatable, dimension(:, :, :) :: Qf, Qfb
    real :: rk, xk

    real, dimension(4) :: Qi, Qb

    real, dimension(2) :: Vf
    real, allocatable, dimension(:) :: Co

! mesh start

    open(newunit=file_unit, file=trim(mesh), action='read')

    ! cells [n1, n2, n3]
    read(file_unit, *)
    read(file_unit, *) line(1), n_cells
    allocate(cells(n_cells, 3))
    do i = 1, n_cells
        read(file_unit, *) line(1), cells(i, 1), cells(i, 2), cells(i, 3)
    end do
    cells = cells + 1

    ! nodes [p1, p2]
    read(file_unit, *) line(1), n_nodes
    allocate(nodes(n_nodes, 2))
    do i = 1, n_nodes
        read(file_unit, *) nodes(i, 1), nodes(i, 2)
    end do

    ! wall faces [n1, n2]
    l1 = [integer ::]
    l2 = [integer ::]
    read(file_unit, *)
    read(file_unit, *)
    read(file_unit, *) line(1), n_wall_faces
    allocate(wall_faces(n_wall_faces, 2))
    do i = 1, n_wall_faces
        read(file_unit, *) line(1), wall_faces(i, 1), wall_faces(i, 2)
    end do
    wall_faces = wall_faces + 1
    l1 = [l1, wall_faces(:, 1)]
    l2 = [l2, wall_faces(:, 2)]

    ! body faces [n1, n2]
    read(file_unit, *)
    read(file_unit, *) line(1), n_body_faces
    allocate(body_faces(n_body_faces, 2))
    do i = 1, n_body_faces
        read(file_unit, *) line(1), body_faces(i, 2), body_faces(i, 1) ! swap
    end do
    body_faces = body_faces + 1
    l1 = [l1, body_faces(:, 1)]
    l2 = [l2, body_faces(:, 2)]

    ! inlet faces [n1, n2]
    read(file_unit, *)
    read(file_unit, *) line(1), n_inlet_faces
    allocate(inlet_faces(n_inlet_faces, 2))
    do i = 1, n_inlet_faces
        read(file_unit, *) line(1), inlet_faces(i, 1), inlet_faces(i, 2)
    end do
    inlet_faces = inlet_faces + 1
    l1 = [l1, inlet_faces(:, 1)]
    l2 = [l2, inlet_faces(:, 2)]

    ! outlet faces [n1, n2]
    read(file_unit, *)
    read(file_unit, *) line(1), n_outlet_faces
    allocate(outlet_faces(n_outlet_faces, 2))
    do i = 1, n_outlet_faces
        read(file_unit, *) line(1), outlet_faces(i, 1), outlet_faces(i, 2)
    end do
    outlet_faces = outlet_faces + 1
    l1 = [l1, outlet_faces(:, 1)]
    l2 = [l2, outlet_faces(:, 2)]

    ! faces [n1, n2]
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        call append_array(l1, l2, n1, n2)
        call append_array(l1, l2, n2, n3)
        call append_array(l1, l2, n3, n1)
    end do
    n_faces = size(l1)
    allocate(faces(n_faces, 2))
    faces(:, 1) = l1
    faces(:, 2) = l2

    ! interior faces [n1, n2]
    n_interior_faces = n_faces - n_wall_faces - n_body_faces - n_inlet_faces - n_outlet_faces
    allocate(interior_faces(n_interior_faces, 2))
    interior_faces(:, 1) = faces(n_faces - n_interior_faces + 1:, 1)
    interior_faces(:, 2) = faces(n_faces - n_interior_faces + 1:, 2)

    start = 1
    finish = n_wall_faces
    m_wall_faces = [(i, i = start, finish)]

    start = start + n_wall_faces
    finish = finish + n_body_faces
    m_body_faces = [(i, i = start, finish)]

    start = start + n_body_faces
    finish = finish + n_inlet_faces
    m_inlet_faces = [(i, i = start, finish)]

    start = start + n_inlet_faces
    finish = finish + n_outlet_faces
    m_outlet_faces = [(i, i = start, finish)]

    start = start + n_outlet_faces
    finish = finish + n_interior_faces
    m_interior_faces = [(i, i = start, finish)]
    
    ! cell faces [f1, f2, f3]
    allocate(cell_faces(n_cells, 3))
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        do j = 1, n_faces
            if (compare_array(l1(j), l2(j), n1, n2)) then
                cell_faces(i, 1) = j
            end if
            if (compare_array(l1(j), l2(j), n2, n3)) then
                cell_faces(i, 2) = j
            end if
            if (compare_array(l1(j), l2(j), n3, n1)) then
                cell_faces(i, 3) = j
            end if
        end do
    end do

    ! face cells [c1, c2]
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

    ! Vc, cc
    allocate(Vc(n_cells))
    allocate(cc(n_cells, 2))
    do i = 1, n_cells
        n1 = cells(i, 1)
        n2 = cells(i, 2)
        n3 = cells(i, 3)
        p1x = nodes(n1, 1)
        p1y = nodes(n1, 2)
        p2x = nodes(n2, 1)
        p2y = nodes(n2, 2)
        p3x = nodes(n3, 1)
        p3y = nodes(n3, 2)

        Vc(i) = .5*abs((p1x*(p2y - p3y) + p2x*(p3y - p1y) + p3x*(p1y - p2y)))
        cc(i, :) = [(p1x + p2x + p3x)/3, (p1y + p2y + p3y)/3]
    end do

    ! Sf, nf, cf, df, ef, gf
    boundary_cells = [integer ::]
    interior_cells = [integer ::]
    allocate(cell_cells(n_cells, 3))
    allocate(cell_adjacent_faces(n_cells, 3))
    allocate(Sf(n_cells, 3))
    allocate(nf(n_cells, 3, 2))
    allocate(cf(n_cells, 3, 2))
    allocate(df(n_cells, 3, 2))
    allocate(ef(n_cells, 3, 2))
    allocate(gf(n_cells, 3))
    do i = 1, n_cells
        do j = 1, 3
            fj = cell_faces(i, j)
            
            if (face_cells(fj, 1) == face_cells(fj, 2)) then
                cj = face_cells(fj, 1)
                cell_cells(i, j) = cj
                sign = 1
            else
                if (face_cells(fj, 1) == i) then
                    cj = face_cells(fj, 2)
                    cell_cells(i, j) = cj
                    sign = -1
                else
                    cj = face_cells(fj, 1)
                    cell_cells(i, j) = cj
                    sign = 1
                end if
            end if
            
            l3 = [1, 2, 3]
            l3 = pack(l3, cell_faces(cj, :) == fj)
            cell_adjacent_faces(i, j) = l3(1)

            n1 = faces(fj, 1)
            n2 = faces(fj, 2)
            p1x = nodes(n1, 1)
            p1y = nodes(n1, 2)
            p2x = nodes(n2, 1)
            p2y = nodes(n2, 2)

            Sf(i, j) = sqrt((p2x - p1x)**2 + (p2y - p1y)**2)
            nf(i, j, :) = sign*[(p2y - p1y)/Sf(i, j), -(p2x - p1x)/Sf(i, j)]
            cf(i, j, :) = [.5*(p1x + p2x), .5*(p1y + p2y)]
            df(i, j, :) = cf(i, j, :) - cc(i, :)
            ef(i, j, :) = cc(cj, :) - cc(i, :)
            gf(i, j) = Vc(i)/(Vc(i) + Vc(cj))
        end do

        if (any(cell_cells(i, :) == i)) then
            boundary_cells = [boundary_cells, i]
        else
            interior_cells = [interior_cells, i]
        end if
    end do
    n_boundary_cells = size(boundary_cells)
    n_interior_cells = size(interior_cells)

! mesh complete

    allocate(Q(n_cells, 4))
    allocate(grad_Q(n_cells, 4, 2))
    allocate(Qf(n_cells, 3, 4))
    allocate(Qfb(n_cells, 3, 4))
    allocate(Co(n_cells))

    Qi = [r, r*u, r*v, .5*r*(u*u + v*v) + p/(g - 1)] ! inlet/initial conditions
    Q(:, 1) = Qi(1)
    Q(:, 2) = Qi(2)
    Q(:, 3) = Qi(3)
    Q(:, 4) = Qi(4)

! loop start

    do n = 1, int(t/dt) + 1

        ! Courant number
        
        do i = 1, n_cells

            Co(i) = 0.

            do j = 1, 3
                
                Vf = [Qf(i, j, 2)/Qf(i, j, 1), Qf(i, j, 3)/Qf(i, j, 1)]
                Co(i) = Co(i) + .5*dt/Vc(i)*abs(dot_product(Vf, nf(i, j, :)))*Sf(i, j)

            end do

        end do

        print *, int(t/dt) + 1, n, minval(Q(:, 1)), maxval(Q(:, 1)), maxval(Co)

        ! Step 1 : Reconstruct

        ! Step 1a : Compute gradients

        do i = 1, n_interior_cells
            
            ci = interior_cells(i)

            grad_Q(ci, :, :) = 0.

            do j = 1, 3

                cj = cell_cells(ci, j)

                do k = 1, 4

                    Qf(ci, j, k) = gf(ci, j)*Q(cj, k) + (1. - gf(ci, j))*Q(ci, k)
                    grad_Q(ci, k, :) = grad_Q(ci, k, :) + Qf(ci, j, k)*Sf(ci, j)*nf(ci, j, :)/Vc(ci)

                end do

            end do

        end do

        do i = 1, n_boundary_cells
            
            ci = boundary_cells(i)

            grad_Q(ci, :, :) = 0.

            do j = 1, 3

                cj = cell_cells(ci, j)
                fj = cell_faces(ci, j)

                if (any(fj == m_inlet_faces)) then

                    do k = 1, 4

                        Qf(ci, j, k) = gf(ci, j)*Qi(k) + (1. - gf(ci, j))*Q(ci, k)
                        grad_Q(ci, k, :) = grad_Q(ci, k, :) + Qf(ci, j, k)*Sf(ci, j)*nf(ci, j, :)/Vc(ci)

                    end do

                else if (any(fj == m_body_faces) .or. any(fj == m_wall_faces)) then

                    Qb = transform_array(Q(ci, :), nf(ci, j, :), 1)
                    Qb(2) = -Qb(2)
                    Qb = transform_array(Qb, nf(ci, j, :), -1)

                    do k = 1, 4

                        Qf(ci, j, k) = gf(ci, j)*Qb(k) + (1. - gf(ci, j))*Q(ci, k)
                        grad_Q(ci, k, :) = grad_Q(ci, k, :) + Qf(ci, j, k)*Sf(ci, j)*nf(ci, j, :)/Vc(ci)

                    end do

                else

                    do k = 1, 4

                        Qf(ci, j, k) = gf(ci, j)*Q(cj, k) + (1. - gf(ci, j))*Q(ci, k)
                        grad_Q(ci, k, :) = grad_Q(ci, k, :) + Qf(ci, j, k)*Sf(ci, j)*nf(ci, j, :)/Vc(ci)

                    end do

                end if
                
            end do

        end do
        
        ! Steps 1b & 1c : Compute limiters & boundary extrapolated values

        do i = 1, n_interior_cells

            ci = interior_cells(i)

            do j = 1, 3

                cj = cell_cells(ci, j)

                do k = 1, 4

                    if (Q(cj, k) - Q(ci, k) /= 0.) then
                        rk = dot_product(grad_Q(ci, k, :), ef(ci, j, :))/(Q(cj, k) - Q(ci, k))
                        if (rk <= 0.) then
                            xk = 0.
                        else if (rk > 0. .and. rk <= .5) then
                            xk = 2*rk
                        else if (rk > .5 .and. rk <= 1.) then
                            xk = 1.
                        else
                            xk = min(rk, 2./(1. + rk), 2.)
                        end if
                    else
                        xk = 0.
                    end if
                    
                    Qf(ci, j, k) = Q(ci, k) + xk*dot_product(grad_Q(ci, k, :), df(ci, j, :))

                end do

            end do

        end do

        do i = 1, n_boundary_cells

            ci = boundary_cells(i)

            do j = 1, 3

                cj = cell_cells(ci, j)
                fj = cell_faces(ci, j)

                if (any(fj == m_inlet_faces)) then
                    
                    do k = 1, 4

                        if (Qi(k) - Q(ci, k) /= 0.) then
                            rk = dot_product(grad_Q(ci, k, :), 2*df(ci, j, :))/(Qi(k) - Q(ci, k))
                            if (rk <= 0.) then
                                xk = 0.
                            else if (rk > 0. .and. rk <= .5) then
                                xk = 2*rk
                            else if (rk > .5 .and. rk <= 1.) then
                                xk = 1.
                            else
                                xk = min(rk, 2./(1. + rk), 2.)
                            end if
                        else
                            xk = 0.
                        end if
                        
                        Qf(ci, j, k) = Q(ci, k) + xk*dot_product(grad_Q(ci, k, :), df(ci, j, :))

                    end do

                else if (any(fj == m_body_faces) .or. any(fj == m_wall_faces)) then

                    Qb = transform_array(Q(ci, :), nf(ci, j, :), 1)
                    Qb(2) = -Qb(2)
                    Qb = transform_array(Qb, nf(ci, j, :), -1)

                    do k = 1, 4

                        if (Qb(k) - Q(ci, k) /= 0.) then
                            rk = dot_product(grad_Q(ci, k, :), 2*df(ci, j, :))/(Qb(k) - Q(ci, k))
                            if (rk <= 0.) then
                                xk = 0.
                            else if (rk > 0. .and. rk <= .5) then
                                xk = 2*rk
                            else if (rk > .5 .and. rk <= 1.) then
                                xk = 1.
                            else
                                xk = min(rk, 2./(1. + rk), 2.)
                            end if
                        else
                            xk = 0.
                        end if
                        
                        Qf(ci, j, k) = Q(ci, k) + xk*dot_product(grad_Q(ci, k, :), df(ci, j, :))

                    end do

                else

                    do k = 1, 4

                        if (Q(cj, k) - Q(ci, k) /= 0.) then
                            rk = dot_product(grad_Q(ci, k, :), ef(ci, j, :))/(Q(cj, k) - Q(ci, k))
                            if (rk <= 0.) then
                                xk = 0.
                            else if (rk > 0. .and. rk <= .5) then
                                xk = 2*rk
                            else if (rk > .5 .and. rk <= 1.) then
                                xk = 1.
                            else
                                xk = min(rk, 2./(1. + rk), 2.)
                            end if
                        else
                            xk = 0.
                        end if
                        
                        Qf(ci, j, k) = Q(ci, k) + xk*dot_product(grad_Q(ci, k, :), df(ci, j, :))

                    end do

                end if

            end do

        end do

        ! Step 2 : Evolve

        do ci = 1, n_cells

            Qfb(ci, :, :) = Qf(ci, :, :)

            do j = 1, 3
                
                Qfh = transform_array(Qf(ci, j, :), nf(ci, j, :), 1)
                Ffh = to_flux(Qfh)

                Qfb(ci, 1, :) = Qfb(ci, 1, :) - .5*dt/Vc(ci)*Sf(ci, j)*transform_array(Ffh, nf(ci, j, :), -1)
                Qfb(ci, 2, :) = Qfb(ci, 2, :) - .5*dt/Vc(ci)*Sf(ci, j)*transform_array(Ffh, nf(ci, j, :), -1)
                Qfb(ci, 3, :) = Qfb(ci, 3, :) - .5*dt/Vc(ci)*Sf(ci, j)*transform_array(Ffh, nf(ci, j, :), -1)
            
            end do

        end do

        ! Step 3 : Solve

        do i = 1, n_interior_cells

            ci = interior_cells(i)
            
            do j = 1, 3
                
                cj = cell_cells(ci, j)
                jj = cell_adjacent_faces(ci, j) ! adjacent cell local face index

                QfL = transform_array(Qfb(ci, j, :), nf(ci, j, :), 1)
                QfR = transform_array(Qfb(cj, jj, :), nf(ci, j, :), 1)
                Ff = get_flux(QfL, QfR)
                
                Q(ci, :) = Q(ci, :) - dt/Vc(ci)*Sf(ci, j)*transform_array(Ff, nf(ci, j, :), -1)
            
            end do

        end do

        do i = 1, n_boundary_cells

            ci = boundary_cells(i)

            do j = 1, 3

                cj = cell_cells(ci, j)
                fj = cell_faces(ci, j)
                jj = cell_adjacent_faces(ci, j) ! adjacent cell local face index

                if (any(fj == m_inlet_faces)) then

                    QfL = transform_array(Qfb(ci, j, :), nf(ci, j, :), 1)
                    QfR = transform_array(Qi(:), nf(ci, j, :), 1)
                    Ff = get_flux(QfL, QfR)

                else if (any(fj == m_body_faces) .or. any(fj == m_wall_faces)) then

                    QfL = transform_array(Qfb(ci, j, :), nf(ci, j, :), 1)
                    QfR = transform_array(Qfb(ci, j, :), nf(ci, j, :), 1)
                    QfR(2) = -QfR(2)
                    Ff = get_flux(QfL, QfR)

                else

                    QfL = transform_array(Qfb(ci, j, :), nf(ci, j, :), 1)
                    QfR = transform_array(Qfb(cj, jj, :), nf(ci, j, :), 1)
                    Ff = get_flux(QfL, QfR)

                end if

                Q(ci, :) = Q(ci, :) - dt/Vc(ci)*Sf(ci, j)*transform_array(Ff, nf(ci, j, :), -1)

            end do

        end do

        if (mod(n, nt) == 0) then

            print *, 'save'

            write(file_name, '(a, i10.10, a)') trim(file), n, '.txt'
            open(100, file=trim(folder)//trim(file_name))
            
            do i = 1, n_cells
                write(100, *) cc(i, 1), cc(i, 2), Q(i, 1)
            end do

        end if

    end do

! loop complete

contains

    subroutine append_array(a1, a2, b1, b2)
    
        integer, allocatable, intent(in out) :: a1(:), a2(:)
        integer, intent(in) :: b1, b2
        
        integer, allocatable :: c1(:), c2(:)
        logical :: in

        allocate(c1, mold=a1)
        allocate(c2, mold=a2)
        c1(:) = b1
        c2(:) = b2

        in = any((c1 == a1) .and. (c2 == a2)) .or. any((c2 == a1) .and. (c1 == a2))
        if (.not. in) then
            a1 = [a1, b1]
            a2 = [a2, b2]
        end if

    end subroutine append_array

    pure function compare_array(a1, a2, b1, b2)
    
        integer, intent(in) :: a1, a2, b1, b2
        logical :: compare_array

        if ((b1 == a1 .and. b2 == a2) .or. (b2 == a1 .and. b1 == a2)) then
            compare_array = .true.
        else
            compare_array = .false.
        end if

    end function compare_array

    pure function transform_array(Q, n, s) result(Qh)
        
        real, intent(in) :: Q(:)
        real, intent(in) :: n(:)
        integer, intent(in) :: s ! setting

        real, dimension(4) :: Qh

        Qh(1) = Q(1)
        Qh(4) = Q(4)
        if (s == 1) then
            Qh(2) = n(1)*Q(2) + n(2)*Q(3)
            Qh(3) = -n(2)*Q(2) + n(1)*Q(3)
        else
            Qh(2) = n(1)*Q(2) - n(2)*Q(3)
            Qh(3) = n(2)*Q(2) + n(1)*Q(3)
        end if

    end function transform_array

    pure function to_flux(Q) result(F)

        real, intent(in) :: Q(:)

        real :: r, u, v, E, p
        real, dimension(4) :: F

        r = Q(1)
        u = Q(2)/r
        v = Q(3)/r
        E = Q(4)
        p = (g - 1)*(E - .5*r*(u*u + v*v))

        F = [r*u, r*u*u + p, r*u*v, u*(E + p)]

    end function to_flux

    pure function get_flux(QfL, QfR) result(Ff)
        
        real, intent(in) :: QfL(:), QfR(:)

        real :: rL, uL, vL, pL, EL
        real :: rR, uR, vR, pR, ER
        real :: aL, aR, sL, sR
        real :: r_bar, a_bar, p_bar, p_star, s_star, wL, wR
        real, dimension(4) :: QL_star, QR_star, FL_star, FR_star
        real, dimension(4) :: QL, QR, FL, FR
        real, dimension(4) :: Ff

        rL = QfL(1)
        uL = QfL(2)/rL
        vL = QfL(3)/rL
        EL = QfL(4)
        pL = (g - 1)*(EL - .5*rL*(uL*uL + vL*vL))

        rR = QfR(1)
        uR = QfR(2)/rR
        vR = QfR(3)/rR
        ER = QfR(4)
        pR = (g - 1)*(ER - .5*rR*(uR*uR + vR*vR))

        ! compute pressure estimate
        aL = sqrt(g*pL/rL)
        aR = sqrt(g*pR/rR)

        r_bar = .5*(rL + rR)
        a_bar = .5*(aL + aR)
        p_bar = .5*(pL + pR) - .5*(uR - uL)*r_bar*a_bar
        p_star = max(0., p_bar)

        ! compute wave speeds
        if (p_star <= pL) then
            wL = 1.
        else
            wL = sqrt(1 + (g + 1)/(2*g)*(p_star/pL - 1))
        end if
        
        if (p_star <= pR) then
            wR = 1.
        else
            wR = sqrt(1 + (g + 1)/(2*g)*(p_star/pR - 1))
        end if

        sL = uL - aL*wL
        sR = uR + aR*wR
        s_star = (pR - pL + rL*uL*(sL - uL) - rR*uR*(sR - uR))/(rL*(sL - uL) - rR*(sR - uR))

        ! compute intercell flux
        QL = (/rL, rL*uL, rL*vL, EL/)
        QR = (/rR, rR*uR, rR*vR, ER/)
        FL = (/rL*uL, rL*uL*uL + pL, rL*uL*vL, uL*(EL + pL)/)
        FR = (/rR*uR, rR*uR*uR + pR, rR*uR*vR, uR*(ER + pR)/)

        QL_star = rL*(sL - uL)/(sL - s_star)*(/1., s_star, vL, EL/rL + (s_star - uL)*(s_star + pL/rL/(sL - uL))/)
        QR_star = rR*(sR - uR)/(sR - s_star)*(/1., s_star, vR, ER/rR + (s_star - uR)*(s_star + pR/rR/(sR - uR))/)
        FL_star = FL + sL*(QL_star - QL)
        FR_star = FR + sR*(QR_star - QR)
    
        if (sL > 0.) then
            Ff = FL
        else if (sL <= 0. .and. s_star >= 0.) then
            Ff = FL_star
        else if (s_star <= 0. .and. sR >= 0.) then
            Ff = FR_star
        else if (sR < 0.) then
            Ff = FR
        end if

    end function get_flux

end program euler
