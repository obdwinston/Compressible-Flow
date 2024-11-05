module mod_utils
    use iso_fortran_env, only: int32, real64, stdout => output_unit
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    implicit none

    private
    public :: &
    get_limiter, &
    get_primitive, &
    get_conserved, &
    get_flux, &
    get_transform, &
    get_inverse_transform, &
    set_field, &
    set_surface, &
    set_elapsed, &
    set_duration, &
    set_data
    
contains

    pure elemental function get_limiter(a, b) result(res)
        real(real64), intent(in) :: a, b
        real(real64) :: res
        integer(int32) :: i
        real(real64) :: numerator, denominator
        if (a*b > 0.) then
            numerator = a*(b**2 + 1e-15) + b*(a**2 + 1e-15)
            denominator = a**2 + b**2 + 2e-15
            res = numerator/denominator
        else
            res = 0.
        end if
    end function get_limiter

    pure function get_primitive(U, g) result(res)
        real(real64), dimension(4), intent(in) :: U
        real(real64), intent(in) :: g
        real(real64), dimension(4) :: res
        real(real64) :: rho, ux, uy, p
        rho = U(1)
        ux = U(2)/rho
        uy = U(3)/rho
        p = (g - 1.)*(U(4) - .5*rho*(ux**2 + uy**2))
        res = [rho, ux, uy, p]
    end function get_primitive
    
    pure function get_conserved(W, g) result(res)
        real(real64), dimension(4), intent(in) :: W
        real(real64), intent(in) :: g
        real(real64), dimension(4) :: res
        res(1) = W(1)
        res(2) = W(1)*W(2)
        res(3) = W(1)*W(3)
        res(4) = .5*W(1)*(W(2)**2 + W(3)**2) + W(4)/(g - 1.)
    end function get_conserved
    
    pure function get_flux(U, g) result(res)
        real(real64), dimension(4), intent(in) :: U
        real(real64), intent(in) :: g
        real(real64), dimension(4) :: res, W
        W = get_primitive(U, g)
        res(1) = W(1)*W(2)
        res(2) = W(1)*W(2)*W(2) + W(4)
        res(3) = W(1)*W(2)*W(3)
        res(4) = W(2)*(U(4) + W(4))
    end function get_flux

    pure function get_transform(U, n) result(res)
        real(real64), dimension(4), intent(in) :: U
        real(real64), dimension(2), intent(in) :: n
        real(real64), dimension(4) :: res
        real(real64), dimension(4, 4) :: T
        T(1, :) = [1., 0., 0., 0.]
        T(2, :) = [0._real64, n(1), n(2), 0._real64]
        T(3, :) = [0._real64, -n(2), n(1), 0._real64]
        T(4, :) = [0., 0., 0., 1.]
        res = matmul(T, U)
    end function get_transform

    pure function get_inverse_transform(U, n) result(res)
        real(real64), dimension(4), intent(in) :: U
        real(real64), dimension(2), intent(in) :: n
        real(real64), dimension(4) :: res
        real(real64), dimension(4, 4) :: T
        T(1, :) = [1., 0., 0., 0.]
        T(2, :) = [0._real64, n(1), -n(2), 0._real64]
        T(3, :) = [0._real64, n(2), n(1), 0._real64]
        T(4, :) = [0., 0., 0., 1.]
        res = matmul(T, U)
    end function get_inverse_transform

    subroutine set_field(msh)
        type(Mesh), intent(in) :: msh
        integer(int32) :: i
        open(10, file='data/cxy.txt') ! field data coordinates
        do i = 1, msh % n_cells
            write(10, *) msh % cells(i) % coordinates(1), &
            msh % cells(i) % coordinates(2)
        end do
        close(10)
    end subroutine set_field

    subroutine set_surface(msh)
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, ni
        open(10, file='data/nxy.txt') ! surface data coordinates
        do i = 1, size(msh % surface_nodes)
            ni = msh % surface_nodes(i)
            write(10, *) msh % nodes(ni) % coordinates(1), &
            msh % nodes(ni) % coordinates(2)
        end do
        close(10)
    end subroutine set_surface
    
    subroutine set_elapsed(ts, te, cs, ce, cr)
        real(real64), intent(in) :: ts, te
        integer(int32), intent(in) :: cs, ce, cr
        write(stdout, *) 'Elapsed time per iteration (cpu_time):', te - ts
        write(stdout, *) 'Elapsed time per iteration (system_clock):', &
        real(ce - cs, real64)/real(cr, real64)
    end subroutine set_elapsed
    
    subroutine set_duration(nt, dt, ts, te, n)
        integer(int32), intent(in) :: nt, n
        real(real64), intent(in) :: dt, ts, te
        integer(int32) :: h, m, s
        real(real64) :: t
        t = (te - ts)*(nt - n)
        h = int(t)/3600
        m = mod(int(t)/60, 60)
        s = mod(int(t), 60)
        write(stdout, *) 'Iteration', n, 'of', nt
        write(stdout, *) 'Timestep:', dt
        write(stdout, *) 'Estimated time remaining (h:m:s):', h, m, s
        write(stdout, *)
    end subroutine set_duration

    subroutine set_data(Uc, Un, msh, g, n)
        real(real64), dimension(:, :), intent(in) :: Uc, Un
        type(Mesh), intent(in) :: msh
        real(real64), intent(in) :: g
        integer(int32), intent(in) :: n
        call set_cell_data(Uc, msh, g, n)
        call set_node_data(Un, msh, g, n)
    end subroutine set_data

    subroutine set_cell_data(Uc, msh, g, n)
        real(real64), dimension(:, :), intent(in) :: Uc
        type(Mesh), intent(in) :: msh
        real(real64), intent(in) :: g
        integer(int32), intent(in) :: n
        character(50) :: file
        integer(int32) :: i
        write(file, '(a, i10.10, a)') 'Wc_', n, '.txt'
        open(10, file='data/'//trim(file))
        do i = 1, msh % n_cells
            write(10, *) get_primitive(Uc(i, :), g)
        end do
        close(10)
    end subroutine set_cell_data

    subroutine set_node_data(Un, msh, g, n)
        real(real64), dimension(:, :), intent(in) :: Un
        type(Mesh), intent(in) :: msh
        real(real64), intent(in) :: g
        integer(int32), intent(in) :: n
        character(50) :: file
        integer(int32) :: i, ni
        write(file, '(a, i10.10, a)') 'Wn_', n, '.txt'
        open(10, file='data/'//trim(file))
        do i = 1, size(msh % surface_nodes)
            ni = msh % surface_nodes(i)
            write(10, *) get_primitive(Un(ni, :), g)
        end do
        close(10)
    end subroutine set_node_data

end module mod_utils
