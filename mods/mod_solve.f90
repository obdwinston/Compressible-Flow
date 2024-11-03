module mod_solve
    use iso_fortran_env, only: int32, real64
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    use mod_flux, only: get_HLLC_flux
    use mod_utils, only: &
    get_limiter, &
    get_primitive, &
    get_conserved, &
    get_transform, &
    get_inverse_transform
    implicit none

    private
    public :: &
    set_initial_conditions, &
    set_timestep, &
    set_node_values, &
    set_face_values, &
    set_face_fluxes, &
    set_cell_values
    
contains

    subroutine set_initial_conditions(U0, Uc, msh, config)
        real(real64), dimension(4), intent(in out) :: U0
        real(real64), dimension(:, :), intent(in out) :: Uc
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        real(real64) :: M, a, g
        M = config % M
        a = config % a
        g = config % g
        U0 = [1._real64, M*cosd(a), M*sind(a), M**2/2. + 1./(g*(g - 1.))]
        Uc = spread(U0, 1, msh % n_cells)
    end subroutine set_initial_conditions

    subroutine set_timestep(config, Uc, msh)
        type(Configuration), intent(in out) :: config
        real(real64), dimension(:, :), intent(in) :: Uc
        type(Mesh), intent(in) :: msh
        integer(int32) :: i
        real(real64) :: dt, a, u
        real(real64), dimension(4) :: W
        if (config % CFL /= 0.) then
            config % dt = 1.
            do concurrent(i = 1:msh % n_cells)
                W = get_primitive(Uc(i, :), config % g)
                a = sqrt(config % g*W(4)/W(1))
                u = sqrt(W(2)**2 + W(3)**2)
                dt = config % CFL*msh % cells(i) % cell_volume/(u + a)
                if (config % dt > dt) config % dt = dt
            end do
        end if
    end subroutine
    
    subroutine set_node_values(Un, Uc, msh)
        real(real64), dimension(:, :), intent(in out) :: Un
        real(real64), dimension(:, :), intent(in) :: Uc
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, j, cj
        real(real64) :: wnj
        real(real64), dimension(4) :: Unj
        do concurrent(i = 1:msh % n_nodes)
            Unj = [0., 0., 0., 0.]
            do concurrent(j = 1:msh % nodes(i) % n_node_cells)
                cj = msh % nodes(i) % node_cells(j)
                wnj = msh % nodes(i) % node_cell_weights(j)
                Unj = Unj + wnj*Uc(cj, :)
            end do
            Un(i, :) = Unj
        end do
    end subroutine set_node_values

    subroutine set_face_values(UL, Un, Uc, msh)
        real(real64), dimension(:, :, :), intent(in out) :: UL
        real(real64), dimension(:, :), intent(in) :: Un, Uc
        type(Mesh), intent(in) :: msh
        integer(int32) :: i, j, fj, na, nb, nc
        real(real64) :: rf, rn
        real(real64), dimension(4) :: dUf, dUn
        do concurrent(i = 1:msh % n_cells)
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                na = msh % faces(fj) % face_nodes(1)
                nb = msh % faces(fj) % face_nodes(2)
                nc = msh % cells(i) % cell_opposite_nodes(j)
                rf = msh % cells(i) % cell_face_distances(j)
                rn = msh % cells(i) % cell_node_distances(j)
                dUf = (.5*(Un(na, :) + Un(nb, :)) - Uc(i, :))/rf
                dUn = (Uc(i, :) - Un(nc, :))/rn
                UL(i, j, :) = Uc(i, :) + rf*get_limiter(dUf, dUn)
            end do
        end do
    end subroutine set_face_values

    subroutine set_face_fluxes(F, UL, U0, msh, config)
        real(real64), dimension(:, :), intent(in out) :: F
        real(real64), dimension(:, :, :), intent(in) :: UL
        real(real64), dimension(4), intent(in) :: U0
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32) :: i, fi, fj, ci, cj
        real(real64), dimension(2) :: nf
        real(real64), dimension(4) :: ULf, URf, UL_hat, UR_hat
        do concurrent(i = 1:msh % n_faces)
            ci = msh % faces(i) % face_cells(1)
            cj = msh % faces(i) % face_cells(2)
            fi = msh % faces(i) % face_faces(1) ! cell local index
            fj = msh % faces(i) % face_faces(2) ! cell local index
            nf = msh % faces(i) % face_normal
            ULf = UL(ci, fi, :)
            UL_hat = get_transform(ULf, nf)
            ! interior faces
            if (trim(msh % faces(i) % face_type) == 'INTERIOR') then
                URf = UL(cj, fj, :)
                UR_hat = get_transform(URf, nf)
            ! slipwall faces
            else if (trim(msh % faces(i) % face_type) == 'SLIPWALL') then
                UR_hat = [UL_hat(1), -UL_hat(2), -UL_hat(3), UL_hat(4)]
            ! freestream faces
            else if (trim(msh % faces(i) % face_type) == 'FREESTREAM') then
                UR_hat = get_transform(U0, nf)
            end if
            F(i, :) = get_HLLC_flux(UL_hat, UR_hat, config % g)
        end do
    end subroutine set_face_fluxes

    subroutine set_cell_values(Uc, Uc0, F, msh, config, k)
        real(real64), dimension(:, :), intent(in out) :: Uc, Uc0
        real(real64), dimension(:, :), intent(in) :: F
        type(Mesh), intent(in) :: msh
        type(Configuration), intent(in) :: config
        integer(int32), intent(in) :: k        
        integer(int32) :: i, j, fj
        real(real64) :: Vc, Sf, sign
        real(real64), dimension(2) :: nf
        real(real64), dimension(4) :: Ff
        Uc0(:, :) = Uc(:, :)
        do concurrent(i = 1:msh % n_cells)
            Vc = msh % cells(i) % cell_volume
            Ff = [0., 0., 0., 0.]
            do concurrent(j = 1:3)
                fj = msh % cells(i) % cell_faces(j)
                Sf = msh % faces(fj) % face_area
                nf = msh % faces(fj) % face_normal
                sign = msh % cells(i) % cell_face_signs(j)
                Ff = Ff + sign*get_inverse_transform(F(fj, :), nf)*Sf
            end do
            Uc(i, :) = Uc0(i, :) - config % RK(k)*config % dt/Vc*Ff
        end do
    end subroutine set_cell_values

end module mod_solve
