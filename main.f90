program main
    use iso_fortran_env, only: int32, real64, stdout => output_unit
    use mod_mesh, only: Mesh
    use mod_config, only: Configuration
    use mod_solve, only: &
    set_initial_conditions, &
    set_timestep, &
    set_node_values, &
    set_face_values, &
    set_face_fluxes, &
    set_cell_values
    use mod_utils, only: &
    set_field, &
    set_surface, &
    set_elapsed, &
    set_duration, &
    set_data
    implicit none

    type(Mesh) :: msh
    type(Configuration) :: config
    integer(int32) :: n, k, cs, ce, cr
    real(real64) :: ts, te
    real(real64), dimension(4) :: U0
    real(real64), dimension(:, :), allocatable :: Un, Uc, Uc0, F
    real(real64), dimension(:, :, :), allocatable :: UL
    
    msh = Mesh()
    call msh % load_mesh('mesh/mesh.txt')
    call set_field(msh)
    call set_surface(msh)
    config = Configuration('config.txt')
    
    allocate(Un(msh % n_nodes, 4))
    allocate(UL(msh % n_cells, 3, 4))
    allocate(Uc(msh % n_cells, 4))
    allocate(Uc0(msh % n_cells, 4))
    allocate(F(msh % n_faces, 4))

    call set_initial_conditions(U0, Uc, msh, config)
    call set_timestep(config, Uc, msh)

    write(stdout, *) 'Running simulation...'
    time_step: do n = 1, config % nt
        call cpu_time(ts)
        call system_clock(cs, cr)

        call set_timestep(config, Uc, msh)
        
        time_stage: do k = 1, size(config % RK)
            call set_node_values(Un, Uc, msh)
            call set_face_values(UL, Un, Uc, msh)
            call set_face_fluxes(F, UL, U0, msh, config)
            call set_cell_values(Uc, Uc0, F, msh, config, k)
        end do time_stage

        call cpu_time(te)
        call system_clock(ce)

        if (mod(n, config % nw) == 0) then
            call set_elapsed(ts, te, cs, ce, cr)
            call set_duration(config % nt, config % dt, ts, te, n)
            call set_data(Uc, Un, msh, config % g, n)
        end if
    end do time_step
    write(stdout, *) 'Done!'
end program main
