module mod_config
    use iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: Configuration

    type :: Configuration
        character(20) :: config_file
        real(real64) :: M
        real(real64) :: a
        real(real64) :: g
        real(real64) :: dt
        integer(int32) :: nt
        integer(int32) :: nw
        real(real64) :: CFL
        integer(int32) :: nRK
        real(real64), dimension(:), allocatable :: RK
    end type Configuration

    interface Configuration
        module procedure :: constructor
    end interface Configuration
    
contains

    type(Configuration) function constructor(config_file) result(res)
        character(*), intent(in) :: config_file
        character(20) :: line
        res % config_file = config_file
        open(10, file=trim(config_file), action='read')
        read(10, *) line, res % M
        read(10, *) line, res % a
        read(10, *) line, res % g
        read(10, *) line, res % dt
        read(10, *) line, res % nt
        read(10, *) line, res % nw
        read(10, *) line, res % CFL
        read(10, *) line, res % nRK
        allocate(res % RK(res % nRK))
        if (res % nRK == 4) then
            res % RK = [.0833, .2069, .4265, 1.]
        else if (res % nRK == 5) then
            res % RK = [.0533, .1263, .2375, .4414, 1.]
        else
            error stop 'Error: invalid Runge-Kutta stages'
        end if
        close(10)
    end function constructor
    
end module mod_config
