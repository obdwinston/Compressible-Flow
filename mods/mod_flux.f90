module mod_flux
    use iso_fortran_env, only: int32, real64
    use mod_utils, only: get_primitive, get_flux
    implicit none

    private
    public :: get_HLLC_flux
    
contains

    pure function get_HLLC_flux(UL, UR, g) result(res)
        real(real64), dimension(4), intent(in) :: UL, UR ! transformed
        real(real64), intent(in) :: g
        real(real64), dimension(4) :: res
        real(real64) :: pS, sS, qL, qR, sL, sR
        real(real64), dimension(4) :: WL, WR, FL, FR, FSL, FSR
        ! estimate pressure
        WL = get_primitive(UL, g)
        WR = get_primitive(UR, g)
        pS = get_pressure_estimate(WL, WR, g)
        ! estimate wave speed
        qL = get_speed_coefficient(WL(4), pS, g)
        qR = get_speed_coefficient(WR(4), pS, g)
        sL = WL(2) - qL*sqrt(g*WL(4)/WL(1))
        sR = WR(2) + qR*sqrt(g*WR(4)/WR(1))
        sS = get_wave_speed_estimate(WL, WR, sL, sR, pS)
        ! compute HLLC flux
        FL = get_flux(UL, g)
        FR = get_flux(UR, g)
        FSL = get_intermediate_flux(FL, UL, WL, sL, sS)
        FSR = get_intermediate_flux(FR, UR, WR, sR, sS)
        if (sL >= 0.) then
            res = FL
        else if (sL <= 0. .and. sS >= 0.) then
            res = FSL
        else if (sS <= 0. .and. sR >= 0.) then
            res = FSR
        else if (sR <= 0.) then
            res = FR
        else
            error stop 'Error: invalid HLLC flux'
        end if
    end function get_HLLC_flux

    pure function get_pressure_estimate(WL, WR, g) result(res)
        real(real64), dimension(4), intent(in) :: WL, WR ! transformed
        real(real64), intent(in) :: g
        real(real64) :: res, rho_bar, a_bar, p_bar
        rho_bar = .5*(WL(1) + WR(1))
        a_bar = .5*(sqrt(g*WL(4)/WL(1)) + sqrt(g*WR(4)/WR(1)))
        p_bar = .5*((WL(4) + WR(4)) - rho_bar*a_bar*(WR(2) - WL(2)))
        res = max(0., p_bar)
    end function get_pressure_estimate

    pure function get_speed_coefficient(p, pS, g) result(res)
        real(real64), intent(in) :: p, pS, g
        real(real64) :: res
        if (pS > p) then
            res = sqrt(1. + (g + 1.)/2./g*(pS/p - 1.))
        else
            res = 1.
        end if
    end function get_speed_coefficient
    
    pure function get_wave_speed_estimate(WL, WR, sL, sR, pS) result(res)
        real(real64), dimension(4), intent(in) :: WL, WR ! transformed
        real(real64), intent(in) :: sL, sR, pS
        real(real64) :: res, numerator, denominator
        numerator = WR(4) - WL(4) + WL(1)*WL(2)*(sL - WL(2)) - WR(1)*WR(2)*(sR - WR(2))
        denominator = WL(1)*(sL - WL(2)) - WR(1)*(sR - WR(2))
        if (abs(numerator) < 1e-15) then ! IMPORTANT
            res = 0.
        else
            res = numerator/denominator
        end if
    end function get_wave_speed_estimate
    
    pure function get_intermediate_flux(F, U, W, s, sS) result(res)
        real(real64), dimension(4), intent(in) :: F, U, W ! transformed
        real(real64), intent(in) :: s, sS
        real(real64) :: scalar
        real(real64), dimension(4) :: res, vector
        scalar = W(1)*(s - W(2))/(s - sS)
        vector(1) = 1.
        vector(2) = sS
        vector(3) = W(3)
        vector(4) = U(4)/W(1) + (sS - W(2))*(sS + W(4)/W(1)/(s - W(2)))
        res = F + s*(scalar*vector - U)
    end function get_intermediate_flux
    
end module mod_flux
