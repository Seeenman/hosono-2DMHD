module eos

    use simulation, only: sim_smallPressure

    implicit none

contains

    function eos_presIdealGas(dens, eint, gama) result(pres)
        implicit none
        real, intent(in) :: dens, eint, gama
        real :: pres

        pres = (gama-1)*dens*eint
        if (pres < sim_smallPressure) then
            print*, "-----------------------------------------------------------------------"
            print*, "DEBUG"
            print*, "Negative pressure in eos_presIdealGas"
            print*, "pressure value of ", pres
            print*, "pressure set to", sim_smallPressure
            print*, "GUBED"
            print*, "-----------------------------------------------------------------------"
            pres = sim_smallPressure ! preserve positivity in a crude manner
        end if

    end function eos_presIdealGas

    function eos_eintIdealGas(pres, dens, gama) result(eint)
        implicit none
        real, intent(in) :: pres, dens, gama
        real :: eint

        eint = pres/(gama-1.0)/dens

    end function eos_eintIdealGas

end module eos
