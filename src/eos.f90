module eos

    use simulation, only: sim_smallPressure

    implicit none

contains

    function eos_presIdealGas(dens, eint, gama) result(pres)
        implicit none
        real, intent(in) :: dens, eint, gama
        real :: pres

        pres = (gama-1)*dens*eint ! basically checking if the internal energy we were given was negative
        if (pres < sim_smallPressure) then
            print*, "-----------------------------------------------------------------------"
            print*, "DEBUG"
            print*, "Negative pressure in eos_presIdealGas"
            print*, "pressure value of ", pres
            print*, "pressure set to", sim_smallPressure
            print*, "GUBED"
            print*, "-----------------------------------------------------------------------"
            pres = sim_smallPressure ! floor the pressure. 
        end if

    end function eos_presIdealGas

    pure function eos_eintIdealGas(pres, dens, gama) result(eint)
        implicit none
        real, intent(in) :: pres, dens, gama
        real :: eint

        eint = pres/(gama-1.0)/dens

    end function eos_eintIdealGas

end module eos
