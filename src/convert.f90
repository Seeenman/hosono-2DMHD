module convert

    ! contains functions to convert 
    !   primitive variables -> flux
    !   primitive variables -> conservative variables
    !   conservative variables -> primitive variables
    ! 

    use definitions
    use eos, only: eos_presIdealGas
    use simulation, only: sim_gamma, sim_smallEnergy, sim_smallDensity, sim_forceHydro

    implicit none

contains

    pure function convert_prim2flux(V, dir) result(Flux)
        implicit none
        real, intent(in) :: V(nPrimVars)
        integer, intent(in) :: dir
        real :: Flux(nConsVars)
        ! local variables
        real :: E, Bp ! energy and magnetic pressure
        real :: Ptot ! total pressure
        real :: rhoEint ! internal energy
        integer :: vel_dirN, vel_dirT1, vel_dirT2, mag_dirN, mag_dirT1, mag_dirT2

        ! magnetic pressure (same as magnetic energy density)
        Bp = 0.0
        if (.not. sim_forceHydro) then
            Bp = 0.5*(V(magx_var)**2 + V(magy_var)**2 + V(magz_var)**2)
        end if

        ! we need to recompute the internal energy here and can't always trust what is stored 
        ! on V(eint_var). This is because when this function is called on the left and right
        ! Riemann states inside of riemannsolver.F90, sometimes only the values
        ! for of DENS, VELX, VELY, VELZ, MAGX, MAGY, MAGZ, PRES are correct for
        ! the Riemann State. This is true, for example, when using PLM.
        rhoEint = V(pres_var)/(sim_gamma-1.0) 

        ! compute energy
        E = rhoEint + V(dens_var)*(0.5*(V(velx_var)**2 + V(vely_var)**2 + V(velz_var)**2)) + Bp

        ! compute total pressure
        Ptot = V(pres_var) + Bp

        ! set indexing variables based on direction
        if (dir==xdir) then
            vel_dirN  = velx_var
            vel_dirT1 = vely_var
            vel_dirT2 = velz_var
            mag_dirN  = magx_var
            mag_dirT1 = magy_var
            mag_dirT2 = magz_var
        else if (dir==ydir) then
            vel_dirN  = vely_var
            vel_dirT1 = velx_var
            vel_dirT2 = velz_var
            mag_dirN  = magy_var
            mag_dirT1 = magx_var
            mag_dirT2 = magz_var
        end if

        Flux = 0.0
        Flux(dens_var)  = V(dens_var)*V(vel_dirN)
        Flux(vel_dirN)  = V(dens_var)*V(vel_dirN)**2 + Ptot
        Flux(vel_dirT1) = V(dens_var)*V(vel_dirN)*V(vel_dirT1)
        Flux(vel_dirT2) = V(dens_var)*V(vel_dirN)*V(vel_dirT2)
        if (.not. sim_forceHydro) then
            Flux(vel_dirN)  = Flux(vel_dirN) - V(mag_dirN)**2
            Flux(vel_dirT1) = Flux(vel_dirT1) - V(mag_dirT1)*V(mag_dirN)
            Flux(vel_dirT2) = Flux(vel_dirT2) - V(mag_dirT2)*V(mag_dirN)
            Flux(mag_dirT1) = V(vel_dirN)*V(mag_dirT1) - V(vel_dirT1)*V(mag_dirN)
            Flux(mag_dirT2) = V(vel_dirN)*V(mag_dirT2) - V(vel_dirT2)*V(mag_dirN)
        end if
        Flux(ener_var)  = (E + Ptot)*V(vel_dirN) &
            - V(mag_dirN)*(V(velx_var)*V(magx_var) + V(vely_var)*V(magy_var) + V(velz_var)*V(magz_var))
                          

    end function convert_prim2flux

    pure function convert_prim2cons(V) result(U)
        implicit none
        real, intent(in) :: V(nPrimVars)
        real :: U(nConsVars)
        ! local variables
        real :: E, Bp ! energy and magnetic pressure
        real :: rhoEint ! internal energy

        ! magnetic pressure (same as magnetic energy density)
        if (sim_forceHydro) then
            Bp = 0.0
        else
            Bp = 0.5*(V(magx_var)**2 + V(magy_var)**2 + V(magz_var)**2)
        end if

        ! we need to recompute the internal energy here and can't always trust what is stored 
        ! on V(eint_var). This is because when this function is called on the left and right
        ! Riemann states inside of riemannsolver.F90, sometimes only the values
        ! for of DENS, VELX, VELY, VELZ, MAGX, MAGY, MAGZ, PRES are correct for
        ! the Riemann State. This is true, for example, when using PLM.
        rhoEint = V(pres_var)/(sim_gamma-1.0) 

        ! compute energy
        E = rhoEint + V(dens_var)*(0.5*(V(velx_var)**2 + V(vely_var)**2 + V(velz_var)**2)) + Bp

        ! convert to conservative variables
        U(dens_var) = V(dens_var)
        U(momx_var) = V(dens_var)*V(velx_var)
        U(momy_var) = V(dens_var)*V(vely_var)
        U(momz_var) = V(dens_var)*V(velz_var)
        if (sim_forceHydro) then
            U(magx_var:magz_var) = 0.0
        else
            U(magx_var:magz_var) = V(magx_var:magz_var)
        end if
        U(ener_var) = E

    end function convert_prim2cons

    function convert_cons2prim(U) result(V)
        implicit none
        real, intent(in) :: U(nConsVars)
        real :: V(nPrimVars)
        ! local variables
        real :: density, velx, vely, velz, internalEnergy, pressure
        real :: Bp ! magnetic pressure (magnetic energy density)

        ! magnetic pressure (same as magnetic energy density)
        if (sim_forceHydro) then
            Bp = 0.0
        else
            Bp = 0.5*(U(magx_var)**2 + U(magy_var)**2 + U(magz_var)**2) 
        end if
        

        density = U(dens_var)
        if (density < sim_smallDensity) then
            print*, "-----------------------------------------------------------------------"
            print*, "DEBUG"
            print*, "Negative density in convert_cons2prim"
            print*, "density value of ", density
            print*, "density set to", sim_smallDensity
            print*, "GUBED"
            print*, "-----------------------------------------------------------------------"
            density = sim_smallDensity ! preserve positivity in a crude manner
        end if
        velx = U(momx_var)/density
        vely = U(momy_var)/density
        velz = U(momz_var)/density
        internalEnergy = (U(ener_var) - Bp)/density - 0.5*(velx**2 + vely**2 + velz**2) 
        if (internalEnergy < sim_smallEnergy) then
            print*, "-----------------------------------------------------------------------"
            print*, "DEBUG"
            print*, "Negative internal energy in convert_cons2prim"
            print*, "internal energy value of ", internalEnergy
            print*, "internal energy set to", sim_smallEnergy
            print*, "GUBED"
            print*, "-----------------------------------------------------------------------"
            internalEnergy = sim_smallEnergy ! preserve positivity in a crude manner
        end if
        pressure = eos_presIdealGas(density, internalEnergy, sim_gamma)

        V(dens_var) = density
        V(velx_var) = velx
        V(vely_var) = vely
        V(velz_var) = velz
        if (sim_forceHydro) then
            V(magx_var:magz_var) = 0.0
        else
            V(magx_var:magz_var) = U(magx_var:magz_var)
        end if
        V(pres_var) = pressure
        V(eint_var) = internalEnergy

    end function convert_cons2prim

end module convert
