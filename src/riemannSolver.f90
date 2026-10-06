module riemannSolver
    ! contains Riemann solver subroutines for MHD
    ! Includes:
    ! - HLLC for MHD


    use definitions, only: nPrimVars, nConsVars, xdir, ydir, &
        dens_var, momx_var, momy_var, momz_var, ener_var, &
        velx_var, vely_var, velz_var, pres_var, magx_var, magy_var, magz_var
    use simulation, only: sim_riemannSolver
    use eigen, only: eigen_valsFromPrim
    use convert, only: convert_cons2flux, convert_cons2prim
    use GP, only: GP_nQuadrature, GP_quadratureWeights, GP_nQuadratureMax

    implicit none

    private

    public :: riemannSolver_getFaceFlux

contains

    ! pure function riemannSolver_getFaceFlux(uL, uR, dir, rr) result(flux)
    function riemannSolver_getFaceFlux(uL, uR, dir, rr) result(flux)
        ! get high order face averaged flux from weighted sum
        ! of pointwise fluxes at quadrature points.
        ! if dir==xdir then returns F, if dir==ydir then returns G
        implicit none
        ! subroutine arguments
        integer, intent(in) :: rr ! GP radius used at this face
        real, dimension(nConsVars, GP_nQuadratureMax), intent(in) :: uL, uR
        integer, intent(in) :: dir
        real :: flux(nConsVars)
        ! local variables
        integer :: i
        integer :: nq ! number of quadrature points

        if (rr==0) then ! FOG
            flux=riemannSolver_getSingleFlux(uL(:,1), uR(:,1), dir)
        else ! high order GP reconstruction
            nq = GP_nQuadrature(rr)
            flux = 0.0
            do i=1,nq
                flux = flux + GP_quadratureWeights(i, rr) * riemannSolver_getSingleFlux(uL(:,i), uR(:,i), dir)
            end do
        end if

    end function riemannSolver_getFaceFlux

    ! pure function riemannSolver_getSingleFlux(uL, uR, dir) result(flux)
    function riemannSolver_getSingleFlux(uL, uR, dir) result(flux)
        implicit none
        ! subroutine arguments
        real, intent(in) :: uL(nConsVars), uR(nConsVars)
        integer, intent(in) :: dir
        real :: flux(nConsVars)

        ! sim_riemannSolver is checked in simulation_init, so the else
        ! branch should be unreachable. error stop (unlike stop and write)
        ! is allowed inside of a pure procedure.
        if ((sim_riemannSolver=="hllc") .or. (sim_riemannSolver=="hllcmhd")) then
            flux = riemannSolver_fromConsHllcMHD(uL, uR, dir)
        else
            error stop "riemannSolver_getSingleFlux: unrecognized sim_riemannSolver "//trim(sim_riemannSolver)
        end if

    end function riemannSolver_getSingleFlux

    ! pure function riemannSolver_fromConsHllcMHD(uL, uR, dir) result(flux)
    function riemannSolver_fromConsHllcMHD(uL, uR, dir) result(flux)
        ! funnction:    riemannSolver_fromConsHllcMHD
        ! Author:       Sean Riedel
        ! purpose:      Given left and right Riemann states 
        !               uL and uR (in terms of the conservative variables), 
        !               computes the HLLC numerical flux function for MHD.
        ! 
        ! Inputs:       - uL (real) the left Riemann state vector of conservative variables
        !               - uR (real) the right Riemann state vector of conservative variables
        !               - dir (integer) the direction (x or y or z) for which we are finding the flux
        !               
        ! Outputs:      - flux (real) the numerical flux function 
        ! ------------------------------------------------------------
        ! Shengtai Li,
        ! Algorithm from
        ! An HLLC Riemann solver for magneto-hydrodynamics,
        ! Journal of Computational Physics,
        ! Volume 203, Issue 1,
        ! 2005,
        ! Pages 344-357,
        ! ISSN 0021-9991,
        ! https://doi.org/10.1016/j.jcp.2004.08.020.
        ! (https://www.sciencedirect.com/science/article/pii/S0021999104003857)
        
        implicit none
        ! subroutine arguments
        real, intent(in) :: uL(nConsVars), uR(nConsVars)
        integer, intent(in) :: dir
        real :: flux(nConsVars)
        ! local variables
        ! eigenvalues of the 7 MHD waves, ordered from the fast left (1)
        ! to the fast right (7) magnetoacoustic wave. See eigen_valsFromPrim
        real :: eigL(7), eigR(7) ! eiegenvalues at left and right states
        real :: eigA(7) ! eigenvalues from arithmetic mean of left and right states
        real :: sL, sR ! fastest signal velocities in left and right direction
        real :: sStr, pStr
        real, dimension(nPrimVars) :: vL, vR
        real, dimension(nConsVars) :: uHLL, uStr
        real, dimension(nConsVars) :: fL, fR
        real, dimension(3) :: magHLL, velHLL
        real :: magDotVelHLL, magDotVelL, magDotVelR
        real :: velNL, velNR, velT1L, velT1R, velT2L, velT2R
        real :: magNL, magNR, magT1L, magT1R, magT2L, magT2R
        real :: magNHLL, magT1HLL, magT2HLL
        real :: pL, pR, pTotL, pTotR, rhoL, rhoR, dL, dR
        real :: momNL, momNR
        integer :: vel_dirN, vel_dirT1, vel_dirT2
        integer :: mom_dirN, mom_dirT1, mom_dirT2
        integer :: mag_dirN, mag_dirT1, mag_dirT2

        ! compute left and right Riemann state in terms of primitive variables
        vL = convert_cons2prim(uL)
        vR = convert_cons2prim(uR)

        ! set indexing variables based on direction
        if (dir==xdir) then
            vel_dirN  = velx_var; mag_dirN  = magx_var; mom_dirN  = momx_var
            vel_dirT1 = vely_var; mag_dirT1 = magy_var; mom_dirT1 = momy_var
            vel_dirT2 = velz_var; mag_dirT2 = magz_var; mom_dirT2 = momz_var
        else if (dir==ydir) then
            vel_dirN  = vely_var; mag_dirN  = magy_var; mom_dirN  = momy_var
            vel_dirT1 = velx_var; mag_dirT1 = magx_var; mom_dirT1 = momx_var
            vel_dirT2 = velz_var; mag_dirT2 = magz_var; mom_dirT2 = momz_var
        end if

        ! get eigenvalues 
        eigL = eigen_valsFromPrim(vL, dir)
        eigR = eigen_valsFromPrim(vR, dir)
        eigA = eigen_valsFromPrim((vL+vR)/2, dir)

        ! find fastest signal velocities. Equation 4.5a-b in Einfeldt et at. 1991
        sL = min(eigL(1), eigA(1)) 
        sR = max(eigA(7), eigR(7))

        ! compute F(uL), F(uR)
        fL = convert_cons2flux(uL, dir)
        fR = convert_cons2flux(uR, dir)

        ! compute velocity annd magnetic field HLL states
        uHLL = (sR*uR - sL*uL + fL - fR)/(sR-sL)
        velHLL = uHLL(momx_var:momz_var)/uHLL(dens_var)
        magHLL = uHLL(magx_var:magz_var)

        ! get dot product of magnetic field and velocity field
        ! in the left, right and, HLL states
        magDotVelHLL = DOT_PRODUCT(magHLL, velHLL)
        magDotVelL = DOT_PRODUCT(vL(magx_var:magz_var), vL(velx_var:velz_var))
        magDotVelR = DOT_PRODUCT(vR(magx_var:magz_var), vR(velx_var:velz_var))

        ! get normal and tangential components of magnetic field and velocity and momentum
        ! in left, right, and HLL states
        velNL  = vL(vel_dirN);  velNR  = vR(vel_dirN)
        velT1L = vL(vel_dirT1); velT1R = vR(vel_dirT1)
        velT2L = vL(vel_dirT2); velT2R = vR(vel_dirT2)
        momNL  = uL(mom_dirN);  momNR  = uR(mom_dirN)
        ! momT1L = uL(mom_dirT1); momT1R = uR(mom_dirT1)
        ! momT2L = uL(mom_dirT2); momT2R = uR(mom_dirT2)
        magNL  = vL(mag_dirN) ; magNR  = vR(mag_dirN) ; magNHLL  = uHLL(mag_dirN)
        magT1L = vL(mag_dirT1); magT1R = vR(mag_dirT1); magT1HLL = uHLL(mag_dirT1)
        magT2L = vL(mag_dirT2); magT2R = vR(mag_dirT2); magT2HLL = uHLL(mag_dirT2)

        ! get other primitive variables at left and right states
        rhoL = vL(dens_var); rhoR = vR(dens_var)
        pL   = vL(pres_var); pR   = vR(pres_var)

        ! total pressure
        pTotL = pL + 0.5*DOT_PRODUCT(vL(magx_var:magz_var), vL(magx_var:magz_var))
        pTotR = pR + 0.5*DOT_PRODUCT(vR(magx_var:magz_var), vR(magx_var:magz_var))

        ! compute S_K - vel_{N,K} for K=L,R
        dL = sL - velNL; dR = sR - velNR

        ! find wavespeed estimate in star region
        sStr = (momNR*dR - momNL*dL + pTotL - pTotR - magNL**2 + magNR**2)/&
               (rhoR*dR - rhoL*dL)

        ! pressure estimate in star region. Note that this can be computed using either
        ! K=L or K=R as the result is the same either way.
        pStr = rhoL*dL*(sStr-velNL) + pTotL - magNL**2 + magNHLL**2 

        ! final choice of HLLC flux
        if (0.0 < sL) then
            flux = fL
        else if ((sL <= 0.0) .and. (0.0 < sStr)) then
            ! approximate solution in left star region
            uStr(dens_var)  = rhoL*dL/(sL-sStr)
            uStr(mom_dirN)  = uStr(dens_var)*sStr
            uStr(mom_dirT1) = uStr(dens_var)*velT1L - (magNHLL*magT1HLL - magNL*magT1L)/(sL-sStr)
            uStr(mom_dirT2) = uStr(dens_var)*velT2L - (magNHLL*magT2HLL - magNL*magT2L)/(sL-sStr)
            uStr(mag_dirN)  = magNHLL
            uStr(mag_dirT1) = magT1HLL
            uStr(mag_dirT2) = magT2HLL
            uStr(ener_var)  = uL(ener_var)*dL + pStr*sStr - pTotL*velNL - magNHLL*magDotVelHLL + magNL*magDotVelL 
            uStr(ener_var)  = uStr(ener_var)/(sL-sStr) 
            ! flux
            flux = fL + sL*(uStr - uL)
        else if ((sStr <= 0.0) .and. (0.0 <= sR)) then
            ! approximate solution in left star region
            uStr(dens_var)  = rhoR*dR/(sR-sStr)
            uStr(mom_dirN)  = uStr(dens_var)*sStr
            uStr(mom_dirT1) = uStr(dens_var)*velT1R - (magNHLL*magT1HLL - magNR*magT1R)/(sR-sStr)
            uStr(mom_dirT2) = uStr(dens_var)*velT2R - (magNHLL*magT2HLL - magNR*magT2R)/(sR-sStr)
            uStr(mag_dirN)  = magNHLL
            uStr(mag_dirT1) = magT1HLL
            uStr(mag_dirT2) = magT2HLL
            uStr(ener_var)  = uR(ener_var)*dR + pStr*sStr - pTotR*velNR - magNHLL*magDotVelHLL + magNR*magDotVelR 
            uStr(ener_var)  = uStr(ener_var)/(sR-sStr)
            ! flux
            flux = fR + sR*(uStr - uR)
        else
            flux = fR
        end if

    end function riemannSolver_fromConsHllcMHD

end module riemannSolver
