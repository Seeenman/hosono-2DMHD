module riemannsolver
    ! contains Riemann solver subroutines for MHD
    ! Includes:
    ! - HLLC for MHD


    use definitions, only: ndim, nPrimVars, nConsVars, &
        velx_var, vely_var, velz_var, magx_var, magy_var, magz_var
    use simulation, only: sim_riemannSolver
    use eigen, only: eigen_valsFromPrim
    use convert, only: convert_cons2flux, convert_cons2prim
    use GP, only: GP_quadratureWeights, GP_nQuadratureMax

    implicit none

    private

    public :: riemannsolver_getFaceFlux

contains

    pure function riemannsolver_getFaceFlux(uL, uR, dir, rr) result(flux)
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

        nq = rr + 1

        flux = 0.0
        do i=1,nq
            flux = flux + GP_quadratureWeights(i, rr) * riemannsolver_getSingleFlux(uL(:,i), uR(:,i), dir)
        end do

    end function riemannsolver_getFaceFlux

    pure function riemannsolver_getSingleFlux(uL, uR, dir) result(flux)
        implicit none
        ! subroutine arguments
        real, intent(in) :: uL(nConsVars), uR(nConsVars)
        integer, intent(in) :: dir
        real :: flux(nConsVars)

        if ((sim_riemannSolver=="hllc") .or. (sim_riemannSolver=="hllcmhd")) then
            flux = riemannsolver_hllcMHD(uL, uR, dir)
        else if (sim_riemannSolver=="roe") then
            write(*,*) "==================================================================="
            write(*,*) "the roe solver is not implemented yet"
            write(*,*) "==================================================================="
            stop
        else
            write(*,*) "==================================================================="
            write(*,*) "Unrecognized choice of riemann solver: ", trim(sim_riemannSolver)
            write(*,*) "Please check that, in your parameter file, the value of"
            write(*,*) "sim_riemannSolver is set to one of hll, roe, hllc"
            write(*,*) "==================================================================="
            stop
        end if

    end function riemannsolver_getSingleFlux

    pure function riemannsolver_fromConsHllcMHD(uL, uR, dir) result(flux)
        ! funnction:    riemannsolver_fromConsHllcMHD
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
        real :: eigL(NUMB_WAVE), eigR(NUMB_WAVE) ! eiegenvalues at left and right states
        real :: eigA(NUMB_WAVE) ! eigenvalues from arithmetic mean of left and right states
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
        integer :: vel_dirN, vel_dirT1, vel_dirT2
        integer :: mom_dirN!, mom_dirT1, mom_dirT2
        integer :: mag_dirN, mag_dirT1, mag_dirT2

        ! compute left and right Riemann state in terms of primitive variables
        vL = convert_cons2prim(uL)
        vR = convert_cons2prim(uR)

        ! set indexing variables based on direction
        if (dir==XDIR) then
            vel_dirN  = velx_var; mag_dirN  = magx_var; mom_dirN  = momx_var
            vel_dirT1 = vely_var; mag_dirT1 = magy_var; mom_dirT1 = momy_var
            vel_dirT2 = velz_var; mag_dirT2 = magz_var; mom_dirT2 = momz_var
        else if (dir==YDIR) then
            vel_dirN  = vely_var; mag_dirN  = magy_var; mom_dirN  = momy_var
            vel_dirT1 = velx_var; mag_dirT1 = magx_var; mom_dirT1 = momx_var
            vel_dirT2 = velz_var; mag_dirT2 = magz_var; mom_dirT2 = momz_var
        end if

        ! get eigenvalues 
        eigL = eigen_valsFromPrim(vL, dir)
        eigR = eigen_valsFromPrim(vR, dir)
        eigA = eigen_valsFromPrim((vL+vR)/2, dir)

        ! find fastest signal velocities. Equation 4.5a-b in Einfeldt et at. 1991
        sL = min(eigL(WAVE_FASTLEFT), eigA(WAVE_FASTLEFT)) 
        sR = max(eigA(WAVE_FASTRGHT), eigR(WAVE_FASTRGHT))

        ! compute F(uL), F(uR)
        fL = convert_cons2flux(vL, dir)
        fR = convert_cons2flux(vR, dir)

        ! compute velocity annd magnetic field HLL states
        uHLL = (sR*uR - sL*uL + fL - fR)/(sR-sL)
        velHLL = uHLL(MOMX_VAR:MOMZ_VAR)/uHLL(DENS_VAR)
        magHLL = uHLL(MAGX_VAR:MAGZ_VAR)

        ! get dot product of magnetic field and velocity field
        ! in the left, right and, HLL states
        magDotVelHLL = DOT_PRODUCT(magHLL, velHLL)
        magDotVelL = DOT_PRODUCT(vL(MAGX_VAR:MAGZ_VAR), vL(VELX_VAR:VELZ_VAR))
        magDotVelR = DOT_PRODUCT(vR(MAGX_VAR:MAGZ_VAR), vR(VELX_VAR:VELZ_VAR))

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
        rhoL = vL(DENS_VAR); rhoR = vR(DENS_VAR)
        pL   = vL(PRES_VAR); pR   = vR(PRES_VAR)

        ! total pressure
        pTotL = pL + 0.5*DOT_PRODUCT(vL(MAGX_VAR:MAGZ_VAR), vL(MAGX_VAR:MAGZ_VAR))
        pTotR = pR + 0.5*DOT_PRODUCT(vR(MAGX_VAR:MAGZ_VAR), vR(MAGX_VAR:MAGZ_VAR))

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
            uStr(DENS_VAR)  = rhoL*dL/(sL-sStr)
            uStr(vel_dirN)  = uStr(DENS_VAR)*sStr
            uStr(vel_dirT1) = uStr(DENS_VAR)*velT1L - (magNHLL*magT1HLL - magNL*magT1L)/(sL-sStr)
            uStr(vel_dirT2) = uStr(DENS_VAR)*velT2L - (magNHLL*magT2HLL - magNL*magT2L)/(sL-sStr)
            uStr(mag_dirN)  = magNHLL
            uStr(mag_dirT1) = magT1HLL
            uStr(mag_dirT2) = magT2HLL
            uStr(ENER_VAR)  = uL(ENER_VAR)*dL + pStr*sStr - pTotL*velNL - magNHLL*magDotVelHLL + magNL*magDotVelL 
            uStr(ENER_VAR)  = uStr(ENER_VAR)/(sL-sStr) 
            ! flux
            flux = fL + sL*(uStr - uL)
        else if ((sStr <= 0.0) .and. (0.0 <= sR)) then
            ! approximate solution in left star region
            uStr(DENS_VAR)  = rhoR*dR/(sR-sStr)
            uStr(vel_dirN)  = uStr(DENS_VAR)*sStr
            uStr(vel_dirT1) = uStr(DENS_VAR)*velT1R - (magNHLL*magT1HLL - magNR*magT1R)/(sR-sStr)
            uStr(vel_dirT2) = uStr(DENS_VAR)*velT2R - (magNHLL*magT2HLL - magNR*magT2R)/(sR-sStr)
            uStr(mag_dirN)  = magNHLL
            uStr(mag_dirT1) = magT1HLL
            uStr(mag_dirT2) = magT2HLL
            uStr(ENER_VAR)  = uR(ENER_VAR)*dR + pStr*sStr - pTotR*velNR - magNHLL*magDotVelHLL + magNR*magDotVelR 
            uStr(ENER_VAR)  = uStr(ENER_VAR)/(sR-sStr)
            ! flux
            flux = fR + sR*(uStr - uR)
        else
            flux = fR
        end if

    end function riemannsolver_fromConsHllcMHD

end module riemannsolver
