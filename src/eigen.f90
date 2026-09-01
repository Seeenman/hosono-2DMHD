module eigen

    ! contains routines that are 
    ! related to the eigenstructure of the
    ! local Riemann problem for the 3D Euler Equations
    ! in the x, y, and z split cases


    use definitions
    use simulation, only: sim_smallB, sim_smallDensity, sim_smallPressure

    implicit none

    private

    public :: eigen_valsFromPrim, eigen_vecsFromPrim, eigen_c_fs

contains

    subroutine eigen_vecsFromPrim(V, dir, primorcons, leigenvecs, reigenvecs)
        implicit none
        real, intent(in) :: V(nPrimVars)
        integer, intent(in) :: dir
        character(len=*), intent(in) :: primorcons
        real, intent(out) :: leigenvecs(nConsVars, nConsVars), reigenvecs(nConsVars, nConsVars)
        ! local variables
        integer :: vel_dirN, vel_dirT1, vel_dirT2, mag_dirN, mag_dirT1, mag_dirT2
        real :: rho, velN, velT1, velT2, B_N, B_T1, B_T2, BB, p, gam, k, &
                c_s, c_f, beta_T1, beta_T2, alpha_f, alpha_s, a, sgnBn
        real, allocatable, dimension(:,:) :: dUdV, dVdU

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
        ! get primitive variables
        rho   = max(V(dens_var), sim_smallDensity)
        velN  = V(vel_dirN)  ! normal velocity component
        velT1 = V(vel_dirT1) ! first tangential velocity component
        velT2 = V(vel_dirT2) ! second tangential velocity component
        B_N   = V(mag_dirN)  ! normal magnetic field component
        B_T1  = V(mag_dirT1) ! first tangential magnetic field component
        B_T2  = V(mag_dirT2) ! second tangential magnetic field component
        p     = max(V(pres_var), sim_smallPressure)
        gam   = V(gama_var)
        BB    = B_N**2 + B_T1**2 + B_T2**2 ! magnitude of magnetic field squared

        ! find wave speeds
        c_f = eigen_c_fs(gam, rho, p, BB, B_N, +1.0)
        c_s = eigen_c_fs(gam, rho, p, BB, B_N, -1.0)
        ! DEBUG
        ! print *, "c_f", c_f
        ! print *, "c_s", c_s
        ! print *, "c_s^2+c_f^2", c_s**2 + c_f**2
        ! print *, "c_s^2-c_f^2", c_s**2 - c_f**2
        ! GUBED
        a = SQRT(gam*p/rho) ! thermal sound speed
        alpha_f = eigen_alpha_fs(c_f, c_s, a, BB, B_N, +1.0)
        alpha_s = eigen_alpha_fs(c_f, c_s, a, BB, B_N, -1.0)

        ! compute normalized tangential magnetic field components
        beta_T1 = eigen_beta_T(B_T1, B_T2)
        beta_T2 = eigen_beta_T(B_T2, B_T1)

        ! get the sign of B_N
        sgnBn = SIGN(1.0, B_N)

        ! DEBUG
        ! print *, "+++inside of eigen+++"
        ! print *, alpha_s**2 + alpha_f**2, "should be 1"
        ! print *, alpha_s**2 * c_s**2 + alpha_f**2 * c_f**2 - a**2, "should be 0"
        ! print *, alpha_s * alpha_f - a*sqrt((B_T1**2 + B_T2**2)/rho)/(c_f**2 - c_s**2), "should be 0"
        ! print *, "^^^inside of eigen^^^"
        ! GUBED

        !! ----------------------
        !! right eigenvectors
        !! ----------------------
        reigenvecs = 0.0
        ! fast left magnetoacoustic wave
        reigenvecs(dens_var,  WAVE_FASTLEFT) = alpha_f*rho
        reigenvecs(vel_dirN,  WAVE_FASTLEFT) = -alpha_f*c_f
        reigenvecs(vel_dirT1, WAVE_FASTLEFT) = alpha_s*c_s*beta_T1*sgnBn
        reigenvecs(vel_dirT2, WAVE_FASTLEFT) = alpha_s*c_s*beta_T2*sgnBn
        reigenvecs(mag_dirT1, WAVE_FASTLEFT) = alpha_s*SQRT(rho)*a*beta_T1
        reigenvecs(mag_dirT2, WAVE_FASTLEFT) = alpha_s*SQRT(rho)*a*beta_T2
        reigenvecs(pres_var,  WAVE_FASTLEFT) = alpha_f*rho*a**2
        ! Alfven wave left
        reigenvecs(vel_dirT1, WAVE_ALFNLEFT) = -beta_T2
        reigenvecs(vel_dirT2, WAVE_ALFNLEFT) = beta_T1
        reigenvecs(mag_dirT1, WAVE_ALFNLEFT) = -SQRT(rho)*beta_T2*sgnBn
        reigenvecs(mag_dirT2, WAVE_ALFNLEFT) = SQRT(rho)*beta_T1*sgnBn
        ! slow left magnetoacoustic wave
        reigenvecs(dens_var,  WAVE_SLOWLEFT) = alpha_s*rho
        reigenvecs(vel_dirN,  WAVE_SLOWLEFT) = -alpha_s*c_s
        reigenvecs(vel_dirT1, WAVE_SLOWLEFT) = -alpha_f*c_f*beta_T1*sgnBn
        reigenvecs(vel_dirT2, WAVE_SLOWLEFT) = -alpha_f*c_f*beta_T2*sgnBn
        reigenvecs(mag_dirT1, WAVE_SLOWLEFT) = -alpha_f*SQRT(rho)*a*beta_T1
        reigenvecs(mag_dirT2, WAVE_SLOWLEFT) = -alpha_f*SQRT(rho)*a*beta_T2
        reigenvecs(pres_var,  WAVE_SLOWLEFT) = alpha_s*rho*a**2
        ! middle entropy wave
        reigenvecs(dens_var,  WAVE_ENTROPY)  = 1.0
        ! slow right magnetoacoustic wave
        reigenvecs(dens_var,  WAVE_SLOWRGHT) = alpha_s*rho
        reigenvecs(vel_dirN,  WAVE_SLOWRGHT) = alpha_s*c_s
        reigenvecs(vel_dirT1, WAVE_SLOWRGHT) = alpha_f*c_f*beta_T1*sgnBn
        reigenvecs(vel_dirT2, WAVE_SLOWRGHT) = alpha_f*c_f*beta_T2*sgnBn
        reigenvecs(mag_dirT1, WAVE_SLOWRGHT) = -alpha_f*SQRT(rho)*a*beta_T1
        reigenvecs(mag_dirT2, WAVE_SLOWRGHT) = -alpha_f*SQRT(rho)*a*beta_T2
        reigenvecs(pres_var,  WAVE_SLOWRGHT) = alpha_s*rho*a**2
        ! Alfven wave right
        reigenvecs(vel_dirT1, WAVE_ALFNRGHT) = beta_T2
        reigenvecs(vel_dirT2, WAVE_ALFNRGHT) = -beta_T1
        reigenvecs(mag_dirT1, WAVE_ALFNRGHT) = -SQRT(rho)*beta_T2*sgnBn
        reigenvecs(mag_dirT2, WAVE_ALFNRGHT) = SQRT(rho)*beta_T1*sgnBn
        ! fast right magnetoacoustic wave
        reigenvecs(dens_var,  WAVE_FASTRGHT) = alpha_f*rho
        reigenvecs(vel_dirN,  WAVE_FASTRGHT) = alpha_f*c_f
        reigenvecs(vel_dirT1, WAVE_FASTRGHT) = -alpha_s*c_s*beta_T1*sgnBn
        reigenvecs(vel_dirT2, WAVE_FASTRGHT) = -alpha_s*c_s*beta_T2*sgnBn
        reigenvecs(mag_dirT1, WAVE_FASTRGHT) = alpha_s*SQRT(rho)*a*beta_T1
        reigenvecs(mag_dirT2, WAVE_FASTRGHT) = alpha_s*SQRT(rho)*a*beta_T2
        reigenvecs(pres_var,  WAVE_FASTRGHT) = alpha_f*rho*a**2

        !! ----------------------
        !! left eigenvectors
        !! ----------------------
        leigenvecs = 0.0
        ! note that we store the left eigenvectors in columns not rows of eigenvecs matrix
        ! fast left magnetoacoustic wave
        leigenvecs(vel_dirN,  WAVE_FASTLEFT) = 1/(2*a**2)*(-alpha_f*c_f)
        leigenvecs(vel_dirT1, WAVE_FASTLEFT) = 1/(2*a**2)*(alpha_s*c_s*beta_T1*sgnBn)
        leigenvecs(vel_dirT2, WAVE_FASTLEFT) = 1/(2*a**2)*(alpha_s*c_s*beta_T2*sgnBn)
        leigenvecs(mag_dirT1, WAVE_FASTLEFT) = 1/(2*a**2)*(alpha_s*a*beta_T1/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_FASTLEFT) = 1/(2*a**2)*(alpha_s*a*beta_T2/SQRT(rho))
        leigenvecs(pres_var,  WAVE_FASTLEFT) = 1/(2*a**2)*(alpha_f/rho)
        ! Alfven wave left
        leigenvecs(vel_dirT1, WAVE_ALFNLEFT) = 0.5*(-beta_T2)
        leigenvecs(vel_dirT2, WAVE_ALFNLEFT) = 0.5*(beta_T1)
        leigenvecs(mag_dirT1, WAVE_ALFNLEFT) = 0.5*(-beta_T2*sgnBn/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_ALFNLEFT) = 0.5*(beta_T1*sgnBn/SQRT(rho))
        ! slow left magnetoacoustic wave
        leigenvecs(vel_dirN,  WAVE_SLOWLEFT) = 1/(2*a**2)*(-alpha_s*c_s)
        leigenvecs(vel_dirT1, WAVE_SLOWLEFT) = 1/(2*a**2)*(-alpha_f*c_f*beta_T1*sgnBn)
        leigenvecs(vel_dirT2, WAVE_SLOWLEFT) = 1/(2*a**2)*(-alpha_f*c_f*beta_T2*sgnBn)
        leigenvecs(mag_dirT1, WAVE_SLOWLEFT) = 1/(2*a**2)*(-alpha_f*a*beta_T1/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_SLOWLEFT) = 1/(2*a**2)*(-alpha_f*a*beta_T2/SQRT(rho))
        leigenvecs(pres_var,  WAVE_SLOWLEFT) = 1/(2*a**2)*(alpha_s/rho)
        ! middle entropy wave
        leigenvecs(dens_var,  WAVE_ENTROPY)  = 1.0
        leigenvecs(pres_var,  WAVE_ENTROPY)  = -1/a**2
        ! slow right magnetoacoustic wave
        leigenvecs(vel_dirN,  WAVE_SLOWRGHT) = 1/(2*a**2)*(alpha_s*c_s)
        leigenvecs(vel_dirT1, WAVE_SLOWRGHT) = 1/(2*a**2)*(alpha_f*c_f*beta_T1*sgnBn)
        leigenvecs(vel_dirT2, WAVE_SLOWRGHT) = 1/(2*a**2)*(alpha_f*c_f*beta_T2*sgnBn)
        leigenvecs(mag_dirT1, WAVE_SLOWRGHT) = 1/(2*a**2)*(-alpha_f*a*beta_T1/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_SLOWRGHT) = 1/(2*a**2)*(-alpha_f*a*beta_T2/SQRT(rho))
        leigenvecs(pres_var,  WAVE_SLOWRGHT) = 1/(2*a**2)*(alpha_s/rho)
        ! Alfven wave right
        leigenvecs(vel_dirT1, WAVE_ALFNRGHT) = 0.5*(beta_T2)
        leigenvecs(vel_dirT2, WAVE_ALFNRGHT) = 0.5*(-beta_T1)
        leigenvecs(mag_dirT1, WAVE_ALFNRGHT) = 0.5*(-beta_T2*sgnBn/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_ALFNRGHT) = 0.5*(beta_T1*sgnBn/SQRT(rho))
        ! fast right magnetoacoustic wave
        leigenvecs(vel_dirN,  WAVE_FASTRGHT) = 1/(2*a**2)*(alpha_f*c_f)
        leigenvecs(vel_dirT1, WAVE_FASTRGHT) = 1/(2*a**2)*(-alpha_s*c_s*beta_T1*sgnBn)
        leigenvecs(vel_dirT2, WAVE_FASTRGHT) = 1/(2*a**2)*(-alpha_s*c_s*beta_T2*sgnBn)
        leigenvecs(mag_dirT1, WAVE_FASTRGHT) = 1/(2*a**2)*(alpha_s*a*beta_T1/SQRT(rho))
        leigenvecs(mag_dirT2, WAVE_FASTRGHT) = 1/(2*a**2)*(alpha_s*a*beta_T2/SQRT(rho))
        leigenvecs(pres_var,  WAVE_FASTRGHT) = 1/(2*a**2)*(alpha_f/rho)

        if (primorcons=='cons') then
            allocate(dUdV(nConsVars, nConsVars))
            allocate(dVdU(nConsVars, nConsVars))

            ! set k
            k = 1.0-gam

            dUdV = 0.0
            ! first col
            dUdV(dens_var, dens_var) = 1.0
            dUdV(velx_var, dens_var) = V(velx_var)
            dUdV(vely_var, dens_var) = V(vely_var)
            dUdV(velz_var, dens_var) = V(velz_var)
            dUdV(pres_var, dens_var) = 0.5*(velN**2 + velT1**2 + velT2**2)
            ! second col
            dUdV(velx_var, velx_var) = rho
            dUdV(pres_var, velx_var) = rho*V(velx_var)
            ! third col
            dUdV(vely_var, vely_var) = rho
            dUdV(pres_var, vely_var) = rho*V(vely_var)
            ! fourth col
            dUdV(velz_var, velz_var) = rho
            dUdV(pres_var, velz_var) = rho*V(velz_var)
            ! fifth col
            dUdV(magx_var, magx_var) = 1.0
            dUdV(pres_var, magx_var) = V(magx_var)
            ! sixth col
            dUdV(magy_var, magy_var) = 1.0
            dUdV(pres_var, magy_var) = V(magy_var)
            ! seventh col
            dUdV(magz_var, magz_var) = 1.0
            dUdV(pres_var, magz_var) = V(magz_var)
            ! eighth col
            dUdV(pres_var, pres_var) = -1/k

            dVdU = 0.0
            ! first col
            dVdU(dens_var, dens_var) = 1.0
            dVdU(velx_var, dens_var) = -V(velx_var)/rho
            dVdU(vely_var, dens_var) = -V(vely_var)/rho
            dVdU(velz_var, dens_var) = -V(velz_var)/rho
            dVdU(pres_var, dens_var) = -k*dUdV(pres_var, dens_var)
            ! second col
            dVdU(velx_var, velx_var) = 1/rho
            dVdU(pres_var, velx_var) = k*V(velx_var)
            ! third col
            dVdU(vely_var, vely_var) = 1/rho
            dVdU(pres_var, vely_var) = k*V(vely_var)
            ! fourth col
            dVdU(velz_var, velz_var) = 1/rho
            dVdU(pres_var, velz_var) = k*V(velz_var)
            ! fifth col
            dVdU(magx_var, magx_var) = 1.0
            dVdU(pres_var, magx_var) = k*V(magx_var)
            ! sixth col
            dVdU(magy_var, magy_var) = 1.0
            dVdU(pres_var, magy_var) = k*V(magy_var)
            ! seventh col
            dVdU(magz_var, magz_var) = 1.0
            dVdU(pres_var, magz_var) = k*V(magz_var)
            ! eighth col
            dVdU(pres_var, pres_var) = -k

            ! convert primitive eigenvectors to conservative
            reigenvecs = MATMUL(dUdV, reigenvecs)
            leigenvecs = MATMUL(TRANSPOSE(dVdU), leigenvecs)

            deallocate(dUdV)
            deallocate(dVdU)
        end if

    end subroutine eigen_vecsFromPrim

    function eigen_valsFromPrim(V, dir) result(eigenvals)
        ! function:     eigen_valsFromPrim
        ! Author:       Sean Riedel
        ! purpose:      Given the primitive variables at a cell
        !               computes the eigenvalues (wave speeds) for
        !               the split system given by dir
        ! 
        ! Inputs:       - V (real) the primitive variable values at a given cell
        !               - dir (character) the direction (x or y or z) of the split system
        !                 for which we are finding the eigenvalues
        !               
        ! Outputs:      - eigenvals (real) the eigenvalues
        ! ------------------------------------------------------------
        implicit none
        real, intent(in) :: V(numb_var)
        integer, intent(in) :: dir
        real :: eigenvals(NUMB_WAVE)
        ! local variables
        integer :: dirN
        real :: rho, velN, B_N, BB, p, gam, &
                c_s, c_f, c_a

        if (dir==XDIR) then
            dirN  = 0
        else if (dir==YDIR) then
            dirN  = 1
        else if (dir==ZDIR) then
            dirN  = 2
        end if
        rho   = max(V(dens_var), sim_smallDensity)
        velN  = V(velx_var+dirN)
        B_N   = V(magx_var+dirN)
        p     = max(V(pres_var), sim_smallPressure)
        gam   = V(gama_var)
        BB    = V(magx_var)**2 + V(magy_var)**2 + V(magz_var)**2

        c_s = eigen_c_fs(gam, rho, p, BB, B_N, -1.0)
        c_f = eigen_c_fs(gam, rho, p, BB, B_N, +1.0)
        c_a = abs(B_N)/SQRT(rho) ! Alfven velocity TODO: why does B_N have absolute value around it?

        eigenvals(WAVE_FASTLEFT) = velN - c_f
        eigenvals(WAVE_ALFNLEFT) = velN - c_a
        eigenvals(WAVE_SLOWLEFT) = velN - c_s
        eigenvals(WAVE_ENTROPY ) = velN
        eigenvals(WAVE_SLOWRGHT) = velN + c_s
        eigenvals(WAVE_ALFNRGHT) = velN + c_a
        eigenvals(WAVE_FASTRGHT) = velN + c_f

    end function eigen_valsFromPrim

    function eigen_c_fs(gam, rho, p, BB, B_N, pm) result(c_fs)
        ! compute fast (slow) magnetoacoustic speed
        ! pm = +1.0 -> c_f
        ! pm = -1.0 -> c_s
        ! TODO (is the above the correct terminology?)
        implicit none
        real, intent(in) :: gam, rho, p, BB, B_N, pm
        real :: c_fs
        ! local variables

        ! compute c_f (or c_s)
        ! equations 2.18 and 2.19 in Ryu and Jones 1995
        c_fs = (gam*p + BB)/rho
        c_fs = 0.5*(c_fs + pm * SQRT(c_fs**2 - 4*gam*p*B_N**2/rho**2))
        c_fs = SQRT(max(1.0e-14,c_fs))
        
    end function eigen_c_fs

    function eigen_alpha_fs(c_f, c_s, a, BB, B_N, pm) result(alpha)
        ! pm = +1 -> alpha_f
        ! pm = -1 -> alpha_s
        ! a is the thermal sound speed
        implicit none
        real, intent(in) :: c_f, c_s, a, BB, B_N, pm
        real :: alpha
        ! local variables
        real :: BB_T, cf2_minus_cs2

        ! Sum of squares of tangential components of magnetic field
        BB_T = BB - B_N**2 

        ! c_f^2 - c_s^2
        cf2_minus_cs2 = c_f**2 - c_s**2

        ! compute alpha
        if (cf2_minus_cs2 > 1.0e-14) then ! from appendix A3. of Athena: new astrophysical MHD code
            alpha = 0.5*(cf2_minus_cs2) + pm*(a**2 - 0.5*(c_f**2 + c_s**2))
            alpha = alpha/(cf2_minus_cs2)
        else if (pm > 0) then ! alpha_f
            alpha = 1.0
        else if (pm < 0) then ! alpha_s
            alpha = 0.0
        end if
        
        alpha = min(1.0, max(0.0, alpha)) ! square root guard that I don't understand
        alpha = sqrt(alpha)

    end function eigen_alpha_fs

    function eigen_beta_T(B_T1, B_T2) result(beta)
        implicit none
        real, intent(in) :: B_T1, B_T2
        real :: beta
        ! local variables
        real :: BB_T

        ! Sum of squares of tangential components of magnetic field
        BB_T = B_T1**2 + B_T2**2

        ! compute beta
        if (BB_T < sim_smallB**2) then
            beta = SIGN(1.0, B_T1)/SQRT(2.0)
        else
            beta = B_T1/SQRT(max(BB_T,1.0e-30))
        end if
    end function eigen_beta_T

end module eigen
