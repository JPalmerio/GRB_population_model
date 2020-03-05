MODULE f90f
    USE stats_f90
    IMPLICIT None

CONTAINS

    SUBROUTINE integrate_1D(x, y, N, z)
        INTEGER :: N
        REAL(8), DIMENSION(0:N) :: x, y 
        REAL(8) :: x1, x2, y1, y2
        INTEGER :: i
        REAL(8) :: z
        
        !f2py intent(in) :: x, y
        !f2py intent(hide), depend(x) :: N = shape(x, 0)
        !f2py intent(out) :: z 

        x1 = x(0)
        y1 = y(0)
        z = 0.d0
        DO i=1, N
            x2 = x(i)
            y2 = y(i)
            z = z + 0.5*(x2-x1)*(y2+y1)
            x1 = x2
            y1 = y2
        END DO
            
    END SUBROUTINE integrate_1D

    SUBROUTINE integrate_2D(x, y, N, M, z)
        INTEGER :: N, M
        REAL(8), DIMENSION(0:N,0:M) :: x, y 
        REAL(8) :: x1, x2, y1, y2, integrated_y
        INTEGER :: i, j
        REAL(8), DIMENSION(0:N) :: z
        
        !f2py intent(in) :: x, y
        !f2py intent(hide), depend(x) :: N = shape(x, 0)
        !f2py intent(hide), depend(x) :: M = shape(x, 1)
        !f2py intent(out) :: z 

        DO i=0, N
            x1 = x(i,0)
            y1 = y(i,0)
            integrated_y = 0.d0
            DO j=1, M
                x2 = x(i,j)
                y2 = y(i,j)
                integrated_y = integrated_y + 0.5*(x2-x1)*(y2+y1)
                x1 = x2
                y1 = y2
            END DO
            z(i) = integrated_y
        END DO
            
    END SUBROUTINE integrate_2D

    SUBROUTINE integrate_2D_w(x, y, w, N, M, z)
        ! 2D integrations with weights on y
        INTEGER :: N, M
        REAL(8), DIMENSION(0:N,0:M) :: x, y, w
        REAL(8) :: x1, x2, y1, y2, w1, w2, integrated_y
        INTEGER :: i, j
        REAL(8), DIMENSION(0:N) :: z
        
        !f2py intent(in) :: x, y, w
        !f2py intent(hide), depend(x) :: N = shape(x, 0)
        !f2py intent(hide), depend(x) :: M = shape(x, 1)
        !f2py intent(out) :: z 

        DO i=0, N
            x1 = x(i,0)
            y1 = y(i,0)
            w1 = w(i,0)
            integrated_y = 0.d0
            DO j=1, M
                x2 = x(i,j)
                y2 = y(i,j)
                w2 = w(i,j)
                integrated_y = integrated_y + 0.5*(x2-x1)*(w2*y2 + w1*y1)
                x1 = x2
                y1 = y2
                w1 = w2
            END DO
            z(i) = integrated_y
        END DO
            
    END SUBROUTINE integrate_2D_w

    SUBROUTINE det_prob_ECLAIRs(cts_ECL, t90obs, Cvar, offax_corr, omega_ECL, n_sigma, threshold, N, N_off, det_prob)
        INTEGER :: N, N_off, GRB
        REAL(8), DIMENSION(0:N_off,0:N_off) :: offax_corr, omega_ECL
        INTEGER, DIMENSION(0:N_off,0:N_off) :: pixel_det_tot, pixel_det_cts, pixel_det_flnc
        REAL(8), DIMENSION(0:N) :: cts_ECL, t90obs, Cvar
        REAL(8), DIMENSION(0:2,0:N) :: det_prob
        REAL(8) :: threshold, n_sigma, omega_ECL_tot

        !f2py intent(in) :: cts_ECL, t90obs, Cvar, offax_corr, omega_ECL, n_sigma, threshold
        !f2py intent(hide), depend(offax_corr) :: N_off = shape(offax_corr, 0)
        !f2py intent(hide), depend(cts_ECL) :: N = shape(cts_ECL, 0)
        !f2py intent(out) :: det_prob


        omega_ECL_tot = SUM(omega_ECL)

        DO GRB=0, N
            WHERE(offax_corr(0:N_off,0:N_off) >= n_sigma * SQRT(threshold)/cts_ECL(GRB) )
              pixel_det_cts = 1
            ELSEWHERE
              pixel_det_cts = 0
            END WHERE
            det_prob(1,GRB) = SUM(pixel_det_cts * omega_ECL) / omega_ECL_tot
            
            WHERE(offax_corr(0:N_off,0:N_off) >= n_sigma * SQRT(threshold/t90obs(GRB)) / (cts_ECL(GRB) * Cvar(GRB)) )
              pixel_det_flnc = 1
            ELSEWHERE
              pixel_det_flnc = 0
            END WHERE
            det_prob(2,GRB) = SUM(pixel_det_flnc * omega_ECL) / omega_ECL_tot
            
            WHERE(pixel_det_cts(0:N_off,0:N_off) == 1)
                pixel_det_tot = 1
            ELSEWHERE(pixel_det_flnc(0:N_off,0:N_off) == 1)
                pixel_det_tot = 1
            ELSEWHERE
                pixel_det_tot = 0
            END WHERE

            det_prob(0,GRB) = SUM(pixel_det_tot * omega_ECL) / omega_ECL_tot
            ! WRITE(*,*) 'GRB number ', GRB,' det_prob: ', det_prob
        END DO

    END SUBROUTINE det_prob_ECLAIRs

    SUBROUTINE calc_ktild(alpha, beta, N, ktild)
        ! Calculates the normalization factor for the spectral shape for a Band function
        ! Uses the Gamma functions for optimized speed
        INTEGER :: i, N
        REAL(8), DIMENSION(0:N) :: alpha, beta, ktild, x_c
        REAL(8) :: xBx, Gamma, gln

        !f2py intent(in) :: alpha, beta
        !f2py intent(hide), depend(alpha) :: N = shape(alpha, 0)
        !f2py intent(out) :: ktild
        
        x_c = (beta - alpha) / (2.d0 - alpha)
        DO i=0, N
            xBx = 0.d0
            CALL GammaSeries(Gamma, 2.d0-alpha(i), beta(i)-alpha(i), gln)
            xBx = Gamma * EXP(gln) /( (2.d0 - alpha(i))**(2.d0-alpha(i)) ) &             ! integral below x_c
                & + x_c(i)**(2.d0-alpha(i)) * EXP(alpha(i)-beta(i)) / (beta(i) - 2.d0)  ! analytic integral above x_c
            ktild(i) = 1.d0/xBx
        END DO

    END SUBROUTINE calc_ktild

END MODULE