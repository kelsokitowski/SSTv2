module ETDRKcoefficients
  use array_dimensions
  implicit none

  contains

  ! Subroutine to compute ETDRK4 coefficients
  ! Based on the contour integral formulation using M points on the unit circle
  subroutine getETDRKcoefficients(linearOperator, kVals, h, kLength, &
                                   EX1, EX2, Q, f1, f2, f3)
    implicit none

    ! Input parameters
    integer, intent(in) :: kLength
    double precision, intent(in) :: linearOperator(kLength)
    double precision, intent(in) :: kVals(kLength)
    double precision, intent(in) :: h

    ! Output parameters (all 1D arrays of length kLength)
    double precision, intent(out) :: EX1(kLength)
    double precision, intent(out) :: EX2(kLength)
    double precision, intent(out) :: Q(kLength)
    double precision, intent(out) :: f1(kLength)
    double precision, intent(out) :: f2(kLength)
    double precision, intent(out) :: f3(kLength)

    ! Local variables
    integer, parameter :: M = 16  ! number of points for periodic trapezoidal rule
    double precision, parameter :: pi = acos(-1.0d0)
    complex(8) :: r(M)  ! M points on unit circle in complex plane
    complex(8) :: LR(M, kLength)  ! h*L + r for each k and each r
    complex(8) :: hL(kLength)  ! h*L for each k
    complex(8) :: temp_sum
    integer :: kj, mj

    ! Compute h*L for each wavenumber
    do kj = 1, kLength
      hL(kj) = cmplx(h * linearOperator(kj), 0.0d0, kind=8)
    end do

    ! Compute exponentials (surprisingly only two are needed)
    ! EX1 = exp(h*L) and EX2 = exp(h*L/2)
    do kj = 1, kLength
      EX1(kj) = exp(h * linearOperator(kj))
      EX2(kj) = exp(h * linearOperator(kj) / 2.0d0)
    end do

    ! Generate M points on the unit circle in the complex plane (roots of unity)
    ! r = exp(1i*pi*((1:M)-0.5)/M)
    do mj = 1, M
      r(mj) = exp(cmplx(0.0d0, pi * (dble(mj) - 0.5d0) / dble(M), kind=8))
    end do

    ! Generate h*L + r for each k and each r value
    ! Used in periodic trapezoidal rule for contour integration
    do kj = 1, kLength
      do mj = 1, M
        LR(mj, kj) = hL(kj) + r(mj)
      end do
    end do

    ! ETDRK4 coefficients using contour integration
    ! Q = h*real(mean((exp(LR/2)-1)/LR, 1))
    do kj = 1, kLength
      temp_sum = cmplx(0.0d0, 0.0d0, kind=8)
      do mj = 1, M
        temp_sum = temp_sum + (exp(LR(mj, kj) / 2.0d0) - 1.0d0) / LR(mj, kj)
      end do
      Q(kj) = h * real(temp_sum / dble(M))
    end do

    ! f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3, 1))
    do kj = 1, kLength
      temp_sum = cmplx(0.0d0, 0.0d0, kind=8)
      do mj = 1, M
        temp_sum = temp_sum + (-4.0d0 - LR(mj, kj) + &
                   exp(LR(mj, kj)) * (4.0d0 - 3.0d0 * LR(mj, kj) + &
                   LR(mj, kj)**2.0d0)) / LR(mj, kj)**3.0d0
      end do
      f1(kj) = h * real(temp_sum / dble(M))
    end do

    ! f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3, 1))
    do kj = 1, kLength
      temp_sum = cmplx(0.0d0, 0.0d0, kind=8)
      do mj = 1, M
        temp_sum = temp_sum + (2.0d0 + LR(mj, kj) + &
                   exp(LR(mj, kj)) * (-2.0d0 + LR(mj, kj))) / LR(mj, kj)**3.0d0
      end do
      f2(kj) = h * real(temp_sum / dble(M))
    end do

    ! f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3, 1))
    do kj = 1, kLength
      temp_sum = cmplx(0.0d0, 0.0d0, kind=8)
      do mj = 1, M
        temp_sum = temp_sum + (-4.0d0 - 3.0d0 * LR(mj, kj) - LR(mj, kj)**2.0d0 + &
                   exp(LR(mj, kj)) * (4.0d0 - LR(mj, kj))) / LR(mj, kj)**3.0d0
      end do
      f3(kj) = h * real(temp_sum / dble(M))
    end do

  end subroutine getETDRKcoefficients

end module ETDRKcoefficients
