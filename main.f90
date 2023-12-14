PROGRAM main

  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 100
  DOUBLE PRECISION, DIMENSION(n,n) :: A
  DOUBLE PRECISION, DIMENSION(n) :: h, b, x, f_ana
  INTEGER :: info, i 
  INTEGER, DIMENSION(n) :: ipiv
  DOUBLE PRECISION :: R, dx, P, sig

  ! Set up parameters for algorithm
  R = 5.0
  dx = (2*R)/(real(n)-1)

  P = -1
  sig = 1

  ! Create x axis and analytical solution
  DO i = 1, n
    x(i) = (REAL(i)-1)*dx - 5
    f_ana(i) = -(P/sig)*(R**2 - x(i)**2)
  END DO

  ! Compute numerical solution
  A = 0
  b = 0

  DO i = 2, n-1
    A(i, i-1) = 1
    A(i, i) = -2
    A(i, i+1) = 1
  END DO

  A(1,1) = -2; A(1,2) = 1
  A(n,n) = -2; A(n,n-1) = 1

  b = (P*dx**2)/sig

  CALL DGESV(n, 1, A, n, ipiv, b, n, info)

  IF(info.ne.0) THEN
    PRINT*, 'Inversion failed'
  ELSE
    PRINT*, 'Inversion successful'
  END IF

  ! Write to text file
  OPEN(9, file='xf.txt', form='formatted')
  23 FORMAT(3 (ES23.12E3))
  DO i = 1, n
    WRITE(9, 23) x(i), f_ana(i), b(i)
  END DO
  CLOSE(9)

END PROGRAM
