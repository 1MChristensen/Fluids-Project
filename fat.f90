module compute
  use iso_fortran_env
  use DynamicalArrays
  implicit none
  contains
  
  subroutine gen_init(s, c, recip, phi, x, z)
    real(real64), dimension(:) :: c
    real(real64) :: s, recip, phi, x, z
    integer :: n, i
    ! A subroutine to generate a few initial points close to the origin using series
    ! expansion of the coordinates in terms of s

    n = size(c)

    recip = 0
    phi = 0

    do i = 1, n
      recip = recip + c(i)*s**(2*(i-1))
      phi = phi + c(i)*s**(2*i-1)/(2*i-1)
    end do
  
    x = s - s**3/6 + (1/120.0 - c(2)/15)*s**5 -(1/5040.0 - c(2)/126 + c(2)**2/126 + c(3)/35)*s**7&
        + (1/362880.0 - c(2)/3240 + c(2)**2/324 + c(3)/270 - c(2)*c(3)/135 - c(4)/63)*s**9 &
        - (1/39916800.0 - c(2)/166320 + c(2)**2/4752 - c(2)**3/1782 + c(3)/6600 -c(2)*c(3)/330 &
        + c(3)**2/550 - c(4)/462 + c(2)*c(4)/231 + c(5)/99)*s**11
        

    z = s**2/2 - (1/24.0 - c(2)/12)*s**4 + (1/720.0 - c(2)/36 + c(3)/30)*s**6 &
        - (1/40320.0 - c(2)/576 + c(2)**2/144 + c(3)/80 - c(4)/56)*s**8 &
        + (1/3628800.0 - c(2)/21600 + c(2)**2/1080 - c(2)**3/1620 + c(3)/1200 - c(2)*c(3)/150 &
        - c(4)/140 + c(5)/90)*s**10 &
        - (1/479001600.0 - c(2)/1451520 + c(2)**2/25920 - c(2)**2/3888 + c(3)/43200 &
        - c(2)*c(3)/1080 + c(3)*c(2)**2/1080 + c(3)**2/600 - c(3)/2016 + c(2)*c(4)/252 &
        + c(5)/216 - c(6)/132)*s**12

  end subroutine

  subroutine find_diff(arr, d, d_store)
    real(real64), dimension(0:) :: arr
    real(real64), dimension(:), allocatable :: d_arr, diff, diff2, b_store, d_store
    real(real64) :: d
    integer :: i, n, n2, j

    ! A subroutine to calcualte and return the successive difference values.

    ! d_store - store of del dphi values
    d_store = 0
    ! d - the next dphi value
    d = 0
    deallocate(d_store)

    ! arr - current dphi values, n - size of arr
    n = size(arr)

    ! d_arr - the doubled array of arr
    allocate(d_arr(-n+1:n-1))

    do i = -n+1, n-1
      d_arr(i) = arr(abs(i))
    end do

    ! n2 - size of doubled arr
    n2 = size(d_arr)
    
    ! diff - difference between successive value of d_arr
    allocate(diff(n2-1))
  
    allocate(b_store(5))
    
    ! calculate successive diffrences and store them in diff
    do i = 1, size(diff)
      diff(i) = d_arr(-n+i+1) - d_arr(-n+i) 
    end do

    ! set first bottom value to last diff value
    b_store(1) = diff(size(diff))
    
    ! start to iterate for each del value
    n2 = n2 - 1
    n = 2
    do j = 1,4
      allocate(diff2(size(diff)-1))

      do i = 1, size(diff2)
        diff2(i) = diff(i+1) - diff(i)
      end do

      deallocate(diff); allocate(diff(size(diff2)))

      diff = diff2

      deallocate(diff2)
      b_store(n) = diff(size(diff))
      n = n + 1
      n2 = n2 - 1
    end do
    
    d = 2*diff(size(diff)) - diff(size(diff)-1)
    n = size(b_store)

    allocate(d_store(size(b_store)))

    d_store(size(d_store)) = d
    do i = 1, size(d_store)-1
      d = d + b_store(size(b_store) - i)
      d_store(size(d_store)-i) = d
    end do
    
    n = ubound(d_arr, dim=1)

    d = d + d_arr(n)
    
  end subroutine

  subroutine succ_diff(arr, nval, d)
    real(real64), dimension(0:) :: arr
    real(real64), dimension(-size(arr)+1:size(arr)-1) :: d_arr
    real(real64), dimension(5) :: lval, d
    real(real64), dimension(:), allocatable :: diff, tmp_diff
    real(real64) :: nval
    integer :: i,j

    ! Another improved subroutine to calculate the successive difference value 
    ! of a given array

    do i = -size(arr)+1, size(arr)-1
      d_arr(i) = arr(abs(i))
    end do

    allocate(diff(size(d_arr)-1))

    do i = 1, size(d_arr)-1
      diff(i) = d_arr(-size(arr)+i+1) - d_arr(-size(arr)+i)
    end do

    lval(1) = diff(size(diff))

    do j = 2,5
      allocate(tmp_diff(size(diff)-1))

      do i = 1, size(tmp_diff)
        tmp_diff(i) = diff(i+1) - diff(i)
      end do
      
      deallocate(diff); allocate(diff(size(tmp_diff)))

      lval(j) = tmp_diff(size(tmp_diff))

      diff = tmp_diff

      deallocate(tmp_diff)
    end do

    d(1) = nval - arr(size(arr)-1)
    do i = 2, 5
      d(i) = d(i-1) - lval(i-1)
    end do
  end subroutine

  subroutine calc_next(phi_store, d, d_store, x_store, z_store, w, phi)
    real(real64), dimension(:), allocatable :: phi_store, d_store, x_store, z_store
    real(real64), dimension(size(x_store)) :: dx, dz 
    real(real64), dimension(5) :: ndx_arr, ndz_arr
    real(real64) :: d, phi, w, ndx, ndz, x, z
    integer :: i

    ! A subroutine to calculate the next coordinate value given its past values and its predicted d value

    phi = phi_store(size(phi_store)) + d - d_store(1)/2 - d_store(2)/12 - d_store(3)/24 - 19*d_store(4)/720 -3*d_store(5)/160

    ndx = w*cos(phi); ndz = w*sin(phi)

    do i = 1, size(x_store)
      dx(i) = w*cos(phi_store(i))
      dz(i) = w*sin(phi_store(i))
    end do

    call succ_diff(dx, ndx, ndx_arr)
  
    x = x_store(size(x_store)) + ndx - ndx_arr(1)/2 - ndx_arr(2)/12 - ndx_arr(3)/24 - 19*ndx_arr(4)/720 - 3*ndx_arr(5)/160

    call succ_diff(dz, ndz, ndz_arr)

    z = z_store(size(z_store)) + ndz - ndz_arr(1)/2 - ndz_arr(2)/12 - ndz_arr(3)/24 - 19*ndz_arr(4)/720 - 3*ndz_arr(5)/160

    call AddToList(phi_store, phi)
    call AddToList(x_store, x)
    call AddToList(z_store, z)
  end subroutine

  subroutine add_correction(phi_store, x_store, z_store, recip_store, w, beta)
    real(real64), dimension(:) :: phi_store, x_store, z_store, recip_store
    real(real64) :: w, eps, beta, eta
    integer :: n
    
    ! A subroutine to calculate the small error and add it on as a correction.

    ! set n for ease 
    n = size(phi_store)

    eps = 2 + beta*z_store(n) - sin(phi_store(n))/x_store(n) - recip_store(n)

    eta = eps/(1 + 2*(w*95.0/288.0)**2 * (sin(phi_store(n))/x_store(n))**2 + (w*95.0*cos(phi_store(n)))/(288.0*x_store(n))&
          - beta*cos(phi_store(n))*(w*95.0/288.0)**2)

    x_store(n) = x_store(n) - (w*95.0/288.0)**2*eta*sin(phi_store(n))

    z_store(n) = z_store(n) + (w*95.0/288.0)**2*eta*cos(phi_store(n))

    phi_store(n) = phi_store(n) + w*eta*95.0/288.0

    recip_store(n) = recip_store(n) + eta
  end subroutine
end module compute

program main
  use compute
  integer, parameter :: n = 5
  integer :: i
  real(real64) :: beta, rhorec, s, phi, x, z, w, d
  real(real64), dimension(6) :: c
  real(real64), dimension(:), allocatable :: d_store, s_store, recip_store, phi_store, x_store, z_store, dphi
  
  ! This is the beta value for the description of the system (Currently only works for negative beta values!)
  beta = -0.7

  ! Store coefficient values
  c(1) = 1
  c(2) = 3.0*beta/8.0
  c(3) = -beta/48.0 + 5.0*beta**2/192.0
  c(4) = 11.0*beta/5760.0 -11.0*beta**2/1920.0 +7*beta**3/9216.0
  c(5) = beta/8960.0 + 629.0*beta**2/645120.0 -487.0*beta**3/645120.0 + beta**4/81920.0
  c(6) = 233.0*beta/14515200.0 + 7.0*beta**2/77414400.0 - 271.0*beta**4/4300800.0 &
        + 11.0*beta**5/88473600.0

  rhorec=0
  phi = 0
  s = 0
  w = 0.1

  allocate(d_store(n), s_store(n), recip_store(n), phi_store(n), x_store(n), z_store(n), dphi(n))

  ! Generate for s=0-0.4 the values of coodinates using the series expansions.
  do i = 0, 4
    s = w*i
    call gen_init(s, c, rhorec, phi, x, z)
    s_store(i+1) = s
    recip_store(i+1) = rhorec
    phi_store(i+1) = phi
    x_store(i+1) = x
    z_store(i+1) = z
    dphi(i+1) = w * recip_store(i+1)
  end do

  ! For larger values of s, iterate using successive differences method
  do while(phi_store(size(phi_store)) > 0) 
    call find_diff(dphi, d, d_store)

    call AddToList(dphi, d)

    rhorec = d/w

    call AddToList(recip_store, rhorec)
    
    call calc_next(phi_store, d, d_store, x_store, z_store, w, phi)

    call add_correction(phi_store, x_store, z_store, recip_store, w, beta)

    dphi(size(dphi)) = w*recip_store(size(recip_store))
  end do

  ! Print out values
  do i = 1, size(x_store)
    print*, w*(i-1),phi_store(i)*57.2958, x_store(i), z_store(i)
  end do

  z_store = z_store - maxval(z_store)


  ! Save values in an unformatted .txt file
  OPEN(9, file='coords.txt', form='formatted')
  23 FORMAT(4 (ES23.12E3))
  DO i = 1, size(x_store)
    WRITE(9, 23) w*(i-1), phi_store(i), x_store(i), z_store(i)
  END DO
  CLOSE(9)

  deallocate(d_store, s_store, recip_store, phi_store, x_store, z_store, dphi)
end program main
