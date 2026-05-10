!********************************************************************************!
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2023 Sarit Dutta <saridut@gmail.com>                             !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!********************************************************************************!

  Module mapping_m

  Implicit none

  Type mapobj
      Character (len=:), allocatable :: nam    ! Name of the map
      Character (len=:), allocatable :: namrev ! Name of the reverse map
      Real (8)  :: a
      Real (8)  :: b
      Real (8)  :: lmap
      Integer   :: m
  End type mapobj

  Private

  Public :: mapobj

  Public :: map_ab2unity,                &
            map_unity2ab,                &
            map_semi_infinite2unity,     &
            map_unity2semi_infinite,     &
            map_unity2itsine,            &
            map_itsine2unity

  Real(8), parameter  :: pi    = 3.141592653589793d0

  Contains

!******************************************************************************

  Elemental real(8) function map_ab2unity ( a, b, z, k)

! map_ab2unity linearly maps z E [a,b] to x E [-1,1]. More precisely, it defines
! the function f(z): [ 2*z-(b+a) ]/(b-a), where z E [a,b], x = f(z) E [-1,1].  k has to
! be an integer >= 0.  If k = 0, the value of the function f(z) at z = z is
! returned.  If k > 0, then the value of the kth derivative of f(z) w.r.t z at z
! = z is returned.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: b
  Real (8), Intent (In)  :: z
  Integer,  Intent (In)  :: k

  If ( k == 0 ) then
     map_ab2unity = (2.d0*z-b-a)/(b-a)
  Elseif ( k == 1) then
     map_ab2unity = 2.d0/(b-a)
  Elseif ( k > 1 ) then
     map_ab2unity = 0.d0
  End if

  End function map_ab2unity

!******************************************************************************

  Elemental real (8) function map_unity2ab ( a, b, x, k )

! map_unity2ab linearly maps x E [-1,1] to z E [a,b]. More precisely, it defines
! the function f(x): (1/2)*[ (b-a)*x + (b+a) ], where x E [-1,1], z = f(x) E [a,b].  k
! has to be an integer >= 0.  If k = 0, the value of the function f(x) at x = x
! is returned.  If k > 0, then the value of the kth derivative of f(x) w.r.t x
! at x = x is returned.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: b
  Real (8), Intent (In)  :: x
  Integer,  Intent (In)  :: k

  If ( k == 0 ) then
      map_unity2ab = 0.5d0 * ( (b-a)*x + b+a )
  Elseif (k == 1 ) then
      map_unity2ab = 0.5d0*(b-a)
  Elseif (k > 1 ) then
      map_unity2ab = 0.d0
  End if

  End function map_unity2ab

!******************************************************************************

  Elemental real (8) function map_semi_infinite2unity ( a, lmap, z, k )

! map_semi_infinte2unity algebraically maps z E [0,Infinity] to x E [-1,1]. More
! precisely, it defines the function f(z): (z-a-lmap)/(z-a+lmap), where z E
! [0,Infinty], x = f(z) E [-1,1], and lmap is the mapping parameter.  k has to
! be an integer = 0, 1, or 2.  If k = 0, the value of the function f(z) at z = z
! is returned.  If k > 0, then the value of the kth derivative of f(z) w.r.t z
! at z = z is returned.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: lmap
  Real (8), Intent (In)  :: z
  Integer,  Intent (In)  :: k

  If ( k == 0 ) then
    map_semi_infinite2unity = (z-a-lmap) / (z-a+lmap)
  Elseif ( k == 1) then
    map_semi_infinite2unity = 2.d0*lmap/( z-a+lmap )**2
  Elseif ( k == 2 ) then
    map_semi_infinite2unity = -4.d0*lmap /( z - a + lmap )**3
  End if

  End function map_semi_infinite2unity

!******************************************************************************

  Elemental real (8) function map_unity2semi_infinite ( a, lmap, x, k )

! map_unity2semi_infinte algebraically maps x E [-1,1] to z E [0,Infinity]. More
! precisely, it defines the function f(x): a + lmap*(1+x)/(1-x), where x E
! [-1,1], z = f(x) E [0,Infinty].  k has to be an integer = 0, 1, or 2.  If k =
! 0, the value of the function f(x) at x = x is returned.  If k > 0, then the
! value of the kth derivative of f(x) w.r.t x at x = x is returned.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: lmap
  Real (8), Intent (In)  :: x
  Integer,  Intent (In)  :: k

  If ( k == 0 ) then
     map_unity2semi_infinite = a + lmap*(1.d0+x) / (1.d0-x)
  Elseif ( k == 1) then
     map_unity2semi_infinite = 2.d0*lmap/( (1.d0-x)*(1.d0-x) )
  Elseif ( k == 2 ) then
     map_unity2semi_infinite = 4.d0*lmap/(1.d0-x)**3
  End if

  End function map_unity2semi_infinite

!******************************************************************************

  Elemental real (8) function map_unity2itsine ( a, b, m, x, k )

! map_unity2itsine maps x E [-1,1] to z E [a,b] using iterated sines. 
! m is a parameter which controls the depth of iterations.
! The map is defined recursively:
! g_0 (x) = x
! g_i (x) = sin [ (pi/2) g_{i-1} ], i >= 1
! The derivates can also be defines through recursion. 
! The first derivative:
! g'_0 (x) = 1
! g'_i (x) = (pi/2) * cos [ (pi/2)*g_{i-1} ] * g'_{i-1}, i >= 1
! The second derivative:
! g''_0 (x) = 0
! g''_i (x) = (pi/2) * [ cos( (pi/2)*g_{i-1} ) * g''_{i-1} 
!                       - (pi/2) * sin( (pi/2)*g_{i-1} ) * (g'{i-1})^2 ]
!
! k has to be an integer = 0, 1, or 2.  If k = 0, the value of the function at x
! = x is returned.  If k > 0, then the value of the kth derivative of the
! function w.r.t x at x = x is returned.

! Reference: 
! Tang, T. and Trummer, M. R., (1996), Boundary layer resolving
! pseudospectral methods for singular perturbation problems, SIAM J. Sci.
! Comput., Vol 17, pp. 430-438.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: b
  Integer,  Intent (In)  :: m
  Real (8), Intent (In)  :: x
  Integer,  Intent (In)  :: k
  Integer                :: i
  Real (8)               :: hpi
  Real (8)               :: gim1
  Real (8)               :: gi
  Real (8)               :: gim1prime
  Real (8)               :: giprime
  Real (8)               :: gim1dblprime
  Real (8)               :: gidblprime

  hpi = pi/2.d0

  If ( k == 0 ) then
     gim1 = x
     Do i = 1, m
         gi   = sin( hpi*gim1 )
         gim1 = gi
     End do
     map_unity2itsine = map_unity2ab( a, b, gi, 0 ) ! <-- OK
  Elseif ( k == 1 ) then
     gim1 = x
     gim1prime = 1.d0
     Do i = 1, m
         gi      = sin( hpi*gim1 )
         giprime = hpi * cos( hpi*gim1 ) * gim1prime
         gim1      = gi
         gim1prime = giprime
     End do
     map_unity2itsine = giprime ! <-- TO BE FIXED
  Elseif ( k == 2 ) then
     gim1 = x
     gim1prime = 1.d0
     gim1dblprime = 0.d0
     Do i = 1, m
         gi         = sin( hpi*gim1 )
         giprime    = hpi * cos( hpi*gim1 ) * gim1prime
         gidblprime = hpi*( gim1dblprime*cos(hpi*gim1) - hpi*sin(hpi*gim1)*gim1prime**2 )
         gim1         = gi
         gim1prime    = giprime
         gim1dblprime = gidblprime
     End do
     map_unity2itsine = gidblprime ! <-- TO BE FIXED
  End if

  End function map_unity2itsine

!******************************************************************************

  Elemental real (8) function map_itsine2unity ( a, b, m, z, k )

! map_itsine2unity maps z E [a,b] to x E [-1,1] using iterated inverse sines. 
! This is the inverse of map_unity2itsine described above.
! m is a parameter which controls the depth of iterations.
! The map is defined recursively:
! h_0 (y) = y
! h_i (y) = sin [ (pi/2) g_{i-1} ], i >= 1
! The derivates can also be defines through recursion. 
! The first derivative:
! h'_0 (y) = 1
! h'_i (y) = (pi/2) * cos [ (pi/2)*g_{i-1} ] * g'_{i-1}, i >= 1
! The second derivative:
! h''_0 (y) = 0
! h''_i (y) = (pi/2) * [ cos( (pi/2)*g_{i-1} ) * g''_{i-1} 
!                       - (pi/2) * sin( (pi/2)*g_{i-1} ) * (g'{i-1})^2 ]
!
! k has to be an integer = 0, 1, or 2.  If k = 0, the value of the function at x
! = x is returned.  If k > 0, then the value of the kth derivative of the
! function w.r.t x at x = x is returned.

! Reference: 
! Tang, T. and Trummer, M. R., (1996), Boundary layer resolving
! pseudospectral methods for singular perturbation problems, SIAM J. Sci.
! Comput., Vol 17, pp. 430-438.

  Real (8), Intent (In)  :: a
  Real (8), Intent (In)  :: b
  Integer,  Intent (In)  :: m
  Real (8), Intent (In)  :: z
  Integer,  Intent (In)  :: k
  Integer                :: i
  Real (8)               :: y
  Real (8)               :: tpi
  Real (8)               :: him1
  Real (8)               :: hi
  Real (8)               :: him1prime
  Real (8)               :: hiprime
  Real (8)               :: him1dblprime
  Real (8)               :: hidblprime

  tpi = 2.d0/pi

! First map to [-1,1] linearly, then apply inverse iterated sine.
  y = map_ab2unity(a, b, z, 0)

  If ( k == 0 ) then
     him1 = y
     Do i = 1, m
         hi   = tpi * asin( him1 )
         him1 = hi
     End do
     map_itsine2unity = hi  ! <-- OK
  Elseif ( k == 1 ) then
     him1 = y
     him1prime = 1.d0
     Do i = 1, m
         hi      = tpi * asin( him1 )
         hiprime = tpi * him1prime / sqrt( 1.d0 - him1**2)
         him1      = hi
         him1prime = hiprime
     End do
     map_itsine2unity = hiprime ! <-- TO BE FIXED
  Elseif ( k == 2 ) then
     him1 = y
     him1prime = 1.d0
     him1dblprime = 0.d0
     Do i = 1, m
         hi         = tpi * asin( him1 )
         hiprime    = tpi * him1prime / sqrt( 1.d0 - him1**2)
         hidblprime = tpi*( him1dblprime + him1*him1prime**2/(1.d0-him1**2) ) &
                           / sqrt(1.d0-him1**2)
         him1         = hi
         him1prime    = hiprime
         him1dblprime = hidblprime
     End do
     map_itsine2unity = hidblprime ! <-- TO BE FIXED
  End if

  End function map_itsine2unity

!******************************************************************************

  End module mapping_m
