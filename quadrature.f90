module quadrature

    implicit none
    double precision, parameter :: PI = acos(-1.0d0)
    
    public

contains

! Arfken p. 541
! Legendre polynomials
double precision pure recursive function Legendre(n,x) result(returnval)
    integer, intent(in) :: n
    double precision, intent(in) :: x
    double precision :: Pn1, Pn2

    if (n == 0) then
        returnval = 1
    else if (n == 1) then
        returnval = x
    else
        Pn2 = Legendre(n-2,x)
        Pn1 = Legendre(n-1,x)
        returnval = 2*x*Pn1 - Pn2 - (x*Pn1-Pn2)/n
    end if
end function

! Arfken p. 541
! Derivatives of Legendre polynomials
double precision pure recursive function DLegendre(n,x) result(returnval)
    integer, intent(in) :: n
    double precision, intent(in) :: x
    double precision :: DPn1, Pn1

    if (n == 0) then
        returnval = 0
    else if (n == 1) then
        returnval = 1
    else
        DPn1 = DLegendre(n-1,x)
        Pn1 = Legendre(n-1, x)
        returnval = n * Pn1 + x * DPn1
    end if
end function

! Second derivative of Legendre polynomials
double precision pure recursive function DDLegendre(n,x) result(returnval)
    integer, intent(in) :: n
    double precision, intent(in) :: x
    double precision :: DPn1, DDPn1

    if (n == 0) then
        returnval = 0
    else if (n == 1) then
        returnval = 0
    else
        DPn1 = DLegendre(n-1,x)
        DDPn1 = DDLegendre(n-1,x)
        returnval = (n + 1.0) * DPn1 + x* DDPn1
    end if
end function

! Gauss-Legendre quadrature points (Roots of Legendre polynomial of order N)
subroutine GL_Points(pts)
    double precision, intent(inout) :: pts(:)
    double precision, parameter :: accuracy = 1.0d-15

    integer :: N, i, j
    double precision :: a(size(pts)), linsp_start, linsp_end, magdiff, tmp

    N = size(pts)

    ! Initial guess of node positions
    ! Newman p. 523
    linsp_start = 3.0
    linsp_end = 4.0*N-1
    a = [(linsp_start + (linsp_end-linsp_start) * (i-1) / (N-1), i=1, N)]
    a = a / (4*N+2)
    pts = cos(PI*a+1/(8*N*N*tan(a)))

    ! Newton's method root finder
    magdiff =  2.0 * accuracy
    do while (magdiff > accuracy)
        a = pts - [( Legendre(N, pts(i)), i=1, N )] / [( DLegendre(N, pts(j)), j=1, N )]
        magdiff = maxval(abs(a - pts))
        pts = a
    end do

    ! Reverse the ordering so the negative ones come first
    do i = 1, N/2
        tmp = pts(N+1-i)
        pts(N+1-i) = pts(i)
        pts(i) = tmp
    end do

end subroutine GL_Points

! Gauss-Legendre-Lobatto quadrature points (Roots of Legendre polynomial derivative of order N-1,
! plus 1 and -1).
subroutine GLL_Points(pts)
    double precision, intent(inout) :: pts(:)
    double precision, parameter :: accuracy = 1.0d-15

    integer :: N, i, j
    double precision :: a(size(pts)-2) , magdiff
    double precision :: interiorPts(size(pts)-2)

    N = size(pts)
    call GL_Points(interiorPts)

    ! Newton's method root finder
    magdiff = accuracy * 2.0
    do while (magdiff > accuracy)
        a = interiorPts - [( DLegendre(N-1, interiorPts(i)), i=1, N-2 )] / [( DDLegendre(N-1, interiorPts(j)), j=1, N-2 )]
        magdiff = maxval(abs(a - interiorPts ))
        interiorPts = a
    end do

    pts(N) = 1.0
    pts(2:N-1) = interiorPts
    pts(1) = -1.0

end subroutine GLL_Points

subroutine GL_Weights(pts, wghts)
    double precision, intent(in) :: pts(:)
    double precision, intent(inout) :: wghts(:)

    double precision :: dleg(size(pts))
    integer :: N , i
    N = size(pts)
    dleg = [( DLegendre(N, pts(i)), i = 1, N )]

    wghts = 2.0 / ((1-pts**2) * (dleg**2))

end subroutine


subroutine GLL_Weights(pts, wghts)
    double precision, intent(in) :: pts(:)
    double precision, intent(inout) :: wghts(:)

    double precision :: leg(size(pts))
    integer :: N , i
    N = size(pts)
    leg = [( Legendre(N-1, pts(i)), i = 1, N )]

    wghts = 2.0 / ( N*(N-1) * leg**2 )

end subroutine

! value of the i-th Lagrange polynomial on the [-1,1] domain at point x
! gx are the quadrature points
double precision pure function Lagrange(gx, i, x)
    double precision, intent(in) :: gx(:), x
    integer, intent(in) :: i

    double precision :: prod
    integer :: j

    prod = 1.0
    do j = 1, size(gx)
        if(j /= i) prod = prod * (x - gx(j))/(gx(i)-gx(j))
    end do
    Lagrange = prod
end function 

! value of the derivative of the i-th Lagrange polynomial on the [-1,1] domain at point x
! gx are the quadrature points
! Calculation takes advantage of log differentiation, f'/f = log(f)'
double precision pure function DLagrange(gx,i,x)
    integer,intent(in) :: i
    double precision,intent(in) :: gx(:), x

    double precision :: L
    integer :: m
    L = 0.0

    do m = 1,size(gx)
        if (m /= i) then
            L = L + 1.0/(x-gx(m))
        end if
    end do

    DLagrange = L * Lagrange(gx,i,x)
end function

end module
