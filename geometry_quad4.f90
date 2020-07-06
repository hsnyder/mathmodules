module geometry

    implicit none
    public


contains

    subroutine ComputeMapCoefs(vertices, coeffs)
        double precision, intent(in) :: vertices(4,2)
        double precision, intent(out) :: coeffs(4,2)

        integer :: ipiv(4), info

        ! Matrix for interpolation from the reference domain [-1,1]x[-1,1] to the deformed element.
        ! A good reference: https://www.particleincell.com/2012/quad-interpolation/
        double precision ::  M(4,4) 
        
        M (1,:) = [1.0d0, -1.0d0, -1.0d0,  1.0d0]
        M (2,:) = [1.0d0, -1.0d0,  1.0d0, -1.0d0]
        M (3,:) = [1.0d0,  1.0d0,  1.0d0,  1.0d0]
        M (4,:) = [1.0d0,  1.0d0, -1.0d0, -1.0d0]
        
        coeffs = vertices
        
        call dgesv(4, 2, M, 4, ipiv, coeffs, 4, info)

        if(info < 0) then
            print *, info, "-th argument had illegal value"
            stop
        else if (info > 0) then
            print *, info, "singular"
            stop
        end if

    end subroutine

    pure function MapRealToReference(co, p) result(res)
        double precision, intent(in) :: co(4,2), p(2)
        double precision:: res(2)

        double precision :: aa, bb, cc, det
        double precision :: a(4), b(4), x, y
        a = co(:,1)
        b = co(:,2)
        x = p(1)
        y = p(2)

        aa = a(4)*b(3) - a(3)*b(4)
        bb = a(4)*b(1) -a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + x*b(4) - y*a(4)
        cc = a(2)*b(1) -a(1)*b(2) + x*b(2) - y*a(2)
        det = sqrt(bb*bb - 4*aa*cc)
        res(2) = (-bb+det)/(2*aa)
        res(1) = (x-a(1)-a(3)*res(2))/(a(2)+a(4)*res(2))
        
    end function MapRealToReference

    pure function dxdr (co, p) result(res)
        double precision, intent(in) :: co(4,2) , p(2) 
        double precision  :: res, r, s 
        r = p(1)
        s = p(2)
        res = (co(2,1) + co(4,1)*s)
    end function

    pure function dyds (co, p) result(res)
        double precision, intent(in) :: co(4,2) , p(2) 
        double precision  :: res, r, s 
        r = p(1)
        s = p(2)
        res = (co(3,2) + co(4,2)*r)
    end function

    pure function dxds (co, p) result(res)
        double precision, intent(in) :: co(4,2) , p(2) 
        double precision  :: res, r, s 
        r = p(1)
        s = p(2)
        res = (co(3,1) + co(4,1)*r)
    end function

    pure function dydr (co, p) result(res)
        double precision, intent(in) :: co(4,2) , p(2) 
        double precision  :: res, r, s 
        r = p(1)
        s = p(2)
        res = (co(2,2) + co(4,2)*s)
    end function


        

    ! Jacobian determinant
    ! p=(r,s) are the x,y vars in the reference domain [-1,1]x[-1,1]
    ! co are the coefficients from ComputeMapCoefs
    pure function Jac(co, p)
        double precision, intent(in) :: co(4,2) , p(2)
        double precision :: Jac
        Jac = dxdr(co,p) * dyds(co,p)  -  dydr(co,p) * dxds(co,p)
    end function

    pure function IsLeft(p0,p1,p2) ! returns > 0 if p2 is left of the line through p0,p1, <0 if it's to the right, 0 if it's on the line
        double precision, intent(in) :: p0(2), p1(2), p2(2)
        double precision :: IsLeft
        IsLeft = (p1(1) - p0(1)) * (p2(2) - p0(2)) - (p2(1) - p0(1)) * (p1(2) - p0(2)) 
    end function

    pure function PtInQuad(pt, quad, V) ! quad is a slice from EtoV, i.e. it's integers containing indices into V
        logical :: PtInQuad
        double precision , intent(in) :: V(:,:) , pt(2)
        integer , intent(in) :: quad(:)   

        integer :: i, wn ! iterator, winding number
        integer, parameter :: n = 4

        ! See http://geomalgorithms.com/a03-_inclusion.html

        wn = 0
        do i = 1, n
            if ( V(2,quad(i)) <= pt(2) ) then   
                if ( V(2, quad(mod(i,n)+1)) > pt(2) ) then
                    if (IsLeft( V(:,quad(i)) , V(:,quad(mod(i,n)+1)), pt) > 0) wn = wn + 1
                end if
            else 
                if ( V(2, quad(mod(i,n)+1)) <= pt(2) ) then
                    if (IsLeft( V(:,quad(i)) , V(:,quad(mod(i,n)+1)), pt) < 0) wn = wn - 1
                end if
            end if
        end do

        PtInQuad = .true.
        if(wn == 0) PtInQuad = .false.
        
        
    end function PtInQuad


    pure function GetMatingElementCoord(face, Np, idx) 
        integer, intent(in) :: face, idx, Np
        integer :: i, GetMatingElementCoord(2)

        ! the face on the corresponding quad is oriented the other way. 
        i = Np+1-idx

        if (face == 1) then 
            GetMatingElementCoord = [1,i]
        else if (face == 2) then 
            GetMatingElementCoord = [i,Np]
        else if (face == 3) then 
            GetMatingElementCoord = [Np,i]
        else if (face == 4) then 
            GetMatingElementCoord = [i,1]
        end if
    end function GetMatingElementCoord 

    subroutine ComputeNormals(Np, co, JacDet, gx, normalpx, normalpy, normalmx, normalmy)

        integer, intent(in) :: Np
        double precision, intent(in) :: JacDet(Np,Np)
        double precision, intent(in)  :: co(:,:), gx(Np)
        double precision,intent(out) :: normalpx(2,Np), normalpy(2,Np), normalmx(2,Np), normalmy(2,Np)

        double precision :: G(2,2), Gtmp(2,2), Gdet

        integer :: i,j

        ! Compute Y positive normals
        do i = 1,Np

            j = Np  
            G(1,:) = [dxdr(co(:,:), [gx(i),gx(j)]), dxds(co(:,:), [gx(i),gx(j)])]
            G(2,:) = [dydr(co(:,:), [gx(i),gx(j)]), dyds(co(:,:), [gx(i),gx(j)])]
            
            Gdet = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            Gtmp = -G
            Gtmp(1,1) = G(2,2)
            Gtmp(2,2) = G(1,1)
            Gtmp = Gtmp/Gdet
            
            G = transpose(Gtmp) 

            normalpy(:,i) = JacDet(i,j) * matmul(G,[ 0, 1])
        end do

        ! Compute Y negative normals
        do i = 1,Np

            j = 1   
            G(1,:) = [dxdr(co(:,:), [gx(i),gx(j)]), dxds(co(:,:), [gx(i),gx(j)])]
            G(2,:) = [dydr(co(:,:), [gx(i),gx(j)]), dyds(co(:,:), [gx(i),gx(j)])]
            
            Gdet = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            Gtmp = -G
            Gtmp(1,1) = G(2,2)
            Gtmp(2,2) = G(1,1)
            Gtmp = Gtmp/Gdet
            
            G = transpose(Gtmp) 

            normalmy(:,i) = JacDet(i,j) * matmul(G,[0, -1])
        end do

        ! Compute X positive normals
        do j = 1,Np

            i = Np  
            G(1,:) = [dxdr(co(:,:), [gx(i),gx(j)]), dxds(co(:,:), [gx(i),gx(j)])]
            G(2,:) = [dydr(co(:,:), [gx(i),gx(j)]), dyds(co(:,:), [gx(i),gx(j)])]
            
            Gdet = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            Gtmp = -G
            Gtmp(1,1) = G(2,2)
            Gtmp(2,2) = G(1,1)
            Gtmp = Gtmp/Gdet
            
            G = transpose(Gtmp) 

            normalpx(:,j) = JacDet(i,j) * matmul(G,[ 1, 0])
        end do

        ! Compute X negative normals
        do j = 1,Np

            i = 1   
            G(1,:) = [dxdr(co(:,:), [gx(i),gx(j)]), dxds(co(:,:), [gx(i),gx(j)])]
            G(2,:) = [dydr(co(:,:), [gx(i),gx(j)]), dyds(co(:,:), [gx(i),gx(j)])]
            
            Gdet = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            Gtmp = -G
            Gtmp(1,1) = G(2,2)
            Gtmp(2,2) = G(1,1)
            Gtmp = Gtmp/Gdet
            
            G = transpose(Gtmp) 

            normalmx(:,j) = JacDet(i,j) * matmul(G,[ -1, 0])
        end do

    end subroutine  


    subroutine CalcJacDet(Np,co,gx,JacDet)
        integer, intent(in) :: Np
        double precision, intent(in)  :: co(:,:), gx(Np)
        double precision, intent(out) :: JacDet(Np,Np)
        integer :: i, j
        ! Construct Jacobian (dx/dr dy/ds - dx/ds dy/dr) for reference domain components r,s
        ! This isn't an actual jacobian, it's a matrix, to be multiplied elementwise into M
        ! where J_km is the value at each node in the local element for the Jacobian determinant. 
        do j = 1, Np
            JacDet(:,j) = [(  Jac(co(:,:), [gx(i), gx(j)])  , i=1, Np  )]
        end do
    end subroutine 

    

end module
