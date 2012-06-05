! -*- f90 -*-

module bspline
contains

function FindSpan(n,p,uu,U) result (span)
  implicit none
  integer(kind=4), intent(in) :: n, p
  real   (kind=8), intent(in) :: uu, U(0:n+p+1)
  integer(kind=4)             :: span
  integer(kind=4) low, high
  if (uu >= U(n+1)) then
     span = n
     return
  end if
  if (uu <= U(p)) then
     span = p
     return
  end if
  low  = p
  high = n+1
  span = (low + high) / 2
  do while (uu < U(span) .or. uu >= U(span+1))
     if (uu < U(span)) then
        high = span
     else
        low  = span
     end if
     span = (low + high) / 2
  end do
end function FindSpan

function FindMult(i,uu,p,U) result (mult)
  implicit none
  integer(kind=4), intent(in)  :: i, p
  real   (kind=8), intent(in)  :: uu, U(0:i+p+1)
  integer(kind=4)              :: mult
  integer(kind=4) :: j
  mult = 0
  do j = -p, p+1
     if (uu == U(i+j)) mult = mult + 1
  end do
end function FindMult

subroutine FindSpanMult(n,p,uu,U,k,s)
  implicit none
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: uu, U(0:n+p+1)
  integer(kind=4), intent(out) :: k, s
  k = FindSpan(n,p,uu,U)
  s = FindMult(k,uu,p,U)
end subroutine FindSpanMult

subroutine BasisFuns(i,uu,p,U,N)
  implicit none
  integer(kind=4), intent(in) :: i, p
  real   (kind=8), intent(in) :: uu, U(0:i+p)
  real   (kind=8), intent(out):: N(0:p)
  integer(kind=4) :: j, r
  real   (kind=8) :: left(p), right(p), saved, temp
  N(0) = 1.0
  do j = 1, p
     left(j)  = uu - U(i+1-j)
     right(j) = U(i+j) - uu
     saved = 0.0
     do r = 0, j-1
        temp = N(r) / (right(r+1) + left(j-r))
        N(r) = saved + right(r+1) * temp
        saved = left(j-r) * temp
     end do
     N(j) = saved
  end do
end subroutine BasisFuns

subroutine DersBasisFuns(i,uu,p,n,U,ders)
  implicit none
  integer(kind=4), intent(in) :: i, p, n
  real   (kind=8), intent(in) :: uu, U(0:i+p)
  real   (kind=8), intent(out):: ders(0:p,0:n)
  integer(kind=4) :: j, k, r, s1, s2, rk, pk, j1, j2
  real   (kind=8) :: saved, temp, d
  real   (kind=8) :: left(p), right(p)
  real   (kind=8) :: ndu(0:p,0:p), a(0:1,0:p)
  ndu(0,0) = 1.0
  do j = 1, p
     left(j)  = uu - U(i+1-j)
     right(j) = U(i+j) - uu
     saved = 0.0
     do r = 0, j-1
        ndu(j,r) = right(r+1) + left(j-r)
        temp = ndu(r,j-1) / ndu(j,r)
        ndu(r,j) = saved + right(r+1) * temp
        saved = left(j-r) * temp
     end do
     ndu(j,j) = saved
  end do
  ders(:,0) = ndu(:,p)
  do r = 0, p
     s1 = 0; s2 = 1;
     a(0,0) = 1.0
     do k = 1, n
        d = 0.0
        rk = r-k; pk = p-k;
        if (r >= k) then
           a(s2,0) = a(s1,0) / ndu(pk+1,rk)
           d =  a(s2,0) * ndu(rk,pk)
        end if
        if (rk > -1) then
           j1 = 1
        else
           j1 = -rk
        end if
        if (r-1 <= pk) then
           j2 = k-1
        else
           j2 = p-r
        end if
        do j = j1, j2
           a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk+j)
           d =  d + a(s2,j) * ndu(rk+j,pk)
        end do
        if (r <= pk) then
           a(s2,k) = - a(s1,k-1) / ndu(pk+1,r)
           d =  d + a(s2,k) * ndu(r,pk)
        end if
        ders(r,k) = d
        j = s1; s1 = s2; s2 = j;
     end do
  end do
  r = p
  do k = 1, n
     ders(:,k) = ders(:,k) * r
     r = r * (p-k)
  end do
end subroutine DersBasisFuns

subroutine CurvePoint(d,n,p,U,Pw,uu,C)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  real   (kind=8), intent(out) :: C(d)
  integer(kind=4) :: j, span
  real   (kind=8) :: basis(0:p)
  span = FindSpan(n,p,uu,U)
  call BasisFuns(span,uu,p,U,basis)
  C = 0.0
  do j = 0, p
     C = C + basis(j)*Pw(:,span-p+j)
  end do
end subroutine CurvePoint

subroutine SurfacePoint(d,n,p,U,m,q,V,Pw,uu,vv,S)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  integer(kind=4), intent(in)  :: m, q
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: V(0:m+q+1)
  real   (kind=8), intent(in)  :: Pw(d,0:m,0:n)
  real   (kind=8), intent(in)  :: uu, vv
  real   (kind=8), intent(out) :: S(d)
  integer(kind=4) :: uj, vj, uspan, vspan
  real   (kind=8) :: ubasis(0:p), vbasis(0:q)
  uspan = FindSpan(n,p,uu,U)
  call BasisFuns(uspan,uu,p,U,ubasis)
  vspan = FindSpan(m,q,vv,V)
  call BasisFuns(vspan,vv,q,V,vbasis)
  S = 0.0
  do uj = 0, p
     do vj = 0, q
        S = S + ubasis(uj)*vbasis(vj)*Pw(:,vspan-q+vj,uspan-p+uj)
     end do
  end do
end subroutine SurfacePoint

subroutine CurvePntByCornerCut(d,n,p,U,Pw,x,Cw)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: x
  real   (kind=8), intent(out) :: Cw(d)
  integer(kind=4) :: i, j, k, s, r
  real   (kind=8) :: uu, alpha, Rw(d,0:p)
  if (x <= U(p)) then
     uu = U(p)
     k = p
     s = FindMult(p,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,0)
        return
     end if
  elseif (x >= U(n+1)) then
     uu = U(n+1)
     k = n+1
     s = FindMult(n,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,n)
        return
     end if
  else
     uu = x
     k = FindSpan(n,p,uu,U)
     s = FindMult(k,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,k-p)
        return
     end if
  end if
  r = p-s
  do i = 0, r
     Rw(:,i) = Pw(:,k-p+i)
  end do
  do j = 1, r
     do i = 0, r-j
        alpha = (uu-U(k-p+j+i))/(U(i+k+1)-U(k-p+j+i))
        Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
     end do
  end do
  Cw(:) = Rw(:,0)
end subroutine CurvePntByCornerCut

subroutine InsertKnot(d,n,p,U,Pw,uu,k,s,r,V,Qw)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: k, s, r
  real   (kind=8), intent(out) :: V(0:n+p+1+r)
  real   (kind=8), intent(out) :: Qw(d,0:n+r)
  integer(kind=4) :: i, j, idx
  real   (kind=8) :: alpha, Rw(d,0:p)
  ! Load new knot vector
  forall (i = 0:k) V(i) = U(i)
  forall (i = 1:r) V(k+i) = uu
  forall (i = k+1:n+p+1) V(i+r) = U(i)
  ! Save unaltered control points
  forall (i = 0:k-p) Qw(:,i)   = Pw(:,i)
  forall (i = k-s:n) Qw(:,i+r) = Pw(:,i)
  forall (i = 0:p-s) Rw(:,i)   = Pw(:,k-p+i)
  ! Insert the knot r times
  do j = 1, r
     idx = k-p+j
     do i = 0, p-j-s
        alpha = (uu-U(idx+i))/(U(i+k+1)-U(idx+i))
        Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
     end do
     Qw(:,idx) = Rw(:,0)
     Qw(:,k+r-j-s) = Rw(:,p-j-s)
  end do
  ! Load remaining control points
  idx = k-p+r
  do i = idx+1, k-s-1
     Qw(:,i) = Rw(:,i-idx)
  end do
end subroutine InsertKnot

subroutine RemoveKnot(d,n,p,U,Pw,uu,r,s,num,t,TOL)
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  real   (kind=8), intent(in)    :: uu
  integer(kind=4), intent(in)    :: r, s, num
  integer(kind=4), intent(out)   :: t
  real   (kind=8), intent(in)    :: TOL

  integer(kind=4) :: m,ord,fout,last,first,off
  integer(kind=4) :: i,j,ii,jj,k
  logical         :: remflag
  real   (kind=8) :: temp(d,0:2*p)
  real   (kind=8) :: alfi,alfj

  m = n + p + 1
  ord = p + 1
  fout = (2*r-s-p)/2
  first = r - p
  last  = r - s
  do t = 0,num-1
     off = first - 1
     temp(:,0) = Pw(:,off)
     temp(:,last+1-off) = Pw(:,last+1)
     i = first; ii = 1
     j = last;  jj = last - off
     remflag = .false.
     do while (j-i > t)
        alfi = (uu-U(i))/(U(i+ord+t)-U(i))
        alfj = (uu-U(j-t))/(U(j+ord)-U(j-t))
        temp(:,ii) = (Pw(:,i)-(1.0-alfi)*temp(:,ii-1))/alfi
        temp(:,jj) = (Pw(:,j)-alfj*temp(:,jj+1))/(1.0-alfj)
        i = i + 1; ii = ii + 1
        j = j - 1; jj = jj - 1
     end do
     if (j-i < t) then
        if (Distance(d,temp(:,ii-1),temp(:,jj+1)) <= TOL) then
           remflag = .true.
        end if
     else
        alfi = (uu-U(i))/(U(i+ord+t)-U(i))
        if (Distance(d,Pw(:,i),alfi*temp(:,ii+t+1)+(1-alfi)*temp(:,ii-1)) <= TOL) then
           remflag = .true.
        end if
     end if
     if (remflag .eqv. .false.) then
        exit ! break out of the for loop
     else
        i = first
        j = last
        do while (j-i > t)
           Pw(:,i) = temp(:,i-off)
           Pw(:,j) = temp(:,j-off)
           i = i + 1
           j = j - 1
        end do
     end if
     first = first - 1
     last  = last  + 1
  end do
  if (t == 0) return
  do k = r+1,m
     U(k-t) = U(k)
  end do
  j = fout
  i = j
  do k = 1,t-1
     if (mod(k,2) == 1) then
        i = i + 1
     else
        j = j - 1
     end if
  end do
  do k = i+1,n
     Pw(:,j) = Pw(:,k)
     j = j + 1
  enddo
contains
  function Distance(d,P1,P2) result (dist)
    implicit none
    integer(kind=4), intent(in) :: d
    real   (kind=8), intent(in) :: P1(d),P2(d)
    integer(kind=4) :: i
    real   (kind=8) :: dist
    dist = 0.0
    do i = 1,d
       dist = dist + (P1(i)-P2(i))*(P1(i)-P2(i))
    end do
    dist = sqrt(dist)
  end function Distance
end subroutine RemoveKnot

subroutine ClampKnot(d,n,p,U,Pw)
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  integer(kind=4) :: k, s
  ! Clamp at left end
  k = p
  s = FindMult(p,U(p),p,U)
  call KntIns(d,n,p,U,Pw,k,s)
  U(0:p-1) = U(p)
  ! Clamp at right end
  k = n+1
  s = FindMult(n,U(n+1),p,U)
  call KntIns(d,n,p,U,Pw,k,s)
  U(n+2:n+p+1) = U(n+1)
contains
  subroutine KntIns(d,n,p,U,Pw,k,s)
      implicit none
      integer(kind=4), intent(in)    :: d
      integer(kind=4), intent(in)    :: n, p
      real   (kind=8), intent(in)    :: U(0:n+p+1)
      real   (kind=8), intent(inout) :: Pw(d,0:n)
      integer(kind=4), intent(in)    :: k, s
      integer(kind=4) :: r, i, j, idx
      real   (kind=8) :: uu, alpha, Rw(d,0:p), Qw(d,0:2*p)
      if (s >= p) return
      uu = U(k)
      r = p-s
      Qw(:,0) = Pw(:,k-p)
      Rw(:,0:p-s) = Pw(:,k-p:k-s)
      do j = 1, r
         idx = k-p+j
         do i = 0, p-j-s
            alpha = (uu-U(idx+i))/(U(i+k+1)-U(idx+i))
            Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
         end do
         Qw(:,j) = Rw(:,0)
         Qw(:,p-j-s+r) = Rw(:,p-j-s)
      end do
      if (k == p) then ! left end
         Pw(:,0:r-1) = Qw(:,r:r+r-1)
      else             ! right end
         Pw(:,n-r+1:n) = Qw(:,p-r:p-1)
      end if
    end subroutine KntIns
end subroutine ClampKnot

subroutine UnclampKnot(d,n,p,U,Pw)
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  integer(kind=4) :: i, j, k
  real   (kind=8) :: alpha
  do i = 0, p-2 ! Unclamp at left end
     U(p-i-1) = U(p-i) - (U(n-i+1)-U(n-i))
     k = p-1
     do j = i, 0, -1
        alpha = (U(p)-U(k))/(U(p+j+1)-U(k))
        Pw(:,j) = (Pw(:,j)-alpha*Pw(:,j+1))/(1-alpha)
        k = k-1
     end do
  end do
  U(0) = U(1) - (U(n-p+2)-U(n-p+1)) ! Set first knot
  do i = 0, p-2  ! Unclamp at right end
     U(n+i+2) = U(n+i+1) + (U(p+i+1)-U(p+i))
     do j = i, 0, -1
        alpha = (U(n+1)-U(n-j))/(U(n-j+i+2)-U(n-j))
        Pw(:,n-j) = (Pw(:,n-j)-(1-alpha)*Pw(:,n-j-1))/alpha
     end do
  end do
  U(n+p+1) = U(n+p) + (U(2*p)-U(2*p-1)) ! Set last knot
end subroutine UnclampKnot

subroutine RefineKnotVector(d,n,p,U,Pw,r,X,Ubar,Qw)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Ubar(0:n+r+1+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n+r+1)
  integer(kind=4) :: m, a, b
  integer(kind=4) :: i, j, k, l
  integer(kind=4) :: idx
  real   (kind=8) :: alpha
  if (r < 0) then
     Ubar = U
     Qw = Pw
     return
  end if
  m = n + p + 1
  a = FindSpan(n,p,X(0),U)
  b = FindSpan(n,p,X(r),U)
  b = b + 1
  forall (j = 0:a-p) Qw(:,j)     = Pw(:,j)
  forall (j = b-1:n) Qw(:,j+r+1) = Pw(:,j)
  forall (j =   0:a) Ubar(j)     = U(j)
  forall (j = b+p:m) Ubar(j+r+1) = U(j)
  i = b + p - 1
  k = b + p + r
  do j = r, 0, -1
     do while (X(j) <= U(i) .and. i > a)
        Qw(:,k-p-1) = Pw(:,i-p-1)
        Ubar(k) = U(i)
        k = k - 1
        i = i - 1
     end do
     Qw(:,k-p-1) = Qw(:,k-p)
     do l = 1, p
        idx = k - p + l
        alpha = Ubar(k+l) - X(j)
        if (abs(alpha) == 0.0) then
           Qw(:,idx-1) = Qw(:,idx)
        else
           alpha = alpha / (Ubar(k+l) - U(i-p+l))
           Qw(:,idx-1) = alpha*Qw(:,idx-1) + (1-alpha)*Qw(:,idx)
        end if
     end do
     Ubar(k) = X(j)
     k = k-1
  end do
end subroutine RefineKnotVector

subroutine DegreeElevate(d,n,p,U,Pw,t,nh,Uh,Qw)
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: t
  integer(kind=4), intent(in)  :: nh
  real   (kind=8), intent(out) :: Uh(0:nh+p+t+1)
  real   (kind=8), intent(out) :: Qw(d,0:nh)

  integer(kind=4) :: i, j, k, kj, tr, a, b
  integer(kind=4) :: m, ph, kind, cind, first, last
  integer(kind=4) :: r, oldr, s, mul, lbz, rbz

  real   (kind=8) :: bezalfs(0:p+t,0:p)
  real   (kind=8) :: bpts(d,0:p), ebpts(d,0:p+t), nextbpts(d,0:p-2)
  real   (kind=8) :: alfs(0:p-2), ua, ub, alf, bet, gam, den
  if (t < 1) then
     Uh = U
     Qw = Pw
     return
  end if
  m = n + p + 1
  ph = p + t
  ! Bezier coefficients
  bezalfs(0,0)  = 1.0
  bezalfs(ph,p) = 1.0
  do i = 1, ph/2
     do j = max(0,i-t), min(p,i)
        bezalfs(i,j) = Bin(p,j)*Bin(t,i-j)*(1.0d+0/Bin(ph,i))
     end do
  end do
  do i = ph/2+1, ph-1
     do j = max(0,i-t), min(p,i)
        bezalfs(i,j) = bezalfs(ph-i,p-j)
     end do
  end do
  kind = ph+1
  cind = 1
  r = -1
  a = p
  b = p+1
  ua = U(a)
  Uh(0:ph) = ua
  Qw(:,0) = Pw(:,0)
  bpts = Pw(:,0:p)
  do while (b < m)
     i = b
     do while (b < m)
        if (U(b) /= U(b+1)) exit
        b = b + 1
     end do
     mul = b - i + 1
     oldr = r
     r = p - mul
     ub = U(b)
     if (oldr > 0) then
        lbz = (oldr+2)/2
     else
        lbz = 1
     end if
     if (r > 0) then
        rbz = ph - (r+1)/2
     else
        rbz = ph
     end if
     ! insert knots
     if (r > 0) then
        do k = p, mul+1, -1
           alfs(k-mul-1) = (ub-ua)/(U(a+k)-ua)
        end do
        do j = 1, r
           s = mul + j
           do k = p, s, -1
              bpts(:,k) = alfs(k-s)  * bpts(:,k) + &
                     (1.0-alfs(k-s)) * bpts(:,k-1)
           end do
           nextbpts(:,r-j) = bpts(:,p)
        end do
     end if
     ! degree elevate
     do i = lbz, ph
        ebpts(:,i) = 0.0
        do j = max(0,i-t), min(p,i)
           ebpts(:,i) = ebpts(:,i) + bezalfs(i,j)*bpts(:,j)
        end do
     end do
     ! remove knots
     if (oldr > 1) then
        first = kind-2
        last = kind
        den = ub-ua
        bet = (ub-Uh(kind-1))/den
        do tr = 1, oldr-1
           i = first
           j = last
           kj = j-kind+1
           do while (j-i > tr)
              if (i < cind) then
                 alf = (ub-Uh(i))/(ua-Uh(i))
                 Qw(:,i) = alf*Qw(:,i) + (1.0-alf)*alf*Qw(:,i-1)
              end if
              if (j >= lbz) then
                 if (j-tr <= kind-ph+oldr) then
                    gam = (ub-Uh(j-tr))/den
                    ebpts(:,kj) = gam*ebpts(:,kj) + (1.0-gam)*ebpts(:,kj+1)
                 else
                    ebpts(:,kj) = bet*ebpts(:,kj) + (1.0-bet)*ebpts(:,kj+1)
                 end if
              end if
              i = i+1
              j = j-1
              kj = kj-1
           end do
           first = first-1
           last = last+1
        end do
     end if
     !
     if (a /= p) then
        do i = 0, ph-oldr-1
           Uh(kind) = ua
           kind = kind+1
        end do
     end if
     do j = lbz, rbz
        Qw(:, cind) = ebpts(:,j)
        cind = cind+1
     end do
     !
     if (b < m) then
        bpts(:,0:r-1) = nextbpts(:,0:r-1)
        bpts(:,r:p) = Pw(:,b-p+r:b)
        a = b
        b = b+1
        ua = ub
     else
        Uh(kind:kind+ph) = ub
     end if
  end do
contains
  pure function Bin(n,k) result (C)
    implicit none
    integer(kind=4), intent(in) :: n, k
    integer(kind=4) :: i, C
    C = 1
    do i = 0, min(k,n-k) - 1
       C = C * (n - i)
       C = C / (i + 1)
    end do
  end function Bin
end subroutine DegreeElevate

end module bspline

!
! ----------
!

module BSp
contains

subroutine FindSpan(p,m,U,uu,span)
  use bspline, FindS => FindSpan
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(out) :: span
  span = FindS(m-(p+1),p,uu,U)
end subroutine FindSpan

subroutine FindMult(p,m,U,uu,span,mult)
  use bspline, FindM => FindMult
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(in)  :: span
  integer(kind=4), intent(out) :: mult
  integer(kind=4) :: k
  if (span >= 0) then
     k = span
  else
     k = FindSpan(m-(p+1),p,uu,U)
  end if
  mult = FindM(k,uu,p,U)
end subroutine FindMult

subroutine FindSpanMult(p,m,U,uu,k,s)
  use bspline, FindSM => FindSpanMult
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(out) :: k, s
  call FindSM(m-(p+1),p,uu,U,k,s)
end subroutine FindSpanMult

subroutine EvalBasisFuns(p,m,U,uu,span,N)
  use bspline
  implicit none
  integer(kind=4), intent(in) :: p, m, span
  real   (kind=8), intent(in) :: U(0:m), uu
  real   (kind=8), intent(out):: N(0:p)
  integer(kind=4) :: i
  if (span >= 0) then
     i = span
  else
     i = FindSpan(m-(p+1),p,uu,U)
  end if
  call BasisFuns(i,uu,p,U,N)
end subroutine EvalBasisFuns

subroutine EvalBasisFunsDers(p,m,U,uu,d,span,dN)
  use bspline
  implicit none
  integer(kind=4), intent(in) :: p, m, d, span
  real   (kind=8), intent(in) :: U(0:m), uu
  real   (kind=8), intent(out):: dN(0:p,0:d)
  integer(kind=4) :: i
  if (span >= 0) then
     i = span
  else
     i = FindSpan(m-(p+1),p,uu,U)
  end if
  call DersBasisFuns(i,uu,p,d,U,dN)
end subroutine EvalBasisFunsDers

subroutine InsertKnot(d,n,p,U,Pw,uu,r,V,Qw)
  use bspline, InsKnt => InsertKnot
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(out) :: V(0:n+p+1+r)
  real   (kind=8), intent(out) :: Qw(d,0:n+r)
  integer(kind=4) :: k, s
  if (r == 0) then
     V = U; Qw = Pw; return
  end if
  call FindSpanMult(n,p,uu,U,k,s)
  call InsKnt(d,n,p,U,Pw,uu,k,s,r,V,Qw)
end subroutine InsertKnot

subroutine RemoveKnot(d,n,p,U,Pw,uu,r,t,V,Qw,TOL)
  use bspline, RemKnt => RemoveKnot
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: r
  integer(kind=4), intent(out) :: t
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  real   (kind=8), intent(in)  :: TOL
  integer(kind=4) :: k, s
  t = 0
  V = U
  Qw = Pw
  if (r == 0) return
  if (uu <= U(p)) return
  if (uu >= U(n+1)) return
  call FindSpanMult(n,p,uu,U,k,s)
  call RemKnt(d,n,p,V,Qw,uu,k,s,r,t,TOL)
end subroutine RemoveKnot

subroutine Clamp(d,n,p,U,Pw,V,Qw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  V  = U
  Qw = Pw
  call ClampKnot(d,n,p,V,Qw)
end subroutine Clamp

subroutine Unclamp(d,n,p,U,Pw,V,Qw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  V  = U
  Qw = Pw
  call UnclampKnot(d,n,p,V,Qw)
end subroutine Unclamp

subroutine RefineKnotVector(d,n,p,U,Pw,r,X,Ubar,Qw)
  use bspline, RefKnt => RefineKnotVector
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Ubar(0:n+r+1+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n+r+1)
  call RefKnt(d,n,p,U,Pw,r,X,Ubar,Qw)
end subroutine RefineKnotVector

subroutine DegreeElevate(d,n,p,U,Pw,t,nh,Uh,Qw)
  use bspline, DegElev => DegreeElevate
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: t
  integer(kind=4), intent(in)  :: nh
  real   (kind=8), intent(out) :: Uh(0:nh+p+t+1)
  real   (kind=8), intent(out) :: Qw(d,0:nh)
  call DegElev(d,n,p,U,Pw,t,nh,Uh,Qw)
end subroutine DegreeElevate

subroutine Extract(d,n,p,U,Pw,x,Cw)
  use bspline, CornerCut => CurvePntByCornerCut
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: x
  real   (kind=8), intent(out) :: Cw(d)
  call CornerCut(d,n,p,U,Pw,x,Cw)
end subroutine Extract

end module BSp

!
! ----------
!

module IGA
contains

subroutine KnotVector(E,p,C,Ui,Uf,periodic,m,U)
  implicit none
  integer(kind=4), intent(in)  :: E,p,C
  real   (kind=8), intent(in)  :: Ui,Uf
  logical(kind=4), intent(in)  :: periodic
  integer(kind=4), intent(in)  :: m
  real   (kind=8), intent(out) :: U(0:m)
  integer(kind=4) :: i, j, k
  real   (kind=8) :: dU
  dU = (Uf-Ui)/E
  k = p+1
  do i = 1, (E-1)
     do j = 1, (p-C)
        U(k) = Ui + i*dU
        k = k+1
     end do
  end do
  do k = 0, p
     U(k)   = Ui
     U(m-k) = Uf
  end do
  if (periodic) then
     do k = 0, C
        U(k)     = Ui - Uf + U(m-C-(p+1)+k)
        U(m-C+k) = Uf - Ui + U(p+1+k)
     end do
  end if
end subroutine KnotVector

subroutine GaussRule(q,X,W)
  implicit none
  integer(kind=4), intent(in)  :: q
  real   (kind=8), intent(out) :: X(0:q-1)
  real   (kind=8), intent(out) :: W(0:q-1)
  select case (q)
  case (1) ! p = 1
     X(0) = 0.0
     W(0) = 2.0
  case (2) ! p = 3
     X(0) = -0.5773502691896257645091487805019576 ! 1/sqrt(3)
     X(1) = -X(0)
     W(0) =  1.0
     W(1) =  W(0)
  case (3) ! p = 5
     X(0) = -0.7745966692414833770358530799564799 ! sqrt(3/5)
     X(1) =  0.0
     X(2) = -X(0)
     W(0) =  0.5555555555555555555555555555555556 ! 5/9
     W(1) =  0.8888888888888888888888888888888889 ! 8/9
     W(2) =  W(0)
  case (4) ! p = 7
     X(0) = -0.8611363115940525752239464888928094 ! sqrt((3+2*sqrt(6/5))/7)
     X(1) = -0.3399810435848562648026657591032448 ! sqrt((3-2*sqrt(6/5))/7)
     X(2) = -X(1)
     X(3) = -X(0)
     W(0) =  0.3478548451374538573730639492219994 ! (18-sqrt(30))/36
     W(1) =  0.6521451548625461426269360507780006 ! (18+sqrt(30))/36
     W(2) =  W(1)
     W(3) =  W(0)
  case (5) ! p = 9
     X(0) = -0.9061798459386639927976268782993929 ! 1/3*sqrt(5+2*sqrt(10/7))
     X(1) = -0.5384693101056830910363144207002086 ! 1/3*sqrt(5-2*sqrt(10/7))
     X(2) =  0.0
     X(3) = -X(1)
     X(4) = -X(0)
     W(0) =  0.2369268850561890875142640407199173 ! (322-13*sqrt(70))/900
     W(1) =  0.4786286704993664680412915148356382 ! (322+13*sqrt(70))/900
     W(2) =  0.5688888888888888888888888888888889 ! 128/225
     W(3) =  W(1)
     W(4) =  W(0)
  case (6)
     X(0) = -0.9324695142031520
     X(1) = -0.6612093864662645
     X(2) = -0.2386191860831969
     X(3) = -X(2)
     X(4) = -X(1)
     X(5) = -X(0)
     W(0) =  0.1713244923791703
     W(1) =  0.3607615730481386
     W(2) =  0.4679139345726911
     W(3) =  W(2)
     W(4) =  W(1)
     W(5) =  W(0)
  case (7)
     X(0) = -0.9491079123427585
     X(1) = -0.7415311855993944
     X(2) = -0.4058451513773972
     X(3) =  0.0
     X(4) = -X(2)
     X(5) = -X(1)
     X(6) = -X(0)
     W(0) =  0.1294849661688697
     W(1) =  0.2797053914892767
     W(2) =  0.3818300505051189
     W(3) =  0.4179591836734694
     W(4) =  W(2)
     W(5) =  W(1)
     W(6) =  W(0)
  case (8)
     X(0) = -0.9602898564975362
     X(1) = -0.7966664774136267
     X(2) = -0.5255324099163290
     X(3) = -0.1834346424956498
     X(4) = -X(3)
     X(5) = -X(2)
     X(6) = -X(1)
     X(7) = -X(0)
     W(0) =  0.1012285362903763
     W(1) =  0.2223810344533745
     W(2) =  0.3137066458778873
     W(3) =  0.3626837833783620
     W(4) =  W(3)
     W(5) =  W(2)
     W(6) =  W(1)
     W(7) =  W(0)
  case (9)
     X(0) = -0.9681602395076261
     X(1) = -0.8360311073266358
     X(2) = -0.6133714327005904
     X(3) = -0.3242534234038089
     X(4) =  0.0
     X(5) = -X(3)
     X(6) = -X(2)
     X(7) = -X(1)
     X(8) = -X(0)
     W(0) =  0.0812743883615744
     W(1) =  0.1806481606948574
     W(2) =  0.2606106964029354
     W(3) =  0.3123470770400029
     W(4) =  0.3302393550012598
     W(5) =  W(3)
     W(6) =  W(2)
     W(7) =  W(1)
     W(8) =  W(0)
  case (10)
     X(0) = -0.973906528517172
     X(1) = -0.865063366688985
     X(2) = -0.679409568299024
     X(3) = -0.433395394129247
     X(4) = -0.148874338981631
     X(5) = -X(4)
     X(6) = -X(3)
     X(7) = -X(2)
     X(8) = -X(1)
     X(9) = -X(0)
     W(0) =  0.066671344308688
     W(1) =  0.149451349150581
     W(2) =  0.219086362515982
     W(3) =  0.269266719309996
     W(4) =  0.295524224714753
     W(5) =  W(4)
     W(6) =  W(3)
     W(7) =  W(2)
     W(8) =  W(1)
     W(9) =  W(0)
  case default
     X = 0.0
     W = 0.0
  end select
end subroutine GaussRule

subroutine BasisData(p,m,U,d,q,r,O,J,W,X,N)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m)
  integer(kind=4), intent(in)  :: d, q, r
  integer(kind=4), intent(out) :: O(r)
  real   (kind=8), intent(out) :: J(r)
  real   (kind=8), intent(out) :: W(q)
  real   (kind=8), intent(out) :: X(q,r)
  real   (kind=8), intent(out) :: N(0:d,0:p,q,r)

  integer(kind=4) i, iq, ir
  real   (kind=8) uu, Xg(q)
  real   (kind=8) basis(0:p,0:d)

  ir = 1
  do i = p, m-p
     if (U(i) /= U(i+1)) then
        O(ir) = i - p
        ir = ir + 1
     end if
  end do
  call GaussRule(q,Xg,W)
  do ir = 1, r
     i = O(ir) + p
     J(ir) = (U(i+1)-U(i))/2.0
     X(:,ir) = (Xg + 1.0) * J(ir) + U(i)
     do iq = 1, q
        uu = X(iq,ir)
        call DersBasisFuns(i,uu,p,d,U,basis)
        N(:,:,iq,ir) = transpose(basis)
     end do
  end do
end subroutine BasisData

end module IGA

!
! ----------
!

module Crv
contains

subroutine Evaluate(d,n,p,U,Pw,r,X,Cw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Cw(d,0:r)
  integer(kind=4) :: i, j, span
  real   (kind=8) :: basis(0:p), C(d)
  !
  do i = 0, r
     span = FindSpan(n,p,X(i),U)
     call BasisFuns(span,X(i),p,U,basis)
     !
     C = 0.0
     do j = 0, p
        C = C + basis(j)*Pw(:,span-p+j)
     end do
     Cw(:,i) = C
     !
  end do
  !
end subroutine Evaluate

subroutine RefineKnotVector(d,n,p,U,Pw,r,X,Ubar,Qw)
  use bspline, RefKnt => RefineKnotVector
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Ubar(0:n+r+1+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n+r+1)
  call RefKnt(d,n,p,U,Pw,r,X,Ubar,Qw)
end subroutine RefineKnotVector

subroutine DegreeElevate(d,n,p,U,Pw,t,nh,Uh,Qw)
  use bspline, DegElev => DegreeElevate
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: t
  integer(kind=4), intent(in)  :: nh
  real   (kind=8), intent(out) :: Uh(0:nh+p+t+1)
  real   (kind=8), intent(out) :: Qw(d,0:nh)
  call DegElev(d,n,p,U,Pw,t,nh,Uh,Qw)
end subroutine DegreeElevate

end module Crv

!
! ----------
!

module Srf
contains

subroutine Evaluate(d,nx,px,Ux,ny,py,Uy,Pw,rx,X,ry,Y,Cw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny
  integer(kind=4), intent(in)  :: px, py
  integer(kind=4), intent(in)  :: rx, ry
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Pw(d,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry)
  real   (kind=8), intent(out) :: Cw(d,0:ry,0:rx)
  integer(kind=4) :: ix, jx, iy, jy, ox, oy
  integer(kind=4) :: spanx(0:rx), spany(0:ry)
  real   (kind=8) :: basisx(0:px,0:rx), basisy(0:py,0:ry)
  real   (kind=8) :: M, C(d)
  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call BasisFuns(spanx(ix),X(ix),px,Ux,basisx(:,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call BasisFuns(spany(iy),Y(iy),py,Uy,basisy(:,iy))
  end do
  !
  do ix = 0, rx
     ox = spanx(ix) - px
     do iy = 0, ry
        oy = spany(iy) - py
        ! ---
        C = 0.0
        do jx = 0, px
           do jy = 0, py
              M = basisx(jx,ix) * basisy(jy,iy)
              C = C + M * Pw(:,oy+jy,ox+jx)
           end do
        end do
        Cw(:,iy,ix) = C
        ! ---
     end do
  end do
  !
end subroutine Evaluate

subroutine RefineKnotVector(d,n,p,U,m,q,V,Pw,r,X,s,Y,Ubar,Vbar,Qw)
  use bspline, RefKnt => RefineKnotVector
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p, m, q
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: V(0:m+q+1)
  real   (kind=8), intent(in)  :: Pw(d,0:m,0:n)
  integer(kind=4), intent(in)  :: r, s
  real   (kind=8), intent(in)  :: X(0:r), Y(0:s)
  real   (kind=8), intent(out) :: Ubar(0:n+p+1+r+1)
  real   (kind=8), intent(out) :: Vbar(0:m+q+1+s+1)
  real   (kind=8), intent(out) :: Qw(d,0:m+s+1,0:n+r+1)
  real   (kind=8) :: Q1(d*(m+1),0:n+r+1)
  real   (kind=8) :: Aw(d,0:n+r+1,0:m)
  real   (kind=8) :: Q2(d*(n+r+1+1),0:m+s+1)
  !
  if (r >= 0) then
     call RefKnt(size(Q1,1),n,p,U,Pw,r,X,Ubar,Q1)
     Aw = reshape(Q1,shape(Aw),order=(/1,3,2/))
  else
     Ubar = U
     Aw = reshape(Pw,shape(Aw),order=(/1,3,2/))
  end if
  !
  if (s >= 0) then
     call RefKnt(size(Q2,1),m,q,V,Aw,s,Y,Vbar,Q2)
     Qw = reshape(Q2,shape(Qw),order=(/1,3,2/))
  else
     Vbar = V
     Qw = reshape(Aw,shape(Qw),order=(/1,3,2/))
  end if
  !
end subroutine RefineKnotVector

subroutine DegreeElevate(d,n,p,U,m,q,V,Pw,r,s,nh,Uh,mh,Vh,Qw)
  use bspline, DegElev => DegreeElevate
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, m
  integer(kind=4), intent(in)  :: p, q
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: V(0:m+q+1)
  real   (kind=8), intent(in)  :: Pw(d,0:m,0:n)
  integer(kind=4), intent(in)  :: r, s
  integer(kind=4), intent(in)  :: nh, mh
  real   (kind=8), intent(out) :: Uh(0:nh+p+r+1)
  real   (kind=8), intent(out) :: Vh(0:mh+q+s+1)
  real   (kind=8), intent(out) :: Qw(d,0:mh,0:nh)
  real   (kind=8) :: Q1(d*(m+1),0:nh)
  real   (kind=8) :: Aw(d,0:nh,0:m)
  real   (kind=8) :: Q2(d*(nh+1),0:mh)
  !
  if (r > 0) then
     call DegElev(size(Q1,1),n,p,U,Pw,r,nh,Uh,Q1)
     Aw = reshape(Q1,shape(Aw),order=(/1,3,2/))
  else
     Uh = U
     Aw = reshape(Pw,shape(Aw),order=(/1,3,2/))
  end if
  !
  if (s > 0) then
     call DegElev(size(Q2,1),m,q,V,Aw,s,mh,Vh,Q2)
     Qw = reshape(Q2,shape(Qw),order=(/1,3,2/))
  else
     Vh = V
     Qw = reshape(Aw,shape(Qw),order=(/1,3,2/))
  end if
  !
end subroutine DegreeElevate

subroutine Extract(d,nx,px,Ux,ny,py,Uy,Pw,ii,uu,n,p,U,Cw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny
  integer(kind=4), intent(in)  :: px, py
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Pw(d,0:ny,0:nx)
  integer(kind=4), intent(in)  :: ii
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: n
  integer(kind=4), intent(in)  :: p
  real   (kind=8), intent(out) :: U(0:n+p+1)
  real   (kind=8), intent(out) :: Cw(d,0:n)
  real   (kind=8)  :: Pw1(d,0:nx,0:ny)
  !
  if (ii == 0) then
     call CurvePntByCornerCut(d*(ny+1),nx,px,Ux,Pw,uu,Cw)
     U = Uy
     return
  end if
  !
  if (ii == 1) then
     Pw1 = reshape(Pw,shape(Pw1),order=(/1,3,2/))
     call CurvePntByCornerCut(d*(nx+1),ny,py,Uy,Pw1,uu,Cw)
     U = Ux
     return
  end if
  !
end subroutine Extract

end module Srf

!
! ----------
!

module Vol
contains

subroutine Evaluate(d,nx,px,Ux,ny,py,Uy,nz,pz,Uz,Pw,rx,X,ry,Y,rz,Z,Cw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny, nz
  integer(kind=4), intent(in)  :: px, py, pz
  integer(kind=4), intent(in)  :: rx, ry, rz
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nz,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry), Z(0:rz)
  real   (kind=8), intent(out) :: Cw(d,0:rz,0:ry,0:rx)
  integer(kind=4) :: ix, jx, iy, jy, iz, jz, ox, oy, oz
  integer(kind=4) :: spanx(0:rx), spany(0:ry), spanz(0:rz)
  real   (kind=8) :: basisx(0:px,0:rx), basisy(0:py,0:ry), basisz(0:pz,0:rz)
  real   (kind=8) :: M, C(d)
  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call BasisFuns(spanx(ix),X(ix),px,Ux,basisx(:,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call BasisFuns(spany(iy),Y(iy),py,Uy,basisy(:,iy))
  end do
  do iz = 0, rz
     spanz(iz) = FindSpan(nz,pz,Z(iz),Uz)
     call BasisFuns(spanz(iz),Z(iz),pz,Uz,basisz(:,iz))
  end do
  !
  do ix = 0, rx
     ox = spanx(ix) - px
     do iy = 0, ry
        oy = spany(iy) - py
        do iz = 0, rz
           oz = spanz(iz) - pz
           ! ---
           C = 0.0
           do jx = 0, px
              do jy = 0, py
                 do jz = 0, pz
                    M = basisx(jx,ix) * basisy(jy,iy) * basisz(jz,iz)
                    C = C + M * Pw(:,oz+jz,oy+jy,ox+jx)
                 end do
              end do
           end do
           Cw(:,iz,iy,ix) = C
           ! ---
        end do
     end do
  end do
  !
end subroutine Evaluate

subroutine RefineKnotVector(d,nx,px,Ux,ny,py,Uy,nz,pz,Uz,Pw,rx,X,ry,Y,rz,Z,Vx,Vy,Vz,Qw)
  use bspline, RefKnt => RefineKnotVector
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny, nz
  integer(kind=4), intent(in)  :: px, py, pz
  integer(kind=4), intent(in)  :: rx, ry, rz
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nz,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry), Z(0:rz)
  real   (kind=8), intent(out) :: Vx(0:nx+px+1+rx+1)
  real   (kind=8), intent(out) :: Vy(0:ny+py+1+ry+1)
  real   (kind=8), intent(out) :: Vz(0:nz+pz+1+rz+1)
  real   (kind=8), intent(out) :: Qw(d,0:nz+rz+1,0:ny+ry+1,0:nx+rx+1)
  real   (kind=8) :: Q1(d*(nz+1)*(ny+1),0:nx+rx+1)
  real   (kind=8) :: Aw(d,0:nz,0:nx+rx+1,0:ny)
  real   (kind=8) :: Q2(d*(nz+1)*(nx+rx+1+1),0:ny+ry+1)
  real   (kind=8) :: Bw(d,0:ny+ry+1,0:nx+rx+1,0:nz)
  real   (kind=8) :: Q3(d*(ny+ry+1+1)*(nx+rx+1+1),0:nz+rz+1)
  !
  if (rx >= 0) then
     call RefKnt(size(Q1,1),nx,px,Ux,Pw,rx,X,Vx,Q1)
     Aw = reshape(Q1,shape(Aw),order=(/1,2,4,3/))
  else
     Vx = Ux
     Aw = reshape(Pw,shape(Aw),order=(/1,2,4,3/))
  end if
  !
  if (ry >= 0) then
     call RefKnt(size(Q2,1),ny,py,Uy,Aw,ry,Y,Vy,Q2)
     Bw = reshape(Q2,shape(Bw),order=(/1,4,3,2/))
  else
     Vy = Uy
     Bw = reshape(Aw,shape(Bw),order=(/1,4,3,2/))
  end if
  !
  if (rz >= 0) then
     call RefKnt(size(Q3,1),nz,pz,Uz,Bw,rz,Z,Vz,Q3)
     Qw = reshape(Q3,shape(Qw),order=(/1,3,4,2/))
  else
     Vz = Uz
     Qw = reshape(Bw,shape(Qw),order=(/1,3,4,2/))
  end if
  !
end subroutine RefineKnotVector

subroutine DegreeElevate(d,nx,px,Ux,ny,py,Uy,nz,pz,Uz,Pw,tx,ty,tz,ox,Vx,oy,Vy,oz,Vz,Qw)
  use bspline, DegElev => DegreeElevate
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny, nz
  integer(kind=4), intent(in)  :: px, py, pz
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nz,0:ny,0:nx)
  integer(kind=4), intent(in)  :: tx, ty, tz
  integer(kind=4), intent(in)  :: ox, oy, oz
  real   (kind=8), intent(out) :: Vx(0:ox+px+tx+1)
  real   (kind=8), intent(out) :: Vy(0:oy+py+ty+1)
  real   (kind=8), intent(out) :: Vz(0:oz+pz+tz+1)
  real   (kind=8), intent(out) :: Qw(d,0:oz,0:oy,0:ox)
  !
  real   (kind=8) :: Q1(d*(nz+1)*(ny+1),0:ox)
  real   (kind=8) :: Aw(d,0:nz,0:ox,0:ny)
  real   (kind=8) :: Q2(d*(nz+1)*(ox+1),0:oy)
  real   (kind=8) :: Bw(d,0:oy,0:ox,0:nz)
  real   (kind=8) :: Q3(d*(oy+1)*(ox+1),0:oz)
  !
  if (tx > 0) then
     call DegElev(size(Q1,1),nx,px,Ux,Pw,tx,ox,Vx,Q1)
     Aw = reshape(Q1,shape(Aw),order=(/1,2,4,3/))
  else
     Vx = Ux
     Aw = reshape(Pw,shape(Aw),order=(/1,2,4,3/))
  end if
  !
  if (ty > 0) then
     call DegElev(size(Q2,1),ny,py,Uy,Aw,ty,oy,Vy,Q2)
     Bw = reshape(Q2,shape(Bw),order=(/1,4,3,2/))
  else
     Vy = Uy
     Bw = reshape(Aw,shape(Bw),order=(/1,4,3,2/))
  end if
  !
  if (tz > 0) then
     call DegElev(size(Q3,1),nz,pz,Uz,Bw,tz,oz,Vz,Q3)
     Qw = reshape(Q3,shape(Qw),order=(/1,3,4,2/))
  else
     Vz = Uz
     Qw = reshape(Bw,shape(Qw),order=(/1,3,4,2/))
  end if
  !
end subroutine DegreeElevate

subroutine Extract(d,nx,px,Ux,ny,py,Uy,nz,pz,Uz,Pw,ii,uu,n0,p0,U0,n1,p1,U1,Cw)
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny, nz
  integer(kind=4), intent(in)  :: px, py, pz
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nz,0:ny,0:nx)
  integer(kind=4), intent(in)  :: ii
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: n0, n1
  integer(kind=4), intent(in)  :: p0, p1
  real   (kind=8), intent(out) :: U0(0:n0+p0+1)
  real   (kind=8), intent(out) :: U1(0:n1+p1+1)
  real   (kind=8), intent(out) :: Cw(d,0:n1,0:n0)
  integer(kind=4) :: d0, d1, d2
  real   (kind=8) :: Pw1(d,0:nz,0:nx,0:ny)
  real   (kind=8) :: Pw2(d,0:ny,0:nx,0:nz)
  !
  if (ii == 0) then
     d0 = d*(ny+1)*(nz+1)
     call CurvePntByCornerCut(d0,nx,px,Ux,Pw,uu,Cw)
     U0 = Uy
     U1 = Uz
     return
  end if
  !
  if (ii == 1) then
     d1 = d*(nx+1)*(nz+1)
     Pw1 = reshape(Pw,shape(Pw1),order=(/1,2,4,3/))
     call CurvePntByCornerCut(d1,ny,py,Uy,Pw1,uu,Cw)
     U0 = Ux
     U1 = Uz
     return
  end if
  !
  if (ii == 2) then
     d2 = d*(nx+1)*(ny+1)
     Pw2 = reshape(Pw,shape(Pw2),order=(/1,4,2,3/))
     call CurvePntByCornerCut(d2,nz,pz,Uz,Pw2,uu,Cw)
     U0 = Ux
     U1 = Uy
     return
  end if
  !
end subroutine Extract

end module Vol
