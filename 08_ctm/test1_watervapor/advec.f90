subroutine tran3d_tend(conc,Nvar,Nz,Ny,Nx,u,v,w,t,lev,lat,lon,Re,tend)
   implicit none
   ! interface
   integer, intent(in) :: Nvar,Nz,Ny,Nx
   real   , intent(in) :: Re
   real   , intent(in) :: lev(Nz), lat(Ny), lon(Nx)
   real, dimension(Nx,Ny,Nz), intent(in) :: u,v,w,t
   real, dimension(Nx,Ny,Nz,Nvar), intent(in ) :: conc
   real, dimension(Nx,Ny,Nz,Nvar), intent(out) :: tend
   ! local
   real, dimension(Nx,Ny,Nz)   :: rho        ! air density
   real, dimension(Nx+1,Ny,Nz) :: fu         ! half-grid fluxes
   real, dimension(Nx,Ny+1,Nz) :: fv
   real, dimension(Nx,Ny,Nz+1) :: fw
   real :: latrad(Ny), lonrad(Nx)
   real, dimension(Ny) :: Re_coslat
   real :: dlon,dlat,dlonh,dlath
   real :: p_up,p_dw
   real, parameter :: Rdryair = 287.058      ! gas constant for dry air [J/(kg K)]
   integer :: ivar, i,j,k
  
   logical :: addW = .True.
   real, parameter :: PI = 3.141592653589793

   ! check info
   !print *, "min/max of u wind",       minval(u), maxval(u)
   !print *, "min/max of v wind",       minval(v), maxval(v)
   !print *, "min/max of w wind",       minval(w), maxval(w)
   !print *, "min/max of temperature",  minval(t), maxval(t)

   ! air density
   rho = 1.0
   do k = 1, Nz
      rho(:,:,k) = 100.0*lev(k) / (Rdryair*t(:,:,k))
      !print *, "min/max of rho(k)", k, minval(rho(:,:,k)), maxval(rho(:,:,k))
   end do

   ! degree to radius
   latrad = lat / 180.0 * PI
   lonrad = lon / 180.0 * PI
   Re_coslat = Re * cos(latrad)

   dlon  = lonrad(2)-lonrad(1)
   dlat  = latrad(2)-latrad(1)
   dlonh = 0.5*dlon
   dlath = 0.5*dlat

   do ivar = 1, Nvar
   !===== TENDENCY =====
      !--- PART 1: U
      do j = 1, Ny
         do k = 1, Nz
            call flux1d_full2half(Nx,u(1:Nx,j,k),conc(1:Nx,j,k,ivar),fu(1:Nx+1,j,k),0)
         end do
         tend(:,j,:,ivar) = tend(:,j,:,ivar) - ( fu(2:Nx+1,j,:) - fu(1:Nx,j,:) ) / ( Re_coslat(j) * dlon ) 
      end do

      !--- PART 2: V
      do i = 1, Nx
         do k = 1, Nz
            call flux1d_full2half(Ny,v(i,1:Ny,k),conc(i,1:Ny,k,ivar),fv(i,1:Ny+1,k),1)
         end do
         do j = 1, Ny
            tend(i,j,:,ivar) = tend(i,j,:,ivar) - &
                               ( fv(i,j+1,:)*cos(latrad(j)+dlath) - fv(i,j,:)*cos(latrad(j)-dlath) ) / &
                               ( Re_coslat(j) * dlat )
         end do
      end do

      !--- PART 3: W
      if (addW) then
      do i = 1, Nx
         do j = 1, Ny
            !call flux1d_full2half(Nz,w(i,j,:),conc(i,j,1:Nz,ivar)*rho(i,j,1:Nz),fw(i,j,1:Nz+1),1)
            call flux1d_full2half(Nz,w(i,j,:),conc(i,j,1:Nz,ivar),fw(i,j,1:Nz+1),1)
         end do
      end do
      do k = 2, Nz-1
         p_up = 0.5*(lev(k)+lev(k+1))*100
         p_dw = 0.5*(lev(k)+lev(k-1))*100
         tend(:,:,k,ivar) = tend(:,:,k,ivar) - ( fw(:,:,k+1) - fw(:,:,k) ) / ( p_dw - p_up ) * rho(:,:,k)
      end do
      ! surface-layer
      p_up = 0.5*(lev(1)+lev(2))*100
      p_dw = lev(1)*100
      tend(:,:,1,ivar) = tend(:,:,1,ivar) - (fw(:,:,2)-fw(:,:,1)) / (p_dw-p_up) * rho(:,:,1)
      ! top-layer
      p_up = lev(Nz)*100
      p_dw = 0.5*(lev(Nz)+lev(Nz-1))*100
      tend(:,:,Nz,ivar) = tend(:,:,Nz,ivar) - (fw(:,:,Nz+1)-fw(:,:,Nz)) / (p_dw-p_up) * rho(:,:,Nz)
      end if
   !===== TENDENCY =====
   end do

end subroutine

subroutine flux1d_full2half(N,u,c,f,option)
   implicit none
   integer, intent(in) :: N
   real   , intent(in) :: u(N), c(N)
   real   , intent(out):: f(N+1)
   integer, intent(in) :: option

   ! local
   real, dimension(N) :: uc
   real :: umid
   integer :: i

   uc = u * c

   ! option=0 : u : 最两边是cyclic边条件
   ! option=1 : v : 最两边半格点通量=0
   ! option=2 : w : 最两边半个点通量只局限在全格点上 也就是半层边条件

   ! leftmost ending
   if (option.eq.0) then
      umid = 0.5*(u(1)+u(N))
      if (umid.gt.0.0) then
         f(1) = uc(N)
      else
         f(1) = uc(1)
      end if
   elseif (option.eq.1) then
      f(1) = 0.0
   elseif (option.eq.2) then
      f(1) = 0.0
   end if

   ! inner area
   do i = 2, N
      umid = 0.5*(u(i)+u(i-1))
      if (umid.gt.0.0) then
         f(i) = uc(i-1)
      else
         f(i) = uc(i)
      end if
   end do

   ! rightmost ending
   if (option.eq.0) then
      f(N+1) = f(1)
   elseif (option.eq.1) then
      f(N+1) = 0.0
   elseif (option.eq.2) then
      f(N+1) = 0.0
   end if
   
end subroutine
