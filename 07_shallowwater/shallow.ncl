;==========================================================;
; Shallow Water Model                                      ;
;----------------------------------------------------------;
; Xinyu Wen, Peking Univ, May 2017                         ;
; xwen@pku.edu.cn                                          ;
;==========================================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


undef("cal_dNdx")
function cal_dNdx(x:numeric,d:float)
begin
	return(center_finite_diff_n(x,d,True,-999,1))
end

undef("cal_dNdy")
function cal_dNdy(x:numeric,d:float)
begin
	return(center_finite_diff_n(x,d,False,-999,0))
end

undef("add_gaussian_peak")
function add_gaussian_peak(N:numeric,x:integer,y:integer,r:float,amp:float)
local ii,jj
begin
	dims    = dimsizes(N)
	ii      = conform_dims((/dims(0),dims(1)/),ispan(1,dims(1),1),1)
	jj      = conform_dims((/dims(0),dims(1)/),ispan(1,dims(0),1),0)
	N       = N + amp*exp(-1.0*((ii-x)^2+(jj-y)^2)/(2.0*r*r))
	return(N)
end

undef("tend")
function tend(x:numeric,f:numeric,g:float,dx:float,dy:float)
local u,v,h,ux,uy,vx,vy,hx,hy,hux,hvy,tend_u,tend_v,tend_h
begin
	u 	= x(0,:,:)
	v 	= x(1,:,:)
	h 	= x(2,:,:)

	ux	= cal_dNdx(u,dx)
	uy	= cal_dNdy(u,dy)

	vx	= cal_dNdx(v,dx)
	vy	= cal_dNdy(v,dy)

	hx	= cal_dNdx(h,dx)
	hy	= cal_dNdy(h,dy)

	hux	= cal_dNdx(h*u,dx)
	hvy	= cal_dNdy(h*v,dy)

	tend_u	= -u*ux-v*uy-g*hx+f*v
	tend_v	= -u*vx-v*vy-g*hy-f*u
	tend_h	= -(hux+hvy)

	return((/tend_u,tend_v,tend_h/))
end

undef("output")
procedure output(fn:string,x:numeric,lat:numeric,lon:numeric,ex:integer,status:integer)
local dims,u,v,h,fo,uu,vv,hh,nt
begin
	u	= x(0,:,:)
	v	= x(1,:,:)
	h	= x(2,:,:)

	if (status.eq.0) then
		; creat a new file
		system("rm -f "+fn)
		fo	= addfile(fn,"c")
		filedimdef(fo,"time",-1,True)

		dims	= dimsizes(x)

		uu	= new((/1,dims(1),dims(2)/),float,-999)
		vv	= new((/1,dims(1),dims(2)/),float,-999)
		hh	= new((/1,dims(1),dims(2)/),float,-999)

		uu(0,:,:)	= u
		vv(0,:,:)	= v
		hh(0,:,:)	= h

		uu!0	= "time"
		vv!0	= "time"
		hh!0	= "time"

		uu!1	= "lat"
		vv!1	= "lat"
		hh!1	= "lat"

		uu&lat	= lat
		vv&lat	= lat
		hh&lat	= lat

		uu!2	= "lon"
		vv!2	= "lon"
		hh!2	= "lon"

		uu&lon	= lon
		vv&lon	= lon
		hh&lon	= lon

		fo->u	= uu
		fo->v	= vv
		fo->h	= hh 
		fo->exp_id	= ex
		fo@history	= "XW, PKU, May 2015"
	else
		; append to an exist file
		fo	= addfile(fn,"w")
		fdims	= getfiledimsizes(fo)
		nt	= fdims(0)
		fo->u(nt,:,:)	= (/ u /)
		fo->v(nt,:,:)	= (/ v /)
		fo->h(nt,:,:)	= (/ h /)
	end if
	delete(fo)
end

undef("shapiro2")
function shapiro2(x:numeric)
local dims,Ny,Nx,xe,i1,i2,j1,j2,shapiro2
begin
	dims	= dimsizes(x)
	Ny	= dims(0)
	Nx	= dims(1)

	xe	= new((/Ny+2,Nx+2/),float,-999)
	xe(1:Ny,1:Nx)	= x
	xe(0,1:Nx)	= x(Ny-1,:)
	xe(Ny+1,1:Nx)	= x(0,:)
	xe(1:Ny,0)	= x(:,Nx-1)
	xe(1:Ny,Nx+1)	= x(:,0)

	i1	= 1
	i2	= Nx
	j1	= 1
	j2	= Ny

	shapiro2 = 1.0 * ( xe(j1:j2,i1+1:i2+1) + \\
                      xe(j1:j2,i1-1:i2-1) + \\
                      xe(j1+1:j2+1,i1:i2) + \\
                      xe(j1-1:j2-1,i1:i2) ) + 4.0 * xe(j1:j2,i1:i2)
	return(shapiro2/8.0)
end


;================
;===   MAIN   === 
;================
begin
	;----- Experiment ID -----
	exp_id	= 20170511	; such as YYYYMMDD

	;----- Physics constants -----
	Re	= 6.371e6	; m
	g	= 9.8		; m/ss
	Omega	= 7.3e-5	; Rad/s
	Hscale	= 4000.0	; averaged depth of shallow water
	Pi	= 3.1415926

	;----- Try-Space -----
	;Re	= Re * 10
	;g	= g * 5
	Omega	= Omega * 10
	;Hscale	= Hscale * 2

	;----- Model parameters -----
	Ntime	= 1000		; total number of time steps
	Nx	= 360		;   0 to 360E
	Ny	= 120		; -60 to 60N

	dt	= 200.0		; seconds
	dx	= 1		; degree
	dy	= 1 		; degree

	lon		= fspan(0.5,  359.5, Nx/dx)
	lon@long_name	= "longitude"
	lon@units	= "degrees_east"

	lat		= fspan(-59.5, 59.5, Ny/dy)
	lat@long_name	= "latitude"
	lat@units	= "degrees_north"

	dxx	= Re*(dx/180.0*Pi)		; grid space in m
	dyy	= Re*(dy/180.0*Pi)

	f	= new((/Ny,Nx/),float,-999)	; coriolis parameter f
	do i = 0, Nx-1
		f(:,i) = 2.0*Omega*sin(lat/180.0*Pi)
	end do

	fname	= "output."+exp_id+".nc"	; output file name

	;----- Core variables -----
	u	= new((/Ny,Nx/),float,-999)
	v	= new((/Ny,Nx/),float,-999)
	h	= new((/Ny,Nx/),float,-999)

	;----- Initial Condition -----
	u	= 0.0
	v	= 0.0
	h	= Hscale

	; in case you want to set any special IC
	;h	= add_gaussian_peak(h,Nx/2,Ny/2-25, 5.0,100.0)
	;h	= add_gaussian_peak(h,Nx/2,Ny/2   , 5.0,100.0)
	;h	= add_gaussian_peak(h,Nx/2,Ny/2+25, 5.0,100.0)

	;----- Checking CFL condition -----
	cfl_lim = dxx/sqrt(g*Hscale)
	if (dt.gt.cfl_lim) then
		print("dt(current="+dt+") should be less than "+cfl_lim)
		exit
	end if

	;----- Wrap-up -----
	exp_id!0       = "id"
	exp_id@g       = g
	exp_id@Re      = Re
	exp_id@Omega   = Omega
	exp_id@Hscale  = Hscale
	exp_id@Ntime   = Ntime
	exp_id@dt      = dt
	exp_id@dx      = dx
	exp_id@dxx     = dxx
	exp_id@dy      = dy
	exp_id@dyy     = dyy
	exp_id@cfl_limit= cfl_lim


	;-------------------
	;---  Main Loop  ---
	;-------------------
	t0 = (/ u, v, h /)
	output(fname,t0,lat,lon,exp_id,0)
	print("Let's Go Go Go!!!")

	do i = 1, Ntime
		; BC forcing
		RandomBC = False
		if (RandomBC) then
		low     = 0.0
		high    = 3.0
		Nstorm  = floattoint(random_uniform(low,high,1))
		if (Nstorm.gt.0) then
			stormx  = floattoint(random_uniform(1,Nx,Nstorm))
			stormy  = floattoint(random_uniform(1,Ny,Nstorm))
			do k    = 0, Nstorm-1
				t0(2,:,:)       = add_gaussian_peak(t0(2,:,:),stormx(k),stormy(k),2,10)
			end do
			delete(stormx)
			delete(stormy)
		end if
		end if

      WesternlyBurst = True 
      if (WesternlyBurst) then
         forcingstep = 100
         forcingamp  = -10.0
         forcinglon  = 180
         if (i.lt.forcingstep) then
            ; North of Eq
            t0(0,Ny/2+4,forcinglon) = t0(0,Ny/2+4,forcinglon) + (0.2 * forcingamp - t0(0,Ny/2+4,forcinglon))*0.5
            t0(0,Ny/2+3,forcinglon) = t0(0,Ny/2+3,forcinglon) + (0.4 * forcingamp - t0(0,Ny/2+3,forcinglon))*0.5
            t0(0,Ny/2+2,forcinglon) = t0(0,Ny/2+2,forcinglon) + (0.6 * forcingamp - t0(0,Ny/2+2,forcinglon))*0.5
            t0(0,Ny/2+1,forcinglon) = t0(0,Ny/2+1,forcinglon) + (0.8 * forcingamp - t0(0,Ny/2+1,forcinglon))*0.5
            t0(0,Ny/2  ,forcinglon) = t0(0,Ny/2  ,forcinglon) + (1.0 * forcingamp - t0(0,Ny/2  ,forcinglon))*0.5
            ; South of Eq
            t0(0,Ny/2-1,forcinglon) = t0(0,Ny/2-1,forcinglon) + (1.0 * forcingamp - t0(0,Ny/2-1,forcinglon))*0.5
            t0(0,Ny/2-2,forcinglon) = t0(0,Ny/2-2,forcinglon) + (0.8 * forcingamp - t0(0,Ny/2-2,forcinglon))*0.5
            t0(0,Ny/2-3,forcinglon) = t0(0,Ny/2-3,forcinglon) + (0.6 * forcingamp - t0(0,Ny/2-3,forcinglon))*0.5
            t0(0,Ny/2-4,forcinglon) = t0(0,Ny/2-4,forcinglon) + (0.4 * forcingamp - t0(0,Ny/2-4,forcinglon))*0.5
            t0(0,Ny/2-5,forcinglon) = t0(0,Ny/2-5,forcinglon) + (0.2 * forcingamp - t0(0,Ny/2-5,forcinglon))*0.5
         end if
      end if

		; Runge-Kutta 4-order
		q1 = dt*tend(t0,		         f,g,dxx,dyy)
		q2 = dt*tend(t0+1.0/2.0*q1,	f,g,dxx,dyy)
		q3 = dt*tend(t0+1.0/2.0*q2,	f,g,dxx,dyy)
		q4 = dt*tend(t0+q3,		      f,g,dxx,dyy)
		t0 = t0+1.0/6.0*(q1+2.0*q2+2.0*q3+q4)

		; exception throw
		if (any(ismissing(t0))) then
			uavg	= avg(t0(0,:,:))
			vavg	= avg(t0(1,:,:))
			havg	= avg(t0(2,:,:))
			t0(0,:,:)	= where(ismissing(t0(0,:,:)),uavg,t0(0,:,:))
			t0(1,:,:)	= where(ismissing(t0(1,:,:)),vavg,t0(1,:,:))
			t0(2,:,:)	= where(ismissing(t0(2,:,:)),vavg,t0(2,:,:))
			print("------------------")
			print(" Missing Found!!! ")
			print("   make a patch   ")
			print("------------------")
		end if

		; non-linear filter
      filter_on         = True
		filter_interval	= 5
		if ((mod(i,filter_interval).eq.0).and.(filter_on)) then
			t0(0,:,:) = shapiro2(t0(0,:,:))	; applying filter on u
			t0(1,:,:) = shapiro2(t0(1,:,:))  ; applying filter on v
			t0(2,:,:) = shapiro2(t0(2,:,:))  ; applying filter on h
			print("... performing Shapiro2 nonlinear filter on U and V ...")
		end if

		; print diagnostic info
		print("-----------------")
		print("time-step = "+i)
		printMinMax(t0(0,:,:),0)	; min & max of u
		printMinMax(t0(1,:,:),0)	; min & max of v
		printMinMax(t0(2,:,:),0)	; min & max of h

		; output
      output_interval   = 1
      if (mod(i,output_interval).eq.0) then
		   output(fname,t0,lat,lon,exp_id,1)
      end if
	end do

	print("----- Successfully Done!! -----")
end
