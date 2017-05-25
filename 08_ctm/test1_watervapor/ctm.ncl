;==========================================================;
; Chemistry Transport Model                                ;
;----------------------------------------------------------;
; Xinyu Wen, Peking Univ, May 2017                         ;
; xwen@pku.edu.cn                                          ;
;==========================================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

external FORTRAN "./advec.so"

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


undef("output")
procedure output(fn:string,x:numeric,lev:numeric,lat:numeric,lon:numeric,ex:integer,status:integer)
local dims,u,v,h,fo,uu,vv,hh,nt
begin
	dims	= dimsizes(x)
	if (status.eq.0) then
		; creat a new file
		system("rm -f "+fn)
		fo	= addfile(fn,"c")
		filedimdef(fo,"time",-1,True)

		xx	= new((/1,dims(0),dims(1),dims(2),dims(3)/),float,-999)
      xx(0,:,:,:,:) = (/ x /)

		xx!0	= "time"
		xx!1	= "tracer"
      xx!2  = "lev"
      xx!3  = "lat"
      xx!4  = "lon"

      xx&lev   = lev
		xx&lat	= lat
		xx&lon	= lon

      xx@long_name = "tracer concentration (mixing ratio)"
      xx@units     = "g/kg"

		fo->conc	= xx
		fo->exp_id	= ex
      fo@model    = "PKU-CTM (Chemistry Transport Model)"
		fo@history	= "X. Wen, Peking Univ., May 2017"
	else
		; append to an exist file
		fo	= addfile(fn,"w")
		fdims	= getfiledimsizes(fo)
		nt	= fdims(0)
		fo->conc(nt,:,:,:,:)	= (/ x /)
	end if
	delete(fo)
end


;================
;===   MAIN   === 
;================
begin
	;----- Experiment ID -----
	exp_id	= 20170521	; such as YYYYMMDD

	;----- Physics constants -----
	Re	= 6.371e6	      ; m
	g	= 9.8		         ; m/ss
	PI	= 3.1415926

	;----- U/V/W from Model -----
   fun   = "uwnd.2010.nc"
   fvn   = "vwnd.2010.nc"
   fwn   = "omega.2010.nc"
   ftn   = "air.2010.nc"

   fu    = addfile(fun,"r")
   fv    = addfile(fvn,"r")
   fw    = addfile(fwn,"r")
   ft    = addfile(ftn,"r")

	Nx	   = dimsizes(fu->lon)     ; total number of longitidues : 0 to 360E
	Ny	   = dimsizes(fu->lat)-1   ; total number of latitudes   : -90 to 90N  (NEED AVERAGE TO REDUCE LATITUDES)
   Nz    = dimsizes(fu->level)   ; total number of levels      : 17 standard levels in NCEP2
   Nt    = dimsizes(fu->time)    ; total number of timesteps   : 1460 for 365 days w/ 6-hourly data
   Ntracer     = 1               ; how many tracers to be simulated

   lon   = fu->lon
   lat2  = fu->lat(::-1)
   lat   = (lat2(1:Ny)+lat2(0:Ny-1))/2.0
   lev   = fu->level

   dtgcm = 21600                 ; GCM timestep in seconds
	dt	   = 600	                  ; CTM timestep in seconds
   Nsmallstep  = dtgcm/dt        ; here = 36
   oneday= 86400.0               ; seconds in ONE day

	fname	= "output."+exp_id+".nc"	; output file name

   ;----- emission -----
   emission = new((/Nz,Ny,Nx/),float,-999)
   emission = 0.0
   do k = 1, 1 ; just add at surface
      do j = 1, Ny
         emission(k-1,j-1,:)  = 5.0 * cos(2.0* lat(j-1)/180.0*PI) + 5.0
      end do
   end do

   switch_tauemi  = True
   tauemi         = new((/Ntracer/),float,-999)
   tauemi         = (/ oneday*20 /)

   ;----- sink -----
   switch_tau  = True
   tau         = new((/Ntracer/),float,-999)
   tau         = (/ oneday*10 /)

	;----- init core variables -----
	conc        = new((/Ntracer,Nz,Ny,Nx/),float,-999)
	conc        = 0.0

   do ivar = 1, Ntracer
      conc(ivar-1,:,:,:) = conc(ivar-1,:,:,:) + emission
   end do

	;----- Wrap-up -----
	exp_id!0       = "id"
	exp_id@g       = g
	exp_id@Re      = Re

	output(fname,conc,lev,lat,lon,exp_id,0)
	print("Let's Go Go Go!!!")


	;-------------------
	;---  Main Loop  ---
	;-------------------
	do i = 1, Nt
		; BC forcing
      u2 = short2flt(fu->uwnd (i-1,:,::-1,:))
      u  = (u2(:,1:Ny,:)+u2(:,0:Ny-1,:))/2.0
      v2 = short2flt(fv->vwnd (i-1,:,::-1,:))
      v  = (v2(:,1:Ny,:)+v2(:,0:Ny-1,:))/2.0
      w2 = short2flt(fw->omega(i-1,:,::-1,:))
      w  = (w2(:,1:Ny,:)+w2(:,0:Ny-1,:))/2.0
      t2 = short2flt(ft->air(i-1,:,::-1,:))
      t  = (t2(:,1:Ny,:)+t2(:,0:Ny-1,:))/2.0

      do j = 1, Nsmallstep
		   ; tend
         tend  = 0.0*conc
		   FORTRAN::tran3d_tend(conc,Ntracer,Nz,Ny,Nx,u,v,w,t,lev,lat,lon,Re,tend)

         ; filter 3-sigma tend
         do ivar = 1, Ntracer
            tend_mean            = avg(tend(ivar-1,:,:,:))
            tend_3stddev         = stddev(tend(ivar-1,:,:,:))*3.0
            tend_anomaly         = tend(ivar-1,:,:,:)-tend_mean
            tend(ivar-1,:,:,:)   = where(abs(tend_anomaly).gt.tend_3stddev, 0.0, tend(ivar-1,:,:,:))
         end do

         ; source and sink
         do ivar = 1, Ntracer
            ; source : add emission
            if (switch_tauemi) then
               tend(ivar-1,:,:,:) = tend(ivar-1,:,:,:) + emission/tauemi(ivar-1)
            end if

            ; sink   : removal with Tau (lift-time)
            if (switch_tau) then
               tend(ivar-1,:,:,:) = tend(ivar-1,:,:,:) - conc(ivar-1,:,:,:)/tau(ivar-1)
            end if
         end do

         ; polar filter
         tend(:,:,0,:)     = 0.0
         tend(:,:,Ny-1,:)  = 0.0

         ; Euler forward
         conc  = conc + dt*tend

	      ; print diagnostic info
         str_i = sprinti("%0.4i",i)
         str_j = sprinti("%0.2i",j)
         str_min = sprintf("%8.2e",min(tend))
         str_max = sprintf("%8.2e",max(tend))
	      print(str_i+" / "+str_j+"     Min(tend)="+str_min+"   Max(tend)="+str_max)	; min & max of tend
      end do

		; output
		output(fname,conc,lev,lat,lon,exp_id,1)
      print("===== writting outputs =====")

	end do

	print("----- Successfully Done!! -----")
end
