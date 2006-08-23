; FILE pltath.pro
;
; PURPOSE: contains IDL procedures to read and make plots from Athena dumps and
;   outputs.  Contains the following procedures:
;
; PRO readbin,filename:        reads .bin file 'filename'
; PRO sod_plot,filename:       plots analytic Sod solution over numerical
; PRO four_plot,filename:      plots d,P,Vx,P/d
; PRO nine_plot,filename,flag: use flag=0 for points, flag=1 for line
; PRO flines,nlev              plot 2D field lines
;
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

;-------------------------------------------------------------------------------
; Procedure READBIN: Reads ATHENA binary dumps
;
PRO readbin,filename
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
; Read number of zones and variables
;
ndata=LONARR(4)
openr,1,filename
readu,1,ndata
nx=ndata[0]
ny=ndata[1]
nz=ndata[2]
nvar=ndata[3]
;
; Read (gamma-1) and isothermal sound speed
;
dat=fltarr(2)
readu,1,dat
gamm1=dat[0]
isocs=dat[1]
;
; Read time,dt
;
readu,1,dat
time=dat[0]
dt=dat[1]
;
; Read grid coordinates
;
x=fltarr(nx)
readu,1,x
y=fltarr(ny)
readu,1,y
z=fltarr(nz)
readu,1,z
;
; Read data.
; nvar=4 means isothermal hydro.  nvar=5 means adiabatic hydro
; nvar=7 means isothermal MHD.    nvar=8 means adiabatic mhd
;
d =fltarr(nx,ny,nz)
e =fltarr(nx,ny,nz)
vx=fltarr(nx,ny,nz)
vy=fltarr(nx,ny,nz)
vz=fltarr(nx,ny,nz)
bx=fltarr(nx,ny,nz)
by=fltarr(nx,ny,nz)
bz=fltarr(nx,ny,nz)
readu,1,d
readu,1,vx
readu,1,vy
readu,1,vz
IF nvar eq 5 OR nvar eq 8 THEN readu,1,e
IF nvar eq 7 OR nvar eq 8 THEN BEGIN
  readu,1,bx
  readu,1,by
  readu,1,bz
ENDIF
vx=vx/d
vy=vy/d
vz=vz/d
IF gamm1 NE 0 THEN p=gamm1*(e-0.5*d*(vx^2+vy^2+vz^2)-0.5*(bx^2+by^2+bz^2))
IF gamm1 EQ 0 THEN p=isocs*isocs*d
;
close,1
END

;-------------------------------------------------------------------------------
; Procedure FOUR_PLOT: plots d,P,Vx,P/d
;
PRO four_plot,filename
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
readbin,filename
!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='P/D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
!P.MULTI=0
END

;-------------------------------------------------------------------------------
; Procedure NINE_PLOT: plots d,P,E,Vx,Vy,Vz,By,Bz,Phi - same plot as in RJ
;   Use flag=0 to plot points, flag=1 to plot line
;
PRO nine_plot,filename,flag
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
readbin,filename
!P.MULTI=[0,3,3,0,0]
dmin=MIN(d)
dmax=MAX(d)
IF (flag EQ 0) THEN plot,x,d,psym=6,symsize=.4,YTITLE='D',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,d,YTITLE='D',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(p)
dmax=MAX(p)
IF (flag EQ 0) THEN plot,x,p,psym=6,symsize=.4,YTITLE='P',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,p,YTITLE='P',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(e)
dmax=MAX(e)
IF (flag EQ 0) THEN plot,x,e,psym=6,symsize=.4,YTITLE='E',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,e,YTITLE='E',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vx)
dmax=MAX(vx)
IF (flag EQ 0) THEN plot,x,vx,psym=6,symsize=.4,YTITLE='Vx',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vx,YTITLE='Vx',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vy)
dmax=MAX(vy)
IF (flag EQ 0) THEN plot,x,vy,psym=6,symsize=.4,YTITLE='Vy',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vy,YTITLE='Vy',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vz)
dmax=MAX(vz)
IF (flag EQ 0) THEN plot,x,vz,psym=6,symsize=.4,YTITLE='Vz',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vz,YTITLE='Vz',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(by)
dmax=MAX(by)
IF (flag EQ 0) THEN plot,x,by,psym=6,symsize=.4,YTITLE='By',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,by,YTITLE='By',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(bz)
dmax=MAX(bz)
IF (flag EQ 0) THEN plot,x,bz,psym=6,symsize=.4,YTITLE='Bz',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,bz,YTITLE='Bz',XTITLE='X',YRANGE=[dmin,dmax]
phi = 180*(atan(bz/(by+1.0e-10))/3.1415927)
dmin=MIN(phi)
dmax=MAX(phi)
IF (flag EQ 0) THEN plot,x,phi,psym=6,symsize=.4,YTITLE='PHI',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,phi,YTITLE='PHI',XTITLE='X',YRANGE=[dmin,dmax]
!P.MULTI=0
END

;------------------------------------------------------------------------------
; Procedure SOD_PLOT: plots analytic solution for Sod shocktube over numerical
;
PRO sod_plot,filename
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
readbin,filename
vs = 1.7522
vc = 0.92745
vf = 0.07027
vh = 1.1832
xint = x[0] + (x[nx-1]-x[0])/2.0
xs = xint+vs*time
xc = xint+vc*time
xf = xint-vf*time
xh = xint-vh*time
dsod=FLTARR(500)
vsod=FLTARR(500)
esod=FLTARR(500)
psod=FLTARR(500)
xsod=FLTARR(500)
xsod[0] = x[0]
FOR I=1,499 DO xsod[I]=xsod[I-1] + (x[nx-1]-x[0])/500.
FOR I=0,499 DO BEGIN
  IF xsod[I] GT xs THEN BEGIN
    dsod[I] = 0.125
    psod[I] = 0.1
    vsod[I] = 0.0
    esod[I] = 2.0
  ENDIF
  IF xsod[I] GT xc AND xsod[I] LT xs THEN BEGIN
    dsod[I] = 0.26557
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 2.8535
  ENDIF
  IF xsod[I] GT xf AND xsod[I] LT xc THEN BEGIN
    dsod[I] = 0.42632
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 1.7776
  ENDIF
  IF xsod[I] GT xh AND xsod[I] LT xf THEN BEGIN
    vsod[I] = 0.92745*(xsod[I]-xh)/(xf-xh)
    dsod[I] = 0.42632*(1.0+0.20046*(0.92745-vsod[I]))^5
    psod[I] = 0.30313*(1.0+0.20046*(0.92745-vsod[I]))^7
    esod[I] = psod[I]/(0.4*dsod[I])
  ENDIF
  IF xsod[I] LT xh THEN BEGIN
    dsod[I] = 1.0
    psod[I] = 1.0
    vsod[I] = 0.0
    esod[I] = 2.5
  ENDIF
ENDFOR
;
!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,dsod
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,psod
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,vsod
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='P/D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,psod/dsod
!P.MULTI=0
END
;------------------------------------------------------------------------------
; Procedure FLINES:  2D plot of field lines
;
PRO flines,nlev
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
vecpot=fltarr(nx,ny)
dx = x/nx
dy = y/ny
vecpot[0,0] = 0.0
FOR J=1,ny-1 DO vecpot[0,J] = vecpot[0,j-1] + dy*bx[0,j]
FOR I=1,nx-1 DO vecpot[I,*] = vecpot[i-1,*] - dx*by[i,*]
contour,vecpot,nlevels=nlev,/ISOTROPIC
END
