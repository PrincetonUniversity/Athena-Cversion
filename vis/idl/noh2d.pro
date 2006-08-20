; makes 1-D scatter plot of d vs r for spherical Noh shock test (noh.c)
;
PRO noh2d_plot,filename
COMMON SHARE1,nx,ny,nz,nvar
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

readbin,filename
 
denr=fltarr(nx*ny)
r=fltarr(nx*ny)

for i=0,nx-1 DO BEGIN
  for j=0,ny-1 DO BEGIN
    r[i*ny + j] = sqrt(float(i)*float(i) + float(j)*float(j))
    denr[i*ny + j] = d[i,j]
  ENDFOR
ENDFOR
r=r/nx

plot,r,denr,XTITLE='R',YTITLE='D',psym=3

END
