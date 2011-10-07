#!/usr/bin/env python
# Plots data in binary dumps.
# Usage: python plot.py nl-shwave.0024.bin

from numpy import *
import sys
from matplotlib import use
use('Agg')
from matplotlib import pyplot

#
# Problem parameters
#
qshear = 1.5
Omega = 1.
Ly = 4.
tswing = 16./3. # Time at which the wave swings from leading to trailing

#
# Read binary file
#
try:
  file = open(sys.argv[1],'rb')
except:
  print 'Usage: ./read.py <binary_dump>'
  raise SystemExit

file.seek(0,2)
eof = file.tell()
file.seek(0,0)

coordsys = fromfile(file,dtype=int32,count=1)[0]

nx,ny,nz = fromfile(file,dtype=int32,count=7)[:3]

gamma1,cs = fromfile(file,dtype=float32,count=2)

t,dt = fromfile(file,dtype=float32,count=2)

x = fromfile(file,dtype=float32,count=nx)
y = fromfile(file,dtype=float32,count=ny)
z = fromfile(file,dtype=float32,count=nz)

shape = (ny,nx)
count = prod(shape)

rho = fromfile(file,dtype=float32,count=count).reshape(shape)
rux = fromfile(file,dtype=float32,count=count).reshape(shape)
ruy = fromfile(file,dtype=float32,count=count).reshape(shape)
ruz = fromfile(file,dtype=float32,count=count).reshape(shape)

if file.tell() != eof: print 'Error: Too few bytes read.'

file.close()

#
# Transform to a shearing coordinate frame in which the radial wave number
# is zero at all times
#
# This is done most easily by Fourier transforming each variable along y,
# then multiplying by an exponential phase factor of the form
#   -q*Omega*(t-tswing)*x,
# which effectively reverts the tilting of the wave fronts due to advection by
# the linear shear, and finally Fourier transforming back along y.
#
ky = 2*pi*arange(ny/2+1)/Ly
phase = -qshear*Omega*(t-tswing)*outer(ky,x)

rho_s = fft.irfft(exp(1j*phase)*fft.rfft(rho,axis=0),axis=0)
rux_s = fft.irfft(exp(1j*phase)*fft.rfft(rux,axis=0),axis=0)
ruy_s = fft.irfft(exp(1j*phase)*fft.rfft(ruy,axis=0),axis=0)
ruz_s = fft.irfft(exp(1j*phase)*fft.rfft(ruz,axis=0),axis=0)

#
# Plot density field
#
fig1,ax1 = pyplot.subplots(1,2,num=1)

ax1[0].imshow(rho,origin='lower')
ax1[1].imshow(rho_s,origin='lower')

fig1.savefig('fig1.png',dpi=fig1.get_dpi())

fig2,ax2 = pyplot.subplots(1,num=2)

ax2.plot(rho_s.mean(axis=1))

fig2.savefig('fig2.png',dpi=fig2.get_dpi())
