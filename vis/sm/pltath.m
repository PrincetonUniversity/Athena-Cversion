# FILE pltath.m
#
# PURPOSE: contains supermongo (sm) macros to read and make plots of Athena 
#   .tab dumps.
#-------------------------------------------------------------------------------
# MACRO adb_hyd_read: Reads adiabatic hydro .tab dumps
adb_hyd_read	1	#
			data $1
			lines 7 100000
			read {x 2 d 3 M1 4 M2 5 M3 6 p 7 e 8}
			set vx=M1/d
			set vy=M2/d
			set vz=M3/d
			set ke=0.5*(M1*M1 + M2*M2 + M3*M3)/d

#-------------------------------------------------------------------------------
# MACRO adb_mhd_read: Reads adiabatic hydro .tab dumps
adb_mhd_read	1	#
			data $1
			lines 7 100000
			read {x 2 d 3 M1 4 M2 5 M3 6 p 7 e 8 b1c 9 b2c 10 b3c 11}
			set vx=M1/d
			set vy=M2/d
			set vz=M3/d
			set ke=0.5*(M1*M1 + M2*M2 + M3*M3)/d
			set me=0.5*(b1c*b1c + b2c*b2c + b3c*b3c)

#-------------------------------------------------------------------------------
# MACRO four_plot: Plots d,P,P/d,Vx
four_plot	0	#
			erase
			window 2 2 1 2
			limits x d
			xlabel X
			ylabel D
			box
			connect x d
			window 2 2 2 2
			limits x p
			xlabel X
			ylabel P
			box
			connect x p
			window 2 2 1 1
			limits x vx
			xlabel X
			ylabel Vx
			box
			connect x vx
			window 2 2 2 1
			limits x (p/d)
			xlabel X
			ylabel (P/D)
			box
			connect x (p/d)

#-------------------------------------------------------------------------------
# MACRO nine_plot: Plots d,P,P/d,Vx,Vy,Vz,By,Bz,atan(Bz/By)
nine_plot	0	#
			erase
			window 3 3 1 3
			limits x d
			xlabel X
			ylabel D
			box
			connect x d
			window 3 3 2 3
			limits x p
			xlabel X
			ylabel P
			box
			connect x p
			window 3 3 3 3
			limits x (p/d)
			xlabel X
			ylabel (P/D)
			box
			connect x (p/d)
			window 3 3 1 2
			limits x vx
			xlabel X
			ylabel Vx
			box
			connect x vx
			window 3 3 2 2
			limits x vy
			xlabel X
			ylabel Vy
			box
			connect x vy
			window 3 3 3 2
			limits x vz
			xlabel X
			ylabel Vz
			box
			connect x vz
			window 3 3 1 1
			limits x b2c
			xlabel X
			ylabel B2
			box
			connect x b2c
			window 3 3 2 1
			limits x b3c
			xlabel X
			ylabel B3
			box
			connect x b3c
			window 3 3 3 1
			set phi=atan(b3c/b2c)
			limits x phi
			xlabel X
			ylabel PHI
			box
			connect x phi
