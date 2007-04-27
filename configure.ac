#-------------------------------------------------------------------------------
# configure.ac: Athena configuration script template file.  Run 'autoconf'
# to generate a new 'configure' script from this file.
#
# If new configure options are added to this file, they should also be
# added in show_config() in src/show_config.c, and to the diagnostic message
# output at the end of this file.
#
# This template adds the following options to the 'configure' script:
#
# PHYSICS "packages":
#   --with-problem='name'                 (use prob/name.c as problem generator)
#   --with-gas=[hydro,mhd]                              (hydro or mhd algorithm)
#   --with-eos=[isothermal, adiabatic]                       (equation of state)
#   --with-nscalars=n                               (number of advected scalars)
#
# ALGORITHM "packages":
#   --with-order=[1,2,3]                       (order of spatial reconstruction)
#   --with-flux=[roe,hlle,hllc,hlld,force,exact]                 (flux function)
#   --with-integrator=[ctu,vl]                (3D unsplit integration algorithm)
#
# ALGORITHM "features":
#   --enable-single                                 (double or single precision)
#   --enable-debug                                       (compile in debug mode)
#   --enable-ghost                      (write out ghost cells in outputs/dumps)
#   --enable-mpi                                          (parallelize with MPI)
#   --enable-h-correction              (turn on H-correction in multidimensions)

#-------------------------------------------------------------------------------
# generic things

AC_INIT(src/athena.h)
AC_CONFIG_HEADER(src/config.h)
AC_PROG_CC
AC_PROG_INSTALL

#-------------------------------------------------------------------------------
# PACKAGE: name of problem generator from src/prob/name.c 
#   --with-problem=name (default is linear_wave1d)

AC_SUBST(PROBLEM)
AC_ARG_WITH(problem, [  --with-problem=PROB     Problem generator from src/prob/name.c],
		with_problem=$withval, with_problem="linear_wave1d")
PROBLEM=$with_problem
if test -e src/prob/$PROBLEM.c; then
   rm -f src/problem.c
   ln -s prob/$PROBLEM.c src/problem.c
else
  ls src/prob
  AC_MSG_ERROR([Invalid problem filename, valid names listed above])
fi

#-------------------------------------------------------------------------------
# PACKAGE: gas mode (hydro or mhd)
#  --with-gas=hydro (default is mhd)
#   mode also sets the default flux handling, since that depends on whether the
#   problem is hydro or mhd (see option --with-flux= below)

AC_SUBST(GAS)
AC_ARG_WITH(gas, [  --with-gas=GAS          Gas properties (mhd,hydro)],
		with_gas=$withval, with_gas="mhd")
if test   "$with_gas" = "hydro"; then
  GAS="HYDRO"
  flux_default="hllc"
elif test "$with_gas" = "mhd";   then
  GAS="MHD"
  flux_default="roe"
else
  AC_MSG_ERROR([expected --with-gas=mhd or hydro])
fi

#-------------------------------------------------------------------------------
# PACKAGE: equation of state
#   --with-eos=[adiabatic,isothermal] (default is adiabatic)

AC_SUBST(EOS)
AC_ARG_WITH(eos, [  --with-eos=EOS          Equation of state (adiabatic,isothermal)],
		with_eos=$withval, with_eos=adiabatic)
if test   "$with_eos" = "adiabatic";  then
  EOS="ADIABATIC"
elif test "$with_eos" = "isothermal"; then
  EOS="ISOTHERMAL"
else
  AC_MSG_ERROR([expected --with-eos=adiabatic or isothermal])
fi

#-------------------------------------------------------------------------------
# PACKAGE: nscalars
#  --with-nscalars=n (n is any integer, default is 0)
#   number of passively advected scalars

AC_SUBST(NSCALARS)
AC_ARG_WITH(nscalars, [  --with-nscalars=n (default is 0)],
                nscalars=$withval, nscalars=0)
if test   "$nscalars" = "0"; then
  NSCALARS="0"
elif test "$nscalars" != "0"; then
  NSCALARS=$nscalars
else
  AC_MSG_ERROR([expected --with-nscalars=n])
fi

#-------------------------------------------------------------------------------
# PACKAGE: set the order of the spatial reconstruction
#   --with-order=[1,2,3] (default is 2)

AC_SUBST(ORDER)
AC_SUBST(ACCURACY)

AC_ARG_WITH(order, [  --with-order=ORDER      Order of Accuracy (1,2,3)],
		with_order=$withval, with_order="2")
ORDER=$with_order
if test   "$with_order" = "1"; then
  ACCURACY="FIRST_ORDER"
elif test "$with_order" = "2"; then
  ACCURACY="SECOND_ORDER"
elif test "$with_order" = "3"; then
  ACCURACY="THIRD_ORDER"
else
  AC_MSG_ERROR([expected --with-order=1, 2 or 3])
fi

#-------------------------------------------------------------------------------
# PACKAGE: flux function
#   --with-flux=[roe,hlle,hllc,hlld,force,exact]
#   (default is roe if gas=mhd, hllc if gas=hydro)

AC_SUBST(FLUX_NAME)
AC_SUBST(FLUX_DEF)
AC_ARG_WITH(flux, [  --with-flux=FLUX_NAME   Flux function (roe, hllc (hydro only), hlld, hlle, force, exact (hydro only))],
		with_flux=$withval, with_flux=$flux_default)
FLUX_NAME=$with_flux
if test   "$with_flux" = "roe"; then
  FLUX_DEF="ROE_FLUX"
elif test "$with_flux" = "hlle"; then
  FLUX_DEF="HLLE_FLUX"
elif test "$with_flux" = "hllc"; then
  FLUX_DEF="HLLC_FLUX"
elif test "$with_flux" = "hlld"; then
  FLUX_DEF="HLLD_FLUX"
elif test "$with_flux" = "force"; then
  FLUX_DEF="FORCE_FLUX"
elif test "$with_flux" = "exact"; then
  FLUX_DEF="EXACT_FLUX"
else
  AC_MSG_ERROR([Invalid flux function, valid types are: roe, hlle, hllc hlld, force, exact])
fi

#-------------------------------------------------------------------------------
# PACKAGE: set the 3D unsplit integration algorithm
#   --with-integrator=[ctu,vl] (default is ctu)

AC_SUBST(INTEGRATOR)
AC_SUBST(THREED_INT_DEF)

AC_ARG_WITH(integrator, [  --with-integrator=INTEGRATOR      Which unsplit 3D integrator (ctu, vl)],
		with_integrator=$withval, with_integrator="ctu")
INTEGRATOR=$with_integrator
if test   "$with_integrator" = "ctu"; then
  THREED_INT_DEF="THREED_CTU"
elif test "$with_integrator" = "vl"; then
  THREED_INT_DEF="THREED_VL"
else
  AC_MSG_ERROR([expected --with-integrator=ctu or vl])
fi

#-------------------------------------------------------------------------------
# FEATURE: precision of floating point arithmetic
#   --enable-single (default is double)

AC_SUBST(PRECISION)
AC_ARG_ENABLE(single, [  --enable-single         single-precision (default is double)], 
		      ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  PRECISION="SINGLE_PREC"
else
  PRECISION="DOUBLE_PREC"
fi

#-------------------------------------------------------------------------------
# FEATURE: compile in debugging mode; turns out diagnostic outputs in code,
#   and uses compiler options: -g -Wall -W -ansi -pedantic
#   --enable-debug (default is debugging off)

AC_SUBST(DEBUG_MODE)
AC_SUBST(COMPILER_OPTS)
AC_ARG_ENABLE(debug, [  --enable-debug          debug (default is optimized)], 
		      ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  COMPILER_OPTS="-g -Wall -W -ansi -pedantic"
  DEBUG_MODE="DEBUG"
else
  COMPILER_OPTS="-O3"
  DEBUG_MODE="OPTIMIZE"
fi

#-------------------------------------------------------------------------------
# FEATURE: write ghost cells in outputs/dumps
#   --enable-ghost

AC_SUBST(WRITE_GHOST_MODE)
AC_ARG_ENABLE(ghost, [  --enable-ghost          write ghost zones], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  WRITE_GHOST_MODE=WRITE_GHOST_CELLS 
else
  WRITE_GHOST_MODE=NO_WRITE_GHOST_CELLS
fi

#-------------------------------------------------------------------------------
# FEATURE: parallelize with MPI, --enable-mpi (default is no MPI)

AC_SUBST(MPI_MODE)

AC_ARG_ENABLE(mpi, [  --enable-mpi           enable MPI parellelization], ok=$enableval, ok=no)
if test "$ok" = "no"; then
  MPI_MODE="MPI_SERIAL"
else
  MPI_MODE="MPI_PARALLEL"
fi

#-------------------------------------------------------------------------------
# FEATURE: turn on H-correction in multidimensional integrators
#   --enable-h-correction
# Note that H-correction only works with Roe flux.  Prints error message if
# H-correction is enabled with any other flux.

AC_SUBST(H_CORRECTION_MODE)
AC_ARG_ENABLE(h-correction, [  --enable-h-correction          turn on H-correction], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
if test "$with_flux" = "roe"; then
  H_CORRECTION_MODE=H_CORRECTION
else
  AC_MSG_ERROR([H-correction only works with Roe flux])
fi
else
  H_CORRECTION_MODE=NO_H_CORRECTION
fi

#-------------------------------------------------------------------------------
# check for various library functions

AC_PROG_GCC_TRADITIONAL
AC_CHECK_FUNCS(strdup)
AC_C_BIGENDIAN

#-------------------------------------------------------------------------------
# date and time of configure, becomes a macro in defs.h

A_CONFIGURE_DATE="`date`"
AC_SUBST(A_CONFIGURE_DATE)

#-------------------------------------------------------------------------------
# write final diagnostic output

AC_SUBST(WARNING1)
AC_SUBST(WARNING2)
WARNING1="WARNING! This file has been automatically generated by configure."
WARNING2="Any changes to it will be overwritten the next time configure is run."

AC_CONFIG_FILES(src/Makefile src/defs.h)
AC_OUTPUT

echo ""
echo "Your athena distribution has now been configured:"
echo ""
echo "Problem:                 $PROBLEM"
echo "Gas properties:          $GAS"
echo "Equation of State:       $EOS"
echo "Advected scalar fields:  $NSCALARS"
echo "Spatial Order:           $ORDER ($ACCURACY)"
echo "Flux:                    $FLUX_NAME"
echo "3D unsplit integrator:   $INTEGRATOR"
echo "Precision:               $PRECISION"
echo "Compilation:"
echo "  OPT:                   $COMPILER_OPTS"
echo "Output modes:"
echo "  Ghost Cells:           $WRITE_GHOST_MODE"
echo "Parallel modes:"
echo "  MPI mode:              $MPI_MODE"
echo "H-correction:            $H_CORRECTION_MODE"
