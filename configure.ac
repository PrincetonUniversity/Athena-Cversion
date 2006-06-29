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
#
# ALGORITHM "packages":
#   --with-order=[1,2,3]                       (order of spatial reconstruction)
#   --with-flux=[roe,hlle,hllc,hlld,force]                       (flux function)
#   --with-integrator=[ctu,vl]                (3D unsplit integration algorithm)
#
# ALGORITHM "features":
#   --enable-single                                 (double or single precision)
#   --enable-debug                                       (compile in debug mode)
#   --enable-ghost                      (write out ghost cells in outputs/dumps)
#   --enable-mpi                                          (parallelize with MPI)

#-------------------------------------------------------------------------------
# generic things

AC_INIT(src/athena.h)
AC_CONFIG_HEADER(src/config.h)
AC_PROG_CC
AC_PROG_INSTALL

#-------------------------------------------------------------------------------
# PACKAGE: name of problem generator from src/prob/name.c 
#   --with-problem=name (default is shkset)

AC_SUBST(PROBLEM)
AC_ARG_WITH(problem, [  --with-problem=PROB     Problem generator from src/prob/name.c],
		with_problem=$withval, with_problem="shkset")
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

rm -f src/lr_states.c
ln -s lr_states_prim$ORDER.c src/lr_states.c

#-------------------------------------------------------------------------------
# PACKAGE: flux function
#   --with-flux=[roe,hlle,hllc,hlld,force]
#   (default is roe if gas=mhd, hllc if gas=hydro)

AC_SUBST(FLUX_TYPE)
AC_ARG_WITH(flux, [  --with-flux=FLUX_TYPE   Flux function (roe, hllc (hydro only), hlld, hlle, force)],
		with_flux=$withval, with_flux=$flux_default)
FLUX_TYPE=$with_flux
if test   "$with_flux" = "roe"; then
  FLUX_TYPE="roe"
elif test "$with_flux" = "hlle"; then
  FLUX_TYPE="hlle"
elif test "$with_flux" = "hllc"; then
  FLUX_TYPE="hllc"
elif test "$with_flux" = "hlld"; then
  FLUX_TYPE="hlld"
elif test "$with_flux" = "force"; then
  FLUX_TYPE="force"
else
  AC_MSG_ERROR([Invalid flux function, valid types are: roe, hlle, hllc hlld, or force])
fi

#-------------------------------------------------------------------------------
# PACKAGE: set the 3D unsplit integration algorithm
#   --with-integrator=[ctu,vl] (default is ctu)

AC_SUBST(INTEGRATOR)
AC_SUBST(THREE_DIM_INT)

AC_ARG_WITH(integrator, [  --with-integrator=INTEGRATOR      Which unsplit 3D integrator (ctu, vl)],
		with_integrator=$withval, with_integrator="ctu")
INTEGRATOR=$with_integrator
if test   "$with_integrator" = "ctu"; then
  THREE_DIM_INT="CTU"
elif test "$with_integrator" = "vl"; then
  THREE_DIM_INT="VL"
else
  AC_MSG_ERROR([expected --with-integrator=ctu or vl])
fi

rm -f src/integrate_3d.c
ln -s integrate_3d_$INTEGRATOR.c src/integrate_3d.c

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
# PACKAGE: set MPI library path, --with-mpi-prefix='path'

AC_SUBST(MPI_MODE)
AC_SUBST(MPIINC)
AC_SUBST(MPILIB)

AC_ARG_ENABLE(mpi, [  --enable-mpi           enable MPI parellelization], ok=$enableval, ok=no)
if test "$ok" = "no"; then
  MPI_MODE="MPI_SERIAL"
  MPIINC=""
  MPILIB=""
else
  AC_MSG_CHECKING(whether we have MPI library)
  AC_ARG_WITH(mpi-prefix,
        [  --with-mpi-prefix=PFX   Prefix where MPI was installed (/usr/local)],
        mpi_prefix="$withval",
        mpi_prefix="/usr/local")

  mpi=$mpi_prefix
  for dir in $mpi/include 
    do
        if test -f $dir/mpi.h; then
                MPIINC="-I$dir"
                break;
        fi
    done
  for dir in $mpi/lib
    do
        if test -f $dir/libmpich.a; then
                MPILIB="-L$dir -lmpich"
                break;
        fi
    done
                                                                              
# MPI library was found, enable MPI
  if test "X$MPIINC" != "X"; then
    AC_DEFINE(HAVE_MPI)
    AC_MSG_RESULT([$MPIINC $MPILIB])
    MPI_MODE="MPI_PARALLEL"

# MPI library was not found, disable MPI
  else
    MPI_MODE="MPI_SERIAL"
    MPIINC=""
    MPILIB=""
    AC_MSG_RESULT(no)
    AC_MSG_WARN(disabling MPI support, or use --with-mpi-prefix=)
  fi
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
echo "Spatial Order:           $ORDER ($ACCURACY)"
echo "Flux:                    $FLUX_TYPE"
echo "3D unsplit integrator:   $THREE_DIM_INT"
echo "Precision:               $PRECISION"
echo "Compilation:"
echo "  OPT:                   $COMPILER_OPTS"
echo "Output modes:"
echo "  Ghost Cells:           $WRITE_GHOST_MODE"
echo "Parallel modes:"
echo "  MPI mode:              $MPI_MODE"
