---
title: The Configure Step
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Configuring

After installing Athena, the next step is to configure the
code for a specific problem and algorithm.  First, the `configure` shell
script must be created using the [GNU autoconf](http://www.gnu.org/software/autoconf) toolkit:

	% cd athena
	% autoconf


The `configure` script is used to
enable or disable *features* in the code, and to
choose between different optional *packages*.  An example of an algorithmic
feature is static mesh refinement (SMR) -- it is either enabled or
disabled.  An example of an algorithmic package is the choice of Riemann solver
-- there are many solvers avialable, and which one is used is controlled by configure.
At the source code level, features and packages are controlled using C precompiler macros.
The `configure` script provides
a powerful way of setting these macros through a command line interface
that does not require the user to edit any special files.

The `configure` script is also used to set compiler
and linker options and flags using environment variables, and to query the
system automatically to see that a C-compiler, linker, and any external
libraries necessary for compilation are installed and accessible.

The `configure` script uses the following syntax

	configure [--with-package=choice] [--enable-feature]

Valid options for the *features* and *packages* implemented in Athena are given in the
following tables.

### Optional physics (package=*choice*) controlled by configure

| Package | Choice | Comments |
|---------|--------|----------|
| problem | file-name | use file-name in directory /src/prob for initial conditions |
| gas    | hydro  | create code for hydrodynamics |
|       | mhd    | **(default)** create code for MHD |
| eos    | adiabatic | **(default)** use adiabatic equation of state |
|       | isothermal | use isothermal equation of state |
| nscalars | # | add # passively advected scalars **(default is 0)** |
| gravity | none      | **(default)** no self-gravity |
|         | fft       | enable self-gravity using FFTs |
|         | multigrid | enable self-gravity using multigrid |
| particles    | none | **(default)** no particles |
|              | feedback | include momentum feedback on gas |
|              | passive  | passive Lagrangian particles |
| coord        | cartesian | **(default)** |
|              | cylindrical | |


### Optional physics *features* controlled by configure

| Feature     | Default | Comments |
|-------------|---------|----------|
| special-relativity | disabled | both relativistic hydro and MHD are implemented |
| conduction         | disabled | both isotropic and anisotropic thermal conduction implemented |
| resistivity        | disabled | Ohmic dissipation, the Hall effect, and ambipolar diffusion all implemented |
| viscosity          | disabled | both isotropic (Navier-Stokes) and anisotropic (Braginskii) viscosity implemented |


### Optional algorithms (package=*choice*) controlled by configure

| Package | Choice | Comments |
|---------|--------|----------|
| flux | roe | **(default)** Roe's linearized Riemann solver |
|      | exact | exact nonlinear Riemann solver (hydro only) |
|      | force | FORCE flux |
|      | hlle | HLLE Riemann solver |
|      | hllc | HLLC Riemann solver (hydro only) |
|      | hlld | HLLD Riemann solver (MHD only) |
|      | two-shock | two-shock approximate solver (hydro only) |
| order | 1 | first-order reconstruction |
|     | 2 | **(default)** second-order reconstruction with limiting in the characteristic variables |
|      | 3 | third-order reconstruction with limiting in the characteristic variables |
|      | 2p | second-order reconstruction with limiting in the primitive variables |
|      | 3p | third-order reconstruction with limiting in the primitive variables |
| integrator | ctu | **default)** corner transport upwind (CTU) unsplit integrator in 3D |
|           | vl | van Leer unsplit integrator in 3D |
| cflags  | opt | **(default)** optimization flags |
|         | debug | add gdb debugger flags |
|         | profile | add gprof flags |


### Optional algorithm *features* controlled by configure

| Feature     | Default | Comments |
|-------------|---------|----------|
| fargo        | disabled | add orbital advection for shearing-box |
| fft         | disabled | compile and link with FFTW block decomposition |
| fofc         | disabled | add first-order flux-correction to VL integrator |
| ghost       | disabled | causes ghost zones to be written during output |
| h-correction | disabled | H-correction to fix carbuncle problem |
| hllallwave  | disabled | add "all-wave" interpolation with HLL Riemann solvers |
| mpi         | disabled | parallelization using MPI library |
| shearing-box | disabled | add source terms for shearing box |
| single      | disabled | computations performed in single precision (default is double) | 
| smr          | disabled | use static mesh refinement |

The `configure` script should be run in the root directory
`./athena3.1`.  Running `configure --help`
gives more information, including a list of all optional features
and packages.

## Examples

To configure Athena to run an isothermal hydrodynamical shocktube
problem initialized with the function `shkset1d.c` using Roe fluxes and third-order interpolation, use:

	configure --with-problem=shkset1d --with-eos=isothermal --with-gas=hydro --with-order=3


To run the linear wave test problem in 3D adiabatic
MHD using HLLD fluxes, the van Leer integrator, and second-order interpolation in double precision 
parallelized with MPI, use the following:

	configure --with-flux=hlld --with-problem=linear_wave3d --with-integrator=vl --enable-mpi


When `configure` runs, it creates a custom versions of `Makefile`, `Makeptions`, and `defs.h` for the
problem using the `Makefile.in`, `Makeoptions.in`, and `defs.h.in` files
as templates.  After successful execution,
`configure` will echo the options that have been set (including all
the default values).

Note the some configure options are mutually exclusive; for example special relativity cannot be run with the CTU
integrator.  **The precompiler should (but may not always) generate an error in such cases.**

