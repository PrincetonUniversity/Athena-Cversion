# Top level makefile for Athena.  The first target (help) documents the rest.
DIRS = bin

help:
	@echo "all:       create ($(DIRS)) subdirectory and compile"
	@echo "dirs:      create ($(DIRS)) subdirectory"
	@echo "compile:	  compile the code"
	@echo "clean:     clean /src subdirectory"
	@echo "test:      run a MHD benchmark"
	@echo "test-all:  configure, compile, and run a test suite"

#-------------------------------------------------------------------------------
#  target all:
all:    dirs compile

#-------------------------------------------------------------------------------
#  target dirs:
dirs:
	-@for i in $(DIRS) ; do \
	(if [ -d $$i ]; \
	then \
	    echo DIR $$i exists; \
	else \
	    echo DIR $$i created; \
	    mkdir $$i; \
	fi); done

#-------------------------------------------------------------------------------
#  target compile:
compile:
	(cd src; $(MAKE) compile)

#-------------------------------------------------------------------------------
#  target clean:
clean:
	(cd src; $(MAKE) clean)

#-------------------------------------------------------------------------------
# test: checks that default configuration runs successfully.  Reports error
# compared to fiducial solution, and a nice speed benchmark.  Requires the
# following steps: > configure
#                  > make all
#                  > make test
test:
	(cd tst/1D-mhd; ./run.brio+wu)

#-------------------------------------------------------------------------------
# test-all: Runs the battery of 1D and 2D tests in the tst/test-all script.
# This script automatically configures, compiles, and runs code for each test.
# So, only step necessary is: > make test-all
test-all:
	(cd tst; ./test-all)
