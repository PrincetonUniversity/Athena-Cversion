# Top level makefile for Athena.  The first target (help) documents the rest.
DIRS = bin

help:
	@echo "all:       create ($(DIRS)) subdirectory and compile"
	@echo "dirs:      create ($(DIRS)) subdirectory"
	@echo "compile:	  compile the code"
	@echo "clean:     clean /src subdirectory"
	@echo "test:      run a MHD benchmark"

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
# compared to fiducial solution, and a speed benchmark.  Requires the
# following steps: > configure
#                  > make all
#                  > make test
test:
	(cd tst/1D-mhd; ./run.test)
