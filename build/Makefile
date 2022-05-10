#!/usr/bin/env make

RTE_DIR    = ../rte
RRTMGP_DIR = ../rrtmgp
RTE_KERNEL_INTERFACE_DIR    = $(RTE_DIR)/kernel-interfaces
RRTMGP_KERNEL_INTERFACE_DIR = $(RRTMGP_DIR)/kernel-interfaces
#
# Choose kernels depending on environment variable
#
ifeq ($(RTE_KERNELS), openacc)
  RTE_KERNEL_DIR    = ../rte/kernels-openacc
  RRTMGP_KERNEL_DIR = ../rrtmgp/kernels-openacc
endif
#
# For some kernel file there is no separate OpenACC implementation
#
RTE_KERNEL_DIR    += ../rte/kernels-fortran
RRTMGP_KERNEL_DIR += ../rrtmgp/kernels-fortran
VPATH  = $(RTE_DIR):$(RTE_KERNEL_INTERFACE_DIR):$(RTE_KERNEL_DIR)
VPATH += $(RRTMGP_DIR):$(RRTMGP_KERNEL_INTERFACE_DIR):$(RRTMGP_KERNEL_DIR)
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# Make all the libraries though we'll only use the interface + kernels
all:           librte.a                           librrtmgp.a
separate-libs: librtekernels.a librteinterfacef.a librrtmgpkernels.a librrtmgpinterfacef.a

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

librte.a: $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)
	ar -rvs librte.a $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)

librtekernels.a: $(RTE_FORTRAN_KERNELS)
	ar -rvs librtekernels.a $(RTE_FORTRAN_KERNELS)

librteinterfacef.a: $(RTE_FORTRAN_INTERFACE)
	ar -rvs librteinterfacef.a $(RTE_FORTRAN_INTERFACE)

librrtmgp.a: $(RRTMGP_FORTRAN_KERNELS)  $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgp.a $(RRTMGP_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_INTERFACE)

librrtmgpkernels.a: $(RRTMGP_FORTRAN_KERNELS)
	ar -rvs librtekernels.a $(RRTMGP_FORTRAN_KERNELS)

librrtmgpinterfacef.a: $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librteinterfacef.a $(RRTMGP_FORTRAN_INTERFACE)

clean:
	rm -f *.optrpt *.mod *.smod *.o lib*.a
