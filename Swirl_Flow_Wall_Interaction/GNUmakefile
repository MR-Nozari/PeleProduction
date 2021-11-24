TOP = /scratch/b/bsavard/mrn7575/PeleProduction_26Oct_Feature_DirichletEB
PELELM_HOME ?= $(TOP)/Submodules/PeleLM
PELE_PHYSICS_HOME ?= $(TOP)/Submodules/PelePhysics
IAMR_HOME ?= $(TOP)/Submodules/IAMR
AMREX_HOME ?= $(TOP)/Submodules/amrex
AMREX_HYDRO_HOME ?= $(TOP)/Submodules/AMReX-Hydro

#
# Build configuration
#
DIM             = 3
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
USE_EB          = TRUE
TINY_PROFILE    = FALSE

#######################
USE_SUNDIALS_PP = TRUE
USE_KLU_PP = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_Model     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
#Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Constant

Blocs   := .

F90EXE_sources += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
