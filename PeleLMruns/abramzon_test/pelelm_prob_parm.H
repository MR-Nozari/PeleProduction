#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real Y_N2;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real Y_O2;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real vel;

}

#endif
