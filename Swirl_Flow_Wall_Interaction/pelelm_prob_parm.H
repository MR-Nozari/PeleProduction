#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
   amrex::Real P_mean = 101325.0;
   amrex::Real splitx = 0.01;
   amrex::Real midtanh = 0.001;
   amrex::Real widthtanh = 0.001;
   amrex::Real Zst = 0.055;
   amrex::Real T_in = 300.0;
   amrex::Real U_b = 10;
   amrex::Real U_c = 0.5;
   amrex::Real x0 = 0.005;
   amrex::Real D = 0.0032;   
   amrex::Real H = 0.0001;   
      
   int bathID = -1;
   int fuelID = -1;
   int oxidID = -1;
};
#endif
