#include <pelelm_prob.H>

namespace ProbParm
{
AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
AMREX_GPU_DEVICE_MANAGED  amrex::Real T0 = 1500.;
AMREX_GPU_DEVICE_MANAGED  amrex::Real Y_O2 = 0.233;
AMREX_GPU_DEVICE_MANAGED  amrex::Real Y_N2 = 0.767;
AMREX_GPU_DEVICE_MANAGED  amrex::Real vel = 0.;
} // namespace ProbParm

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("ref_p", ProbParm::P_mean);
        pp.query("ref_T", ProbParm::T0);
        pp.query("init_vel", ProbParm::vel);
        pp.query("init_N2", ProbParm::Y_N2);
        pp.query("init_O2", ProbParm::Y_O2);
    }
}
