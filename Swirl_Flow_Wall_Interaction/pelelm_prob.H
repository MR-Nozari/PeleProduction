#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <PelePhysics.H>
#include <EOS.H>
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 PmfData const *pmf_data)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real L_y = prob_hi[1] - prob_lo[1];


    amrex::Real massfrac_eq[NUM_SPECIES] = {0.0};

    amrex::Real y_lo = 0.0001;
    amrex::Real Gauss_maxT = 1800.0;
    amrex::Real Gauss_T_width = 0.0024;
    amrex::Real Gauss_Spec_width = 0.0004;


    amrex::Real T_sigma = 0.002;
    amrex::Real Spec_sigma = 0.0001;
    amrex::Real T_max_exponential = 2000;


  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};


  amrex::Real Ro_swirler = 0.0162;
  amrex::Real Ri_swirler = 0.0101;
  amrex::Real R = Ri_swirler/Ro_swirler;
  amrex::Real x0 = 0.0;
  amrex::Real y0 = 0.0;
  amrex::Real H = 0.00001;
  amrex::Real Do_swirler = 0.0324;
  amrex::Real Di_swirler = 0.0202;
  amrex::Real D_coflow = Do_swirler;
  amrex::Real Ro_coflow = prob_hi[0] ;
  amrex::Real Ri_coflow = Ro_swirler;
  amrex::Real c1 = -0.0067;
  amrex::Real c2 = 3.39e-7;
  amrex::Real c3 = -21.081;
  amrex::Real A = 0.0096;
  amrex::Real B = -3.869e-11;
  amrex::Real F_theta = -7.1884e5;


  massfrac_air[prob_parm.oxidID] = 0.233;  
  massfrac_air[prob_parm.bathID] = 0.767;

   
    
    state(i,j,k,DEF_Temp) = 300;

    for (int n = 0; n < NUM_SPECIES; ++n) {  
        massfrac[n] = massfrac_air[n];
    }     
       
    state(i,j,k,Xvel) = 0;
    state(i,j,k,Yvel) = 0;
    state(i,j,k,Zvel) = 0;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  PmfData const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};  
  amrex::Real massfrac[NUM_SPECIES] = {0.0}; 

  amrex::Real Ro_swirler = 0.0162;
  amrex::Real Ri_swirler = 0.0101;
  amrex::Real R = Ri_swirler/Ro_swirler;
  amrex::Real D_jet = 0.02;
  amrex::Real x0 = 0.0;
  amrex::Real y0 = 0.0;
  amrex::Real H = 0.00001;
  amrex::Real Do_swirler = 0.0324;
  amrex::Real Di_swirler = 0.0202;
  amrex::Real D_coflow = Do_swirler;
  amrex::Real Ro_coflow = prob_hi[0] ;
  amrex::Real Ri_coflow = Ro_swirler;
  amrex::Real c1 = -0.0067;
  amrex::Real c2 = 3.39e-7;
  amrex::Real c3 = -21.081;
  amrex::Real A = 0.0096;
  amrex::Real B = -3.869e-11;
  amrex::Real F_theta = -7.1884e5;

 
  massfrac_air[prob_parm.oxidID] = 0.233;
  massfrac_air[prob_parm.bathID] = 0.767;
 
    
    amrex::Real vel_swirler_x = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*((x[1]/(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))),2))/3)+(c1*(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))+((c2)/((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))))));

 

amrex::Real vel_swirler_y = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*(x[0]/(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))),2))/3)+(c1*(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))+((c2)/((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))));
 
 s_ext[Xvel] = -vel_swirler_x;
 s_ext[Yvel] = vel_swirler_y;

    amrex::Real vel_swirler_z = 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*prob_parm.U_b*(c3)*((std::pow(((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))/(Ro_swirler)),2))-(((std::pow(R,2)-1)/(std::log(R)))*std::log((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))/Ro_swirler))-1);

   amrex::Real vel_coflow_z = prob_parm.U_c * 0.5 * (1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_coflow/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_coflow/2))/(H))));

amrex::Real vel_jet = prob_parm.U_b*( 1 - (0.5*(1+(std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_jet/2))/(H)))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_jet/2))/(H))))));

s_ext[Zvel] = vel_swirler_z + vel_coflow_z + vel_jet;
   

    s_ext[DEF_Temp] = prob_parm.T_in;
      


    for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = massfrac_air[n];

     }

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    
    eos.PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
    s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
    }
  

}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(
  const amrex::Real xEBface[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm)
{

  // Setup a fixed EB temperature

  s_ext[DEF_Temp] = 300.0;
                        
 
 
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];
}
#endif
