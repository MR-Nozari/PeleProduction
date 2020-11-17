
#include <AMReX_Particles.H>
#include "SprayParticles.H"
//#include <PeleLM.H>
//#include "pelelm_prob.H"

using namespace amrex;

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  return false;
} 

void
SprayParticleContainer::InitSprayParticles()
{
}
