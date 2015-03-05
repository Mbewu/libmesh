
#include "particles3d.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// boundary ids must have the inflow boundary id followed  by the outlfow 
// boundary ids in whatever order actually, must just match up
// these values are just magnitude values
void Particles3D::init_particles () 
{
	Point new_point(0.5,0.5,0.);
	Particle new_particle(*es,new_point)
	particles.push_back(new_particle);
}

// move particles over one timestep
void Particles3D::move_particles () 
{
	// it's as easy as this
	for(unsigned int i=0; i<particles.size(); i++)
		particles[i].move();

}


