#include "MK_particle.h"

MK_particle::MK_particle()
{
	this->coord[0] = 0.0;
	this->coord[1] = 0.0;
	this->coord[2] = 0.0;

	this->Vel[0] = 0.0;
	this->Vel[1] = 0.0;
	this->Vel[2] = 0.0;

	this->mu = 0.0;
	this->KSI = 0.0;
	this->I_do = 0.0;
	this->sort = 0;
	this->cel = nullptr;
}

void MK_particle::AddVel(const double& a, const double& b, const double& c)
{
	this->Vel[0] = a;
	this->Vel[1] = b;
	this->Vel[2] = c;
}

void MK_particle::AddVel(const Eigen::Vector3d& a)
{
	this->Vel[0] = a[0];
	this->Vel[1] = a[1];
	this->Vel[2] = a[2];
}


void MK_particle::Addcoord(const double& a, const double& b, const double& c)
{
	this->coord[0] = a;
	this->coord[1] = b;
	this->coord[2] = c;
}

void MK_particle::Addcoord(const Eigen::Vector3d& a)
{
	this->coord[0] = a[0];
	this->coord[1] = a[1];
	this->coord[2] = a[2];
}

void MK_particle::Move(const Eigen::Vector3d& a)
{
	this->coord[0] += a[0];
	this->coord[1] += a[1];
	this->coord[2] += a[2];
}

double MK_particle::Vel_norm(void)
{
	return norm2(this->Vel[0], this->Vel[1], this->Vel[2]);
}
