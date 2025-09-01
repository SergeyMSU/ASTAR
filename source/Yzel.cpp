#include "Yzel.h"


Yzel::Yzel()
{
	this->coord[0][0] = this->coord[0][1] = this->coord[0][2] = 0.0;
	this->coord[1][0] = this->coord[1][1] = this->coord[1][2] = 0.0;
	//this->luch = nullptr;
	this->number = 0;
	this->is_inner = false;
	this->type = Type_yzel::Us;
}

Yzel::Yzel(const double& a, const double& b, const double& c)
{
	this->coord[0][0] = a;
	this->coord[0][1] = b;
	this->coord[0][2] = c;

	this->coord[1][0] = a;
	this->coord[1][1] = b;
	this->coord[1][2] = c;

	///this->luch = nullptr;
}

double Yzel::func_R(int i_time)
{
	return sqrt(kvv(this->coord[i_time][0], this->coord[i_time][1], this->coord[i_time][2]));
}

double Yzel::func_Ryz(int i_time)
{
	return sqrt(kvv(0.0, this->coord[i_time][1], this->coord[i_time][2]));
}
