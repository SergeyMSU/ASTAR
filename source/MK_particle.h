#pragma once
#include "Header.h"

class MK_particle
{
public:
	double Vel[3];    // Скорость частицы
	double coord[3];  // Положение частицы
	short int sort;         // Сорт частицы  1, 2, 3, 4
	double mu;              // Вес частицы
	double KSI;              // Для вычисления длины свободного провбега
	double I_do;             // Для розыгрыша длины свободного пробега
	Cell* cel;              // Ячейка, в которой находится частица

	MK_particle();

	void AddVel(const double& a, const double& b, const double& c);
	void AddVel(const Eigen::Vector3d& a);
	void Addcoord(const double& a, const double& b, const double& c);
	void Addcoord(const Eigen::Vector3d& a);

	void Move(const Eigen::Vector3d& a); // Передвигает центр на вектор a

	double Vel_norm(void);

};

