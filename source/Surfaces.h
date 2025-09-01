#pragma once
#include "Header.h"


class Surfaces // ‘айл считывани€ сетки из записанных файлов
// сделан дл€ старого формата (из программы MIK на фортране)
{
public:
	string name;
	map<string, int> parameters;

	vector<double> phi_angle;

	// ƒалее перва€ координата массивов - цилинтрический угол phi (так как вс€ сетка 
	// в старой программе получена цилиндрическим вращением

	// √оловна€ область - все поверхности радиальны
	boost::multi_array<double, 2> the_angle;
	boost::multi_array<double, 2> TS_radial;
	boost::multi_array<double, 2> HP_radial;
	boost::multi_array<double, 2> BS_radial;

	// TS в задней части - она радиальна
	boost::multi_array<double, 2> the_angle2; 
	boost::multi_array<double, 2> TS_radial2;

	// HP в задней части - она цилиндрическа€
	boost::multi_array<double, 2> x_cilindr;  // x < 0
	boost::multi_array<double, 2> HP_cilindr;

	~Surfaces();

	void Read_old(string nam);

	double Get_TS(const double& phi, const double& the);
	double Get_HP(const double& phi, const double& the, int met);
	double Get_BS(const double& phi, const double& the);

	// ¬изуализируем поверхность

	void Print_TS(); // печатаем точки TS
	void Print_HP(); // печатаем точки TS
	void Print_BS(); // печатаем точки TS
};

