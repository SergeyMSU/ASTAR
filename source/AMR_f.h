#pragma once
#include "Header.h"
#include <vector>

class AMR_f
{
public:
	double xL;
	double xR;
	double yL;
	double yR;
	double zL;
	double zR;

	double procent_signif = 0.3;  // 0.3;
	double procent_devide = 1.0;  // 1.0;

	array<double, 3> Vn;
	array<double, 3> Vt;
	array<double, 3> Vm;

	unsigned int xn;
	unsigned int yn;
	unsigned int zn;

	// Параметры для определения надо ли сгущать функцию
	// Интеграллы от функции распределения
	double Sf;
	double Sfu;
	double Sfux;
	double Sfuu;

	unordered_map<string, double> parameters;  // дополнительные МК-параметры

	double SpotokV = 0.0;  // Поток через данную грань функции распределения
	// этот поток уже умножен на площадь грани (чтоб удобно его сразу брать)

	AMR_f* AMR_self;
	mutex mut;

	boost::multi_array<AMR_cell*, 3> cells;

	void Culk_SpotokV(const double& Squ);
	// Считает поток функции распределения через границу
	// Также вычисляет поток в родительских (разделённых) ячейках

	void Get_random_velosity(AMR_f* AMR, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens);

	AMR_f();
	AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);

	void AMR_resize(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);

	void Get_real_koordinate(const double& x, const double& y, 
		const double& z, double& Vx, double& Vy, double& Vz);
	// По локальным координатам получает реальные глобальные координаты

	void Get_lokal_koordinate(const double& Vx, const double& Vy, 
		const double& Vz, double& x, double& y, double& z);

	void Add_particle(const double& Vx, const double& Vy,
		const double& Vz, const double& mu);

	void Normir_velocity_volume(const double& squ);
	// Нормировка в конце расчёта монте-карло

	void Set_bazis(void);
	// По нормали определяет два других вектора

	AMR_cell* find_cell(const double& x, const double& y, const double& z);
	// Ищет ячейку (указатель на неё) по координате

	void Get_all_cells(std::vector<AMR_cell*>& cells); 
	// Получить список действительных ячеек (т.е. если ячейка разделена, она не
	// включаются, а включаются её дети и т.д.).

	void Fill_maxwel_inf(const double& Vinf); 
	// Заполнить максевеллом на бесконечности
	// И провести процедуру мельчения

	double Integrate_Maxwell_V(const double& xL, const double& xR,
		const double& yL, const double& yR,
		const double& zL, const double& zR, const double& Vinf);

	void Fill_test(void);
	// Заполнить ячейки максвеллом

	void Fill_null(void);
	// Заполняем значищие ячейки нулём

	unsigned int Refine(void);
	unsigned int de_Refine(void);

	void Save(string namef);
	void Read(string namef);

	unsigned int Size(void);


	void Print_info(void);

	void Print_all_center_Tecplot(AMR_f* AMR, const string& name = "_");
	void Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d);
	// плоскость  a x + b y + c z + d = 0

	void Print_all_sosed_Tecplot(AMR_f* AMR);

	void Print_1D_Tecplot(AMR_f* AMR, const double& VV);

	void Delete(void);
	// Удаляет сетку, очищает память и т.д.
};

