#pragma once
#include "Header.h"

class AMR_cell
{
public:
	double f = 0.0;
	double Spotok = 0.0; // Это поток в ячейке не умноженный на грань!

	AMR_cell* I_self;            // Указатель на себя

	uint8_t level = 0;
	AMR_cell* parent = nullptr;    // Ячейка - родитель
	uint8_t nx = 0;           // Номер данной ячейки в ячейке-родителе
	uint8_t ny = 0;
	uint8_t nz = 0;

	struct Flags {
		unsigned is_divided : 1;     // 1 бит  // Разделена ли ячейка
		unsigned is_signif : 1;      // 1 бит // Сущестывенная ячейка, та, которую можно делить, если надо 
	// определяется по процентру плотности к всему объёму
		unsigned need_devide_x : 1;  // 1 бит  // Нужно ли её делить вдоль x?
		unsigned need_devide_y : 1;  // 1 бит  // Нужно ли её делить вдоль y?
		unsigned need_devide_z : 1;  // 1 бит  // Нужно ли её делить вдоль z?
	} flags;  // Размер: 1 байт (вместо 5!)

	boost::multi_array<AMR_cell*, 3> cells;  // Ячейки - дети


	AMR_cell();

	double Get_SpotokV(void);
	void Get_Moment(AMR_f* AMR, double & m, double& mu, double& mux, double& muu);
	void Get_f(AMR_f* AMR, double& S);

	void divide(AMR_f* AMR, unsigned short int n1, unsigned short int n2, unsigned short int n3); // Разделить ячейку

	AMR_cell* find_cell(const double& x, const double& y, const double& z, const double& xL, 
		const double& xR, const double& yL, const double& yR, const double& zL, const double& zR);
	// Ищет ячейку по её соседям

	AMR_cell* get_sosed(AMR_f* AMR, short int nn);
	// nn = 0, 1, 2, 3, 4, 5
	//     по х вперёд - назад, по y ...

	void Print_info(void);

	void Get_random_velosity_in_cell(AMR_f* AMR, const double& ksi, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens);

	void Get_index(std::vector<std::array<unsigned int, 3>>& numbers);
	// Получить индекс ячейки

	void Get_Center(AMR_f* AMR, std::array<double, 3>& center); // Получить центр ячейки (даже если она разбита)
	void Get_Center(AMR_f* AMR, std::array<double, 3>& center, std::array<double, 3>& razmer); // Получить центр ячейки (даже если она разбита)

	void Get_Centers(AMR_f* AMR, std::vector<std::array<double, 3>>& centers); // Получить центры ячейки (включая центры подъечеек)
	void Get_all_cells(std::vector< AMR_cell*>& cells); // Получить список действительных ячеек (неразделённых)

	void Slice_plane(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d, std::vector < std::vector<std::array<double, 3>>>& points);
	// Разрезать ячейку плоскостью

	void Save_cell(std::ofstream& out);
	void Read_cell(std::ifstream& in);

	void Delete(void); 
	// Удаляет текущую ячейку и все её соседние
};

