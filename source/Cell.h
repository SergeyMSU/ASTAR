#pragma once
#include"Header.h"

// Ячейкам задаются зоны в функции New_append_surfaces
enum class Type_cell {
	none,   // 0
	Zone_1,  // 1   
	Zone_2,  // 2  
	Zone_3,  // 3   
	Zone_4  // 4   
};


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	//vector<Edge*> edges;          // Вектор рёбер нужен только для правильного вычисления ротеров
	

	int number = 0;                // номера начинаются с единицы
	bool is_inner = false;         // Является ли ячейка внутренней (которые считаются отдельно)
	bool is_TVD = true;         // Можно ли делать TVD в этой ячейке
	// Фактически равно false только для первых ячеек  вблизи нуля

	mutex mut;
	short int is_need = 0;          // параметр, который можно исползовать для разных нужд

	Type_cell type = Type_cell::none;  // по умолчанию создаём обычный узел

	short unsigned int MK_zone_r = 0;
	short unsigned int MK_zone_phi = 0;
	short unsigned int MK_zone = 0;

	std::array < unordered_map<string, double>, 2> parameters;   // Параметры в ячейке (тут могут быть значения
	// плазменных полей, значения дивергенций и т.д.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - маячок для переноса и определения HP

	unordered_map<string, double> geo_parameters;  // Геометрические параметры для 
	// Удобства и ускорения расчёта Монте-Карло
	// l_size  характерный размер ячейки

	unordered_map<string, Eigen::VectorXd> interpolate_alpha;

	// Источники рождения пикапов
	vector<double> pui_Sm;
	vector<double> pui_Sp;


	void Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par);
	void Get_IDW_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par, Phys_param* phys_param);

	Cell* Get_Sosed(Gran* gr); // Получить соседа через данную грань
	// возвращает nullptr, если соседа через данную грань нет, т.е. грань граничная
	// Для работы функции номера ячеек должны быть актуальны

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	// Как предыдущая, но сравнение идёт по указателю, не требует номеров ячеек

	bool Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next);
	// Принадлежит ли точка ячейке

	void Set_Cell_Geo_for_MK(void);


	double center[2][3];           // Центр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];


	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // Расстояние от центра до начала координат


	void Tecplot_print_cell(void);

	void MK_Add_particle(MK_particle& P, const double& time);
	void MK_Add_pui_source(const double& wr, const double& nu_ex, const double& mu,
		const double& time, Phys_param* phys_param);

	void MK_calc_Sm(Phys_param* phys_param);
	void MK_normir_Moments(Phys_param* phys_param);

};

