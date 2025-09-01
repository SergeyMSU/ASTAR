#pragma once
#include "Header.h"

// Следующая классификация для постановки граничных условий
enum class Type_Gran {
	Us,    // 0    обычная грань
	Inner_Hard,      // 1       граничная грань на внутренней сфере с жёсткими граничными условиями
	Outer_Hard,      // 2       граничная грань на внешней границе с жёсткими граничными условиями
	Outer_Soft      // 3       граничная грань на внешней границе с мягкими граничными условиями
};

enum class Type_Gran_surf {
	Us,    // 0    обычная грань
	TS,      // 1       граничная грань на внутренней сфере с жёсткими граничными условиями
	HP,      // 2       граничная грань на внешней границе с жёсткими граничными условиями
	BS      // 3       граничная грань на внешней границе с мягкими граничными условиями
};

class Gran
{
public:
	vector<Yzel*> yzels;  // Узлы грани (узлы должны быть расположены по кругу)
	vector<Cell*> cells;  // две ячейки грани (нормаль должна смотреть от первой ко второй
	// Расположение узлов по кругу необходимо для того, чтобы любые выбранные подряд две точки 
	// образовывали ребро грани
	vector<Cell*> cells_TVD;  // Ячейки для сноса ТВД процедуры на грань

	vector<Edge*> edges;


	vector<Gran*> grans_surf;  // Грани-соседи, принадлежащие поверхностям разрыва
	// Этот массив актуален только для граней, принадлежащих поверхностям
	// Он заполняется в функции 


	bool work1 = false;            // Рабочий маячок, используется для разных задачь (перед использованием надо задать значения)
	int number = 0;                // номера начинаются с единицы
	Type_Gran type = Type_Gran::Us;  // по умолчанию создаём обычный
	Type_Gran_surf type2 = Type_Gran_surf::Us;  // по умолчанию создаём обычный

	vector<short int> MK_type;

	unordered_map<string, double> parameters;   // Параметры на грани (тут могут быть значения
	// плазменных полей (большие величины и т.д.), могут быть значения потоков

	unordered_map<string, double> geo_parameters;

	short int Get_method(); // Каким методом Laks/HLL/HLLC/HLLD считаем распад через грань

	double normal[2][3];           // Нормаль грани (также в предыдущий и следуюoий момент времени)
	double center[2][3];           // Центр грани (также в предыдущий и следуюoий момент времени)
	double area[2];           // Площадь грани (также в предыдущий и следуюoий момент времени)

	// Для Монте-Карло
	vector<array<AMR_f*, 2>> AMR;
	// Первая направлена по нормали грани, вторая против нормали
	double MK_Potok;  // Суммарный поток по всем сортам у грани
	unsigned short int N_particle = 0; // Число частиц, попавших в грань
	mutex mut;
	
	void Read_AMR(short int ni, short int nH, bool need_refine = false);

	bool Have_zone_number(short int z);
	// Проверяет есть ли в векторе MK_type зона с номером z?

	void Culc_measure(unsigned short int st_time);
	// вычисляет normal, center, area
	// для корректной работы определения нормали центры ячеек должны быть актуальны

	double culc_velosity(short int now1, const double& time);
	double func_R(unsigned short int i_time); // Расстояние от центра до начала координат

	// Функция сравнения по набору узлов
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	// Для Монте-Карло
	void Get_Random_pozition(Eigen::Vector3d& poz, Sensor* Sens);
	// Получить случайную позоцию на грани

	void Set_Gran_Geo_for_MK(void);
	// Заполняем необходимые параметры для МК
	// x_min, x_max, y_min, y_max, z_min, z_max - для быстрого пересечения траекторий с гранью
	//

	bool Luch_iz_cross_approx(const Eigen::Vector3d& R, const Eigen::Vector3d& V);
	// Пересекает ли луч данную грань (приблизительно, по параллелепипеду вокруг грани).

	bool Luch_crossing(const Eigen::Vector3d& R, const Eigen::Vector3d& V, double& time);

	bool rayTriangleIntersect(
		const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
		const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
		double& t);

	Gran();
};

