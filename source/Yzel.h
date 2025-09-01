#pragma once
#include "Header.h"
#include <vector>

class Gran;

enum class Type_yzel{
	TS,    // 0    узел на TS
	HP, // 1       узел на HP
	BS, // 2       узел на BS
	Us, // 3       обычный узел (не идентифицирован), 
	//такие есть на границе между зонами (если там нет поверхности разрыва)
	// точно есть в местах где можно продолжть HP и BS
	Zone_1,  // 4   узел до TS
	Zone_2,  // 5   узел между TS и HP
	Zone_3,  // 6   узел между HP и BS
	Zone_4  // 7   узел после BS
};

class Yzel
{
public:
	double coord[2][3];   // [¬ременной слой координат, три координаты]
	//Luch* luch;           // Ћуч на котором расположен узел (каждый узел может лежать только на одном луче)
	int number = 0;                // номера начинаютс€ с единицы
	vector<Gran*> grans;
	Type_yzel type = Type_yzel::Us;  // по умолчанию создаЄм обычный узел
	bool is_inner = false;

	double velocity[3];
	int num_velocity = 0;
	mutex mut;

	unordered_map<string, Yzel*> Yzel_sosed_sglag;
	// ”злы - соседи дл€ сглаживани€ поверхности
	// AA33  AA3  (AA)  AA1  AA11 - вдоль HP  направление в апвинд
	// A22  AA2  (AA)  AA4  AA44  - вращение HP
	std::unordered_set<Yzel*> Yzel_sosed_sglag2;

	unordered_map<string, double> parameters;   // ѕараметры в узле (дл€ интерпол€ции надо)

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // –ассто€ние до начала координат
	double func_Ryz(int i_time); // –ассто€ние до начала координат

	friend double Yzel_distance(Yzel* A, Yzel* B, int time); // рассто€ние между двум€ узлами в пространстве
	friend double Yzel_distance_x(Yzel* A, Yzel* B, int time); // рассто€ние между двум€ узлами в пространстве без учЄта х координаты
};
