#pragma once

#include"Header.h"

using namespace std;


class Luch
{
public:
	string type;        // "A_Luch", "A2_Luch", .....
	static Geo_param* geo;             // Задаётся при инициализации класса Setka
	map<string, double> parameters;    // Параметры луча (могут быть свои для каждого типа)
	// "the" - угол в сферической ск от оси х - ось набегающего потока - зенитный угол
	// "phi" - угол в сферической ск - азимутальный угол
	//Yzel* start;            
	//Yzel* end;              
	vector<Yzel*> Yzels;
	vector<Yzel*> Yzels_opor;
	//vector<Luch*> Luch_soseds;    // Лучи-соседи (по ним легко искать соседние точки)

	void dvigenie(int i_time);

	Yzel* get_yzel_near_opor(int num_opor, int shift);
	// Взять узел возле такого-то опорного с таким-то сдвигом
};
