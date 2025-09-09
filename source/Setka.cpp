#include "Setka.h"
#include <algorithm>

#define  macros1(n, x, y, z) A->yzels[n]->coord[0][0] = x; \
	A->yzels[n]->coord[0][1] = y; \
	A->yzels[n]->coord[0][2] = z


// Этот макрос заполняет Gran_TS   Gran_HP    Gran_BS
// Они заполняются после того, как узлам распределили их Type_yzel::
#define  macros2(name) for (auto& i : this->All_Gran)\
	{ \
		b1 = true; \
		for (auto& j : i->yzels) \
		{ \
			if (j->type != Type_yzel::name) \
			{ \
				b1 = false; \
				break; \
			} \
		} \
		if (b1 == true) \
		{ \
			this->Gran_##name.push_back(i); \
		} \
	}


#define  macros3(name, proc) if ((100.0 - d2 * 100.0 / d1) > proc) \
		{\
			k++;\
			izmen = true;\
			if (i->parameters.find(#name) != i->parameters.end())\
			{\
				i->parameters[#name] *= (1.0 + procent / 100.0);\
			}\
			else\
			{\
				i->parameters[#name] = this->geo->name * (1.0 + procent / 100.0);\
			}\
		}\
		else if ((100.0 - d2 * 100.0 / d1) < -proc)\
		{\
			k++;\
			izmen = true;\
			if (i->parameters.find(#name) != i->parameters.end())\
			{\
				i->parameters[#name] *= (1.0 - procent / 100.0);\
			}\
			else\
			{\
				i->parameters[#name] = this->geo->name * (1.0 - procent / 100.0);\
			}\
		}

// Этот макрос добавляет грани-соседи к Gran_TS   Gran_HP    Gran_BS
#define  macros4(name) for (auto& i : this->Gran_##name)\
	{\
	auto A = i->cells[0];\
	for (auto& j : A->grans)\
	{\
		if (j->number == i->number) continue;\
		auto B = A->Get_Sosed(j);\
		for (auto& k : B->grans)\
		{\
			if (k->type2 == Type_Gran_surf::name && k->number != i->number)\
			{\
				i->grans_surf.push_back(k);\
			}\
		}\
	}\
	}

using namespace std;

Setka::Setka(bool real)
{
	this->Surf1 = nullptr;
	this->geo = new Geo_param();
	this->phys_param = new Phys_param();
	Luch::geo = this->geo;

	if (real == false) return;

	// ! Это число зависит от сетки триангуляции круга (сколько там вращений по углу она подразумевает)
	//cout << "AAAAAAAAAAAAAAAAAAAAA " << endl;
	this->All_name_luch["A_Luch"] = &this->A_Luch;
	this->name_luch.reserve(1);
	this->name_luch.push_back("A_Luch");

	this->Cell_Center = new Cell;

	this->New_initial();                     // Начальное создание узлов и ячеек

	this->New_connect();                     // Начальное создание граней и связывание их с ячейками
	// Также добавляет соседей для каждой ячеки


	for (auto& i : this->All_Cell)
	{
		i->type = Type_cell::Zone_1;
	}

	// Переворачиваем сетку
	for (auto& i : this->All_Yzel)
	{
		double ccc = i->coord[0][2];
		i->coord[0][2] = i->coord[0][0];
		i->coord[0][0] = i->coord[0][1];
		i->coord[0][1] = ccc;

	}

	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}

	return;

	this->New_append_surfaces();
	// Создание граней на TS, HP, BS 
	// связывание граничных граней с соседями

	this->Renumerate();

	// заполняем коррдинаты узлов на другом временном слое
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}

	this->Tecplot_print_all_lush_in_2D();

}

Setka::~Setka()
{
	for (auto& i : this->All_Yzel)
	{
		delete i;
	}
	this->All_Yzel.clear();
	cout << "A1" << endl;
	for (auto& i : this->All_Cell)
	{
		delete i;
	}
	this->All_Cell.clear();

	cout << "A2" << endl;
	for (auto& i : this->All_Gran)
	{
		delete i;
	}
	this->All_Gran.clear();

	cout << "A3" << endl;
	for (auto& i : this->All_Luch)
	{
		delete i;
	}
	this->All_Luch.clear();

	cout << "A4" << endl;
	delete this->Cell_Center;
	cout << "A5" << endl;
	delete this->geo;
	cout << "A6" << endl;

	//delete this->phys_param;
	cout << "A7" << endl;

}

void Setka::Algoritm(short int alg)
{
	// 2 - Монте-Карло

	cout << "Start Algoritm " << alg << endl;

	this->Test_geometr();
	this->Calculating_measure(0);
	this->Calculating_measure(1);

	if (alg == 2)
	{
		// Определим зоны для МК
		this->Set_MK_Zone();

		//Проверим зоны
		this->Tecplot_print_gran_with_condition(0);
		this->Tecplot_print_gran_with_condition(1);
		this->Tecplot_print_gran_with_condition(2);
		this->Tecplot_print_gran_with_condition(3);
		this->Tecplot_print_gran_with_condition(4);
		this->Tecplot_print_gran_with_condition(5);
		this->Tecplot_print_gran_with_condition(6);

		// Надо проверить зону 7 на постоянных полях с перезарядкой (на Максвелле) - не уверен в правильности

		vector<short int> zones_number;

		zones_number.push_back(6);
		zones_number.push_back(4);
		zones_number.push_back(2);
		zones_number.push_back(1);
		zones_number.push_back(3);
		zones_number.push_back(5);
		zones_number.push_back(7);
		zones_number.push_back(5);
		zones_number.push_back(3);
		zones_number.push_back(1);
		zones_number.push_back(2);
		zones_number.push_back(4);

		for (const auto& zone_play : zones_number)
		{
			cout << "Start zone = " << zone_play << endl;
			this->MK_prepare(zone_play);
			this->MK_go(zone_play);
			this->MK_delete(zone_play);
		}
	}
}

Cell* Setka::Find_cell_point(const double& x, const double& y, const double& z, short int now, Cell*& previos)
{
	// now - какие координаты сейчас актуальны (0 по умолчанию)

	if (previos == nullptr)
	{
		previos = this->All_Cell[50];
	}

	Cell* next;
	Cell* next2;
	bool b, b2;
	unsigned int kkk = 0;

	while (true)
	{
		kkk++;
		if (kkk > 10000)
		{
			//cout << "Ne naydeno! A " << endl;
			return nullptr;
		}


		b = previos->Belong_point(x, y, z, now, true, next);

		if (b == true || next == nullptr)
		{
			next2 = next;
			b2 = previos->Belong_point(x, y, z, now, false, next);
			if (b2 == true)
			{
				return previos;
			}
			else
			{
				for (auto& gr : previos->grans)
				{
					auto C = previos->Get_Sosed(gr);
					if (C == nullptr) continue;
					b2 = C->Belong_point(x, y, z, now, false, next);
					if(b2 == true) return C;
				}

				//if(next2 == nullptr) cout << "Ne naydeno! B1    " << x << " " << y << " " << z << endl;
				//if(next2 != nullptr) cout << "Ne naydeno! B2    " << endl;

				return nullptr;
			}
		}
		else
		{
			previos = next;
			continue;
		}

	}
}

void Setka::Renumerate(void)
{
	int kkk = 0;
	for (auto& i : this->All_Yzel)
	{
		i->number = kkk;
		kkk++;
	}

	kkk = 0;
	for (auto& i : this->All_Cell)
	{
		i->number = kkk;
		kkk++;
	}

	kkk = 0;
	for (auto& i : this->All_Gran)
	{
		i->number = kkk;
		kkk++;
	}
}

bool Setka::Test_geometr(void)
{
	cout << endl;
	cout << "Testing" << endl;

	this->Calculating_measure(0);
	this->Calculating_measure(1);

	// НЕ РАБОТАЕТ
	cout << "Test 1:" << endl; // Грань ссылается либо на 1, либо на 2 ячейки
	// нормаль у грани от первой ячейки ко второй
	for (auto& i : this->All_Gran)
	{
		if (i->cells.size() != 1 && i->cells.size() != 2)
		{
			cout << "Failure: i->cells.size()  = " << i->cells.size() << endl;
			exit(-1);
		}

		if (i->cells.size() == 2)
		{
			double n1 = i->normal[0][0];
			double n2 = i->normal[0][1];
			double n3 = i->normal[0][2];

			double a1 = i->cells[1]->center[0][0] - i->cells[0]->center[0][0];
			double a2 = i->cells[1]->center[0][1] - i->cells[0]->center[0][1];
			double a3 = i->cells[1]->center[0][2] - i->cells[0]->center[0][2];

			/*if (scalarProductFast(n1, n2, n3, a1, a2, a3) < 0.0)
			{
				cout << "Failure: 5678675634   = " << scalarProductFast(n1, n2, n3, a1, a2, a3) << endl;
				cout << i->cells[0]->center[0][0] << " " << i->cells[0]->center[0][1] << " " <<
					i->cells[0]->center[0][2] << endl;
				cout << i->cells[1]->center[0][0] << " " << i->cells[1]->center[0][1] << " " <<
					i->cells[1]->center[0][2] << endl;
				cout << n1 << " " << n2 << " " << n3 << endl;
				cout << i->center[0][0] << " " << i->center[0][1] << " " << i->center[0][2] << endl;
				exit(-1);
			}*/
		}


	}
	cout << "Success" << endl;


	cout << "Test 2:" << endl;  // У ячейки 6 граней и 8 узлов
	for (auto& i : this->All_Cell)
	{
		if (i->grans.size() != 6)
		{
			cout << "Failure: i->grans.size()  = " << i->grans.size() << endl;
			exit(-1);
		}

		if (i->yzels.size() != 8)
		{
			cout << "Failure: i->yzels.size()  = " << i->yzels.size() << endl;
			exit(-1);
		}
	}
	cout << "Success" << endl;

	cout << "Test 5:" << endl; // Все нормали грани единичные и в ячейках выполнен геометрический закон суммы нормалей
	for (auto& i : this->All_Cell)
	{
		double V1 = 0.0;
		double V2 = 0.0;
		double V3 = 0.0;

		for (auto& j : i->grans)
		{
			double n1 = j->normal[0][0];
			double n2 = j->normal[0][1];
			double n3 = j->normal[0][2];
			if (fabs(norm2(n1, n2, n3) - 1.0) > 0.00001)
			{
				cout << "Error 5654354876" << endl;
				whach(norm2(n1, n2, n3));
				exit(-1);
			}

			if (j->cells[0] != i)
			{
				n1 *= -1;
				n2 *= -1;
				n3 *= -1;
			}
			V1 += n1 * j->area[0];
			V2 += n2 * j->area[0];
			V3 += n3 * j->area[0];
		}

		if (false)//(norm2(V1, V2, V3)/i->volume[0] > 0.1)
		{
			cout << "Error 8768978654" << endl;
			whach(V1);
			whach(V2);
			whach(V3);
			whach(i->volume[0]);
			whach(i->center[0][0]);
			whach(i->center[0][1]);
			whach(i->center[0][2]);
			whach(norm2(V1, V2, V3) / i->volume[0]);
			exit(-1);
		}
	}
	cout << "Success" << endl;


	cout << "Test 6:" << endl; // У граней все узлы разные и их 4 штуки!
	for (auto& i : this->All_Gran)
	{
		if (i->yzels[0]->number == i->yzels[1]->number ||
			i->yzels[0]->number == i->yzels[2]->number ||
			i->yzels[0]->number == i->yzels[3]->number ||
			i->yzels[1]->number == i->yzels[2]->number ||
			i->yzels[1]->number == i->yzels[3]->number ||
			i->yzels[2]->number == i->yzels[3]->number)
		{
			cout << "Error 9865134089" << endl;
			whach(i->yzels[0]->number);
			whach(i->yzels[1]->number);
			whach(i->yzels[2]->number);
			whach(i->yzels[3]->number);
			exit(-1);
		}
	}
	cout << "Success" << endl << endl;


	cout << "Test 6:" << endl; // У граней все узлы есть в ячейках
	for (auto& i : this->All_Gran)
	{
		auto cell = i->cells[0];
		for (const auto& j : i->yzels)
		{
			//TODOOOOOOOO
		}
	}
	cout << "Success" << endl << endl;

	return true;
}

void Setka::Set_luch_parametr()
{
	// Добавляет лучам необходимые параметры 
	// (например, углы the и phi для радиальных лучей A, B, C типов) - это ускорит расчёты, так как 
	// они остаются постоянными для луча, и можно постоянно не пересчитывать значения
	for (auto& i : this->A2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->A_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->B_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->D_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->E_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->H_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->G_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& i : this->C2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->C_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}
}

void Setka::Read_old_surface(string name)
{
	cout << "Start Read_old_surface: " << name << endl;
	if (this->Surf1 != nullptr)
	{
		this->Surf1->~Surfaces();
		delete this->Surf1;
	}

	this->Surf1 = new Surfaces();
	this->Surf1->Read_old(name);
	cout << "End Read_old_surface: " << name << endl;
}

void Setka::Move_to_surf(Surfaces* Surf)
{
	double x, y, z, phi, the, r, rr;
	cout << "Start: Move_to_surf" << endl;
	double R_BS; // положение внешней ударной волны (для движение B E D лучей)
	for (int st = 0; st < 15; st++)
	{
		for (auto& i : this->A_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);

				

				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}


				x = j->Yzels_opor[2]->coord[0][0];
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(x, y, z));

				rr = Surf->Get_HP(phi, the, 0);

				j->Yzels_opor[2]->coord[0][0] *= rr / r;
				j->Yzels_opor[2]->coord[0][1] *= rr / r;
				j->Yzels_opor[2]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "9443563295   errjr" << endl;
				}

				x = j->Yzels_opor[3]->coord[0][0];
				y = j->Yzels_opor[3]->coord[0][1];
				z = j->Yzels_opor[3]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);


				rr = Surf->Get_BS(phi, the);

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "5794671565   errjr" << endl;
				}

				j->Yzels_opor[3]->coord[0][0] *= rr / r;
				j->Yzels_opor[3]->coord[0][1] *= rr / r;
				j->Yzels_opor[3]->coord[0][2] *= rr / r;

			}
		}

		for (auto& j : this->A2_Luch)
		{
			x = j->Yzels_opor[1]->coord[0][0];
			y = j->Yzels_opor[1]->coord[0][1];
			z = j->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);



			rr = Surf->Get_TS(phi, the);

			j->Yzels_opor[1]->coord[0][0] *= rr / r;
			j->Yzels_opor[1]->coord[0][1] *= rr / r;
			j->Yzels_opor[1]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "0989898653   errjr" << endl;
			}


			x = j->Yzels_opor[2]->coord[0][0];
			y = j->Yzels_opor[2]->coord[0][1];
			z = j->Yzels_opor[2]->coord[0][2];
			r = sqrt(kvv(x, y, z));

			rr = Surf->Get_HP(phi, the, 0);

			j->Yzels_opor[2]->coord[0][0] *= rr / r;
			j->Yzels_opor[2]->coord[0][1] *= rr / r;
			j->Yzels_opor[2]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "9443563295   errjr" << endl;
			}

			x = j->Yzels_opor[3]->coord[0][0];
			y = j->Yzels_opor[3]->coord[0][1];
			z = j->Yzels_opor[3]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);


			rr = Surf->Get_BS(phi, the);

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "5794671565   errjr" << endl;
			}

			j->Yzels_opor[3]->coord[0][0] *= rr / r;
			j->Yzels_opor[3]->coord[0][1] *= rr / r;
			j->Yzels_opor[3]->coord[0][2] *= rr / r;

		}

		for (auto& i : this->B_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);

				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}


				x = j->Yzels_opor[2]->coord[0][0];
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];

				rr = Surf->Get_HP(phi, x, 1);
				r = sqrt(kvv(0.0, y, z));
				j->Yzels_opor[2]->coord[0][1] *= rr / r;
				j->Yzels_opor[2]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "6510292073   errjr" << endl;
				}

			}
		}

		for (auto& i : this->C_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);



				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}
			}
		}

		for (auto& j : this->C2_Luch)
		{
			x = j->Yzels_opor[1]->coord[0][0];
			y = j->Yzels_opor[1]->coord[0][1];
			z = j->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);



			rr = Surf->Get_TS(phi, the);

			j->Yzels_opor[1]->coord[0][0] *= rr / r;
			j->Yzels_opor[1]->coord[0][1] *= rr / r;
			j->Yzels_opor[1]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "0989898653   errjr" << endl;
			}
		}

		for (auto& i : this->D_Luch)
		{
			x = i[0]->Yzels_opor[1]->coord[0][0];
			y = i[0]->Yzels_opor[1]->coord[0][1];
			z = i[0]->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(0.0, y, z));
			phi = polar_angle(y, z);
			rr = Surf->Get_HP(phi, x, 1);

			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				//phi = polar_angle(y, z);
				//rr = Surf->Get_HP(phi, x, 1);

				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

			}
		}

		for (auto& i : this->E_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				phi = polar_angle(y, z);
				rr = Surf->Get_HP(phi, x, 1);

				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "6510292073   errjr" << endl;
				}
			}
		}
	
		// Двигаем BS для B E D лучей
		for (int i = 0; i < this->B_Luch.size(); i++)
		{
			R_BS = this->A_Luch[i].back()->Yzels_opor[3]->func_R(0);
			for (auto& j : this->B_Luch[i])
			{
				y = j->Yzels_opor[3]->coord[0][1];
				z = j->Yzels_opor[3]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[3]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[3]->coord[0][2] *= R_BS / r;
			}

			for (auto& j : this->E_Luch[i])
			{
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[2]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[2]->coord[0][2] *= R_BS / r;
			}

			for (auto& j : this->D_Luch[i])
			{
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[2]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[2]->coord[0][2] *= R_BS / r;
			}
		}

		for (auto& i : this->All_Luch)
		{
			i->dvigenie(0);
		}

	}

	cout << "End: Move_to_surf" << endl;
}

void Setka::New_initial()
{
	cout << "---START New_initial---" << endl;

	int n, i, j, n2, n3;
	double x, y;

	// Считываем файл круга
	ifstream fin("SDK1_krug_setka.bin", ios::binary | ios::in);  // SDK1_krug_setka
	if (!fin)
	{
		cout << "Net takogo fajla (fajl setki v krugu)" << endl;
		exit(-1);
	}

	fin.read((char*)&n, sizeof n);  // Версия файла
	cout << "Versiya fajla setki v krugu: " << n << endl;

	// Считываем узлы (сохраняем их в первый слой головных узлов)
	
	this->Krug_Yzel.resize(this->geo->M1);
	this->Krug_Yzel_2.resize(this->geo->M1);

	fin.read((char*)&n, sizeof n); // Число узлов

	for (int i = 0; i < this->Krug_Yzel.size(); i++)
	{
		this->Krug_Yzel[i].resize(n);
	}

	for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
	{
		this->Krug_Yzel_2[i].resize(n);
	}
		
	for (int i = 0; i < n; i++)
	{
		fin.read((char*)&x, sizeof x);
		fin.read((char*)&y, sizeof y);

		for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
		{
			auto yzel = new Yzel(x, y, 0.0);
			yzel->number = -1;
			this->Krug_Yzel[ii][i] = yzel;
			this->All_Yzel.push_back(yzel);
		}

		for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			auto yzel = new Yzel(x, y, 0.0);
			yzel->number = -1;
			this->Krug_Yzel_2[ii][i] = yzel;
			this->All_Yzel.push_back(yzel);
		}
	}

	double Rmax = 1.0;
	double thetamin = 4.0;
	// Считаем радиус круга

	Yzel* AAA = nullptr;
	for (int i = 0; i < this->Krug_Yzel[0].size(); i++)
	{
		if (this->Krug_Yzel[0][i]->func_R(0) > 1.0 - 0.00001)
		{
			AAA = this->Krug_Yzel[0][i];
			break;
		}
	}

	thetamin = polar_angle(AAA->coord[0][0], AAA->coord[0][1]);

	// Ищем номера узлов, которые лежат на границе
	vector<int> num_gr;
	for (int i = 0; i < this->Krug_Yzel[0].size(); i++)
	{
		if (this->Krug_Yzel[0][i]->func_R(0) > 1.0 - 0.00001)
		{
			num_gr.push_back(i);
		}
	}

	if (num_gr.size() != this->geo->Nphi)
	{
		cout << "Error 93it9hgirnig" << endl;
		cout << num_gr.size() << " " << this->geo->Nphi << endl;
		exit(-1);
	}



	// Двигаем координаты круга 1 на сферу
	for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
	{
		for (auto& i : this->Krug_Yzel[ii])
		{
			i->coord[0][0] /= Rmax;
			i->coord[0][1] /= Rmax;
			i->coord[0][2] = i->coord[0][1];
			i->coord[0][1] = i->coord[0][0];
			i->coord[0][0] = 1.0 / tan(this->geo->tetta0);
			double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
			i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
			i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
			i->coord[0][2] *= this->geo->R0 / r * (ii + 1);
		}
	}

	// Вращаем координаты для совпадения с точками второй сетки
	for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
	{
		for (auto& i : this->Krug_Yzel[ii])
		{
			double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
			double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
			i->coord[0][1] = newy;
			i->coord[0][2] = newz;
		}
	}

	// Двигаем координаты круга 2 на сферу
	for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
	{
		for (auto& i : this->Krug_Yzel_2[ii])
		{
			i->coord[0][0] /= Rmax;
			i->coord[0][1] /= Rmax;
			i->coord[0][2] = i->coord[0][1];
			i->coord[0][1] = i->coord[0][0];
			i->coord[0][0] = 1.0 / tan(const_pi - this->geo->tetta0);
			double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
			i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
			i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
			i->coord[0][2] *= this->geo->R0 / r * (ii + 1);

		}
	}

	// Вращаем координаты для совпадения с точками второй сетки
	for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel_2[ii])
			{
				double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
				double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}
		}


	//теперь сортируем точки на границе кругов в правильном порядке
	std::sort(num_gr.begin(), num_gr.end(),
		[this](int a, int b)
		{
			auto yz1 = this->Krug_Yzel[0][a];
			auto yz2 = this->Krug_Yzel[0][b];
			double phi1 = polar_angle(yz1->coord[0][1], yz1->coord[0][2]);
			double phi2 = polar_angle(yz2->coord[0][1], yz2->coord[0][2]);
			return phi1 < phi2; // Собственная логика сравнения
		}
	);


	// Двигаем все точки на нужные радиусы

	for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
	{
		double rr = this->geo->R0 + ii * (this->geo->R1 - this->geo->R0) / (this->Krug_Yzel_2.size() - 1);
		for (auto& i : this->Krug_Yzel_2[ii])
		{
			double rrr = norm2(i->coord[0][0], i->coord[0][1], i->coord[0][2]);
			i->coord[0][0] *= rr / rrr;
			i->coord[0][1] *= rr / rrr;
			i->coord[0][2] *= rr / rrr;
		}
	}

	for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
	{
		double rr = this->geo->R0 + ii * (this->geo->R1 - this->geo->R0) / (this->Krug_Yzel.size() - 1);
		for (auto& i : this->Krug_Yzel[ii])
		{
			double rrr = norm2(i->coord[0][0], i->coord[0][1], i->coord[0][2]);
			i->coord[0][0] *= rr / rrr;
			i->coord[0][1] *= rr / rrr;
			i->coord[0][2] *= rr / rrr;
		}
	}

	// Перенумеровываем все узлы
	if (true)
	{
		int kkk = 0;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// Считываем ячейки (для точек в кругах)
	if (true)
	{
		this->Cell_layer_head.resize(this->Krug_Yzel.size() - 1);
		this->Cell_layer_tail.resize(this->Krug_Yzel_2.size() - 1);

		int nn, ny;
		fin.read((char*)&n, sizeof n); // Число ячеек

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			this->Cell_layer_head[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			this->Cell_layer_tail[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_head[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_tail[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_tail[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}


		for (int i = 0; i < n; i++)
		{
			fin.read((char*)&nn, sizeof nn); // Число узлов в ячейке
			for (int j = 0; j < nn; j++)
			{
				fin.read((char*)&ny, sizeof ny); // Номер узла

				for (int k = 0; k < this->Cell_layer_head.size(); k++)
				{
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k][ny - 1]);
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k + 1][ny - 1]);
				}

				for (int k = 0; k < this->Cell_layer_tail.size(); k++)
				{
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k][ny - 1]);
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k + 1][ny - 1]);
				}
			}

		}

		// Меняем местами узлы, чтобы они шли в правильном порядке
		if (true)
		{
			for (int i = 0; i < this->Cell_layer_head.size(); i++)
			{
				for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
				{
					auto a1 = this->Cell_layer_head[i][j]->yzels[0];
					auto a2 = this->Cell_layer_head[i][j]->yzels[1];
					auto a3 = this->Cell_layer_head[i][j]->yzels[2];
					auto a4 = this->Cell_layer_head[i][j]->yzels[3];
					auto a5 = this->Cell_layer_head[i][j]->yzels[4];
					auto a6 = this->Cell_layer_head[i][j]->yzels[5];
					auto a7 = this->Cell_layer_head[i][j]->yzels[6];
					auto a8 = this->Cell_layer_head[i][j]->yzels[7];

					this->Cell_layer_head[i][j]->yzels[1] = a3;
					this->Cell_layer_head[i][j]->yzels[2] = a5;
					this->Cell_layer_head[i][j]->yzels[3] = a7;
					this->Cell_layer_head[i][j]->yzels[4] = a2;
					this->Cell_layer_head[i][j]->yzels[5] = a4;
					this->Cell_layer_head[i][j]->yzels[6] = a6;
					this->Cell_layer_head[i][j]->yzels[7] = a8;
				}
			}

			for (int i = 0; i < this->Cell_layer_tail.size(); i++)
			{
				for (int j = 0; j < this->Cell_layer_tail[0].size(); j++)
				{
					auto a1 = this->Cell_layer_tail[i][j]->yzels[0];
					auto a2 = this->Cell_layer_tail[i][j]->yzels[1];
					auto a3 = this->Cell_layer_tail[i][j]->yzels[2];
					auto a4 = this->Cell_layer_tail[i][j]->yzels[3];
					auto a5 = this->Cell_layer_tail[i][j]->yzels[4];
					auto a6 = this->Cell_layer_tail[i][j]->yzels[5];
					auto a7 = this->Cell_layer_tail[i][j]->yzels[6];
					auto a8 = this->Cell_layer_tail[i][j]->yzels[7];

					this->Cell_layer_tail[i][j]->yzels[1] = a3;
					this->Cell_layer_tail[i][j]->yzels[2] = a5;
					this->Cell_layer_tail[i][j]->yzels[3] = a7;
					this->Cell_layer_tail[i][j]->yzels[4] = a2;
					this->Cell_layer_tail[i][j]->yzels[5] = a4;
					this->Cell_layer_tail[i][j]->yzels[6] = a6;
					this->Cell_layer_tail[i][j]->yzels[7] = a8;
				}
			}
		}


	}

	// Теперь строим регулярную сетку
	this->Yzel_regular.resize(boost::extents[this->geo->Nphi][this->geo->N1][this->geo->M1]);
	this->Cell_regular.resize(boost::extents[this->geo->Nphi][this->geo->N1 - 1][this->geo->M1 - 1]);

	for (int i = 0; i < this->geo->Nphi; ++i) 
	{
		for (int j = 0; j < this->geo->N1; ++j) 
		{
			for (int k = 0; k < this->geo->M1; ++k)
			{
				double phi = 2.0 * const_pi * i / (this->geo->Nphi);    // от 0 до 2π
				double theta = this->geo->tetta0 + (j + 1) * (const_pi - 2 * this->geo->tetta0) / (this->geo->N1 + 1);        // от 0 до π
				double r = this->geo->R0 + k * (this->geo->R1 - this->geo->R0) / (this->geo->M1 - 1);

				double y = r * cos(phi) * sin(theta);   
				double z = r * sin(phi) * sin(theta);
				double x = r * cos(theta);

				auto yzel = new Yzel(x, y, z);
				yzel->number = -1;
				this->Yzel_regular[i][j][k] = yzel;
				this->All_Yzel.push_back(yzel);
			}
		}
	}

	for (int i = 0; i < this->geo->Nphi; ++i)
	{
		for (int j = 0; j < this->geo->N1 - 1; ++j)
		{
			for (int k = 0; k < this->geo->M1 - 1; ++k)
			{
				auto cell = new Cell();
				this->Cell_regular[i][j][k] = cell;
				this->All_Cell.push_back(cell);

				int i1 = i + 1;
				if (i1 > this->geo->Nphi - 1) i1 = 0;

				cell->yzels.push_back(this->Yzel_regular[i][j][k]);
				cell->yzels.push_back(this->Yzel_regular[i][j + 1][k]);
				cell->yzels.push_back(this->Yzel_regular[i1][j + 1][k]);
				cell->yzels.push_back(this->Yzel_regular[i1][j][k]);

				cell->yzels.push_back(this->Yzel_regular[i][j][k + 1]);
				cell->yzels.push_back(this->Yzel_regular[i][j + 1][k + 1]);
				cell->yzels.push_back(this->Yzel_regular[i1][j + 1][k + 1]);
				cell->yzels.push_back(this->Yzel_regular[i1][j][k + 1]);
			}
		}
	}

	// Теперь создаём ячейки на стыке сеток
	for (int i = 0; i < this->geo->Nphi; ++i)
	{
		for (int k = 0; k < this->geo->M1 - 1; ++k)
		{
			auto cell = new Cell();
			this->All_Cell.push_back(cell);

			int i1 = i + 1;
			if (i1 > this->geo->Nphi - 1) i1 = 0;

			int jj1 = num_gr[i];
			int jj2 = num_gr[i1];

			cell->yzels.push_back(this->Krug_Yzel[k][jj1]);
			cell->yzels.push_back(this->Yzel_regular[i][0][k]);
			cell->yzels.push_back(this->Yzel_regular[i1][0][k]);
			cell->yzels.push_back(this->Krug_Yzel[k][jj2]);

			cell->yzels.push_back(this->Krug_Yzel[k + 1][jj1]);
			cell->yzels.push_back(this->Yzel_regular[i][0][k + 1]);
			cell->yzels.push_back(this->Yzel_regular[i1][0][k + 1]);
			cell->yzels.push_back(this->Krug_Yzel[k + 1][jj2]);
		}
	}

	for (int i = 0; i < this->geo->Nphi; ++i)
	{
		for (int k = 0; k < this->geo->M1 - 1; ++k)
		{
			auto cell = new Cell();
			this->All_Cell.push_back(cell);

			int i1 = i + 1;
			if (i1 > this->geo->Nphi - 1) i1 = 0;

			int jj1 = num_gr[i];
			int jj2 = num_gr[i1];

			
			cell->yzels.push_back(this->Yzel_regular[i][this->geo->N1 - 1][k]);
			cell->yzels.push_back(this->Krug_Yzel_2[k][jj1]);
			cell->yzels.push_back(this->Krug_Yzel_2[k][jj2]);
			cell->yzels.push_back(this->Yzel_regular[i1][this->geo->N1 - 1][k]);
			

			cell->yzels.push_back(this->Yzel_regular[i][this->geo->N1 - 1][k + 1]);
			cell->yzels.push_back(this->Krug_Yzel_2[k + 1][jj1]);
			cell->yzels.push_back(this->Krug_Yzel_2[k + 1][jj2]);
			cell->yzels.push_back(this->Yzel_regular[i1][this->geo->N1 - 1][k + 1]);
		}
	}



	fin.close();
	cout << "---END New_initial---" << endl;
}

struct GranHash {
	size_t operator()(const Gran* g) const {
		// Создаем хэш на основе отсортированных ID узлов
		std::vector<int> ids;
		for (auto yzel : g->yzels) ids.push_back(yzel->number);
		std::sort(ids.begin(), ids.end());
		size_t seed = 0;
		for (int id : ids) {
			seed ^= std::hash<int>()(id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

struct GranEqual {
	bool operator()(const Gran* a, const Gran* b) const {
		// 1. Если количество узлов разное, грани не равны
		if (a->yzels.size() != b->yzels.size()) {
			return false;
		}

		// 2. Собираем номера узлов из обеих граней
		std::vector<int> a_ids, b_ids;
		for (auto* yzel : a->yzels) a_ids.push_back(yzel->number);
		for (auto* yzel : b->yzels) b_ids.push_back(yzel->number);

		// 3. Сортируем для независимости от порядка
		std::sort(a_ids.begin(), a_ids.end());
		std::sort(b_ids.begin(), b_ids.end());

		// 4. Сравниваем отсортированные вектора
		return a_ids == b_ids;
	}
};

void Setka::New_connect()
{
	// План следующий: сначала создаём для каждой ячейки 6 граней
	// Потом находим повторы, удаляем лишние и связываем грани с ячейками
	// Для ускорения удаления и т.д. сначала всё делаем в локально созданном списке, потом перенесём всё в вектор
	this->Renumerate();
	cout << "START: New_connect" << endl;

	//std::unordered_set<Gran*, GranHash> unique_grans;
	std::unordered_set<Gran*, GranHash, GranEqual> unique_grans;
	

	// Создаём все грани
	for (auto& i : this->All_Cell)
	{
		auto G = new Gran();
		G->yzels.push_back(i->yzels[0]); // 1
		G->yzels.push_back(i->yzels[3]); // 4
		G->yzels.push_back(i->yzels[2]); // 3
		G->yzels.push_back(i->yzels[1]); // 2

		auto it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}


		G = new Gran();
		G->yzels.push_back(i->yzels[4]); // 5
		G->yzels.push_back(i->yzels[5]); // 6
		G->yzels.push_back(i->yzels[6]); // 7
		G->yzels.push_back(i->yzels[7]); // 8
		
		it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}


		G = new Gran();
		G->yzels.push_back(i->yzels[0]); // 1
		G->yzels.push_back(i->yzels[1]); // 2
		G->yzels.push_back(i->yzels[5]); // 6
		G->yzels.push_back(i->yzels[4]); // 5
		
		it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}


		G = new Gran();
		G->yzels.push_back(i->yzels[2]); // 3
		G->yzels.push_back(i->yzels[3]); // 4
		G->yzels.push_back(i->yzels[7]); // 8
		G->yzels.push_back(i->yzels[6]); // 7
		
		it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}


		G = new Gran();
		G->yzels.push_back(i->yzels[1]); // 2
		G->yzels.push_back(i->yzels[2]); // 3
		G->yzels.push_back(i->yzels[6]); // 7
		G->yzels.push_back(i->yzels[5]); // 6
		
		it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}


		G = new Gran();
		G->yzels.push_back(i->yzels[0]); // 1
		G->yzels.push_back(i->yzels[4]); // 5
		G->yzels.push_back(i->yzels[7]); // 8
		G->yzels.push_back(i->yzels[3]); // 4
		
		it = unique_grans.find(G);
		if (it == unique_grans.end()) {
			// Дубликата нет, добавляем
			unique_grans.insert(G);
			// Добавляем ячейку к грани
			G->cells.push_back(i);
		}
		else
		{
			// Дубликат найден, добавляем ячейку к существующей грани
			(*it)->cells.push_back(i);
			delete G; // Удаляем временную грань
		}

	}

	// Переносим уникальные грани в All_Gran
	for (auto gran : unique_grans) 
	{
		All_Gran.push_back(gran);
	}

	// Сдесь было наделано много лишних граней

	// Удаляем дубликаты граней

	// нумеруем все грани (для того, чтобы запомнить, какие надо удалить)
	int kkk = 0;
	for (auto& i : All_Gran)
	{
		i->number = kkk;
		kkk++;
	}

	// Пробегаемся по всем граням и связываем узлы с гранями
	for (auto& i : this->All_Gran)
	{
		for (auto& j : i->yzels)
		{
			j->grans.push_back(i);
		}
	}


	// Пробегаемся по всем граням и связываем грани с ячейками
	for (auto& i : this->All_Gran)
	{
		for (auto& j : i->cells)
		{
			j->grans.push_back(i);
		}
	}

	unique_grans.clear();

	cout << "END: New_connect" << endl;

}

void Setka::New_append_surfaces()
{
	// Сначала зададим типы узлам (через опорные точки на лучах проще всего), 
	// потом пробежимся по граням и найдём нужные

	for (auto& i : this->All_Yzel)
	{
		i->type = Type_yzel::Us;
	}

	for (auto& i : this->All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 + 1 + this->geo->M2 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 
					+ 1 + this->geo->M2 + 1 + this->geo->M3 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_3;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_4;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
			i->Yzels_opor[3]->type = Type_yzel::BS;
		}
		else if (i->type == "B_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 + 1 + this->geo->M2 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1
					+ 1 + this->geo->M2 + 1 + this->geo->M3 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_3;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_4;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
		}
		else if (i->type == "C_Luch" || i->type == "C2_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
		}
		else if (i->type == "E_Luch")
		{
			short int kl = 1;
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (i->Yzels[j] == i->Yzels_opor[1]) kl++;
				if (i->Yzels[j] == i->Yzels_opor[2]) kl++;

				if(kl == 1) i->Yzels[j]->type = Type_yzel::Zone_2;
				if(kl == 2) i->Yzels[j]->type = Type_yzel::Zone_3;
				if(kl == 3) i->Yzels[j]->type = Type_yzel::Zone_4;
			}

			i->Yzels_opor[2]->type = Type_yzel::Us;

			i->Yzels_opor[1]->type = Type_yzel::HP;
		}
		else if (i->type == "D_Luch")
		{
			short int kl = 1;
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (i->Yzels[j] == i->Yzels_opor[1]) kl++;
				if (i->Yzels[j] == i->Yzels_opor[2]) kl++;

				if (kl == 1) i->Yzels[j]->type = Type_yzel::Zone_2;
				if (kl == 2) i->Yzels[j]->type = Type_yzel::Zone_3;
				if (kl == 3) i->Yzels[j]->type = Type_yzel::Zone_4;
			}

			i->Yzels_opor[1]->type = Type_yzel::Us;
			i->Yzels_opor[2]->type = Type_yzel::Us;

			if (i->Yzels_opor[1]->coord[0][0] >= this->geo->L6 - 0.00001)
			{
				i->Yzels_opor[1]->type = Type_yzel::HP;
			}
		}
		else if (i->type == "H_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				i->Yzels[j]->type = Type_yzel::Zone_2;
			}
		}
		else if (i->type == "G_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size() - 1; j++)
			{
				i->Yzels[j]->type = Type_yzel::Zone_2;
			}
		}
	}

	bool b1 = false;
	// Тут заполняются Gran_TS, Gran_HP, Gran_BS
	macros2(TS);
	macros2(HP);
	macros2(BS);

	for (auto& i : this->Gran_TS)
	{
		i->type2 = Type_Gran_surf::TS;
	}

	for (auto& i : this->Gran_BS)
	{
		i->type2 = Type_Gran_surf::BS;
	}

	for (auto& i : this->Gran_HP)
	{
		i->type2 = Type_Gran_surf::HP;
	}

	this->Renumerate();

	// Зададим зону ячейкам
	if (true)
	{
		for (auto& Cell : this->All_Cell)
		{
			int z1 = 0;
			int z2 = 0;
			int z3 = 0;
			int z4 = 0;

			for (auto& yz : Cell->yzels)
			{
				if (yz->type == Type_yzel::Zone_1)
				{
					z1++;
				}
				if (yz->type == Type_yzel::Zone_2)
				{
					z2++;
				}
				if (yz->type == Type_yzel::Zone_3)
				{
					z3++;
				}
				if (yz->type == Type_yzel::Zone_4)
				{
					z4++;
				}
			}

			if (z1 > z2 && z1 > z3 && z1 > z4)
			{
				Cell->type = Type_cell::Zone_1;
			}
			else if (z2 > z1 && z2 > z3 && z2 > z4)
			{
				Cell->type = Type_cell::Zone_2;
			}
			else if (z3 > z1 && z3 > z2 && z3 > z4)
			{
				Cell->type = Type_cell::Zone_3;
			}
			else
			{
				Cell->type = Type_cell::Zone_4;
			}
		}
	}


	macros4(TS);
	macros4(HP);
	macros4(BS);

}

void Setka::Write_file_for_FCMHD(void)
{
	this->Renumerate();

	ofstream file("FCMHD_2.bin", ios::binary);
	if (!file.is_open()) 
	{
		cout << "ERROR FCMHD_2.bin" << endl;
		exit(-1);
	}

	int n1 = this->All_Cell.size();
	int n2 = this->All_Gran.size();

	file.write(reinterpret_cast<const char*>(&this->phys_param->T_all), sizeof(double));
	file.write(reinterpret_cast<const char*>(&n1), sizeof(int));
	file.write(reinterpret_cast<const char*>(&n2), sizeof(int));

	size_t host_N_cell = this->All_Cell.size();
	size_t host_N_gran = this->All_Gran.size();
	const vector<string> param_order = { "rho", "Vx", "Vy", "Vz", "p", "Bx", "By", "Bz" };

	// host_Cell_par
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (const string& param : param_order) 
		{
			double value = this->All_Cell[i]->parameters[0].at(param);
			if (i == 0)
			{
				cout << "S1 - " << value << endl;
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// host_Cell_center
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 3; k++) 
		{
			double value = this->All_Cell[i]->center[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// host_Cell_Volume
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			double value = this->All_Cell[i]->volume[k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// host_Cell_gran
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 6; k++)
		{
			int value = this->All_Cell[i]->grans[k]->number + 1;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	// host_Gran_normal
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 3; k++)
		{
			double value = this->All_Gran[i]->normal[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// host_Gran_square
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			double value = this->All_Gran[i]->area[0];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}


	// host_Gran_center
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 3; k++)
		{
			double value = this->All_Gran[i]->center[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// host_Gran_neighbour
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 2; k++)
		{
			Cell* cc = nullptr;
			if(this->All_Gran[i]->cells.size() >= k + 1) cc = this->All_Gran[i]->cells[k];
			int value = 0;
			if (cc != nullptr)
			{
				value = cc->number + 1;
			}

			if (value > host_N_cell)
			{
				cout << "ERROR 38u4rh8fur   " << value << " " << host_N_cell << endl;
				cout << cc->number << endl;
				cout << this->All_Gran[i]->cells.size() << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	// host_Gran_neighbour_TVD
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 2; k++)
		{
			Cell* cc = nullptr;
			if (this->All_Gran[i]->cells_TVD.size() >= k + 1) cc = this->All_Gran[i]->cells_TVD[k];
			int value = 0;
			if (cc != nullptr)
			{
				value = cc->number + 1;
			}

			if (i == 1)
			{
				cout << "S2 - " << value << endl;
			}

			if (value < 0 || value > host_N_cell)
			{
				cout << "ERRORRcv 8u39874fy8eh" << endl;
				cout << value << endl;
				exit(-1);
			}

			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}
	double vv = 121.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	// host_Gran_type
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			auto cc = this->All_Gran[i];
			int value = 0;
			if (cc->type == Type_Gran::Us) value = 1;
			if (cc->type == Type_Gran::Inner_Hard) value = 2;
			if (cc->type == Type_Gran::Outer_Soft) value = 3;

			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	vv = 122.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	const vector<string> param_order2 = { "Prho", "PVx", "PVy", "PVz", "Pp", "PBx", "PBy", "PBz", "PdivB"};

	// host_Gran_POTOK
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (const string& param : param_order2)
		{
			auto cc = this->All_Gran[i];
			double value = 0.0;
			if (cc->type == Type_Gran::Inner_Hard)
			{
				string name;
				// Здесь надо посчитать потоки через граничную грань (так как они всегда одинаковые, их можно посчитать один раз)
				double time = Culc_Gran_Potok(cc, 0, 3, name, 1.0);
				value = cc->parameters[param];
			}

			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	// Проверяющий параметр
	vv = 123.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	file.close();
}

void Setka::Read_file_for_FCMHD(void)
{
	this->Renumerate();

	std::ifstream file("FCMHD_1.10_out.bin", std::ios::binary);
	if (!file.is_open()) {
		throw std::runtime_error("Error opening file: FCMHD_1.8_out.bin");
	}

	int n1 = this->All_Cell.size();
	int n2 = this->All_Gran.size();

	file.read(reinterpret_cast<char*>(&this->phys_param->T_all), sizeof(double));

	size_t host_N_cell = this->All_Cell.size();
	size_t host_N_gran = this->All_Gran.size();
	const vector<string> param_order = { "rho", "Vx", "Vy", "Vz", "p", "Bx", "By", "Bz" };

	// host_Cell_par
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (const string& param : param_order)
		{
			double value;
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			this->All_Cell[i]->parameters[0][param] = value;
		}
	}

	// Проверяющий параметр
	double vv;
	file.read(reinterpret_cast<char*>(&vv), sizeof(double));

	cout << "Proverka (321) = " << vv << endl;

	file.close();
}


void Setka::Calculating_measure(unsigned short int st_time)
{
	// Следующий блок был для тестирования
	/*auto A = this->All_Cell[0];
	macros1(0, 0, 0, 0);
	macros1(1, 1, 0, 0);
	macros1(2, 1, 1, 0);
	macros1(3, 0, 1, 0);
	macros1(4, 0, 0, 1);
	macros1(5, 1, 0, 1);
	macros1(6, 1, 1, 1);
	macros1(7, 0, 1, 1);*/

	// Следующий порядкок вычисления является важным
#pragma omp parallel for
	//for (auto& i : this->All_Cell)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
	{
		this->All_Cell[idx]->Culc_center(st_time);
	}
#pragma omp barrier

#pragma omp parallel for
	//for (auto& i : this->All_Gran)
	for (size_t idx = 0; idx < this->All_Gran.size(); ++idx) 
	{
		this->All_Gran[idx]->Culc_measure(st_time);
	}

#pragma omp barrier

#pragma omp parallel for
	//for (auto& i : this->All_Cell)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
	{
			this->All_Cell[idx]->Culc_volume(st_time);
	}
}

void Setka::auto_set_luch_geo_parameter(int for_new)
{
	// автоматическая настройки сгущений сетки с разных областях
	// for_new = 0 - значит это первый запуск функции и параметры далеки от идеальных
	// в этом случае движение изначально будет большое
	// for_new = 1 - небольшое движение, если это не первый запуск и поверхности уже стоят где надо

	cout << "Start: Izmenenie geo parameters" << endl;
	// Эту функцию совместно с функцией движения сетки (luch.cpp) можно улучшать и дополнять в процессе
	bool izmen = true;
	int k = 0;
	int ii = 0;

	double procent = 50.0;
	if (for_new == 1)
	{
		procent = 1.0;
	}

	while(izmen == true)
	//for(int ii = 0; ii < 50; ii++)
	{
		ii++;

		if (ii > 80)
		{
			procent = 0.1;
		}

		if (for_new == 0) // динамическое изменение коэффициента
		{
			if (ii > 80)
			{
				procent = 0.1;
			}
			else if (ii > 50)
			{
				procent = 1.0;
			}
			else if (ii > 20)
			{
				procent = 3.0;
			}
			else if (ii > 13)
			{
				procent = 10.0;
			}
			else if (ii > 3)
			{
				procent = 20.0;
			}
		}



		izmen = false;
		k = 0;

		for (auto& i : this->All_Luch)
		{

			if (i->type == "A_Luch" || i->type == "A2_Luch")
			{
				// da1
				auto a1 = i->get_yzel_near_opor(1, -1);
				auto a2 = i->Yzels_opor[1];

				auto b1 = i->get_yzel_near_opor(1, this->geo->M11);
				auto b2 = i->get_yzel_near_opor(1, this->geo->M11 + 1);

				double d1 = fabs(a1->func_R(0) - a2->func_R(0));
				double d2 = fabs(b1->func_R(0) - b2->func_R(0));

				macros3(da1, 5.0);

				// da2
				b1 = i->Yzels_opor[2];
				b2 = i->get_yzel_near_opor(2, -1);

				d2 = fabs(b1->func_R(0) - b2->func_R(0));

				macros3(da2, 5.0);

				//da3
				if (i->parameters.find("da3") == i->parameters.end() || i->parameters["da3"] > 0.05)
				{
					b1 = i->Yzels_opor[2];
					b2 = i->get_yzel_near_opor(2, 1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da3, 5.0);
				}

				// da4
				if (i->parameters.find("da4") == i->parameters.end() || i->parameters["da4"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, -1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da4, 5.0);
				}

				//da5
				if (i->parameters.find("da5") == i->parameters.end() || i->parameters["da5"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, 1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da5, 5.0);
				}

			}

			if (i->type == "B_Luch")
			{
				// ba1
				auto a1 = i->get_yzel_near_opor(1, -1);
				auto a2 = i->Yzels_opor[1];

				auto b1 = i->get_yzel_near_opor(1, this->geo->M11);
				auto b2 = i->get_yzel_near_opor(1, this->geo->M11 + 1);

				double d1 = Yzel_distance(a1, a2, 0);
				double d2 = Yzel_distance_x(b1, b2, 0);

				macros3(ba1, 5.0);

				// ba2
				b1 = i->Yzels_opor[2];
				b2 = i->get_yzel_near_opor(2, -1);

				d2 = Yzel_distance(b1, b2, 0);

				macros3(ba2, 3.0);

				// ba3
				if (i->parameters.find("ba3") == i->parameters.end() || i->parameters["ba3"] > 0.05)
				{
					a1 = i->get_yzel_near_opor(1, -1);
					a2 = i->Yzels_opor[1];

					b1 = i->Yzels_opor[2];
					b2 = i->get_yzel_near_opor(2, 1);

					d1 = Yzel_distance(a1, a2, 0);
					d2 = Yzel_distance_x(b1, b2, 0);

					macros3(ba3, 1.0);
				}
				// ba4
				if (i->parameters.find("ba4") == i->parameters.end() || i->parameters["ba4"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, -1);

					d2 = Yzel_distance(b1, b2, 0);

					macros3(ba4, 1.0);
				}
			}
		}

		// для B лучей за TS
		for (int i = 0; i < this->B_Luch.size(); i++)
		{
			auto A = this->A_Luch[i].back();
			auto B = this->B_Luch[i][0];
			auto a1 = A->get_yzel_near_opor(3, 1);
			auto a2 = A->Yzels_opor[3];

			auto b1 = B->Yzels_opor[3];
			auto b2 = B->get_yzel_near_opor(3, 1);

			double d1 = Yzel_distance(a1, a2, 0);
			double d2 = Yzel_distance(b1, b2, 0);

			if ((100.0 - d2 * 100.0 / d1) > 1.0)
			{
				k++;
				izmen = true;
				for (auto& kk : this->B_Luch[i])
				{
					if (kk->parameters.find("ba5") != kk->parameters.end())
					{
						kk->parameters["ba5"] *= (1.0 + procent / 100.0);
					}
					else
					{
						kk->parameters["ba5"] = this->geo->ba5 * (1.0 + procent / 100.0);
					}
				}
			}
			else if((100.0 - d2 * 100.0 / d1) < -1.0)
			{
				k++;
				izmen = true;
				for (auto& kk : this->B_Luch[i])
				{
					if (kk->parameters.find("ba5") != kk->parameters.end())
					{
						kk->parameters["ba5"] *= (1.0 - procent / 100.0);
					}
					else
					{
						kk->parameters["ba5"] = this->geo->ba5 * (1.0 - procent / 100.0);
					}
				}
			}
		}

		// для E лучей
		for (int i = 0; i < this->E_Luch.size(); i++)
		{
			auto A = this->E_Luch[i][0];
			auto B = this->B_Luch[i].back();
			
			// ea1
			auto a1 = A->get_yzel_near_opor(1, 1);
			auto a2 = A->Yzels_opor[1];

			auto b1 = B->Yzels_opor[2];
			auto b2 = B->get_yzel_near_opor(2, 1);

			double d2 = Yzel_distance(a1, a2, 0);
			double d1 = Yzel_distance(b1, b2, 0);

			if (this->E_Luch[i][0]->parameters.find("ea1") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea1"] > 0.01)
			{
				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea1") != kk->parameters.end())
						{
							kk->parameters["ea1"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea1"] = this->geo->ea1 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea1") != kk->parameters.end())
						{
							kk->parameters["ea1"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea1"] = this->geo->ea1 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// ea2

			if (this->E_Luch[i][0]->parameters.find("ea2") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea2"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, -1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[3];
				b2 = B->get_yzel_near_opor(3, -1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea2") != kk->parameters.end())
						{
							kk->parameters["ea2"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea2"] = this->geo->ea2 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea2") != kk->parameters.end())
						{
							kk->parameters["ea2"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea2"] = this->geo->ea2 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// ea3

			if (this->E_Luch[i][0]->parameters.find("ea3") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea3"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, 1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[3];
				b2 = B->get_yzel_near_opor(3, 1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea3") != kk->parameters.end())
						{
							kk->parameters["ea3"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea3"] = this->geo->ea3 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea3") != kk->parameters.end())
						{
							kk->parameters["ea3"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea3"] = this->geo->ea3 * (1.0 - procent / 100.0);
						}
					}
				}
			}
		}

		// для D лучей
		for (int i = 0; i < this->D_Luch.size(); i++)
		{
			auto A = this->D_Luch[i][0];
			auto B = this->E_Luch[i].back();

			// md1
			auto a1 = A->get_yzel_near_opor(1, 1);
			auto a2 = A->Yzels_opor[1];

			auto b1 = B->Yzels_opor[1];
			auto b2 = B->get_yzel_near_opor(1, 1);

			double d2 = Yzel_distance(a1, a2, 0);
			double d1 = Yzel_distance(b1, b2, 0);

			if (this->D_Luch[i][0]->parameters.find("md1") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md1"] > 0.01)
			{
				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md1") != kk->parameters.end())
						{
							kk->parameters["md1"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md1"] = this->geo->md1 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md1") != kk->parameters.end())
						{
							kk->parameters["md1"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md1"] = this->geo->md1 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// md2

			if (this->D_Luch[i][0]->parameters.find("md2") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md2"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, -1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[2];
				b2 = B->get_yzel_near_opor(2, -1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md2") != kk->parameters.end())
						{
							kk->parameters["md2"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md2"] = this->geo->md2 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md2") != kk->parameters.end())
						{
							kk->parameters["md2"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md2"] = this->geo->md2 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// md3

			if (this->D_Luch[i][0]->parameters.find("md3") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md3"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, 1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[2];
				b2 = B->get_yzel_near_opor(2, 1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md3") != kk->parameters.end())
						{
							kk->parameters["md3"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md3"] = this->geo->md3 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md3") != kk->parameters.end())
						{
							kk->parameters["md3"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md3"] = this->geo->md3 * (1.0 - procent / 100.0);
						}
					}
				}
			}
		}

		// для G лучей 
		for (int ik = 0; ik < this->G_Luch.size(); ik++)
		{
				auto A = this->G_Luch[ik][0];
				auto B = this->B_Luch[ik].back();

				auto a1 = A->get_yzel_near_opor(2, -3);
				auto a2 = A->Yzels_opor[2];

				auto b1 = B->Yzels_opor[2];
				auto b2 = B->get_yzel_near_opor(2, -3);

				double d2 = Yzel_distance(a1, a2, 0);
				double d1 = Yzel_distance(b1, b2, 0);

				auto i = this->G_Luch[ik][0];

				if (i->parameters.find("dd2") == i->parameters.end() || i->parameters["dd2"] > 0.01)
				{
					if ((100.0 - d2 * 100.0 / d1) > 0.5)
					{
						k++;
						izmen = true;
						for (auto& kk : this->G_Luch[ik])
						{
							if (kk->parameters.find("dd2") != kk->parameters.end())
							{
								kk->parameters["dd2"] *= (1.0 + procent / 100.0);
							}
							else
							{
								kk->parameters["dd2"] = this->geo->dd2 * (1.0 + procent / 100.0);
							}
						}
					}
					else if ((100.0 - d2 * 100.0 / d1) < -0.5)
					{
						k++;
						izmen = true;
						for (auto& kk : this->G_Luch[ik])
						{
							if (kk->parameters.find("dd2") != kk->parameters.end())
							{
								kk->parameters["dd2"] *= (1.0 - procent / 100.0);
							}
							else
							{
								kk->parameters["dd2"] = this->geo->dd2 * (1.0 - procent / 100.0);
							}
						}
					}
				}
	

		}

		for (auto& i : this->All_Luch)
		{
			i->dvigenie(0);
		}

		//cout << "k = " << k << endl;
	}


	cout << "END: Izmenenie geo parameters" << endl;
}

void Setka::Winslow_method(void)
{
	int N, M;
	double R = 1.0;
	double** x;
	double** y;
	double** z;
	double** x2;
	double** y2;
	double xx, yy;
	double du, dv;
	double dxu, dxv, dyu, dyv, ddxuv, ddyuv, A, B, C, err;

	N = 16;
	M = 16;
	

	// Создаём массивы
	if (true)
	{
		x = new double* [N];
		for (int i = 0; i < N; i++) x[i] = new double[M];

		y = new double* [N];
		for (int i = 0; i < N; i++) y[i] = new double[M];

		x2 = new double* [N];
		for (int i = 0; i < N; i++) x2[i] = new double[M];

		y2 = new double* [N];
		for (int i = 0; i < N; i++) y2[i] = new double[M];

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				x[i][j] = 0.0;
				x2[i][j] = 0.0;
				y[i][j] = 0.0;
				y2[i][j] = 0.0;
			}
		}
	}

	// Инициализация массивов
	if (true)
	{
		for (int i = 0; i < N; i++)
		{
			xx = R * cos(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][0] = xx;
			x2[i][0] = xx;
			y[i][0] = yy;
			y2[i][0] = yy;
		}

		for (int i = 0; i < N; i++)
		{
			xx = R * cos(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][M - 1] = xx;
			x2[i][M - 1] = xx;
			y[i][M - 1] = yy;
			y2[i][M - 1] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			x[N - 1][j] = xx;
			x2[N - 1][j] = xx;
			y[N - 1][j] = yy;
			y2[N - 1][j] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			x[0][j] = xx;
			x2[0][j] = xx;
			y[0][j] = yy;
			y2[0][j] = yy;
		}
	}

	du = 1.0 / N;
	dv = 1.0 / M;

	do
	{
		err = 0.0;
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < M - 1; j++)
			{
				dxu = (x[i + 1][j] - x[i - 1][j]) / (2.0 * du);
				dxv = (x[i][j + 1] - x[i][j - 1]) / (2.0 * dv);
				dyu = (y[i + 1][j] - y[i - 1][j]) / (2.0 * du);
				dyv = (y[i][j + 1] - y[i][j - 1]) / (2.0 * dv);
				ddxuv = (x[i + 1][j + 1] + x[i - 1][j - 1] - x[i + 1][j - 1] - x[i - 1][j + 1]) / (4.0 * du * dv);
				ddyuv = (y[i + 1][j + 1] + y[i - 1][j - 1] - y[i + 1][j - 1] - y[i - 1][j + 1]) / (4.0 * du * dv);
				A = dxv * dxv + dyv * dyv;
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddxuv;
				C = dxu * dxu + dyu * dyu;
				/*if (i == 1 && j == 1)
				{
					cout << const_pi<< " " << dxu << " " << dxv << " " << A << " a " << B << " " << C << endl;
					cout << x[i + 1][j] << " " << x[i - 1][j] << " " << du << endl;
				}*/
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					x2[i][j] = (A * (x[i + 1][j] + x[i - 1][j]) / (du * du) - B + C * (x[i][j + 1] + x[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(x2[i][j] - x[i][j]));
				}
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddyuv;
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					y2[i][j] = (A * (y[i + 1][j] + y[i - 1][j]) / (du * du) - B + C * (y[i][j + 1] + y[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(y2[i][j] - y[i][j]));
				}
			}
		}
		z = x2;
		x2 = x;
		x = z;

		z = y2;
		y2 = y;
		y = z;

		//cout << err << " " << x[1][1] << " " << x2[1][1] << endl;
		//return;
	} while (err > 0.000001);

	for (int i = 0; i < int(N/4); i++)
	{
		for (int j = 0; j < int(M/4); j++)
		{
			x[i][j] = (x[i][j] - x[N - 1 - i][M - 1 - j] + x[i][M - 1 - j] - x[N - 1 - i][j]) / 4.0;
			x[i][M - 1 - j] = x[i][j];
			x[N - 1 - i][j] = -x[i][j];
			x[N - 1 - i][M - 1 - j] = -x[i][j];

			y[i][j] = (y[i][j] - y[N - 1 - i][M - 1 - j] - y[i][M - 1 - j] + y[N - 1 - i][j]) / 4.0;
			y[i][M - 1 - j] = -y[i][j];
			y[N - 1 - i][j] = y[i][j];
			y[N - 1 - i][M - 1 - j] = -y[i][j];
		}
	}

	ofstream fout;
	fout.open("winslow.txt");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//if( i%10 == 0 && j%10 == 0)  
			fout << x[i][j] << " " << y[i][j] << endl;
			/*else if( i == N - 1 && j % 10 == 0)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i % 10 == 0 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i == N - 1 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;*/
		}
	}
	fout.close();

}

void Setka::Tecplot_print_all_yzel_in_3D(string name)
{
	// name - это имя сетки
	int k = 1;
	ofstream fout;
	string name_f = "Tecplot_all_yzel_in_3D_" + name + ".txt";
	

	fout.open(name_f);



	fout << "TITLE = HP  VARIABLES = X, Y, Z"  << endl;

	for (auto& ii : this->Yzel_2D)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_2_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel_2)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_yzel_with_condition()
{
	// name - это имя сетки
	int k = 1;
	ofstream fout;
	string name_f = "Tecplot_all_yzel_with_condition.txt";


	fout.open(name_f);



	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;
	fout << "ZONE T = HP" << endl;

	for (auto& i : this->All_Yzel)
	{
		//if (i->type == Type_yzel::TS)
		//{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		//}
	}

	fout.close();
}

void Setka::Tecplot_print_krug_yzel_in_3D(int num)
{
	// num - какой круг выводим - головной или хвостовой
	int k = 1;
	vector <vector<Yzel*>> VVVV;
	ofstream fout;
	string name_f = "Tecplot_krug_yzel_in_3D_" + to_string(num) +  ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	VVVV = Krug_Yzel;
	if (num == 2) VVVV = Krug_Yzel_2;

	for (auto& ii : VVVV)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_opor_yzel_in_luchs_3D(string name)
{
	ofstream fout;
	string name_f = "Tecplot_opor_yzel_in_luchs_3D_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->All_Luch)
	{
		if (ii->type != name) continue;
		for (auto& i : ii->Yzels_opor)
		//for (auto& i : ii->Yzels)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
	}

	fout.close();
}
 
void Setka::Tecplot_print_all_lush_in_3D(string name)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_3D_" + name + ".txt";

	auto A = *All_name_luch[name];
	int k = 0;
	for (auto& j : A[0])
	{
		k += (j->Yzels.size() - 1);
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;

	
	for(auto& i : A)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i)
		{
			for (int m = 0; m < j->Yzels.size() - 1; m++)
			{
				fout << j->Yzels[m]->coord[0][0] << " " << j->Yzels[m]->coord[0][1] << " " << j->Yzels[m]->coord[0][2] << endl;
				fout << j->Yzels[m + 1]->coord[0][0] << " " << j->Yzels[m + 1]->coord[0][1] << " " << j->Yzels[m + 1]->coord[0][2] << endl;
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_lush_in_2D()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_2D.txt";

	int k = 0;

	for (auto& nn : this->name_luch)
	{
		auto A = *(this->All_name_luch[nn]);
		for (auto& j : A[0])
		{
			k += (j->Yzels.size() - 1);
		}
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;

	int kkk = 1;


	for (int i = 0; i < this->geo->Nphi; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& nn : this->name_luch)
		{
			auto A = *(this->All_name_luch[nn]);
			for (auto& j : A[i])
			{
				for (int m = 0; m < j->Yzels.size() - 1; m++)
				{
					fout << j->Yzels[m]->coord[0][0] << " " << sqrt(kv(j->Yzels[m]->coord[0][1]) + kv(j->Yzels[m]->coord[0][2])) << endl;
					fout << j->Yzels[m + 1]->coord[0][0] << " " << sqrt(kv(j->Yzels[m + 1]->coord[0][1]) + kv(j->Yzels[m + 1]->coord[0][2])) << endl;
				}
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_plane_lush(int plane)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_plane_lush_" + to_string(plane) + ".txt";

	int k = 0;

	for (auto& nn : this->name_luch)
	{
		auto A = *(this->All_name_luch[nn]);
		for (auto& j : A[0])
		{
			k += (j->Yzels.size() - 1);
		}
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;


	for (int i = plane; i <= plane; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG" << endl;
	
		for (auto& nn : this->name_luch)
		{
			auto A = *(this->All_name_luch[nn]);
			for (auto& j : A[i])
			{
				for (int m = 0; m < j->Yzels.size() - 1; m++)
				{
					fout << j->Yzels[m]->coord[0][0] << " " << sqrt(kv(j->Yzels[m]->coord[0][1]) + kv(j->Yzels[m]->coord[0][2])) << endl;
					fout << j->Yzels[m + 1]->coord[0][0] << " " << sqrt(kv(j->Yzels[m + 1]->coord[0][1]) + kv(j->Yzels[m + 1]->coord[0][2])) << endl;
				}
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_plane_surfase(int plane)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_plane_surfase_" + to_string(plane) + ".txt";

	int kTS = 0;
	kTS += this->A_Luch[0].size();
	kTS += this->B_Luch[0].size();
	kTS += this->C_Luch[0].size();

	int kHP = 0;
	kHP += this->A_Luch[0].size();
	kHP += this->B_Luch[0].size();
	kHP += this->E_Luch[0].size();
	kHP += this->D_Luch[0].size();

	int kBS = 0;
	kBS += this->A_Luch[0].size();


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;

	for (int i = plane; i <= plane; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << kTS + kHP + kBS << ", E = " << kTS + kHP + kBS - 3 << ", F=FEPOINT, ET=LINESEG" << endl;

		// TS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->C_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// HP
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->E_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->D_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// BS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[3]->coord[0][0];
			double y = i->Yzels_opor[3]->coord[0][1];
			double z = i->Yzels_opor[3]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}


		for (int m = 1; m < kTS; m++)
		{
			fout << m << " " << m + 1 << endl;
		}

		for (int m = 1; m < kHP; m++)
		{
			fout << kTS + m << " " << kTS + m + 1 << endl;
		}

		for (int m = 1; m < kBS; m++)
		{
			fout << kTS + kHP + m << " " << kTS + kHP + m + 1 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_All_surfase_in_2D()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_All_surfase_in_2D.txt";

	int kTS = 0;
	kTS += this->A_Luch[0].size();
	kTS += this->B_Luch[0].size();
	kTS += this->C_Luch[0].size();

	int kHP = 0;
	kHP += this->A_Luch[0].size();
	kHP += this->B_Luch[0].size();
	kHP += this->E_Luch[0].size();
	kHP += this->D_Luch[0].size();

	int kBS = 0;
	kBS += this->A_Luch[0].size();


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;
	int kkk = 0;

	for (int i = 0; i < this->geo->Nphi; i++)
	{
		kkk++;
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << kTS + kHP + kBS << ", E = " << kTS + kHP + kBS - 3 << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;

		// TS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->C_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// HP
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->E_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->D_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// BS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[3]->coord[0][0];
			double y = i->Yzels_opor[3]->coord[0][1];
			double z = i->Yzels_opor[3]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}


		for (int m = 1; m < kTS; m++)
		{
			fout << m << " " << m + 1 << endl;
		}

		for (int m = 1; m < kHP; m++)
		{
			fout << kTS + m << " " << kTS + m + 1 << endl;
		}

		for (int m = 1; m < kBS; m++)
		{
			fout << kTS + kHP + m << " " << kTS + kHP + m + 1 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_cell_in_3D()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_all_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;


	for (auto& i : this->All_Cell)
	{
		fout << "ZONE T=HP, N = " << i->yzels.size() << ", E = " << i->yzels.size() - 1 << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}

		for (int m = 0; m < i->yzels.size() - 1; m++)
		{
			fout << m + 1 << " " << m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_cut_plane_parameters(const Eigen::Vector3d& A,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2)
{
	Eigen::Vector3d B, C;
	Eigen::Vector3d B1, B2;

	Eigen::Vector3d n;
	n = v1.cross(v2);

	bool a1, a2;

	// Строим уравнение плоскости



	for (auto& cell : this->All_Cell)
	{
		a1 = false;
		a2 = false;
		// Сначала надо понять разрезает ли плоскость ячейку
		for (auto& yz : cell->yzels)
		{
			B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			C = B - A;
			if (n.dot(C) >= 0.0)
			{
				a1 = true;
			}
			else
			{
				a2 = true;
			}
		}

		if (a1 == false || a2 == false) continue; // Ячейка не пересекается плоскостью

		// Далее надо найти все стороны в пересечении
		for (auto& gran : cell->grans)
		{
			int ii;
			for (int i = 0; i < gran->yzels.size(); i++)
			{
				ii = i + 1;
				if (ii >= gran->yzels.size()) ii = 0;
				auto Y1 = gran->yzels[i];
				auto Y2 = gran->yzels[ii];
				a1 = false;
				a2 = false;

				// проверим находятся ли эти узлы по разную сторону от плоскости
				B1 << Y1->coord[0][0], Y1->coord[0][1], Y1->coord[0][2];
				C = B1 - A;
				if (n.dot(C) >= 0.0)
				{
					a1 = true;
				}
				else
				{
					a2 = true;
				}

				B2 << Y1->coord[0][0], Y1->coord[0][1], Y1->coord[0][2];
				C = B2 - A;
				if (n.dot(C) >= 0.0)
				{
					a1 = true;
				}
				else
				{
					a2 = true;
				}

				if (a1 == false || a2 == false) continue;

			}
		}
	}
}

void Setka::Tecplot_print_cell_plane_parameters()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_cell_plane_parameters.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y";
	for (auto& nn : this->phys_param->param_names)
	{
		fout << ", " << nn;
	}
	fout << ", my_zone, Mach";
	fout << endl;

	int kkk = 1;


	for (auto& i : this->Cell_2D[0])
	{
		fout << i->center[0][0] << " " << norm2(0.0, i->center[0][1], i->center[0][2]);
		for (auto& nn : this->phys_param->param_names)
		{
			if (i->parameters[0].find(nn) != i->parameters[0].end())
			{
				if (nn != "Q")
				{
					fout << " " << i->parameters[0][nn];
				}
				else
				{
					fout << " " << i->parameters[0]["Q"]/ i->parameters[0]["rho"];
				}
			}
			else
			{
				fout << " " << 0.0;
			}
		}

		if (i->type == Type_cell::Zone_1)
		{
			fout << " " << 1.0;
		}
		else if (i->type == Type_cell::Zone_2) 
		{
			fout << " " << 2.0;
		}
		else if (i->type == Type_cell::Zone_3)
		{
			fout << " " << 3.0;
		}
		else if (i->type == Type_cell::Zone_4)
		{
			fout << " " << 4.0;
		}
		else
		{
			fout << " " << 0.0;
		}

		fout << " " << norm2(i->parameters[0]["Vx"], i->parameters[0]["Vy"], i->parameters[0]["Vz"])/
			sqrt(this->phys_param->gamma * i->parameters[0]["p"] / i->parameters[0]["rho"]);

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_gran_in_cell()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_all_gran_in_3D.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (int i = 0; i < 1000; i++) // this->All_Gran.size()
	{
		auto C = this->All_Cell[i];

		fout << "ZONE T=HP, N = " << C->grans.size() * 4 << ", E = " << C->grans.size() << ", F=FEPOINT, ET=quadrilateral, SOLUTIONTIME = " << i + 1 << endl;

		for (auto& A : C->grans)
		{
			for (auto& j : A->yzels)
			{
				fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
			}
		}

		for (int k = 0; k < C->grans.size(); k++)
		{
			fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
		}
	}
}

void Setka::Tecplot_print_all_gran_in_surface(string name_surf)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_all_gran_in_surface_" + name_surf + ".txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	vector<Gran*>* C = nullptr;

	if (name_surf == "TS")
	{
		C = &this->Gran_TS;
	}
	else if (name_surf == "HP")
	{
		C = &this->Gran_HP;
	}
	else if (name_surf == "BS")
	{
		C = &this->Gran_BS;
	}
	else
	{
		cout << "Error  9785656643" << endl;
	}


	fout << "ZONE T=HP, N = " << C->size() * 4 << ", E = " << C->size() << ", F=FEPOINT, ET=quadrilateral" <<  endl;

	for (auto& A : *C)
	{
		for (auto& j : A->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}
	}

	for (int k = 0; k < C->size(); k++)
	{
		fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
	}
	
}


void Setka::Tecplot_print_gran_with_condition(short int iii)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_gran_with_condition_" + to_string(iii) + "_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int k = 0; 
	int kk = 0;

	for (auto& i : this->MK_Grans[iii])
	{
		k++;
		kk += i->yzels.size();
	}

	fout << "ZONE T=HP, N = " << kk << ", E = " << k << ", F=FEPOINT, ET=quadrilateral" << endl;

	for (auto& C : this->MK_Grans[iii])
	{
		for (auto& j : C->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}
	}
	
	k = 1;
	for (auto& C : this->MK_Grans[iii])
	{
		for (int i = k; i < k + C->yzels.size(); i++)
		{
			fout << i << " ";
		}
		fout << endl;
		k += C->yzels.size();
	}
}

void Setka::Tecplot_print_1D(Interpol* Int1, const Eigen::Vector3d& Origin, 
	const Eigen::Vector3d& vec, string name, const double& leng)
{
	cout << "Start: Tecplot_Tecplot_print_1D" + name + ".txt" << endl;
	ofstream fout;
	string name_f = "Tecplot_Tecplot_print_1D" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = R, X, Y, Z";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", BB_8Pi, rho_Th, T_Th";
	if (std::find(Int1->param_names.begin(), Int1->param_names.end(), "rho_Pui_1") != Int1->param_names.end())
	{
		fout << ", T_Pui_1";
	}
	if (std::find(Int1->param_names.begin(), Int1->param_names.end(), "rho_Pui_2") != Int1->param_names.end())
	{
		fout << ", T_Pui_2";
	}
	fout << endl;

	unsigned int N = 5000;
	Eigen::Vector3d C;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;

	for (size_t i = 0; i < N; i++)
	{
		C = Origin + i * vec / N * leng;

		//cout << "Int1->Get_param A" << endl;
		fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
		//cout << "Int1->Get_param B" << endl;
		if (fine_int == false) continue;
		prev_cell = next_cell;

		fout << 1.0 * i / N * leng << " " << C(0) << " " << C(1) << " " << C(2);
		for (auto& nam : Int1->param_names)
		{
			if (nam == "Q")
			{
				fout << " " << parameters["Q"] / parameters["rho"];
			}
			else if (nam == "rho_Pui_2")
			{
				if (parameters[nam] < 1e-7)
				{
					fout << " " << 0.0;
				}
				else
				{
					fout << " " << parameters[nam];
				}
			}
			else
			{
				fout << " " << parameters[nam];
			}
		}

		unordered_map<string, double> param;
		short int zone = 1;
		if (parameters["Q"] / parameters["rho"] < 50.0)
		{
			zone = 2;
		}
		else
		{
			zone = 3;
		}
		this->phys_param->Plasma_components(zone, parameters, param);

		if (param["rho_Th"] < 1e-8) param["rho_Th"] = 1e-8;

		fout << " " << kvv(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi);
		fout << " " << param["rho_Th"];
		fout << " " << param["T_Th"];
		if (std::find(Int1->param_names.begin(), Int1->param_names.end(), "rho_Pui_1") != Int1->param_names.end())
		{
			fout << " " << 2.0 * parameters["p_Pui_1"] / parameters["rho_Pui_1"];
		}
		if (std::find(Int1->param_names.begin(), Int1->param_names.end(), "rho_Pui_2") != Int1->param_names.end())
		{
			if (parameters["rho_Pui_2"] < 1e-6)
			{
				fout << " " << 0.0;
			}
			else
			{
				fout << " " << 2.0 * parameters["p_Pui_2"] / parameters["rho_Pui_2"];
			}
			
		}
		fout << endl;

	}

	fout.close();

}

void Setka::Tecplot_print_2D(Interpol* Int1, const double& a, 
	const double& b, const double& c, const double& d, string name, bool razmer, bool moove_system)
{
	// Находим нормаль к плоскости
	cout << "Start: Tecplot_print_2D " << name << endl;
	std::array<double, 3> normal;
	normal[0] = a;
	normal[1] = b;
	normal[2] = c;

	double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	if (length > 0)
	{
		normal[0] /= length;
		normal[1] /= length;
		normal[2] /= length;
	}

	// Выбираем произвольное направление для сортировки (например, ось OX в плоскости)
	std::array<double, 3> reference_dir;
	if (std::abs(normal[0]) > 0.9)
	{ // Если нормаль близка к OX, выбираем OY
		reference_dir = { 0.0, 1.0, 0.0 };
	}
	else
	{
		reference_dir = { 1.0, 0.0, 0.0 };
	}

	// Находим вектор в плоскости, перпендикулярный нормали
	std::array<double, 3>  tangent_dir = {
		reference_dir[1] * normal[2] - reference_dir[2] * normal[1],
		reference_dir[2] * normal[0] - reference_dir[0] * normal[2],
		reference_dir[0] * normal[1] - reference_dir[1] * normal[0]
	};



	std::vector< std::array<short int, 2> > rebro;
	rebro.resize(12);
	rebro[0][0] = 1; rebro[0][1] = 2;
	rebro[1][0] = 2; rebro[1][1] = 3;
	rebro[2][0] = 3; rebro[2][1] = 4;
	rebro[3][0] = 4; rebro[3][1] = 1;
	rebro[4][0] = 5; rebro[4][1] = 6;
	rebro[5][0] = 6; rebro[5][1] = 7;
	rebro[6][0] = 7; rebro[6][1] = 8;
	rebro[7][0] = 8; rebro[7][1] = 5;
	rebro[8][0] = 1; rebro[8][1] = 5;
	rebro[9][0] = 2; rebro[9][1] = 6;
	rebro[10][0] = 3; rebro[10][1] = 7;
	rebro[11][0] = 4; rebro[11][1] = 8;

	Yzel* A1;
	Yzel* A2;

	std::vector < std::vector< std::array<double, 3> >> all_setka;

	bool aa;

	cout << "Tecplot_print_2D: srez" << endl;

	// Находим срез сетки
	for (const auto& cell : this->All_Cell)
	{
		std::vector< std::array<double, 3> > all_point;
		std::array<double, 3> P1;
		std::array<double, 3> P2;
		std::array<double, 3> outIntersection;

		for (size_t i = 0; i < 12; i++) // Пробегаемся по рёбрам
		{
			A1 = cell->yzels[rebro[i][0] - 1];
			A2 = cell->yzels[rebro[i][1] - 1];

			P1 = { A1->coord[0][0], A1->coord[0][1], A1->coord[0][2] };
			P2 = { A2->coord[0][0], A2->coord[0][1], A2->coord[0][2] };

			//cout << "A" << endl;
			aa = findIntersection(P1, P2, a, b, c, d, outIntersection);
			//cout << "B" << endl;
			if (aa == true)
			{
				all_point.push_back(outIntersection);
			}
			else
			{
				continue;
			}

			if (fabs(a * outIntersection[0] + b * outIntersection[1] + c * outIntersection[2] + d) > 0.001)
			{
				cout << "Error  798646575467845456" << endl;
				cout << outIntersection[0] << endl;
				cout << outIntersection[1] << endl;
				cout << outIntersection[2] << endl;
				cout << fabs(a * outIntersection[0] + b * outIntersection[1] + c * outIntersection[2] + d) << endl;
				exit(-1);
			}

			//cout << "C" << endl;
		}

		//cout << "________" << endl;

		if (all_point.size() == 0) continue;
		
		if (all_point.size() < 3)
		{
			cout << "Error 9867531090    " << all_point.size() << endl;
			continue;
		}


		// std::vector< std::array<double, 3> > all_point;
		// Сейчас тут хранятся все найденные точки, которые надо рассортировать по кругу

		std::array<double, 3> centroid = { 0.0, 0.0, 0.0 };
		for (const auto& p : all_point) {
			centroid[0] += p[0];
			centroid[1] += p[1];
			centroid[2] += p[2];
		}
		centroid[0] /= all_point.size();
		centroid[1] /= all_point.size();
		centroid[2] /= all_point.size();

		/*for (const auto& i : all_point)
		{
			cout << i[0] << " " << i[1] << " " << i[2] << endl;
		}
		cout << "___" << endl;*/


		// Сортируем точки по углу относительно tangent_dir
		std::sort(all_point.begin(), all_point.end(), [centroid, &tangent_dir]
		(std::array<double, 3> aa, std::array<double, 3> bb)
			{
				std::array<double, 3> vec_a = { aa[0] - centroid[0], aa[1] - centroid[1], aa[2] - centroid[2] };
				std::array<double, 3> vec_b = { bb[0] - centroid[0], bb[1] - centroid[1], bb[2] - centroid[2] };;

				// Угол между vec_a и tangent_dir
				double dot_a = vec_a[0] * tangent_dir[0] + vec_a[1] * tangent_dir[1] + vec_a[2] * tangent_dir[2];
				double cross_a =
					tangent_dir[1] * vec_a[2] - tangent_dir[2] * vec_a[1] -
					tangent_dir[0] * vec_a[2] + tangent_dir[2] * vec_a[0] +
					tangent_dir[0] * vec_a[1] - tangent_dir[1] * vec_a[0];

				double angle_a = std::atan2(cross_a, dot_a);

				// Угол между vec_b и tangent_dir
				double dot_b = vec_b[0] * tangent_dir[0] + vec_b[1] * tangent_dir[1] + vec_b[2] * tangent_dir[2];
				double cross_b =
					tangent_dir[1] * vec_b[2] - tangent_dir[2] * vec_b[1] -
					tangent_dir[0] * vec_b[2] + tangent_dir[2] * vec_b[0] +
					tangent_dir[0] * vec_b[1] - tangent_dir[1] * vec_b[0];

				double angle_b = std::atan2(cross_b, dot_b);

				return angle_a < angle_b;
			});

		/*for (const auto& i : all_point)
		{
			cout << i[0] << " " << i[1] << " " << i[2] << endl;
		}
		cout << endl;

		exit(-1);*/


		all_setka.push_back(all_point);
	}

	unsigned int NN = 0;
	unsigned int NN2 = 0;

	for (const auto& i : all_setka)
	{
		NN += i.size() + 1;
		NN2 += i.size();
	}

	
	cout << "Tecplot_print_2D: print" << endl;

	// Рисуем саму сетку
	ofstream fout;
	string name_f = "Tecplot_setka_srez_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z" << endl;
	fout << "ZONE T=HP, NODES = " << NN2 << ", ELEMENTS = " << NN2 << ", F = FEPOINT, ET = LINESEG" << endl;

	Eigen::Vector3d C;
	for (const auto& i : all_setka)
	{
		for (const auto& j : i)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			fout << C(0) << " " << C(1) << " " << C(2) << endl;
		}
	}


	size_t all_k1 = 1;
	for (const auto& i : all_setka)
	{
		size_t k1 = i.size();
		for (size_t ii = 0; ii < k1; ii++)
		{
			size_t k2 = ii + 1;
			if (k2 >= k1) k2 = 0;
			fout << all_k1 + ii << " " << all_k1 + k2 << endl;
		}

		all_k1 = all_k1 + k1;
	}

	fout.close();

	// Рисуем поля параметров
	name_f = "Tecplot_Tecplot_print_2D_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", |V|, Mach, BB_8pi, rho_Th, p_Th, T_Th";
	fout << endl;

	fout << "ZONE T=HP, ";

	fout << "NODES = " << NN << ", ELEMENTS = " << NN2 << ", F = FEPOINT, ET = TRIANGLE" << endl;

	Eigen::Vector3d C2;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;

	for (const auto& i : all_setka)
	{
		std::vector< std::array<double, 3> > kj(i);
		C2 << 0.0, 0.0, 0.0;

		for (const auto& j : i)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			C2 += C;
		}
		C2 /= (1.0 * i.size());

		kj.push_back({ C2[0], C2[1], C2[2] });

		bool visible = true;

		for (const auto& j : kj)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];

			Eigen::Vector3d koord = C;
			if (moove_system == false)
			{
				double alpha = this->phys_param->T_all * this->phys_param->Omega;
				Eigen::Vector3d point(0.0, 0.0, this->phys_param->Omega);
				Eigen::Vector3d rotate_velosity = // Угловая скорость (как вектор)
					rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * point;
				koord = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * C;
			}

			fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);

			if (fine_int == true) next_cell = prev_cell;
			if (fine_int == false)
			{
				visible = false;
				break;
			}
		}

		for (const auto& j : kj)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];

			Eigen::Vector3d koord = C;
			Eigen::Vector3d rotate_velosity;
			Eigen::Vector3d point;
			double alpha;
			if (moove_system == false)
			{
				double alpha = this->phys_param->T_all * this->phys_param->Omega;
				Eigen::Vector3d point(0.0, 0.0, this->phys_param->Omega);
				Eigen::Vector3d rotate_velosity = // Угловая скорость (как вектор)
					rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * point;
				koord = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * C;
			}

			fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);

			if (fine_int == true) next_cell = prev_cell;

			if (visible == false)
			{
				fout << 0.0 << " " << 0.0 << " " << 0.0;
			}
			else
			{
				double kk = 1.0;
				if (razmer == true) kk = this->phys_param->Get_razmer("r");
				fout << koord(0) * kk << " " << koord(1) * kk << " " << koord(2) * kk;
			}

			if (moove_system == false)
			{
				Eigen::Vector3d VV(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
				Eigen::Vector3d VV1;

				koord(2) = 0.0;
				VV1 = 
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * VV + point.cross(koord);

				parameters["Vx"] = VV1[0];
				parameters["Vy"] = VV1[1];
				parameters["Vz"] = VV1[2];

			}


			for (auto& nam : Int1->param_names)
			{
				if (fine_int == true && visible == true)
				{
					if (nam != "Q")
					{
						double kk = 1.0;
						if (razmer == true) kk = this->phys_param->Get_razmer(nam);

						if (std::isnan(parameters[nam]) || std::fpclassify(parameters[nam]) == FP_SUBNORMAL)
						{
							fout << " " << 0.0;
						}
						else
						{
							fout << " " << parameters[nam] * kk;
						}
					}
					else
					{
						fout << " " << parameters["Q"] / parameters["rho"];
					}
				}
				else
				{
					fout << " " << 0.0;
				}
			}

			short int zone;
			if (parameters["Q"] / parameters["rho"] < 50)
			{
				if (sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]) /
					sqrt(this->phys_param->gamma * parameters["p"]) > 2)
				{
					zone = 1;
				}
				else
				{
					zone = 2;
				}
			}
			else
			{
				if (sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]) /
					sqrt(this->phys_param->gamma * parameters["p"]) > 1)
				{
					zone = 4;
				}
				else
				{
					zone = 3;
				}
			}


			double rho_Th;
			double p_Th;
			double T_Th;

			unordered_map<string, double> param;

			this->phys_param->Plasma_components(zone, parameters, param);

			rho_Th = param["rho_Th"];
			p_Th = param["p_Th"];
			T_Th = param["T_Th"];

			double kT = 1.0;
			if (razmer == true) kT = this->phys_param->Get_razmer("T");

			double kp = 1.0;
			if (razmer == true) kp = this->phys_param->Get_razmer("p");

			double krho = 1.0;
			if (razmer == true) krho = this->phys_param->Get_razmer("rho");

			double VV = norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
			double Mach = sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]) /
				sqrt(this->phys_param->gamma * parameters["p"]);
			if (parameters["rho"] <= 0.0) Mach = 0.0;

			fout << " " << VV << " " << Mach << " "
				<< norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi) << 
				" " << rho_Th * krho  << " " << p_Th * kp << " " << T_Th * kT;


			fout << endl;
		}

		kj.clear();
	}

	size_t all_k = 1;
	for (const auto& i : all_setka)
	{
		size_t k1 = i.size();
		for (size_t ii = 0; ii < k1; ii++)
		{
			size_t k2 = ii + 1;
			if (k2 >= k1) k2 = 0;
			fout << all_k + k1 << " " << all_k + ii << " " << all_k + k2 << endl;
		}

		all_k = all_k + k1 + 1;
	}

	fout.close();

	// Отдельно находим пересечение плоскости с поверхностями разрыва
	
	for (short int ik = 0; ik < 3; ik++)
	{
		std::vector< std::array<double, 3> > all_point_surf;
		vector<Gran*>* AA = nullptr;
		if (ik == 0) AA = &this->Gran_HP;
		if (ik == 1) AA = &this->Gran_TS;
		if (ik == 2) AA = &this->Gran_BS;

		for (const auto& gr : *AA)
		{
			std::vector< std::array<double, 3> > all_point;
			std::array<double, 3> P1;
			std::array<double, 3> P2;
			std::array<double, 3> outIntersection;

			for (size_t i = 0; i < 4; i++) // Пробегаемся по рёбрам
			{
				size_t ii = i + 1;
				if (ii >= 4) ii = 0;
				A1 = gr->yzels[i];
				A2 = gr->yzels[ii];

				P1 = { A1->coord[0][0], A1->coord[0][1], A1->coord[0][2] };
				P2 = { A2->coord[0][0], A2->coord[0][1], A2->coord[0][2] };

				//cout << "A" << endl;
				aa = findIntersection(P1, P2, a, b, c, d, outIntersection);
				//cout << "B" << endl;
				if (aa == true)
				{
					all_point.push_back(outIntersection);
				}
				else
				{
					continue;
				}
			}


			if (all_point.size() != 2) continue;

			for (auto& ii : all_point)
			{
				all_point_surf.push_back(ii);
			}
		}

		if (ik == 0) name_f = "HP_srez_" + name + ".txt";
		if (ik == 1) name_f = "TS_srez_" + name + ".txt";
		if (ik == 2) name_f = "BS_srez_" + name + ".txt";

		fout.open(name_f);
		fout << "TITLE = HP" << endl;
		fout << "VARIABLES = X, Y" << endl;
		fout << "ZONE T=HP, NODES = " << all_point_surf.size() << ", ELEMENTS = " << all_point_surf.size() / 2 << ", F = FEPOINT, ET = LINESEG" << endl;

		for (auto& ii : all_point_surf)
		{
			fout << ii[0] << " " << ii[1] << endl;
		}

		for (size_t ii = 0; ii < all_point_surf.size() / 2; ii++)
		{
			fout << 2 * ii + 1 << " " << 2 * ii + 2 << endl;
		}

		fout.close();
	}


	cout << "End: Tecplot_print_2D " << name << endl;

}

void Setka::Tecplot_print_2D_dekard(Interpol* Int1, Eigen::Vector3d V1, Eigen::Vector3d V2,
	double L1, double R1, double L2, double R2, string name, bool moove_system)
{
	string name_f = "Tecplot_Tecplot_print_dekard_2D_" + name + ".txt";

	ofstream fout;
	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", |V|, Mach, BB_8pi, V_perenos_x, V_perenos_y, V_perenos_z";
	fout << endl;

	Eigen::Vector3d C2;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;


	Eigen::Vector3d X;

	for (double u = L1; u < R1; u = u + (R1 - L1) / 1000.0)
	{
		for (double v = L2; v < R2; v = v + (R2 - L2) / 1000.0)
		{
			X = u * V1 + v * V2;
			Eigen::Vector3d koord = X;
			Eigen::Vector3d rotate_velosity;
			Eigen::Vector3d point(0.0, 0.0, this->phys_param->Omega);
			double alpha;

			if (moove_system == false)
			{
				alpha = this->phys_param->T_all * this->phys_param->Omega;
				rotate_velosity = // Угловая скорость (как вектор)
					 rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * point;
				koord = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * X;
			}

			fine_int = Int1->Get_param(koord(0), koord(1), koord(2), parameters, prev_cell, next_cell);
			if (fine_int == false) continue;

			fout << X(0) << " " << X(1) << " " << X(2);

			if (moove_system == false)
			{
				Eigen::Vector3d VV(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
				Eigen::Vector3d VV1;

				X(2) = 0.0;
				VV1 =
					 rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * VV + point.cross(X);

				parameters["Vx"] = VV1[0];
				parameters["Vy"] = VV1[1];
				parameters["Vz"] = VV1[2];

				VV << parameters["Bx"], parameters["By"], parameters["Bz"];
				VV1 =
					rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * VV;
				parameters["Bx"] = VV1[0];
				parameters["By"] = VV1[1];
				parameters["Bz"] = VV1[2];
			}


			for (auto& nam : Int1->param_names)
			{
				fout << " " << parameters[nam];
			}

			short int zone = 1;
			double rho_Th;
			double p_Th;
			double T_Th;

			rho_Th = 0.0;
			p_Th = 0.0;
			T_Th = 0.0;

			double kT = 1.0;
			double kp = 1.0;
			double krho = 1.0;

			double VV = norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
			double Mach = sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]) /
				sqrt(this->phys_param->gamma * parameters["p"]);
			if (parameters["rho"] <= 0.0) Mach = 0.0;

			fout << " " << VV << " " << Mach << " "
				<< norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi) <<
				" " << point.cross(X)[0] << " " << point.cross(X)[1] << " " << point.cross(X)[2];


			fout << endl;

		}
	}

	fout.close();

	cout << "End: Tecplot_print_2D " << name << endl;

}


void Setka::Tecplot_print_spherik(Interpol* Int1, double R, int Nphi, int Ntheta, string name, bool moove_system)
{
	string name_f = "Tecplot_Tecplot_print_dekard_2D_" + name + ".txt";

	ofstream fout;
	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = phi, the";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", |V|, Mach, BB_8pi, V_perenos_x, V_perenos_y, V_perenos_z";
	fout << endl;

	Eigen::Vector3d C2;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;


	Eigen::Vector3d X;
	

	for (double phi = 0.0; phi < 2.0 * const_pi; phi = phi + 2.0 * const_pi/ Nphi)
	{
		for (double the = 0.0; the <= const_pi; the = the + const_pi / Ntheta)
		{
			X(0) = R * sin(the) * cos(phi);
			X(1) = R * sin(the) * sin(phi);
			X(2) = R * cos(the);
			Eigen::Vector3d koord = X;
			Eigen::Vector3d rotate_velosity;
			Eigen::Vector3d point(0.0, 0.0, this->phys_param->Omega);
			double alpha;

			if (moove_system == false)
			{
				alpha = this->phys_param->T_all * this->phys_param->Omega;
				rotate_velosity = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * point;
				koord = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * X;
			}

			fine_int = Int1->Get_param(koord(0), koord(1), koord(2), parameters, prev_cell, next_cell);
			if (fine_int == false) continue;

			fout << phi << " " << the;

			if (moove_system == false)
			{
				Eigen::Vector3d VV(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
				Eigen::Vector3d VV1;

				X(2) = 0.0;
				VV1 =
					rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * VV + point.cross(X);

				parameters["Vx"] = VV1[0];
				parameters["Vy"] = VV1[1];
				parameters["Vz"] = VV1[2];

				VV << parameters["Bx"], parameters["By"], parameters["Bz"];
				VV1 =
					rotationMatrixZ(alpha) * rotationMatrixX(this->phys_param->lambda) * VV;
				parameters["Bx"] = VV1[0];
				parameters["By"] = VV1[1];
				parameters["Bz"] = VV1[2];
			}


			for (auto& nam : Int1->param_names)
			{
				fout << " " << parameters[nam];
			}

			short int zone = 1;
			double rho_Th;
			double p_Th;
			double T_Th;

			rho_Th = 0.0;
			p_Th = 0.0;
			T_Th = 0.0;

			double kT = 1.0;
			double kp = 1.0;
			double krho = 1.0;

			double VV = norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]);
			double Mach = sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"]) /
				sqrt(this->phys_param->gamma * parameters["p"]);
			if (parameters["rho"] <= 0.0) Mach = 0.0;

			fout << " " << VV << " " << Mach << " "
				<< norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi) <<
				" " << point.cross(X)[0] << " " << point.cross(X)[1] << " " << point.cross(X)[2];


			fout << endl;

		}
	}

	fout.close();

	cout << "End: Tecplot_print_2D " << name << endl;

}



void Setka::Tecplot_print_2D_sphere(int SS, string name, bool razmer, bool moove_system)
{
	ofstream fout;
	string name_f = "Tecplot_Tecplot_print_2D_" + name + ".txt";
	fout.open(name_f);

	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = phi, the";

	for (auto& nam : this->phys_param->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", |V|, BB_8pi";
	fout << endl;

	fout << "ZONE T=HP, ";

	Eigen::Vector3d C2;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;

	int i, j, k;
	k = SS;
	// Получаем размеры массива
	const size_t* dims = this->Cell_regular.shape();

	for (i = 0; i < dims[0]; ++i)
	{
		for (j = 0; j < dims[1]; ++j)
		{
			auto cell = this->Cell_regular[i][j][k];
			Eigen::Vector3d C(cell->center[0][0], cell->center[0][1], cell->center[0][2]);
			//Eigen::Vector3d C0 = C;
			Eigen::Vector3d koord;
			Eigen::Vector3d rotate_velosity;
			Eigen::Vector3d point(0.0, 0.0, this->phys_param->Omega);
			double alpha;
			if (moove_system == false)
			{
				alpha = this->phys_param->T_all * this->phys_param->Omega;
				rotate_velosity = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * point;
				koord = // Угловая скорость (как вектор)
					rotationMatrixX(-this->phys_param->lambda) * rotationMatrixZ(-alpha) * C;

				C = koord;
			}


			double x = C[0];
			double y = C[1];
			double z = C[2];

			double r_1 = sqrt(x * x + y * y + z * z);
			double the_1 = acos(z / r_1);
			double phi_1 = polar_angle(x, y);

			std::array<double, 5> phi;
			std::array<double, 5> the;

			phi[0] = phi_1;
			the[0] = the_1;

			phi[1] = phi_1 + const_pi * 2;
			the[1] = the_1;

			phi[2] = phi_1 - const_pi * 2;
			the[2] = the_1;

			phi[3] = phi_1;
			the[3] = 2 * const_pi - the_1;

			phi[4] = phi_1;
			the[4] = -the_1;

			Eigen::Vector3d VV1(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);
			if (moove_system == false)
			{
				Eigen::Vector3d VV(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);

				VV1 = rotationMatrixZ(-alpha) *
					rotationMatrixX(-this->phys_param->lambda) * VV + point.cross(C);
			}

			
			for (int ii = 0; ii < 5; ii++)
			{
				fout << phi[ii] << " " << the[ii];

				for (auto& nam : this->phys_param->param_names)
				{
					if (nam == "Vx")
					{
						fout << " " << VV1[0];
					}
					else if (nam == "Vy")
					{
						fout << " " << VV1[1];
					}
					else if (nam == "Vz")
					{
						fout << " " << VV1[2];
					}
					else
					{
						fout << " " << cell->parameters[0][nam];
					}
				}

				double VV = norm2(VV1[0], VV1[1], VV1[2]);
				double BB = kvv(cell->parameters[0]["Bx"], cell->parameters[0]["By"],
					cell->parameters[0]["Bz"]) / 8.0 / const_pi;
				fout << " " << VV << " " << BB << endl;
			}

		}
	}

	
	for (j = 0; j < this->Cell_layer_head[k].size(); ++j)
	{
		auto cell = this->Cell_layer_head[k][j];
		Eigen::Vector3d C(cell->center[0][0], cell->center[0][1], cell->center[0][2]);
		Eigen::Vector3d koord = C;
		Eigen::Vector3d rotate_velosity;
		Eigen::Vector3d point;
		double alpha;
		if (moove_system == false)
		{
			alpha = this->phys_param->T_all * this->phys_param->Omega;
			point << 0.0, 0.0, this->phys_param->Omega;
			rotate_velosity = rotationMatrixX(this->phys_param->lambda) * // Угловая скорость (как вектор)
				rotationMatrixZ(alpha) * point;
			koord = rotationMatrixZ(-alpha) * rotationMatrixX(-this->phys_param->lambda) * C;

			C = koord;
		}


		double x = C[0];
		double y = C[1];
		double z = C[2];

		double r_1 = sqrt(x * x + y * y + z * z);
		double the_1 = acos(z / r_1);
		double phi_1 = polar_angle(x, y);

		std::array<double, 5> phi;
		std::array<double, 5> the;

		phi[0] = phi_1;
		the[0] = the_1;

		phi[1] = phi_1 + const_pi * 2;
		the[1] = the_1;

		phi[2] = phi_1 - const_pi * 2;
		the[2] = the_1;

		phi[3] = phi_1;
		the[3] = const_pi - the_1;

		phi[4] = phi_1;
		the[4] = -the_1;

		Eigen::Vector3d VV1(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);
		if (moove_system == false)
		{
			Eigen::Vector3d VV(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);

			VV1 = rotationMatrixZ(-alpha) *
				rotationMatrixX(-this->phys_param->lambda) * VV + point.cross(C);
		}


		for (int ii = 0; ii < 5; ii++)
		{
			fout << phi[ii] << " " << the[ii];


			for (auto& nam : this->phys_param->param_names)
			{
				if (nam == "Vx")
				{
					fout << " " << VV1[0];
				}
				else if (nam == "Vy")
				{
					fout << " " << VV1[1];
				}
				else if (nam == "Vz")
				{
					fout << " " << VV1[2];
				}
				else
				{
					fout << " " << cell->parameters[0][nam];
				}
			}

			double VV = norm2(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);
			double BB = kvv(cell->parameters[0]["Bx"], cell->parameters[0]["By"],
				cell->parameters[0]["Bz"]) / 8.0 / const_pi;
			fout << " " << VV << " " << BB << endl;
		}

	}

	for (j = 0; j < this->Cell_layer_tail[k].size(); ++j)
	{
		auto cell = this->Cell_layer_tail[k][j];
		Eigen::Vector3d C(cell->center[0][0], cell->center[0][1], cell->center[0][2]);
		Eigen::Vector3d koord = C;
		Eigen::Vector3d rotate_velosity;
		Eigen::Vector3d point;
		double alpha;
		if (moove_system == false)
		{
			alpha = this->phys_param->T_all * this->phys_param->Omega;
			point << 0.0, 0.0, this->phys_param->Omega;
			rotate_velosity = rotationMatrixX(this->phys_param->lambda) * // Угловая скорость (как вектор)
				rotationMatrixZ(alpha) * point;
			koord = rotationMatrixZ(-alpha) * rotationMatrixX(-this->phys_param->lambda) * C;

			C = koord;
		}


		double x = C[0];
		double y = C[1];
		double z = C[2];

		double r_1 = sqrt(x * x + y * y + z * z);
		double the_1 = acos(z / r_1);
		double phi_1 = polar_angle(x, y);

		std::array<double, 5> phi;
		std::array<double, 5> the;

		phi[0] = phi_1;
		the[0] = the_1;

		phi[1] = phi_1 + const_pi * 2;
		the[1] = the_1;

		phi[2] = phi_1 - const_pi * 2;
		the[2] = the_1;

		phi[3] = phi_1;
		the[3] = const_pi - the_1;

		phi[4] = phi_1;
		the[4] = -the_1;

		Eigen::Vector3d VV1(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);
		if (moove_system == false)
		{
			Eigen::Vector3d VV(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);

			VV1 = rotationMatrixZ(-alpha) *
				rotationMatrixX(-this->phys_param->lambda) * VV + point.cross(C);
		}


		for (int ii = 0; ii < 5; ii++)
		{
			fout << phi[ii] << " " << the[ii];



			for (auto& nam : this->phys_param->param_names)
			{
				if (nam == "Vx")
				{
					fout << " " << VV1[0];
				}
				else if (nam == "Vy")
				{
					fout << " " << VV1[1];
				}
				else if (nam == "Vz")
				{
					fout << " " << VV1[2];
				}
				else
				{
					fout << " " << cell->parameters[0][nam];
				}
			}

			double VV = norm2(cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"]);
			double BB = kvv(cell->parameters[0]["Bx"], cell->parameters[0]["By"],
				cell->parameters[0]["Bz"]) / 8.0 / const_pi;
			fout << " " << VV << " " << BB << endl;
		}

	}



	fout.close();


}


void Setka::Tecplot_print_plane_interpolation(Eigen::Vector3d A, Eigen::Vector3d v1, 
	Eigen::Vector3d v2, int l1, int r1, int l2, int r2)
{
	cout << "Start: Tecplot_print_plane_interpolation" << endl;
	ofstream fout;
	string name_f = "Tecplot_print_plane_interpolation.txt";
	Eigen::Vector3d P;
	Cell* previos = nullptr;
	Cell* cel = nullptr;
	unordered_map<string, double> par;

	auto names = this->phys_param->param_names;

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y";

	for (auto& nam : names)
	{
		fout << ", " << nam;
	}

	fout << endl;


	for (int i = l1; i <= r1; i++)
	{
		for (int j = l2; j <= r2; j++)
		{
			//cout << "i - " << i << "    j - " << j << endl;
			par.clear();
			P = A + i * v1 + j * v2;

			if (P.norm() < this->geo->R0 - 0.0000001) continue;
			if(P(0) < this->geo->L7 - 0.000001) continue;
			if(norm2(0.0, P(1), P(2)) > this->geo->R5 + 0.000001) continue;
			if(P.norm() > this->geo->R5 + 0.000001 && P(0) > 0.00000001) continue;


			cel = this->Find_cell_point(P(0), P(1), P(2), 0, previos);
			previos = cel;

			if (cel != nullptr)
			{
				//cel->Get_RBF_interpolation(P(0), P(1), P(2), par);
				cel->Get_IDW_interpolation(P(0), P(1), P(2), par, this->phys_param);

				fout << i * v1.norm() << " " << j * v2.norm();
				for (auto& nam : names)
				{
					fout << " " << par[nam];
				}
				fout << endl;
			}
			else
			{
				//cout << "NO = " << P(0) << " " << P(1) << " " << P(2) << endl;
			}
		}
	}

	cout << "End: Tecplot_print_plane_interpolation" << endl;
}

struct PairHash {
	size_t operator()(const std::pair<Yzel*, Yzel*>& p) const {
		// Комбинируем хеши указателей
		return std::hash<Yzel*>()(p.first) ^ (std::hash<Yzel*>()(p.second) << 1);
	}
};
// Возможны коллизии (например, пары (A, B) и (B, A) дадут одинаковый хеш, 
// если ^ используется без дополнительных преобразований).

void Setka::Edges_create(void)
{
	cout << "Start: Edges_create" << endl;
	this->Renumerate();

	Edge* Edge_;
	std::unordered_map<std::pair<Yzel*, Yzel*>, Edge*, PairHash>  map;


	for (auto& gr: this->All_Gran)
	{
		for (short int i = 0; i < gr->yzels.size(); i++)
		{
			short int ii = i + 1;
			if (ii >= gr->yzels.size()) ii = 0;

			std::pair<Yzel*, Yzel*> p1 = { gr->yzels[i] , gr->yzels[ii] };
			std::pair<Yzel*, Yzel*> p2 = { gr->yzels[ii] , gr->yzels[i] };

			if (map.find(p1) != map.end() || map.find(p2) != map.end())
			{
				Edge* AA;
				if (map.find(p1) != map.end())
				{
					AA = map[p1];
				}
				else
				{
					AA = map[p2];
				}

				AA->grans.push_back(gr);
				gr->edges.push_back(AA);
			}
			else
			{
				Edge_ = new Edge();
				Edge_->grans.push_back(gr);
				Edge_->A = gr->yzels[i];
				Edge_->B = gr->yzels[ii];
				map[p1] = Edge_;
				this->All_Edges.push_back(Edge_);
				gr->edges.push_back(Edge_);
			}
		}
	}

	unsigned int nn = 0;

	// Теперь добавляем ячейки в рёбра
	for (auto& ed : this->All_Edges)
	{
		for (const auto& gr : ed->grans)
		{
			if (std::find(ed->cells.begin(),
				ed->cells.end(), gr->cells[0]) ==
				ed->cells.end())
			{
				ed->cells.push_back(gr->cells[0]);
			}

			if (gr->cells.size() >= 2)
			{
				if (std::find(ed->cells.begin(),
					ed->cells.end(), gr->cells[1]) ==
					ed->cells.end())
				{
					ed->cells.push_back(gr->cells[1]);
				}
			}
		}

		nn += ed->cells.size();
	}

	cout << "V srednem = " << (1.0 * nn) / this->All_Edges.size() << endl;



	cout << "End: Edges_create" << endl;
}