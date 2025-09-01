#include "Setka.h"


bool findSphereIntersectionTime(
	const Eigen::Vector3d& X,  // Положение частицы
	const Eigen::Vector3d& V,  // Скорость частицы
	const double& R,        // Радиус сферы
	double& time) 
{
	const double a = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
	const double b = 2.0 * (X[0] * V[0] + X[1] * V[1] + X[2] * V[2]);
	const double c = (X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) - R * R;

	const double D = b * b - 4 * a * c;

	if (D < 0) 
	{
		return false;  // Нет пересечений
	}

	const double sqrtD = std::sqrt(D);
	const double t1 = (-b - sqrtD) / (2 * a);
	const double t2 = (-b + sqrtD) / (2 * a);

	// Находим минимальное положительное время
	if (t1 >= 0)
	{
		time = t1;
		return true;
	}
	if (t2 >= 0)
	{
		time = t2;
		return true;
	}

	return false;  // Оба времени отрицательные (пересечение было в прошлом)
}

double Velosity_1(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 2.0 * cp / sqrtpi_ + 2.0 * u * u / (3.0 * cp * sqrtpi_) - u * u * u * u / (15.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp / sqrtpi_ + (u + kv(cp) / (2.0 * u)) * erf(u / cp);
	}
}

double Velosity_2(const double& u, const double& cp)  // Считает на совсем скорость, а только её числитель (см. статью)
{
	if (u < 0.00001)
	{
		return (8.0 / 3.0) * kv(cp) * kv(cp) * const_pi * u + 
			(8.0 / 15.0) * kv(cp) * const_pi * u * u * u -
			(4.0 / 105.0) * const_pi * kv(u) * kv(u) * u;
	}
	else
	{
		return  cp * cp * cp * const_pi * (exp(-u * u / kv(cp)) * cp * u * 2.0 * (kv(cp) +
			2.0 * kv(u)) +//
			sqrtpi_ * (4.0 * kv(u) * kv(u) + 
				4.0 * cp * cp * kv(u) - kv(cp) * kv(cp)) * erf(u / cp)) / (4.0 * u * u);
	}
}

double Velosity_3(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 8.0 * cp / (3.0 * sqrtpi_) + 8.0 * u * u / (9.0 * cp * sqrtpi_) - 
			44.0 * u * u * u * u / (135.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp * (5.0 * kv(cp) + 2.0 * kv(u)) / 
			(sqrtpi_ * (3.0 * kv(cp) + 2.0 * kv(u))) +//
			(4.0 * kv(u) * kv(u) + 12.0 * cp * cp * kv(u) + 3.0 * kv(cp) * kv(cp)) * 
			erf(u / cp) / (2.0 * u * (3.0 * kv(cp) + 2.0 * kv(u)));
	}
}

// Генерация случайного числа в диапазоне [-scale, +scale]
double randomNoise(double scale) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(-scale, scale);
	return dist(gen);
}

// Функция, которая отклоняет вектор, сохраняя его длину
void perturbVectorKeepingMagnitude(double& Vx, double& Vy, double& Vz, const double& noiseScale = 0.01) {
	// 1. Вычисляем исходную длину
	const double magnitude = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
	if (magnitude == 0.0) return;  // нулевой вектор нельзя отклонить

	// 2. Нормализуем вектор (делаем единичным)
	const double invMag = 1.0 / magnitude;
	Vx *= invMag;
	Vy *= invMag;
	Vz *= invMag;

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise(noiseScale);
	double noiseY = randomNoise(noiseScale);
	double noiseZ = randomNoise(noiseScale);

	// 4. Делаем шум строго перпендикулярным исходному вектору (чтобы не менять длину)
	const double dot = Vx * noiseX + Vy * noiseY + Vz * noiseZ;
	noiseX -= dot * Vx;
	noiseY -= dot * Vy;
	noiseZ -= dot * Vz;

	// 5. Добавляем шум и нормализуем
	Vx += noiseX;
	Vy += noiseY;
	Vz += noiseZ;

	const double newMag = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
	const double correction = magnitude / newMag;
	Vx *= correction;
	Vy *= correction;
	Vz *= correction;
}

// Функция, которая отклоняет вектор
void perturbVector(double& Vx, double& Vy, double& Vz, double noiseScale = 0.01) {

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise(noiseScale);
	double noiseY = randomNoise(noiseScale);
	double noiseZ = randomNoise(noiseScale);

	// 5. Добавляем шум и нормализуем
	Vx += noiseX;
	Vy += noiseY;
	Vz += noiseZ;
}

void Setka::Set_MK_Zone(void)
{
	cout << "Start Set_MK_Zone" << endl;
	this->Renumerate();
	Eigen::Vector3d Centr;

	this->Cell_Center->MK_zone_r = 1;
	this->Cell_Center->MK_zone_phi = 0;

	// Задаём зону для каждой ячейки
	for (auto& cell : this->All_Cell)
	{
		Centr[0] = cell->center[0][0];
		Centr[1] = cell->center[0][1];
		Centr[2] = cell->center[0][2];

		if (cell->type == Type_cell::Zone_1)
		{
			cell->MK_zone_r = 1;
			cell->MK_zone = 1;
			cell->MK_zone_phi = 0;
		}
		else if (cell->type == Type_cell::Zone_2)
		{
			cell->MK_zone_r = 2;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 2;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 3;
			}
		}
		else if (cell->type == Type_cell::Zone_3)
		{
			cell->MK_zone_r = 3;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 4;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 5;
			}
		}
		else if (cell->type == Type_cell::Zone_4)
		{
			cell->MK_zone_r = 4;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 6;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 7;
			}
		}
		else
		{
			cout << "Error 86743207564" << endl;
			exit(-1);
		}
	}

	this->MK_Grans.resize(7);
	this->MK_Potoks.resize(7);
	for (short int i = 0; i < 7; ++i) this->MK_Potoks[i] = 0.0;

	// 1 зона
	for (auto& gr : this->Gran_TS)
	{
		this->MK_Grans[0].push_back(gr);
		gr->MK_type.push_back(1);
	}

	// 2 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] < -0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP)
		{
			this->MK_Grans[1].push_back(gr);
			gr->MK_type.push_back(2);
		}


		if (gr->type2 == Type_Gran_surf::TS)
		{
			this->MK_Grans[1].push_back(gr);
			gr->MK_type.push_back(2);
		}
		

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2)
			{
				this->MK_Grans[1].push_back(gr);
				gr->MK_type.push_back(2);
			}
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[1].push_back(gr);
					gr->MK_type.push_back(2);
				}
			}
		}
	}

	// 3 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] > 0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP)
		{
			this->MK_Grans[2].push_back(gr);
			gr->MK_type.push_back(3);
		}

		if (gr->type2 == Type_Gran_surf::TS)
		{
			this->MK_Grans[2].push_back(gr);
			gr->MK_type.push_back(3);
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2)
			{
				this->MK_Grans[2].push_back(gr);
				gr->MK_type.push_back(3);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) 
			{
				this->MK_Grans[2].push_back(gr);
				gr->MK_type.push_back(3);
			}

			
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[2].push_back(gr);
					gr->MK_type.push_back(3);
				}
			}
		}
	}

	// 4 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if (Centr[0] > 0)
			{
				this->MK_Grans[3].push_back(gr);
				gr->MK_type.push_back(4);
			}
		}

		if (gr->type2 == Type_Gran_surf::BS)
		{
			this->MK_Grans[3].push_back(gr);
			gr->MK_type.push_back(4);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[3].push_back(gr);
					gr->MK_type.push_back(4);
				}
			}

		}
	}

	// 5 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if (Centr[0] < 0.0)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}


			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[4].push_back(gr);
					gr->MK_type.push_back(5);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0)
				{
					this->MK_Grans[4].push_back(gr);
					gr->MK_type.push_back(5);
				}
			}

		}
	}

	// 6 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::BS)
		{
			this->MK_Grans[5].push_back(gr);
			gr->MK_type.push_back(6);
		}

		if (gr->type == Type_Gran::Outer_Hard)
		{
			if (Centr[0] > 0)
			{
				this->MK_Grans[5].push_back(gr);
				gr->MK_type.push_back(6);
			}
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[5].push_back(gr);
					gr->MK_type.push_back(6);
				}
			}

		}
	}

	// 7 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type != Type_Gran::Us && Centr[0] < 0)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4)
			{
				this->MK_Grans[6].push_back(gr);
				gr->MK_type.push_back(7);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0)
				{
					this->MK_Grans[6].push_back(gr);
					gr->MK_type.push_back(7);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[6].push_back(gr);
					gr->MK_type.push_back(7);
				}
			}

		}


	}

	for (size_t jj = 0; jj < 7; jj++)
	{
		cout << "MK_grans:  " << jj << "   size = " << this->MK_Grans[jj].size() << endl;
		unsigned int N_ = 0;
		for (auto& gr : this->MK_Grans[jj])
		{
			if (gr->type == Type_Gran::Us) N_++;
		}
		cout << "Only active:  " << N_ << endl;

	}

	// Проверяем, что в массивах нет повторов
	if (true)
	{
		for (size_t jj = 0; jj < 7; jj++)
		{
			for (auto& i : this->All_Gran)
			{
				i->work1 = false;
			}

			for (auto& i : this->MK_Grans[jj])
			{
				if (i->work1 == false)
				{
					i->work1 = true;
				}
				else
				{
					cout << "ERROR 6435856408" << endl;
					cout << jj << endl;
					exit(-1);
				}
			}
		}
			

	}





	cout << "END Set_MK_Zone" << endl;
}

void Setka::MK_prepare(short int zone_MK)
{
	cout << "Start MK_prepare   zone_MK = " << zone_MK << endl;
	// zone_MK должно начинаться с единицы
	if (zone_MK == 0)
	{
		cout << "Error 2341238507" << endl;
		exit(-1);
	}

	this->Renumerate();

	// Блок загрузки датчиков случайных чисел
	if (true)
	{
		ifstream fin2;
		fin2.open("rnd_my.txt");
		if (fin2.is_open() == false)
		{
			cout << "ERROR open  rnd_my.txt " << endl;
			exit(-100);
		}
		double a1, b1, c;
		for (int i = 0; i < 1021; i++)
		{
			fin2 >> a1 >> b1 >> c;
			auto s = new Sensor(a1, b1, c);
			this->Sensors.push_back(s);
		}
		fin2.close();
	}

	if (this->MK_Grans.size() < zone_MK || this->MK_Grans[zone_MK - 1].size() == 0)
	{
		cout << "Error 2341963846 " << this->MK_Grans.size() << endl;
		exit(-1);
	}

	// Готовим/загружаем AMR сетку для граней
	if (true)
	{
		cout << "Start: Zagruzka AMR" << endl;
		unsigned int NN1_ = 0;
		unsigned int NN2_ = 0;

		unsigned short int NNall = 0;
		double S = 0.0;

#pragma omp parallel for schedule(dynamic)
		for(size_t ijk = 0; ijk < this->MK_Grans[zone_MK - 1].size(); ijk++)
		//for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			auto gr = this->MK_Grans[zone_MK - 1][ijk];
			gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			gr->MK_Potok = 0.0;
			if (gr->AMR.size() < this->phys_param->num_H)
			{
				gr->AMR.resize(this->phys_param->num_H);
			}

			short int ni = 1;  // Определяем выходящую функцию распределения
			if (gr->cells[0]->MK_zone == zone_MK) ni = 0;

			short int ni2 = 0; // Определяем входящую функцию распределения
			if (gr->cells[0]->MK_zone == zone_MK) ni2 = 1;


			for(short int ii = 0; ii <= 1; ii++)
			{
				for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
				{
					if (ni == ii) gr->Read_AMR(ii, iH, true && this->phys_param->refine_AMR);
					if (ni2 == ii) gr->Read_AMR(ii, iH, false && this->phys_param->refine_AMR);

					if (false)
					{
						if (gr->AMR[iH - 1][ii] == nullptr)
						{
							string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
								to_string(iH) + "_" + to_string(gr->number) + ".bin";
							if (file_exists("data_AMR/" + name_f) && gr->type == Type_Gran::Us)
							{
								//cout << "Exist  " << iH << " " << ii << endl;
								// В этом случае просто считываем AMR - сетку
								gr->AMR[iH - 1][ii] = new AMR_f();
								gr->AMR[iH - 1][ii]->AMR_self = gr->AMR[iH - 1][ii];
								gr->AMR[iH - 1][ii]->Read("data_AMR/" + name_f);


								if (ni == ii && this->phys_param->refine_AMR == true && gr->type == Type_Gran::Us &&
									gr->AMR[iH - 1][ii]->Size() < 3000)
								{
									unsigned short int NN = 1;
									NN = gr->AMR[iH - 1][ii]->Refine();
									NNall += NN;
								}

								if (gr->AMR[iH - 1][ii]->Size() > 5000)
								{
									cout << "_____________________________" << endl;
									cout << "Anomalious function  H = " << iH << "  ii = " << ii << endl;
									cout << gr->center[0][0] << " " << gr->center[0][1] << " " <<
										gr->center[0][2] << " " << endl;
									cout << gr->AMR[iH - 1][ii]->Size() << endl;
									cout << "_____________________________" << endl;
								}
							}
							else
							{
								// В этом случае создаём новую AMR - сетку
								//cout << "ih = " << iH << " " << ii << endl;
								if (gr->type == Type_Gran::Us)
								{
									gr->AMR[iH - 1][ii] = new AMR_f(0.0, 20.0, -20.0, 20.0,
										-20.0, 20.0, 3, 6, 6);
									gr->AMR[iH - 1][ii]->AMR_self = gr->AMR[iH - 1][ii];
								}
								else
								{
									gr->AMR[iH - 1][ii] = new AMR_f(0.0, 20.0, -20.0, 20.0,
										-20.0, 20.0, 1, 1, 1);
									gr->AMR[iH - 1][ii]->AMR_self = gr->AMR[iH - 1][ii];
								}
							}

							// На всякий случай задаём нормаль
							if (ii == 0)
							{
								gr->AMR[iH - 1][ii]->Vn[0] = gr->normal[0][0];
								gr->AMR[iH - 1][ii]->Vn[1] = gr->normal[0][1];
								gr->AMR[iH - 1][ii]->Vn[2] = gr->normal[0][2];
							}
							else
							{
								gr->AMR[iH - 1][ii]->Vn[0] = -gr->normal[0][0];
								gr->AMR[iH - 1][ii]->Vn[1] = -gr->normal[0][1];
								gr->AMR[iH - 1][ii]->Vn[2] = -gr->normal[0][2];
							}
							gr->AMR[iH - 1][ii]->Set_bazis();

							// Заполняем параметры на AMR
							gr->AMR[iH - 1][ii]->parameters["n"] = 0.0;
							gr->AMR[iH - 1][ii]->parameters["nn"] = 0.0;
							gr->AMR[iH - 1][ii]->parameters["Smu"] = 0.0;
						}
					}
				}
			}

			if (gr->type == Type_Gran::Us)
			{
				for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
				{
					short int ni = 1;  // определяем выходящую функцию распределения
					if (gr->cells[0]->MK_zone == zone_MK)
					{
						ni = 0;
					}
					unsigned int njn = gr->AMR[iH - 1][ni]->Size();
					#pragma omp critical (erfgwerwe) 
					{
						NN1_ += njn;
						NN2_++;
					}
				}
			}

			// Считаем сразу поток атомов через грань
			if (true)
			{
				if (gr->type == Type_Gran::Us)
				{
					for (auto& ai : gr->AMR)
					{
						ai[ni2]->Culk_SpotokV(gr->area[0]);

						#pragma omp critical (erfgwerweS) 
						{
							S += ai[ni2]->SpotokV;
						}
						gr->MK_Potok += ai[ni2]->SpotokV;
					}
				}
				else // Вручную посчитаем поток с границы расчётной области
				{
					for (size_t ijk = 0; ijk < this->phys_param->num_H; ijk++)
					{
						gr->AMR[ijk][ni2]->SpotokV = 0.0;
					}

					Eigen::Vector3d n;
					n << -gr->AMR[3][ni2]->Vn[0], -gr->AMR[3][ni2]->Vn[1], -gr->AMR[3][ni2]->Vn[2];
					// Так как нормаль должна быть внешняя к грани
					double sjv = Get_Spotok_inf(n);
					gr->AMR[3][ni2]->SpotokV = sjv * gr->area[0];
					#pragma omp critical (erfgwerweS) 
					{
						S += gr->AMR[3][ni2]->SpotokV;
					}
					gr->MK_Potok += gr->AMR[3][ni2]->SpotokV;
				}
			}

			// Удалим сразу входящие функции распределения
			for (short int iH = 1; iH <= gr->AMR.size(); iH++)
			{
				gr->AMR[iH - 1][ni2]->Delete();
			}

		}
		cout << "Izmelcheno  " << NNall << "  yacheek" << endl;
		cout << "Srednee chislo yacheek v AMR na vixodnix granyax =  " << 1.0 * NN1_/ NN2_ << endl;
		cout << "End: Zagruzka AMR" << endl;
		this->MK_Potoks[zone_MK - 1] = S; // Входящий поток через всю границу зоны
	}

	//pause_seconds(15);

	// Можно удалить файлы граничных граней, так как они только место занимают
	if (true)
	{
		cout << "Start: Ydalenie lishnix failov" << endl;
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			if (gr->type != Type_Gran::Us)
			{
				for (short int ii = 0; ii <= 1; ii++)
				{
					for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
					{
						if (gr->AMR[iH - 1][ii] == nullptr)
						{
							string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
								to_string(iH) + "_" + to_string(gr->number) + ".bin";
							if (file_exists("data_AMR/" + name_f))
							{
								std::filesystem::remove(name_f);
							}
						}
					}
				}
			}
		}
		cout << "END: Ydalenie lishnix failov" << endl;
	}

	// Заполняем граничные условия для AMR сетки (на границе расчётной области)
	// В том случае, если среди граней зоны есть граничные грани, на которых надо задавать максвелл
	// можно не задавать, так как я сделал аналитическое разыгрывание скорости
	if (false)
	{
		//for (auto& gr : this->All_boundary_Gran)
#pragma omp parallel for schedule(dynamic)
		//for (auto& gr : this->MK_Grans[zone_MK - 1])
		for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
		{
			auto& gr = this->MK_Grans[zone_MK - 1][idx];
			if (gr->type == Type_Gran::Outer_Hard || gr->type == Type_Gran::Outer_Soft)
			//if(true)  // тестирование, задаём максвел для всех граней
			//if (gr->type2 == Type_Gran_surf::BS)
			{
				if (gr->AMR.size() != 0)
				{
					short int ni = 0;  // определяем входящую функцию распределения
					if (gr->cells[0]->MK_zone == zone_MK)
					{
						ni = 1;
					}

					/*gr->AMR[3][ni]->Print_all_center_Tecplot(gr->AMR[3][ni], "1");
					gr->AMR[3][ni]->Culk_SpotokV(gr->area[0]);
					cout << "1 POTOK = " << gr->AMR[3][ni]->SpotokV/ gr->area[0] << endl;
					cout << gr->AMR[3][ni]->Vn[0] << " " << gr->AMR[3][ni]->Vn[1] << "  " <<
						gr->AMR[3][ni]->Vn[2] << endl;*/

					// Здесь надо задать граничные условия для четвёртого сорта 
					// а также помельчить сетку, если необходимо
					//cout << "_----------------------_" << endl;
					gr->AMR[3][ni]->Fill_maxwel_inf(this->phys_param->Velosity_inf);
					/*cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					exit(-1);*/



					// <<1>> - внутрь, <<3>> - 4 сорт водорода
					//gr->AMR[3][ni]->Print_all_center_Tecplot(gr->AMR[3][ni], "2");
					//gr->AMR[3][ni]->Culk_SpotokV(gr->area[0]);
					//cout << "2 POTOK = " << gr->AMR[3][ni]->SpotokV/ gr->area[0] << endl;
					//exit(-1);
				}
			}
		}
	}

	// Теперь для каждой функции распределения вычисляем поток через неё
	if (false)
	{
		double S = 0.0;
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			// Нужно вычислять поток только у входящей части функции распределения
			short int ni = 0; // Определяем входящую функцию
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 1;
			}
			gr->MK_Potok = 0.0;

			if (gr->type == Type_Gran::Us)
			{
				for (auto& ai : gr->AMR)
				{
					ai[ni]->Culk_SpotokV(gr->area[0]);

					S += ai[ni]->SpotokV;
					gr->MK_Potok += ai[ni]->SpotokV;
				}
			}
			else // Вручную посчитаем поток с границы расчётной области
			{
				gr->AMR[0][ni]->SpotokV = 0.0;
				gr->AMR[1][ni]->SpotokV = 0.0;
				gr->AMR[2][ni]->SpotokV = 0.0;
				Eigen::Vector3d n;
				n << -gr->AMR[3][ni]->Vn[0], -gr->AMR[3][ni]->Vn[1], -gr->AMR[3][ni]->Vn[2];
				// Так как нормаль должна быть внешняя к грани
				double sjv = Get_Spotok_inf(n);
				gr->AMR[3][ni]->SpotokV = sjv * gr->area[0];
				S += gr->AMR[3][ni]->SpotokV;
				gr->MK_Potok += gr->AMR[3][ni]->SpotokV;

				//cout << " Potok = " << sjv << endl;
				//exit(-1);
			}


		}
		this->MK_Potoks[zone_MK - 1] = S; // Входящий поток через всю границу зоны
	}

	// Считаем необходимые геометрические параметры для МК
	// И добавляем переменные в ячейки
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			i->Set_Cell_Geo_for_MK();
			if (i->MK_zone == zone_MK && this->phys_param->culc_cell_moments == true)
			//if(true)
			{
				for (const auto& nam : this->phys_param->MK_param)
				{
					i->parameters[0][nam] = 0.0;
				}
			}
		}

		for (auto& i : this->All_Gran)
		{
			i->Set_Gran_Geo_for_MK();
		}
	}

	// Обнулим функции распределения, в которые будем накапливать информацию
	if (true)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			// gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			short int ni = 1;
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}
			for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
			{
				auto funk = gr->AMR[iH - 1][ni];

				// Это для того, чтобы помельчить функцию как надо
				/*if (iH == 4)
				{
					funk->Fill_maxwel_inf(this->phys_param->Velosity_inf);
				}*/

				funk->Fill_null();
			}
			
		}
	}

	cout << "END MK_prepare   zone_MK = " << zone_MK << endl;
}

void Setka::MK_delete(short int zone_MK)
{
	cout << "Start MK_delete" << endl;
	// Блок удаления датчиков случайных чисел
	if (true)
	{
		for (auto& i : this->Sensors)
		{
			delete i;
		}
		this->Sensors.clear();
	}

	// Мельчим AMR сетку, если надо
	if (this->phys_param->de_refine_AMR == true)
	{
		cout << "start de_refine_AMR" << endl;
		unsigned int nmnm = 0;
		unsigned int num_ = 0;
		unsigned int sr_num = 0;
#pragma omp parallel for schedule(dynamic)
		//for (auto& gr : this->MK_Grans[zone_MK - 1])
		for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
		{
			auto& gr = this->MK_Grans[zone_MK - 1][idx];
			if (gr->type != Type_Gran::Us) continue;

			short int ni = 1;
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			for (short int iH = 1; iH <= gr->AMR.size(); iH++)
			{
				int iki = gr->AMR[iH - 1][ni]->de_Refine();

				#pragma omp critical (third) 
				{
					nmnm += iki;
					num_++;
					sr_num += gr->AMR[iH - 1][ni]->Size();
				}
					
			}
		}
		cout << "Ydaleno  " << nmnm << "  yacheek" << endl;
		cout << "Srednee chislo yacheek =  " << (1.0 * sr_num) / num_ << endl;
	}

	// Записываем AMR сетку для граней
	if (true)
	{
		// Нужно сохрянять только выходящие грани!
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			short int ni = 1;
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			for (short int ii = 0; ii <= 1; ii++) 
			{
				for (short int iH = 1; iH <= gr->AMR.size(); iH++)
				{
					if (ii == ni)
					{
						string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
							to_string(iH) + "_" + to_string(gr->number) + ".bin";

						if (this->phys_param->save_AMR == true)
						{
							gr->AMR[iH - 1][ii]->Save("data_AMR/" + name_f);
						}
						//cout << "Delete " << ii << " " << iH << " " << 
						//	gr->number << endl;
						gr->AMR[iH - 1][ii]->Delete();
						//cout << "Delete2" << endl;
						delete gr->AMR[iH - 1][ii];
					}
					else
					{
						gr->AMR[iH - 1][ii]->Delete();
						delete gr->AMR[iH - 1][ii];
					}
				}
			}
			gr->AMR.clear();
			gr->AMR.resize(0);
		}

		//this->MK_Grans.clear();  // Не надо удалять, чтобы сразу считать следующую зону
	}

	// Записываем моменты
	if (this->phys_param->culc_cell_moments == true)
	{
		if (file_exists(this->phys_param->MK_file))
		{
			this->Download_cell_MK_parameters(this->phys_param->MK_file, zone_MK);
		}
		this->Save_cell_MK_parameters(this->phys_param->MK_file);
	}

	cout << "END MK_delete" << endl;
}

void Setka::MK_go(short int zone_MK)
{
	auto start = std::chrono::high_resolution_clock::now();
	cout << "Start MK_go " << zone_MK << "   N_on_gran = " << this->phys_param->N_per_gran << endl;
	int N_on_gran = this->phys_param->N_per_gran;   // Сколько запускаем частиц на грань в среднем
	double mu_expect = 0.0;
	mu_expect = this->MK_Potoks[zone_MK - 1] / 
		(1.0 * N_on_gran * this->MK_Grans[zone_MK - 1].size());

	cout << "All potok = " << this->MK_Potoks[zone_MK - 1] << endl;
	//exit(-1);

	unsigned int ALL_N = 0;  // Общее число запущенных в итоге частиц
	unsigned int k1 = 0;
#pragma omp parallel for schedule(dynamic)
	//for (auto& gr : this->MK_Grans[zone_MK - 1])
	for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
	{
		auto& gr = this->MK_Grans[zone_MK - 1][idx];
		Eigen::Vector3d n;
		Eigen::Vector3d t;
		Eigen::Vector3d m;

		#pragma omp critical (first) 
		{
			k1++;
			if (k1 % 1000 == 0)
			{
				cout << "Gran = " << k1 << "    Iz: " << this->MK_Grans[zone_MK - 1].size() << endl;
			}
		}
		// Выбираем конкретный номер датчика случайных чисел
		unsigned int sens_num1 = 2 * omp_get_thread_num();
		unsigned int sens_num2 = 2 * omp_get_thread_num() + 1;

		//this->Sensors[sens_num1]->MakeRandom();
		//this->Sensors[sens_num1]->MakeRandom();
		//this->Sensors[sens_num2]->MakeRandom();

		short int ni = 0; // Номер "входящей" функции распределения
		if (gr->cells[0]->MK_zone == zone_MK)
		{
			ni = 1;
		}
		double full_gran_potok = gr->MK_Potok;


		if (full_gran_potok < 0.0000001 * MK_Potoks[zone_MK - 1]/ this->MK_Grans[zone_MK - 1].size())
		{
			continue;
		}

		// Получаем нормаль для граничных граней (это будет внешняя нормаль)
		if (gr->type != Type_Gran::Us)
		{
			n << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];
			get_bazis(n, t, m);
		}

		// Разыгрываем каждый сорт отдельно
		for (short int nh_ = 0; nh_ < this->phys_param->num_H; ++nh_)
		{
			auto& func = gr->AMR[nh_][ni];

			if (gr->type == Type_Gran::Us)
			{
				gr->Read_AMR(ni, nh_ + 1, false);
				func->Culk_SpotokV(gr->area[0]);
			}

			if (func->SpotokV < 0.0000001 * full_gran_potok)
			{
				// Функция распределния нулевая, можно не разыгрывать
				continue;
			}

			// Расчитываем число запускаемых частиц
			unsigned int N_particle = max(static_cast<int>(func->SpotokV / mu_expect) + 1, 
				min(N_on_gran, 100));
			double mu = func->SpotokV / N_particle; // Вес каждой частицы

			#pragma omp critical (second) 
			{
				ALL_N += N_particle;
			}

			// Запускаем каждую частицу
			for (unsigned int num = 0; num < N_particle; ++num)
			{
				MK_particle P = MK_particle();
				if (ni == 0)
				{
					P.cel = gr->cells[1];               // Ячейка в которой находится частица
				}
				else
				{
					P.cel = gr->cells[0];
				}

				if (P.cel->MK_zone != zone_MK)
				{
					cout << "Error 9767653421" << endl;
					exit(-1);
				}
				
				P.mu = mu;                           // Вес частицы
				P.sort = nh_ + 1;                    // Сорт частицы

				Eigen::Vector3d poz;

				// Находим положение точки на грани
				gr->Get_Random_pozition(poz, this->Sensors[sens_num1]);
				P.Addcoord(poz);

				// Находим скорость частицы
				if (gr->type != Type_Gran::Us)
				{
					// Можно разыгрывать аналитическую функцию распределения на границе, а не табличную
					this->Velosity_initial(this->Sensors[sens_num1], poz, n, t, m);
				}
				else
				{
					func->Get_random_velosity(func, gr->area[0], poz, this->Sensors[sens_num1]);
				}

				P.AddVel(poz);

				// Некоторые проверки разыгрынной скорости частицы
				if (this->regim_otladki)
				{
					if (ni == 0)
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) < 0.0)
						{
							cout << "Error  8765656431" << endl;
							whach((int)gr->type);
							whach(gr->normal[0][0]);
							whach(gr->normal[0][1]);
							whach(gr->normal[0][2]);
							whach(poz(0));
							whach(poz(1));
							whach(poz(2));
							whach(nh_);
							whach(num);
							whach(func->Vn[0]);
							whach(func->Vn[1]);
							whach(func->Vn[2]);
							exit(-1);
						}
					}
					else
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) > 0.0)
						{
							cout << "Error  7411100090" << endl;
							whach(scalarProductFast(poz(0), poz(1), poz(2),
								gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]));
							whach(poz(0));
							whach(poz(1));
							whach(poz(2));
							whach(gr->normal[0][0]);
							whach(gr->normal[0][1]);
							whach(gr->normal[0][2]);
							whach(gr->center[0][0]);
							whach(gr->center[0][1]);
							whach(gr->center[0][2]);
							for (int i = 0; i < gr->MK_type.size(); i++)
							{
								whach(i);
								whach(gr->MK_type[i]);
							}
							whach(gr->cells[0]->MK_zone);
							whach(gr->cells[1]->MK_zone);
							whach(func->Vn[0]);
							whach(func->Vn[1]);
							whach(func->Vn[2]);
							exit(-1);
						}
					}
				}

				// Надо проверить, (1) что точка находится в нужной ячейке и
				// (2) что она будет находиться в ней через время dt
				if (true)
				{
					Cell* previos = P.cel;
					Eigen::Vector3d Center_cell, Move;
					Center_cell << P.cel->center[0][0], P.cel->center[0][1], P.cel->center[0][2];
					double dt = previos->geo_parameters["l_size"] / P.Vel_norm() / 1000.0;

					Cell* ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);

					short unsigned int klk = 0;
					while (ppp == nullptr || ppp != P.cel)
					{
						klk++;
						if (klk > 100)
						{
							cout << "Error 6439011209" << endl;
							cout << P.coord[0] << " " << P.coord[1] <<
								" " << P.coord[2] << endl;
							exit(-1);
						}
						Move[0] = (-P.coord[0] + Center_cell[0]) / 1000.0;
						Move[1] = (-P.coord[1] + Center_cell[1]) / 1000.0;
						Move[2] = (-P.coord[2] + Center_cell[2]) / 1000.0;
						P.Move(Move);
						ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);
					}

					ppp = this->Find_cell_point(P.coord[0] + P.Vel[0] * dt,
						P.coord[1] + P.Vel[1] * dt,
						P.coord[2] + P.Vel[2] * dt,
						0, previos);

					klk = 0;
					while (ppp == nullptr || ppp != P.cel)
					{
						klk++;
						if (klk > 100)
						{
							cout << "Error 3296411221" << endl;
							exit(-1);
						}
						Move[0] = (-P.coord[0] + Center_cell[0]) / 1000.0;
						Move[1] = (-P.coord[1] + Center_cell[1]) / 1000.0;
						Move[2] = (-P.coord[2] + Center_cell[2]) / 1000.0;
						P.Move(Move);
						ppp = this->Find_cell_point(P.coord[0] + P.Vel[0] * dt,
							P.coord[1] + P.Vel[1] * dt,
							P.coord[2] + P.Vel[2] * dt,
							0, previos);
					}
				}

				// Отладочная информация
				if (false)
				{
					cout << "_______________________________" << endl;
					cout << "Zapusk test  " << P.sort << "    " << gr->number << endl;
					whach(P.coord[0]);
					whach(P.coord[1]);
					whach(P.coord[2]);
					whach(P.Vel[0]);
					whach(P.Vel[1]);
					whach(P.Vel[2]);
					cout << "mu = " << P.mu << endl;
					cout << "_______________________________" << endl;
				}

				Cell* previos = P.cel;
				Cell* ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);
				if (ppp != P.cel)
				{
					cout << "Error  8756121199" << endl;
					P.coord[0] = P.cel->center[0][0];
					P.coord[1] = P.cel->center[0][1];
					P.coord[2] = P.cel->center[0][2];
				}


				P.KSI = -log(1.0 - this->Sensors[sens_num1]->MakeRandom());
				P.I_do = 0.0;

				//cout << "FLY" << endl;

				this->MK_fly_immit(P, zone_MK, this->Sensors[sens_num2]); // Запускаем частицу в полёт   // !! Не написана
				//exit(-1);
				//cout << "END" << endl;
			}

			if (gr->type == Type_Gran::Us) func->Delete();
		}
	}

	cout << "**********************************" << endl;
	cout << "Obshee chislo chastic = " << ALL_N << endl;

	// Нормировка функции распределения
	for (auto& gr : this->MK_Grans[zone_MK - 1])
	{
		short int ni = 1; // Номер "выходящей" функции распределения
		if (gr->cells[0]->MK_zone == zone_MK)
		{
			ni = 0;
		}

		for (short int nh_ = 0; nh_ < this->phys_param->num_H; ++nh_)
		{
			auto& func = gr->AMR[nh_][ni];
			func->Normir_velocity_volume(gr->area[0]);
		}
	}

	// Нормировка Моментов в ячейках
	k1 = 0;
	cout << "Start: Normir moment in cells" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
	{
		auto cell = this->All_Cell[idx];
		if (cell->MK_zone != zone_MK) continue;
		#pragma omp critical (first) 
		{
			k1++;
			if (k1 % 10000 == 0)
			{
				cout << "Cells = " << k1 << endl;
			}
		}

		cell->MK_normir_Moments(this->phys_param);
		cell->MK_calc_Sm(this->phys_param);  // Нужно параллелить, так как эта функция долго обрабатывается
	}
	cout << "End: Normir moment in cells" << endl;

	// Выведем одну функцию посмотреть что получилось)
	if (true)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			if (gr->type2 != Type_Gran_surf::BS) continue;

			if( fabs(fabs(gr->normal[0][0]) - 0.936503) > 0.0001) continue;
			if( fabs(fabs(gr->normal[0][1]) - 0.350126) > 0.0001) continue;
			if( fabs(fabs(gr->normal[0][2]) - 0.0193568) > 0.0001) continue;


			short int ni = 1; // Номер "выходящей" функции распределения
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			cout << gr->center[0][0] << " " << gr->center[0][1] << " " <<
				gr->center[0][2] << endl;
			cout << "Normal = " << gr->normal[0][0] << " " << gr->normal[0][1] << " " <<
				gr->normal[0][2] << endl;
			cout << "N_particle = " << gr->N_particle << endl;
			auto& func = gr->AMR[3][ni];
			cout << "n_func = " << func->parameters["n"] << endl;
			cout << "nn_func = " << func->parameters["nn"] << endl;
			cout << "Potok cherez funkcion = " << func->parameters["Smu"] << endl;
			double dF = func->parameters["nn"] / func->parameters["Smu"] * gr->area[0]
				- kv(func->parameters["n"] * gr->area[0] / func->parameters["Smu"]);
			cout << "dF = " << dF << endl;
			cout << "delta = " << sqrt(dF / func->parameters["Smu"]) / func->parameters["n"] << endl;
			func->Print_all_center_Tecplot(func, "test");
			func->Print_1D_Tecplot(func, -2.54327);
			func->Culk_SpotokV(gr->area[0]);
			cout << "Spotok = " << func->SpotokV / gr->area[0] << endl;
			break;
		}
	}


	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "MK all time: " << duration.count() / 1000.0 / 60.0 << " minutes" << std::endl;
}

void Setka::MK_fly_immit(MK_particle& P, short int zone_MK, Sensor* Sens)
{
	/*cout << "______Start_MK_fly_immit___________" << endl;
	whach(P.coord[0]);
	whach(P.coord[1]);
	whach(P.coord[2]);
	cout << "_______________________________" << endl;*/
	//cout << "Start " << endl;
	//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;


	Eigen::Vector3d coord_init;
	Eigen::Vector3d Vel_init;

	coord_init << P.coord[0], P.coord[1], P.coord[2];
	Vel_init << P.Vel[0], P.Vel[1], P.Vel[2];


	unsigned int k_cikl = 0;
	bool vtoroy_shans = false;
	bool vtoroy_shans2 = false;

	// Главный цикл по ячейкам
	// Выйти из него можно только если частица достигнет конца области
	while (true)
	{
		Eigen::Vector3d coord_do;
		coord_do[0] = P.coord[0];
		coord_do[1] = P.coord[1];
		coord_do[2] = P.coord[2];

		//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
		if (P.cel->MK_zone != zone_MK)
		{
			cout << "Error 8767567487" << endl;
			cout << P.cel->MK_zone << "   " << zone_MK << endl;
			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
			cout << k_cikl << endl;
			return;
			exit(-1);
		}

		k_cikl++;
		if (k_cikl > 10000)
		{
			cout << "Error 8675498765" << endl;
			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
			cout << P.sort << " " << P.KSI << " " << P.I_do << endl;
			//exit(-1);
		}

		Cell* Cell_do = P.cel;  // На всякий случай сохраним стартовую ячейку, вдруг надо будет вернуться
		if (Cell_do == nullptr)
		{
			cout << "Error 9865749586" << endl;
			exit(-1);
		}

		double time = 0.0;            // время нахождения частицы в ячейке
		Gran* gran = nullptr;         // Через какую грань ячейка выйдет из ячейки

		//cout << "A1 " << endl;
		//cout << "P.coord = " << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
		// Находим время до выхода частицы из ячейки, а также через какую грань будет выход
		bool b1 = false;
		unsigned short int k1 = 0;
		while (b1 == false)
		{
			k1++;
			if (P.cel == nullptr)
			{
				cout << "Error 12569834678" << endl;
				exit(-1);
			}
			b1 = this->Time_to_vilet(P, time, gran);
			if (P.cel->MK_zone != zone_MK)
			{
				cout << "Error 9871216655" << endl;
				cout << P.cel->MK_zone << "   " << zone_MK << endl;
				cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
				cout << k_cikl << endl;
				exit(-1);
			}

			if (b1 == true && gran == nullptr)
			{
				cout << "Error 7786341271" << endl;
				exit(-1);
			}

			if (b1 == false)
			{
				if (P.cel == nullptr)
				{
					cout << "Error 8675463895" << endl;
					exit(-1);
				}
				Eigen::Vector3d Cell_centerr;
				Cell_centerr << P.cel->center[0][0], P.cel->center[0][1],
					P.cel->center[0][2];

				// Подвинем немного точку к центру ячейки

				for (short int i = 0; i < 3; i++)
				{
					/*if (k1 < 2)
					{
						P.coord[i] += 1e-6 * P.Vel[i];
					}
					else if (k1 < 4)
					{
						P.coord[i] += 1e-5 * P.Vel[i];
					}*/
					if(k1 < 10)
					{
						P.coord[i] = P.coord[i] + (Cell_centerr[i] - P.coord[i]) / 800.0;
					}
					else
					{
						P.coord[i] = Cell_centerr[i];
					}
				}
				// Немного двигаем точку
				//P.coord += 1e-6 * P.Vel;

				auto cepp_prev = P.cel;
				P.cel = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, cepp_prev);
				// Здесь надо проверить, что во время микро-движения точка не
				// вышла в другую ячейку или за пределы расчётной области
				if (P.cel == nullptr)
				{
					//cout << "TUT  1875408695" << endl;
					P.cel = Cell_do;
				}
			}

			if (k1 > 11)
			{
				if (vtoroy_shans2 == false)
				{
					vtoroy_shans2 = true;
					P.cel = Cell_do;
					P.coord[0] = coord_do[0];
					P.coord[1] = coord_do[1];
					P.coord[2] = coord_do[2];
					perturbVectorKeepingMagnitude(P.Vel[0], P.Vel[1], P.Vel[2], 0.1 * norm2(P.Vel[0], P.Vel[1], P.Vel[2]));

					if (P.cel->MK_zone != zone_MK)
					{
						cout << "Error 9856312132" << endl;
					}


					continue;
				}

				cout << "Poteryal D" << endl;
				cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
				return;
				//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				//whach(P.cel->number);
				//exit(-1);
			}
		}

		//cout << "B " << endl;
		/*cout << "D " << time << endl;
		whach(P.coord[0]);
		whach(P.coord[1]);
		whach(P.coord[2]);
		cout << "_______________________________" << endl;*/

		if (gran == nullptr)
		{
			cout << "Error 6438609412" << endl;
			exit(-1);
		}

		vtoroy_shans2 = false;

		// Здесь время до выхода из ячейки определено time
		// Также определено через какую грань это произойдёт  gran

		// далее блок основной программы в ячейке
		// ****************************************************************************

		// Получаем параметры плазмы в ячейке ----------------------------
		double ro = P.cel->parameters[0]["rho"];
		double p = P.cel->parameters[0]["p"];
		double rho_He = P.cel->parameters[0]["rho_He"];
		double cp;// = sqrt(P.cel->parameters[0]["p"] / ro);
		double vx = P.cel->parameters[0]["Vx"];			// Скорости плазмы в ячейке
		double vy = P.cel->parameters[0]["Vy"];
		double vz = P.cel->parameters[0]["Vz"];

		double rho_Th, p_Th;
		
		//Sootnosheniya(ro, p, rho_He, 0.0, 0.0, (int)(P.cel->type),
		//	rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E);

		unordered_map<string, double> param;

		this->phys_param->Plasma_components((int)(P.cel->type), P.cel->parameters[0], param);

		rho_Th = param["rho_Th"];
		p_Th = param["p_Th"];

		if (rho_Th <= 1e-8) rho_Th = 1e-8;
		if (p_Th <= 1e-8/2.0) p_Th = 1e-8/2.0;

		ro = rho_Th;
		cp = sqrt(2.0 * p_Th / rho_Th);

		// ---------------------------------------------------------------

		// Найдём время до перезарядки

		double l = sqrt(kvv(time * P.Vel[0], time * P.Vel[1], time * P.Vel[2]));
		double I = P.I_do;
		double u = sqrt(kvv(P.Vel[0] - vx, P.Vel[1] - vy, P.Vel[2] - vz));
		double u1 = vx - P.Vel[0];
		double u2 = vy - P.Vel[1];
		double u3 = vz - P.Vel[2];
		double skalar = u1 * P.Vel[0] + u2 * P.Vel[1] + u3 * P.Vel[2];
		double Vel_norm = sqrt(kvv(P.Vel[0], P.Vel[1], P.Vel[2]));

		double uz = Velosity_1(u, cp);
		double nu_ex = ro * uz * sigma(uz) / this->phys_param->par_Kn;
		double sig = Vel_norm / nu_ex;
		I += l / sig;


		if (P.cel->MK_zone != zone_MK)
		{
			cout << P.cel->number << " " << Cell_do->number << endl;
			cout << "Error 1654875068" << endl;
			exit(-1);
		}


		if (vtoroy_shans == false)
		{
			if (I < P.KSI)
			{
				P.I_do = I;  // В этом случае перезарядки в ячейке не произошло
				// Здесь записываем необходимые моменты в ячейку ---------------------
				if (this->phys_param->culc_cell_moments == true)
				{
					P.cel->MK_Add_particle(P, time);
				}

				if (this->phys_param->MK_source_S == true)
				{
					P.cel->MK_Add_pui_source(u, nu_ex, P.mu, time, this->phys_param);
				}
				// -------------------------------------------------------------------
			}
			else
			{
				double ksi = (P.KSI - P.I_do) * sig;
				double t_ex = ksi / Vel_norm;
				P.I_do = 0.0;
				for (short int i = 0; i < 3; i++) P.coord[i] += t_ex * P.Vel[i];
				Cell* Cnow = P.cel;
				Cell* CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);

				if (P.cel != CC)
				{
					// Если перезарядка произошла за пределами текущей ячейки
					for (short int i = 0; i < 3; i++)
					{
						P.coord[i] -= t_ex / 1000.0 * P.Vel[i];
					}

					Cnow = CC;
					CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);
					if (P.cel != CC)
					{
						Eigen::Vector3d Cell_center;
						Cell_center << P.cel->center[0][0], P.cel->center[0][1],
							P.cel->center[0][2];
						unsigned short int kklk = 0;
					dchj12:
						kklk++;
						if (kklk > 20)
						{
							cout << "Poteryal C" << endl;
							return;
						}

						// Подвинем немного точку к центру ячейки
						for (short int i = 0; i < 3; i++)
						{
							if (kklk < 18)
							{
								P.coord[i] += (Cell_center[i] - P.coord[i]) / 100.0;
							}
							else
							{
								P.coord[i] = Cell_center[i];
							}
						}

						Cnow = CC;
						CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);

						if (P.cel != CC)
						{
							goto dchj12;
							cout << "Error 8674539765" << endl;
							whach(CC->number);
							whach(P.cel->number);
							whach(P.coord[0]);
							whach(P.coord[1]);
							whach(P.coord[2]);
							whach(P.Vel[0]);
							whach(P.Vel[1]);
							whach(P.Vel[2]);
							whach(t_ex);
							whach(time);
							exit(-1);
						}
					}
				}

				double uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * const_pi * sqrtpi_);
				double uz_E = Velosity_3(u, cp);

				// Здесь записываем необходимые моменты в ячейку ---------------------
				if (this->phys_param->culc_cell_moments == true)
				{
					P.cel->MK_Add_particle(P, t_ex);
				}

				if (this->phys_param->MK_source_S == true)
				{
					P.cel->MK_Add_pui_source(u, nu_ex, P.mu, t_ex, this->phys_param);
				}
				// -------------------------------------------------------------------

				// Разыгрываем новую скорость
				double Ur, Uphi, Uthe;
				double Vr, Vphi, Vthe;
				double Wr, Wthe, Wphi;
				spherical_skorost(P.coord[0], P.coord[1], P.coord[2],
					vx, vy, vz, Ur, Uphi, Uthe);
				spherical_skorost(P.coord[0], P.coord[1], P.coord[2],
					P.Vel[0], P.Vel[1], P.Vel[2], Vr, Vphi, Vthe);
				this->M_K_Change_Velosity(Sens, Ur / cp, Uthe / cp, Uphi / cp,
					Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, cp);
				Wr *= cp;
				Wthe *= cp;
				Wphi *= cp;

				dekard_skorost(P.coord[0], P.coord[1], P.coord[2],
					Wr, Wphi, Wthe, P.Vel[0], P.Vel[1], P.Vel[2]);

				P.sort = (short int)(P.cel->type);
				P.KSI = -log(1.0 - Sens->MakeRandom());
				vtoroy_shans = false;
				if (P.cel->MK_zone != zone_MK)
				{
					cout << P.cel->number << " " << Cell_do->number << endl;
					cout << "Error 9088565453" << endl;
				}
				continue;
				//return this->MK_fly_immit(P, zone_MK, Sens);
			}
		}

		// ****************************************************************************
		// Находим следующую ячейку
		for (short int i = 0; i < 3; i++)
		{
			P.coord[i] += 1.000001 * time * P.Vel[i];
		}

		if (norm2(P.coord[0], P.coord[1], P.coord[2]) < 1.01 * this->geo->R0)
		{
			// Частица попала во внутреннюю сферу, надо, чтобы они пролетели мимо неё
			Eigen::Vector3d X(P.coord[0], P.coord[1], P.coord[2]);
			Eigen::Vector3d V(P.Vel[0], P.Vel[1], P.Vel[2]);
			double time_;
			if (findSphereIntersectionTime(X, V, 1.01 * this->geo->R0,
				time_) == true)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += time_ * P.Vel[i];
				}
			}
			else
			{
				double norm_ = norm2(P.Vel[0], P.Vel[1], P.Vel[2]);
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (2.02 * this->geo->R0) * P.Vel[i] / norm_;
				}
			}
		}

		if (gran->Have_zone_number(zone_MK))
		{
			a1:
			// В этом случае долетели до границы, записываем что надо и выключаем частицу
			short int nn = 1;
			if (gran->cells[0]->MK_zone == zone_MK) nn = 0;
			auto AMR = gran->AMR[P.sort - 1][nn];

			if (this->phys_param->culc_AMR == true && gran->type == Type_Gran::Us)
			{
				// В грани на границе нет смысла ничего записывать
				AMR->Add_particle(P.Vel[0], P.Vel[1], P.Vel[2], P.mu); // мьютексы внутри
			}

			gran->mut.lock(); // Мбютекс для записи в гранб
			gran->N_particle++;
			gran->mut.unlock();

			return;
		}

		Cell* Cell_next = P.cel->Get_Sosed(gran);
		short unsigned int kkk2 = 0;
	vv1:
		kkk2++;
		// точно находим следующую ячейку
		Cell* Cell_next_ = Cell_next;
		P.cel = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cell_next_);

		// Кажется в случае проблем надо просто подтянуть ячейку к центру грани

		if (P.cel == nullptr)
		{
			for (auto& gr : Cell_do->grans)
			{
				if (gr->Have_zone_number(zone_MK))
				{
					gran = gr;
					goto a1;
				}
			}
		}

		//cout << "F " << endl;
		// В этом случае попали в следующую зону, пропустив граничную грань
		if (P.cel != nullptr && P.cel->MK_zone != zone_MK)
		{
			for (auto& gr : P.cel->grans)
			{
				if (gr->Have_zone_number(zone_MK))
				{
					gran = gr;
					goto a1;
				}
			}

			// В этом случае надо либо отключать ячейку (но мы потеряем часть массы)
			// либо запускать ей заново в этой ячейке (тогда наоборот получим лишнюю массу, так
			// как она уже записалась в данную ячейку

			if (vtoroy_shans == false)
			{
				P.coord[0] = coord_do[0] + (Cell_do->center[0][0] - coord_do[0]) / 300.0;
				P.coord[1] = coord_do[1] + (Cell_do->center[0][1] - coord_do[1]) / 300.0;
				P.coord[2] = coord_do[2] + (Cell_do->center[0][2] - coord_do[2]) / 300.0;
				P.cel = Cell_do;
				perturbVectorKeepingMagnitude(P.Vel[0], P.Vel[1], P.Vel[2], 0.01);
				vtoroy_shans = true;

				if (P.cel->MK_zone != zone_MK)
				{
					cout << "Error 1213563589" << endl;
				}

				continue;
			}

			cout << "Poteryal A" << endl;
			return;
		}

		if (P.cel == nullptr)
		{
			// В этом случае точка часто попадает в угол ячеки или на грань

			if (kkk2 < 3)
			{
				Eigen::Vector3d normal;
				normal << gran->normal[0][0], gran->normal[0][1], gran->normal[0][2];
				if (gran->cells[0] != Cell_do)
				{
					normal *= -1.0;
				}
				double l_ = norm2(gran->yzels[0]->coord[0][0] - gran->yzels[1]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[1]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[1]->coord[0][2]);
				l_ = min(l_, norm2(gran->yzels[0]->coord[0][0] - gran->yzels[2]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[2]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[2]->coord[0][2]));
				l_ = min(l_, norm2(gran->yzels[0]->coord[0][0] - gran->yzels[3]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[3]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[3]->coord[0][2]));

				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (gran->center[0][i] - P.coord[i])/200.0 + normal[i] * l_ / 200.0;
				}
				goto vv1;
			}
			if (kkk2 < 4)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += 0.001 * time * P.Vel[i];
				}
				goto vv1;
			}
			else if (kkk2 < 6)
			{
				perturbVector(P.coord[0], P.coord[1], P.coord[2], Cell_do->geo_parameters["l_size"]/100);
				goto vv1;
			}
			else if (kkk2 < 9)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (-Cell_do->center[0][i] + P.coord[i]) / 100.0;
				}
				goto vv1;
			}

			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << "Poteryal B" << endl;
			return;
		}

		//cout << "G " << endl;
		//cout << "2 Soburausi otpravit " << endl;
		//cout << "P.coord = " << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
		vtoroy_shans = false;

		if (P.cel->MK_zone != zone_MK)
		{
			cout << "Error 4562190567" << endl;
		}

		continue;
		//return this->MK_fly_immit(P, zone_MK, Sens);
	}
}

double Setka::Get_Spotok_inf(const Eigen::Vector3d& n)
{
	Eigen::Vector3d V_inf_;
	V_inf_ << this->phys_param->Velosity_inf, 0.0, 0.0;

	double Ux = V_inf_.dot(n);

	return -(-exp(-kv(Ux)) + sqrtpi_ * Ux * erfc(Ux)) / (2.0 * sqrtpi_);
}

void Setka::Velosity_initial(Sensor* s, Eigen::Vector3d& V,
	const Eigen::Vector3d& n, const Eigen::Vector3d& t, 
	const Eigen::Vector3d& m)
{
	Eigen::Vector3d V_inf_;
	Eigen::Vector3d V_inf;
	V_inf_ << this->phys_param->Velosity_inf, 0.0, 0.0;

	V_inf[0] = V_inf_.dot(n);
	V_inf[1] = V_inf_.dot(t);
	V_inf[2] = V_inf_.dot(m);

	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(V_inf[0]) * sqrtpi_ /
		(1.0 + fabs(V_inf[0]) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(const_pi * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + V_inf[0]) / (fabs(V_inf[0]) + fabs(z)) < ksi6 || z > -V_inf[0]);

	double V1 = z + V_inf[0];

	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(ksi2));
	double V2 = V_inf[1] + a * cos(2.0 * const_pi * ksi1);
	double V3 = V_inf[2] + a * sin(2.0 * const_pi * ksi1);

	V = V1 * n + V2 * t + V3 * m;

	return;
}

bool Setka::Time_to_vilet(const MK_particle& P, double& time, Gran*& gran)
{
	Cell* C = P.cel;
	Eigen::Vector3d R, V;

	R[0] = P.coord[0];
	R[1] = P.coord[1];
	R[2] = P.coord[2];
	V[0] = P.Vel[0];
	V[1] = P.Vel[1];
	V[2] = P.Vel[2];

	Gran* gran_min = nullptr;
	double time_min = 1e10;
	double time1;

	bool b1 = false;

	for (const auto& gr : C->grans)
	{
		if (gr == nullptr)
		{
			cout << "Error 9453286475" << endl;
			exit(-1);
		}

		if (gr->Luch_iz_cross_approx(R, V) == true)
		{
			if (gr->Luch_crossing(R, V, time1) == true)
			{
				//if (time1 > 1e-06) continue;

				if (time_min > time1)
				{
					time_min = time1;
					gran_min = gr;
					b1 = true;
				}
			}
		}
	}

	if (b1 == false)
	{
		// Пересечение ни с одной гранью не произошло!
		return false;
	}
	else
	{
		gran = gran_min;
		time = time_min;

		if (gran == nullptr)
		{
			cout << "Error 1112389046" << endl;
			exit(-1);
		}

		return true;
	}
}


void Setka::M_K_Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe,
	const double& Uphi, const double& Vr, const double& Vthe, 
	const double& Vphi, double& Wr, double& Wthe, double& Wphi, const double& cp)
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h;

	double gg = 0.0;
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * const_pi * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * const_pi * ksi5);
			// Более экономичный алгоритм   --  выйгрыша нет вроде от него
			/*do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1.0);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(const_pi * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * const_pi * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * const_pi * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(X, cp) * (X + yy)));
	} while (h < ksi6); 


	Wr = v1;
	Wthe = v2;
	Wphi = v3;


	return;
}