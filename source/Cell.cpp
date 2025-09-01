#include "Cell.h"

void Cell::Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par)
{
	Eigen::Vector3d point;
	Eigen::Vector3d point2;

	point << x, y, z;

	if (this->interpolate_alpha.find("rho") != this->interpolate_alpha.end())
	{
		double result = 0.0;
		for (size_t i = 0; i < this->yzels.size(); ++i) 
		{
			point2 << this->yzels[i]->coord[0][0], this->yzels[i]->coord[0][1],
				this->yzels[i]->coord[0][2];
			result += this->interpolate_alpha["rho"][i] * rbfKernel((point - point2).norm());
		}
		par["rho"] = result;
	}
}

void Cell::Get_IDW_interpolation(const double& x, const double& y, const double& z, 
	unordered_map<string, double>& par, Phys_param* phys_param)
{
	Eigen::Vector3d point;
	Eigen::Vector3d point2;

	point << x, y, z;

	unordered_map<string, double> sum_weights;
	unordered_map<string, double> sum_weighted_values;

	for (auto& nam : phys_param->param_names)
	{
		sum_weights[nam] = 0.0;
		sum_weighted_values[nam] = 0.0;
	}

	//double sum_weights = 0.0;
	//double sum_weighted_values = 0.0;

	for (size_t i = 0; i < this->yzels.size(); ++i)
	{
			
		point2 << this->yzels[i]->coord[0][0], 
			this->yzels[i]->coord[0][1], this->yzels[i]->coord[0][2];
		double dist = 0.0001 + (point - point2).norm();
		double weight = 1.0 / std::pow(dist, 1.0);

		for (auto& nam : phys_param->param_names)
		{
			if (this->yzels[i]->parameters.find(nam) != this->yzels[i]->parameters.end())
			{

				sum_weights[nam] += weight;
				sum_weighted_values[nam] += weight * this->yzels[i]->parameters[nam];
				//cout << this->yzels[i]->parameters["rho"] << endl;
			}
		}
	}


	// Добавим центр ячейки
	point2 << this->center[0][0],
		this->center[0][1], this->center[0][2];
	double dist = 0.0001 + (point - point2).norm();
	double weight = 1.0 / std::pow(dist, 1.0);

	for (auto& nam : phys_param->param_names)
	{
		if (this->parameters[0].find(nam) != this->parameters[0].end())
		{

			sum_weights[nam] += weight;
			sum_weighted_values[nam] += weight * this->parameters[0][nam];
			//cout << this->yzels[i]->parameters["rho"] << endl;
		}
	}

	for (auto& nam : phys_param->param_names)
	{
		if (sum_weights[nam] > 0.00001)
		{
			par[nam] = sum_weighted_values[nam] / sum_weights[nam];
		}
		else
		{
			par[nam] = 0.0;
		}
	}
}

Cell* Cell::Get_Sosed(Gran* gr)
{
	if (gr == nullptr)
	{
		cout << "Error 5434323121" << endl;
		exit(-1);
	}

	if(gr->cells.size() == 1) return nullptr;

	if (gr->cells[0]->number == this->number)
	{
		return gr->cells[1];
	}
	else
	{
		return gr->cells[0];
	}
	return nullptr;
}

double determinant4x4(const Eigen::Vector4d& a, const Eigen::Vector4d& b, 
	const Eigen::Vector4d& c, const Eigen::Vector4d& d) 
{
	Eigen::Matrix4d mat;
	mat << a, b, c, d;
	return mat.determinant();
}

bool isPointInsideTetrahedron(const Eigen::Vector3d& P, const Eigen::Vector3d& A, const Eigen::Vector3d& B, 
	const Eigen::Vector3d& C, const Eigen::Vector3d& D) 
{
	// Добавляем 4-ю координату (1) для работы с 4×4 матрицами
	Eigen::Vector4d A4(A(0), A(1), A(2), 1.0);
	Eigen::Vector4d B4(B(0), B(1), B(2), 1.0);
	Eigen::Vector4d C4(C(0), C(1), C(2), 1.0);
	Eigen::Vector4d D4(D(0), D(1), D(2), 1.0);
	Eigen::Vector4d P4(P(0), P(1), P(2), 1.0);

	double detABCD = determinant4x4(A4, B4, C4, D4);
	if (std::abs(detABCD) < 1e-10) {
		std::cerr << "Degenerate tetrahedron (zero volume)!" << std::endl;
		cout << A4(0) << " " << A4(1) << " " << A4(2) << endl;
		cout << B4(0) << " " << B4(1) << " " << B4(2) << endl;
		cout << C4(0) << " " << C4(1) << " " << C4(2) << endl;
		cout << D4(0) << " " << D4(1) << " " << D4(2) << endl;
		cout << "__________________________________________" << endl;
		return false;
	}

	double u = determinant4x4(P4, B4, C4, D4) / detABCD;
	if (u < -1e-10) return false;
	double v = determinant4x4(A4, P4, C4, D4) / detABCD;
	if (v < -1e-10) return false;
	double w = determinant4x4(A4, B4, P4, D4) / detABCD;
	if (w < -1e-10) return false;
	double t = determinant4x4(A4, B4, C4, P4) / detABCD;
	if (t < -1e-10) return false;

	return (std::fabs(u + v + w + t - 1.0) < 1e-6);
}


bool Cell::Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next)
{
	// fast = false - правильный алгоритм
	// fast = true - быстрый алгоритм (но имеет погрешность)

	Eigen::Vector3d P(x, y, z); // Точка
	Next = nullptr;  // Next - следующая ячейка определяется только если fast == true
	// алгоритм для fast == false не позволяет найти следующую ячейку

	//double xc = this->center[now][0];
	//double yc = this->center[now][1];
	//double zc = this->center[now][2];

	// Честное определение принадлежности
	if (fast == false)
	{
		Eigen::Vector3d A(this->center[now][0], this->center[now][1], this->center[now][2]);
		Eigen::Vector3d B, C, D;
		bool is_in = false;
		for (auto& gr : this->grans)
		{
			//double x1 = gr->center[now][0];
			//double y1 = gr->center[now][1];
			//double z1 = gr->center[now][2];
			if (false) // старый вариант, где грань = 4 треугольника
			{
				B << gr->center[now][0], gr->center[now][1], gr->center[now][2];
				int n_yz = gr->yzels.size();
				for (short int i = 0; i < n_yz; i++)
				{
					int j = i + 1;
					if (j >= n_yz) j = 0;
					C << gr->yzels[i]->coord[now][0], gr->yzels[i]->coord[now][1], gr->yzels[i]->coord[now][2];
					D << gr->yzels[j]->coord[now][0], gr->yzels[j]->coord[now][1], gr->yzels[j]->coord[now][2];
					is_in = isPointInsideTetrahedron(P, A, B, C, D);
					if (is_in == true) return true;
				}
			}
			else // грань = 2 треугольника
			{
				if (gr->yzels.size() != 4)
				{
					cout << "Error 9765486573" << endl;
					exit(-1);
				}

				B << gr->yzels[0]->coord[now][0], gr->yzels[0]->coord[now][1], gr->yzels[0]->coord[now][2];
				C << gr->yzels[1]->coord[now][0], gr->yzels[1]->coord[now][1], gr->yzels[1]->coord[now][2];
				D << gr->yzels[2]->coord[now][0], gr->yzels[2]->coord[now][1], gr->yzels[2]->coord[now][2];
				is_in = isPointInsideTetrahedron(P, A, B, C, D);
				if (is_in == true) return true;

				
				C << gr->yzels[3]->coord[now][0], gr->yzels[3]->coord[now][1], gr->yzels[3]->coord[now][2];
				is_in = isPointInsideTetrahedron(P, A, B, C, D);
				if (is_in == true) return true;

			}
		}

		return false;
	}
	else
	{
		Eigen::Vector3d B, C, D;

		for (auto& gr : this->grans)
		{
			if (gr->cells.size() != 2) continue;

			B << gr->center[now][0], gr->center[now][1], gr->center[now][2];
			C << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];

			if (gr->cells[0]->number == this->number) C = -C;

			D = P - B;

			if (C.dot(D) < -1e-10) 
			{
				Next = this->Get_Sosed(gr);
				return false;
			}
		}
		return true;
	}

	return false;
}

void Cell::Set_Cell_Geo_for_MK(void)
{
	Eigen::Vector3d A;
	Eigen::Vector3d B;
	Eigen::Vector3d C;


	// Считаем l_size   характерный размер ячейки
	if (true)
	{
		double S = 0.0;

		A[0] = this->yzels[0]->coord[0][0];
		A[1] = this->yzels[0]->coord[0][1];
		A[2] = this->yzels[0]->coord[0][2];

		B[0] = this->yzels[0]->coord[1][0];
		B[1] = this->yzels[0]->coord[1][1];
		B[2] = this->yzels[0]->coord[1][2];
		C = B - A;

		S = C.norm();

		B[0] = this->yzels[0]->coord[4][0];
		B[1] = this->yzels[0]->coord[4][1];
		B[2] = this->yzels[0]->coord[4][2];
		C = B - A;
		S = min(S, C.norm());


		B[0] = this->yzels[0]->coord[2][0];
		B[1] = this->yzels[0]->coord[2][1];
		B[2] = this->yzels[0]->coord[2][2];
		C = B - A;
		S = min(S, C.norm());

		this->geo_parameters["l_size"] = S;
	}
}


void Cell::Culc_center(unsigned short int st_time)
{
	double xc, yc, zc;
	int Ny = this->yzels.size();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	// Вычисляем центр грани
	for (auto& i : this->yzels)
	{
		xc += i->coord[st_time][0];
		yc += i->coord[st_time][1];
		zc += i->coord[st_time][2];
	}
	xc /= Ny;
	yc /= Ny;
	zc /= Ny;
	this->center[st_time][0] = xc;
	this->center[st_time][1] = yc;
	this->center[st_time][2] = zc;
}

void Cell::Culc_volume(unsigned short int st_time, unsigned short int method)
{
	// 0 - быстрый вариант (нужны посчанные площади граней и их центры)
	// 1 - медленный вариант (нужны посчитанные только центры граней)
	// для ровных граней оба методы работают одинаково, но чем более кривая грань, 
	// тем большее расхождения получается.
	// медленный вариант более правильный (но это не точно). 
	// Думаю всегда можно выбирать вариант "0"
	if (method == 0)
	{
		double xg, yg, zg, h;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			xg = i->center[st_time][0];
			yg = i->center[st_time][1];
			zg = i->center[st_time][2];
			h = fabs(scalarProductFast(i->normal[st_time][0], i->normal[st_time][1],
				i->normal[st_time][2],
				this->center[st_time][0] - xg, this->center[st_time][1] - yg,
				this->center[st_time][2] - zg));
			V += h * i->area[st_time];
		}
		this->volume[st_time] = V / 3.0; 
		return;
	}
	if (method == 1)
	{
		int j1;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			if (false) // старый вариант где грань = 4 треугольника
			{
				for (int j = 0; j < i->yzels.size(); j++)
				{
					j1 = j + 1;
					if (j1 >= i->yzels.size()) j1 = 0;
					V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
						this->center[st_time][2],
						i->center[st_time][0], i->center[st_time][1], i->center[st_time][2],
						i->yzels[j]->coord[st_time][0], i->yzels[j]->coord[st_time][1],
						i->yzels[j]->coord[st_time][2],
						i->yzels[j1]->coord[st_time][0], i->yzels[j1]->coord[st_time][1],
						i->yzels[j1]->coord[st_time][2]);
				}
			}
			else
			{
				V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
					this->center[st_time][2],
					i->yzels[0]->coord[st_time][0], i->yzels[0]->coord[st_time][1],
					i->yzels[0]->coord[st_time][2],
					i->yzels[1]->coord[st_time][0], i->yzels[1]->coord[st_time][1],
					i->yzels[1]->coord[st_time][2],
					i->yzels[2]->coord[st_time][0], i->yzels[2]->coord[st_time][1],
					i->yzels[2]->coord[st_time][2]);
				V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
					this->center[st_time][2],
					i->yzels[0]->coord[st_time][0], i->yzels[0]->coord[st_time][1],
					i->yzels[0]->coord[st_time][2],
					i->yzels[3]->coord[st_time][0], i->yzels[3]->coord[st_time][1],
					i->yzels[3]->coord[st_time][2],
					i->yzels[2]->coord[st_time][0], i->yzels[2]->coord[st_time][1],
					i->yzels[2]->coord[st_time][2]);
			}
		}
		this->volume[st_time] = V;
		return;
	}
	else
	{
		cout << "Error  0967452122" << endl;
		return;
	}
}

double Cell ::func_R(unsigned short int i_time)
{
	return norm2(this->center[i_time][0], this->center[i_time][1], this->center[i_time][2]);
}

void Cell::Tecplot_print_cell(void)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	fout << "ZONE T=HP, N = " << this->yzels.size() << ", E = " << this->yzels.size() - 1 << ", F=FEPOINT, ET=LINESEG" << endl;

	for (auto& j : this->yzels)
	{
		fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
	}

	for (int m = 0; m < this->yzels.size() - 1; m++)
	{
		fout << m + 1 << " " << m + 2 << endl;
	}

	fout << endl;


	fout.close();
}

void Cell::MK_Add_particle(MK_particle& P, const double& time)
{
		this->mut.lock();

		this->parameters[0]["MK_n_H"] += time * P.mu;

		switch (P.sort)
		{
		case 1:
			this->parameters[0]["MK_n_H1"] += time * P.mu;
			break;
		case 2:
			this->parameters[0]["MK_n_H2"] += time * P.mu;
			break;
		case 3:
			this->parameters[0]["MK_n_H3"] += time * P.mu;
			break;
		case 4:
			this->parameters[0]["MK_n_H4"] += time * P.mu;
			break;
		default:
			break;
		}

		this->mut.unlock();
}

void Cell::MK_Add_pui_source(const double& wr, const double& nu_ex, const double& mu, 
	const double& time, Phys_param* phys_param)
{
	int index = static_cast<int>(wr / phys_param->pui_wR * phys_param->pui_nW);
	if (index < 0)  index = 0;
	if (index >= phys_param->pui_nW)  index = phys_param->pui_nW - 1;
	
	this->mut.lock();

	this->pui_Sm[index] += mu * time;
	this->pui_Sp[index] +=  mu * time * nu_ex;

	this->mut.unlock();
}

void Cell::MK_calc_Sm(Phys_param* phys_param)
{
	vector<double> pui_Sm2(phys_param->pui_nW);
	for (auto& i : pui_Sm2) i = 0.0;
	double dthe = const_pi / 40.0;
	double Vh, ff, the, d;

	for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
	{
		double w = ((ij + 0.5) * phys_param->pui_wR / phys_param->pui_nW);
		for (size_t j = 0; j < phys_param->pui_nW; j++)
		{
			Vh = ((j + 0.5) * phys_param->pui_wR / phys_param->pui_nW);
			ff = this->pui_Sm[j];
			if (ff <= 0.0) continue;
			for (size_t k = 1; k <= 40; k++)
			{
				the = dthe * k;
				d = sqrt(kv(Vh) + kv(2) - 2.0 * w * Vh * cos(the));
				if (d > 0.000000001)
				{
					pui_Sm2[ij] += ff * d * sigma(d) * sin(the) * dthe * 2.0 * const_pi;
				}
			}
		}
		pui_Sm2[ij] = pui_Sm2[ij] / (4.0 * const_pi);
	}

	for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
	{
		this->pui_Sm[ij] = pui_Sm2[ij] / phys_param->par_Kn;
	}

}

void Cell::MK_normir_Moments(Phys_param* phys_param)
{
	if (this->parameters[0].find("MK_n_H") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H"] /= this->volume[0];
	}

	if (this->parameters[0].find("MK_n_H1") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H1"] /= this->volume[0];
	}

	if (this->parameters[0].find("MK_n_H2") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H2"] /= this->volume[0];
	}

	if (this->parameters[0].find("MK_n_H3") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H3"] /= this->volume[0];
	}

	if (this->parameters[0].find("MK_n_H4") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H4"] /= this->volume[0];
	}

	if (phys_param->MK_source_S == true)
	{
		for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
		{
			double pui_w1 = ij * phys_param->pui_wR / phys_param->pui_nW;
			double pui_w2 = (ij + 1) * phys_param->pui_wR / phys_param->pui_nW;
			this->pui_Sm[ij] /= (this->volume[0]);
			this->pui_Sp[ij] /= (this->volume[0] * 4.0 * const_pi * (1.0 / 3.0) * (pow3(pui_w2) - pow3(pui_w1)));
		}
	}
}
