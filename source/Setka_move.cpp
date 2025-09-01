#include "Setka.h"

// Макрос параболического сглаживания HP
#define  macros1() p[0] = AA->coord[now][0]; \
p[1] = norm2(0.0, AA->coord[now][1], AA->coord[now][2]);\
\
p1[0] = AA1->coord[now][0];\
p1[1] = norm2(0.0, AA1->coord[now][1], AA1->coord[now][2]);\
\
p11[0] = AA11->coord[now][0];\
p11[1] = norm2(0.0, AA11->coord[now][1], AA11->coord[now][2]);\
\
p3[0] = AA3->coord[now][0];\
p3[1] = norm2(0.0, AA3->coord[now][1], AA3->coord[now][2]);\
\
p33[0] = AA33->coord[now][0];\
p33[1] = norm2(0.0, AA33->coord[now][1], AA33->coord[now][2]);\
\
solveQuadraticEquation(p3[0], p3[1],\
	p1[0], p1[1], p11[0], p11[1], a, b, c);\
r1 = a * kv(p[0]) + b * p[0] + c;\
\
solveQuadraticEquation(p33[0], p33[1],\
	p3[0], p3[1], p1[0], p1[1], a, b, c);\
r2 = a * kv(p[0]) + b * p[0] + c;\
\
p2[0] = polar_angle(AA2->coord[now][1], AA2->coord[now][2]);\
p2[1] = norm2(0.0, AA2->coord[now][1], AA2->coord[now][2]);\
\
p[0] = polar_angle(AA->coord[now][1], AA->coord[now][2]);\
\
p4[0] = polar_angle(AA4->coord[now][1], AA4->coord[now][2]);\
p4[1] = norm2(0.0, AA4->coord[now][1], AA4->coord[now][2]);\
\
p22[0] = polar_angle(AA22->coord[now][1], AA22->coord[now][2]);\
p22[1] = norm2(0.0, AA22->coord[now][1], AA22->coord[now][2]);\
\
p44[0] = polar_angle(AA44->coord[now][1], AA44->coord[now][2]);\
p44[1] = norm2(0.0, AA44->coord[now][1], AA44->coord[now][2]);\
\
solveQuadraticEquation(p2[0], p2[1],\
	p4[0], p4[1], p44[0], p44[1], a, b, c);\
r3 = a * kv(p[0]) + b * p[0] + c;\
\
solveQuadraticEquation(p4[0], p4[1],\
	p2[0], p2[1], p22[0], p22[1], a, b, c);\
r4 = a * kv(p[0]) + b * p[0] + c;\
\
A << AA->coord[now][0], AA->coord[now][1], AA->coord[now][2];\
\
B = A;\
B[1] *= r1 / p[1];\
B[2] *= r1 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_along * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r2 / p[1];\
B[2] *= r2 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_along * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r3 / p[1];\
B[2] *= r3 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_angle * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r4 / p[1];\
B[2] *= r4 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_angle * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();


#define  macros2(p, AA) p[0] = polar_angle(AA->coord[now][0],\
	norm2(0.0, AA->coord[now][1], AA->coord[now][2]));\
	if (p[0] > const_pi) p[0] = p[0] - 2 * const_pi;\
	p[1] = norm2(AA->coord[now][0], AA->coord[now][1], AA->coord[now][2]);



void Setka::Culc_Velocity_surface(short int now, const double& time, short int metod)
{
	double dsr, dsc, dsl;

	int now2 = (now + 1) % 2;

	for (auto& yz : this->All_Yzel)
	{
		yz->num_velocity = 0;
		yz->velocity[0] = 0.0;
		yz->velocity[1] = 0.0;
		yz->velocity[2] = 0.0;
	}

	// Считаем движение TS
	if (this->phys_param->move_TS == true)
	{
		#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
		{
			auto gr = this->Gran_TS[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];


			unordered_map<string, double> par_left;
			unordered_map<string, double> par_right;

			Snos_on_Gran(gr, par_left, par_right, now);

			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_TS";

			//A->parameters[now]
			qqq1[0] = par_left["rho"];
			qqq1[1] = par_left["Vx"];
			qqq1[2] = par_left["Vy"];
			qqq1[3] = par_left["Vz"];
			qqq1[4] = par_left["p"];
			qqq1[5] = par_left["Bx"];
			qqq1[6] = par_left["By"];
			qqq1[7] = par_left["Bz"];

			qqq2[0] = par_right["rho"];
			qqq2[1] = par_right["Vx"];
			qqq2[2] = par_right["Vy"];
			qqq2[3] = par_right["Vz"];
			qqq2[4] = par_right["p"];
			qqq2[5] = par_right["Bx"];
			qqq2[6] = par_right["By"];
			qqq2[7] = par_right["Bz"];

			double w = 0.0;


			this->phys_param->chlld(gr->Get_method(), gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += this->phys_param->velocity_TS * dsl * gr->normal[now][0];
				yz->velocity[1] += this->phys_param->velocity_TS * dsl * gr->normal[now][1];
				yz->velocity[2] += this->phys_param->velocity_TS * dsl * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}

	// Считаем движение HP
	if (this->phys_param->move_HP == true)
	{
#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
		{
			auto gr = this->Gran_HP[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];
			short int metod_ = gr->Get_method();

			if (gr->type2 != Type_Gran_surf::HP) cout << "Error 7823467276345679264978234" << endl;


			if (this->phys_param->null_bn_on_HP == true)
			{
				// Обнуляем bn в ячейках
				Eigen::Vector3d nnn;
				nnn << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];
				double dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
					A->parameters[now]["Bx"], A->parameters[now]["By"],
					A->parameters[now]["Bz"]);
				A->parameters[now]["Bx"] -= dtpr * nnn[0];
				A->parameters[now]["By"] -= dtpr * nnn[1];
				A->parameters[now]["Bz"] -= dtpr * nnn[2];

				dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
					B->parameters[now]["Bx"], B->parameters[now]["By"],
					B->parameters[now]["Bz"]);
				B->parameters[now]["Bx"] -= dtpr * nnn[0];
				B->parameters[now]["By"] -= dtpr * nnn[1];
				B->parameters[now]["Bz"] -= dtpr * nnn[2];
			}

			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_HP";

			// Получаем парамытры для снесённых значений
			unordered_map<string, double> par_left, par_right;
			this->Snos_on_Gran(gr, par_left, par_right, now);

			if (false) // Параметры в центрах ячеек
			{
				qqq1[0] = A->parameters[now]["rho"];
				qqq1[1] = A->parameters[now]["Vx"];
				qqq1[2] = A->parameters[now]["Vy"];
				qqq1[3] = A->parameters[now]["Vz"];
				qqq1[4] = A->parameters[now]["p"];
				qqq1[5] = A->parameters[now]["Bx"];
				qqq1[6] = A->parameters[now]["By"];
				qqq1[7] = A->parameters[now]["Bz"];

				qqq2[0] = B->parameters[now]["rho"];
				qqq2[1] = B->parameters[now]["Vx"];
				qqq2[2] = B->parameters[now]["Vy"];
				qqq2[3] = B->parameters[now]["Vz"];
				qqq2[4] = B->parameters[now]["p"];
				qqq2[5] = B->parameters[now]["Bx"];
				qqq2[6] = B->parameters[now]["By"];
				qqq2[7] = B->parameters[now]["Bz"];
			}
			else
			{
				qqq1[0] = par_left["rho"];
				qqq1[1] = par_left["Vx"];
				qqq1[2] = par_left["Vy"];
				qqq1[3] = par_left["Vz"];
				qqq1[4] = par_left["p"];
				qqq1[5] = par_left["Bx"];
				qqq1[6] = par_left["By"];
				qqq1[7] = par_left["Bz"];

				qqq2[0] = par_right["rho"];
				qqq2[1] = par_right["Vx"];
				qqq2[2] = par_right["Vy"];
				qqq2[3] = par_right["Vz"];
				qqq2[4] = par_right["p"];
				qqq2[5] = par_right["Bx"];
				qqq2[6] = par_right["By"];
				qqq2[7] = par_right["Bz"];

				// Для снесённых значений надо также обнулить Bn
				if (this->phys_param->null_bn_on_HP == true)
				{
					// Обнуляем bn в ячейках
					Eigen::Vector3d nnn;
					nnn << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];
					double dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
						qqq1[5], qqq1[6],
						qqq1[7]);
					qqq1[5] -= dtpr * nnn[0];
					qqq1[6] -= dtpr * nnn[1];
					qqq1[7] -= dtpr * nnn[2];

					dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
						qqq2[5], qqq2[6], qqq2[7]);
					qqq2[5] -= dtpr * nnn[0];
					qqq2[6] -= dtpr * nnn[1];
					qqq2[7] -= dtpr * nnn[2];
				}
			}

			double w = 0.0;

			// Записываем магнитное давление в обычное
			// И удаляем магнитные поля
			if (this->phys_param->bn_in_p_on_HP == true)
			{
				metod_ = 2;

				qqq1[4] += kvv(qqq1[5], qqq1[6], qqq1[7]) / (8.0 * const_pi);
				qqq1[5] = qqq1[6] = qqq1[7] = 0.0;

				qqq2[4] += kvv(qqq2[5], qqq2[6], qqq2[7]) / (8.0 * const_pi);
				qqq2[5] = qqq2[6] = qqq2[7] = 0.0;
			}

			if (this->phys_param->bn_in_p_on_HP == true)
			{
				//cout << "A" << endl;
				std::vector<double> n(3);
				n[0] = gr->normal[now][0];
				n[1] = gr->normal[now][1];
				n[2] = gr->normal[now][2];

				this->phys_param->Godunov_Solver_Alexashov(qqq1, qqq2,//
					n, qqq, dsl, dsr, dsc, w);
				qqq[5] = qqq[6] = qqq[7] = 0.0;
			}
			else
			{
				//cout << "B" << endl;
				this->phys_param->chlld(metod_, gr->normal[now][0], gr->normal[now][1],
					gr->normal[now][2],
					w, qqq1, qqq2, qqq, false, 3,
					konvect_left, konvect_right, konvect, dsr, dsc, dsl,
					Option);
			}

			if (std::isnan(dsc) || std::fpclassify(dsc) == FP_SUBNORMAL)
			{
				cout << "Error  5436542867" << endl;
				whach(dsc);
				whach(dsl);
				whach(dsr);
				whach(qqq1[0]);
				whach(qqq2[0]);
				whach(qqq[0]);
				whach(qqq[1]);
				whach(qqq[2]);
				whach(qqq[3]);
				whach(qqq[4]);
				exit(-1);
			}

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += this->phys_param->velocity_HP * dsc * gr->normal[now][0];
				yz->velocity[1] += this->phys_param->velocity_HP * dsc * gr->normal[now][1];
				yz->velocity[2] += this->phys_param->velocity_HP * dsc * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}

	// Считаем движение BS
	if (this->phys_param->move_BS == true)
	{
#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_BS.size(); i_step++)
		{
			auto gr = this->Gran_BS[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];


			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_BS";

			// Получаем параметры для снесённых значений
			unordered_map<string, double> par_left, par_right;
			this->Snos_on_Gran(gr, par_left, par_right, now);

			if (false) // Параметры в центрах ячеек
			{
				qqq1[0] = A->parameters[now]["rho"];
				qqq1[1] = A->parameters[now]["Vx"];
				qqq1[2] = A->parameters[now]["Vy"];
				qqq1[3] = A->parameters[now]["Vz"];
				qqq1[4] = A->parameters[now]["p"];
				qqq1[5] = A->parameters[now]["Bx"];
				qqq1[6] = A->parameters[now]["By"];
				qqq1[7] = A->parameters[now]["Bz"];

				qqq2[0] = B->parameters[now]["rho"];
				qqq2[1] = B->parameters[now]["Vx"];
				qqq2[2] = B->parameters[now]["Vy"];
				qqq2[3] = B->parameters[now]["Vz"];
				qqq2[4] = B->parameters[now]["p"];
				qqq2[5] = B->parameters[now]["Bx"];
				qqq2[6] = B->parameters[now]["By"];
				qqq2[7] = B->parameters[now]["Bz"];
			}
			else
			{
				qqq1[0] = par_left["rho"];
				qqq1[1] = par_left["Vx"];
				qqq1[2] = par_left["Vy"];
				qqq1[3] = par_left["Vz"];
				qqq1[4] = par_left["p"];
				qqq1[5] = par_left["Bx"];
				qqq1[6] = par_left["By"];
				qqq1[7] = par_left["Bz"];

				qqq2[0] = par_right["rho"];
				qqq2[1] = par_right["Vx"];
				qqq2[2] = par_right["Vy"];
				qqq2[3] = par_right["Vz"];
				qqq2[4] = par_right["p"];
				qqq2[5] = par_right["Bx"];
				qqq2[6] = par_right["By"];
				qqq2[7] = par_right["Bz"];
			}


			double w = 0.0;

			this->phys_param->chlld(gr->Get_method(), gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += this->phys_param->velocity_BS * dsr * gr->normal[now][0];
				yz->velocity[1] += this->phys_param->velocity_BS * dsr * gr->normal[now][1];
				yz->velocity[2] += this->phys_param->velocity_BS * dsr * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}


	// Сглаживание поверхностей (уже используя новые координаты узлов на поверхности)

	if (this->phys_param->sglag_TS == true)
	{
		//double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
		//	norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
		//		this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

		Eigen::Vector3d A, B, B2, V;

		// Сглаживание в головной области х > 0
		// Здесть просто Лаплас в сферических СК. или декартовых
		// Лаплас в декартовых работает очень плохо и "сплющивает" поверхность
		if (true) // Это старый вариант сглаживания, сейчас работает другой
		{
			if (true)
			{   // Лаплас в сферических новых
#pragma omp parallel for private(A, B, V, B2) schedule(dynamic)
				for(size_t ili = 0; ili < this->All_Yzel.size(); ili++)
				{
					auto yz = this->All_Yzel[ili];
					if (yz->type != Type_yzel::TS) continue;
					A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];

					//double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// Включаем радиально-сферическое сглаживание только в головной зоне
					//if (phi > phi_a + const_pi / 180.0) continue;

					B << yz->parameters["xc"], yz->parameters["yc"], yz->parameters["zc"];

					double rr = (A - B).norm();
					double r = 0.0;
					for (auto& j : yz->Yzel_sosed_sglag2)
					{
						B2 << j->coord[now][0], j->coord[now][1], j->coord[now][2];
						r += (B2 - B).norm();
					}
					r /= yz->Yzel_sosed_sglag2.size();

					double k = (r / rr);
					B2 = k * (A - B) + B;

					V = this->phys_param->velocity_TS * 
						this->phys_param->sglag_TS_k_sphere * (B2 - A) / time;

					yz->mut.lock();
					yz->velocity[0] += V(0);
					yz->velocity[1] += V(1);
					yz->velocity[2] += V(2);
					yz->num_velocity++;
					yz->mut.unlock();

				}
			}
		}

		if (false)
		{
			double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
				norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
					this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

			#pragma omp parallel for
			for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
			{
				Eigen::Vector3d A, B, V;
				auto gr = this->Gran_TS[i_step];
				A << gr->center[now][0], gr->center[now][1], gr->center[now][2];
				double rr = A.norm();
				double r = 0.0;
				for (auto& j : gr->grans_surf)
				{
					r += norm2(j->center[now][0], j->center[now][1], j->center[now][2]);
				}
				r /= gr->grans_surf.size();

				B = A * r / rr;

				double pk = this->phys_param->sglag_TS_k;
				double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
				if (phi > this->geo->tetta1)
				{
					pk = this->phys_param->sglag_TS_k_sphere_tail;
				}
				else if (phi < this->geo->tetta0)
				{
					pk = this->phys_param->sglag_TS_k_sphere_head;
				}

				V = this->phys_param->velocity_TS * pk * (B - A) / time;

				for (auto& yz : gr->yzels)
				{
					yz->mut.lock();
					yz->velocity[0] += V(0);
					yz->velocity[1] += V(1);
					yz->velocity[2] += V(2);
					yz->num_velocity++;
					yz->mut.unlock();
				}
			}
		}
		
	}


	if (this->phys_param->sglag_HP == true)
	{
		vector<vector <vector<Luch*>>> VVV;
		vector<int> VVV_int;

		VVV.push_back(this->D_Luch);
		VVV.push_back(this->E_Luch);
		VVV.push_back(this->B_Luch);
		VVV_int.push_back(1);
		VVV_int.push_back(1);
		VVV_int.push_back(2);

		double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0], 
			norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1], 
				this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

		Eigen::Vector3d A, B, B2, V;

		// Сглаживание в головной области х > 0
		// Здесть просто Лаплас в сферических СК. или декартовых
		// Лаплас в декартовых работает очень плохо и "сплющивает" поверхность
		if (true) // Это старый вариант сглаживания, сейчас работает другой
		{
			if (true)
			{   // Сглаживание к локальной сфере
#pragma omp parallel for private(A, B, V, B2) schedule(dynamic)
				for (size_t ili = 0; ili < this->All_Yzel.size(); ili++)
				{
					auto yz = this->All_Yzel[ili];
					if (yz->type != Type_yzel::HP) continue;
					A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
					
					double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// Включаем радиально-сферическое сглаживание только в головной зоне
					if (phi > phi_a + const_pi / 180.0) continue;

					B << yz->parameters["xc"], yz->parameters["yc"], yz->parameters["zc"];
						
					double rr = (A - B).norm();
					double r = 0.0;
					for (auto& j : yz->Yzel_sosed_sglag2)
					{
						B2 << j->coord[now][0], j->coord[now][1], j->coord[now][2];
						r += (B2 - B).norm();
					}
					r /= yz->Yzel_sosed_sglag2.size();

					double k = (r / rr);
					B2 = k * (A - B) + B;

					V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_k_sphere * (B2 - A) / time;
					//V = this->phys_param->sglag_HP_k_sphere * (B - A) / time;

					yz->mut.lock();
					yz->velocity[0] += V(0);
					yz->velocity[1] += V(1);
					yz->velocity[2] += V(2);
					yz->num_velocity++;
					yz->mut.unlock();

				}
			}
			if (false)
			{   // Лаплас в сферических
				#pragma omp parallel for private(A, B, V)
				for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
				{
					auto gr = this->Gran_HP[i_step];
					A << gr->center[now][0], gr->center[now][1], gr->center[now][2];

					double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// Включаем радиально-сферическое сглаживание только в головной зоне
					if (phi > phi_a + const_pi / 180.0) continue;


					double rr = A.norm();
					double r = 0.0;
					for (auto& j : gr->grans_surf)
					{
						r += norm2(j->center[now][0], j->center[now][1], j->center[now][2]);
					}
					r /= gr->grans_surf.size();

					B = A * r / rr;

					V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_k_sphere * (B - A) / time;

					for (auto& yz : gr->yzels)
					{
						yz->mut.lock();
						yz->velocity[0] += V(0);
						yz->velocity[1] += V(1);
						yz->velocity[2] += V(2);
						yz->num_velocity++;
						yz->mut.unlock();
					}
				}
			}
			if (false)
			{
				// Лаплас в декартовых
				#pragma omp parallel for private(A, B, V)
				for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
				{
					auto gr = this->Gran_HP[i_step];
					A << gr->center[now][0], gr->center[now][1], gr->center[now][2];
					
					double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// Включаем радиально-сферическое сглаживание только в головной зоне
					if (phi > phi_a + const_pi / 180.0) continue;
					B = 2 * A;

					for (auto& j : gr->grans_surf)
					{
						B[0] += j->center[now][0];
						B[1] += j->center[now][1];
						B[2] += j->center[now][2];
					}
					B /= (gr->grans_surf.size() + 2.0);


					V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_k_sphere * (B - A) / time;

					for (auto& yz : gr->yzels)
					{
						yz->mut.lock();
						yz->velocity[0] += V(0);
						yz->velocity[1] += V(1);
						yz->velocity[2] += V(2);
						yz->num_velocity++;
						yz->mut.unlock();
					}
				}
			}

		}


		// Сглаживание к центру локальных сфер для HP в хвосте
		for (size_t ii = 0; ii < VVV.size(); ii++)
		{
			int nn = VVV[ii].size();
			int mm = VVV[ii][0].size();

#pragma omp parallel for private(A, B, V) schedule(dynamic)
			for (size_t i = 0; i < nn; i++)
			{
				if (i > 3 && i < nn - 3) continue;

				for (size_t j = 0; j < mm; j++)
				{
					auto AA = VVV[ii][i][j]->Yzels_opor[VVV_int[ii]];

					A << AA->coord[now][0], AA->coord[now][1], AA->coord[now][2];

					Eigen::Vector3d center;
					center[0] = A[0];
					center[1] = AA->parameters["xcc"];
					center[2] = AA->parameters["ycc"];
						
					
					double rr = (A - center).norm();
					double r = 0.0;
					for (auto& j : AA->Yzel_sosed_sglag2)
					{
						B << A[0], j->coord[now][1], j->coord[now][2];
						r += (center - B).norm();
					}
					r /= AA->Yzel_sosed_sglag2.size();

					double k = (r / rr);
					B = k * (A - center) + center;

					V = this->phys_param->velocity_HP * 
						this->phys_param->sglag_HP_k_angle * (B - A) / time;

					AA->mut.lock();
					AA->velocity[0] += V(0);
					AA->velocity[1] += V(1);
					AA->velocity[2] += V(2);
					AA->num_velocity++;
					AA->mut.unlock();
				}
			}
		}



		Yzel* AA, * AA1, * AA11, * AA2, * AA22, * AA3, * AA33, * AA4, * AA44;
		// AA33  AA3  AA  AA1  AA11 - вдоль HP  направление в апвинд
		// A22  AA2  AA  AA4  AA44  - вращение HP
		std::array<double, 2> p, p1, p11, p3, p33, p2, p22, p4, p44;
		double a, b, c;
		double r1, r2, r3, r4;

		// Сглаживание в головной области (новое - сплайнами)
		if (false)
		{
#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
			//for (auto& yz : this->Yzels_HP_sglag)
			for (size_t idx = 0; idx < this->Yzels_HP_sglag.size(); ++idx)
			{
				auto& yz = this->Yzels_HP_sglag[idx];
				AA = yz;
				AA1 = yz->Yzel_sosed_sglag["AA1"];
				AA11 = yz->Yzel_sosed_sglag["AA11"];
				AA2 = yz->Yzel_sosed_sglag["AA3"];
				AA22 = yz->Yzel_sosed_sglag["AA33"];
				AA3 = yz->Yzel_sosed_sglag["AA2"];
				AA33 = yz->Yzel_sosed_sglag["AA22"];
				AA4 = yz->Yzel_sosed_sglag["AA4"];
				AA44 = yz->Yzel_sosed_sglag["AA44"];


				macros2(p, AA);
				macros2(p1, AA1);
				macros2(p11, AA11);
				macros2(p3, AA3);
				macros2(p33, AA33);

				solveQuadraticEquation(p3[0], p3[1],
					p1[0], p1[1], p11[0], p11[1], a, b, c);
				r1 = a * kv(p[0]) + b * p[0] + c;

				solveQuadraticEquation(p33[0], p33[1],
					p3[0], p3[1], p1[0], p1[1], a, b, c);
				r2 = a * kv(p[0]) + b * p[0] + c;


				p2[0] = polar_angle(AA2->coord[now][1], AA2->coord[now][2]);
				p2[1] = norm2(AA2->coord[now][0], AA2->coord[now][1], AA2->coord[now][2]);

				p[0] = polar_angle(AA->coord[now][1], AA->coord[now][2]);

				p4[0] = polar_angle(AA4->coord[now][1], AA4->coord[now][2]);
				p4[1] = norm2(AA4->coord[now][0], AA4->coord[now][1], AA4->coord[now][2]);

				p22[0] = polar_angle(AA22->coord[now][1], AA22->coord[now][2]);
				p22[1] = norm2(AA22->coord[now][0], AA22->coord[now][1], AA22->coord[now][2]);

				p44[0] = polar_angle(AA44->coord[now][1], AA44->coord[now][2]);
				p44[1] = norm2(AA44->coord[now][0], AA44->coord[now][1], AA44->coord[now][2]);

				solveQuadraticEquation(p2[0], p2[1],
					p4[0], p4[1], p44[0], p44[1], a, b, c);
				r3 = a * kv(p[0]) + b * p[0] + c;

				solveQuadraticEquation(p4[0], p4[1],
					p2[0], p2[1], p22[0], p22[1], a, b, c);
				r4 = a * kv(p[0]) + b * p[0] + c;

				A << AA->coord[now][0], AA->coord[now][1], AA->coord[now][2];

				B = A;
				B[0] *= r1 / p[1];
				B[1] *= r1 / p[1];
				B[2] *= r1 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time;

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r2 / p[1];
				B[1] *= r2 / p[1];
				B[2] *= r2 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r3 / p[1];
				B[1] *= r3 / p[1];
				B[2] *= r3 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r4 / p[1];
				B[1] *= r4 / p[1];
				B[2] *= r4 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

			}
		}


		// Параболическое сглаживание на A лучах
		int nn = this->A_Luch.size();
		int mm = this->A_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 2; j < mm; j++)
			{
				AA = this->A_Luch[i][j]->Yzels_opor[2];
				AA2 = this->A_Luch[im][j]->Yzels_opor[2];
				AA22 = this->A_Luch[imm][j]->Yzels_opor[2];
				AA4 = this->A_Luch[ip][j]->Yzels_opor[2];
				AA44 = this->A_Luch[ipp][j]->Yzels_opor[2];

				if (j == mm - 1)
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][0]->Yzels_opor[2];
					AA33 = this->B_Luch[i][1]->Yzels_opor[2];
				}
				else if (j == mm - 2)
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->A_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][0]->Yzels_opor[2];
				}
				else
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->A_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->A_Luch[i][j + 2]->Yzels_opor[2];
				}

				macros1();
			}
		}


		// Параболическое сглаживание на B лучах
		nn = this->B_Luch.size();
		mm = this->B_Luch[0].size();
#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->B_Luch[i][j]->Yzels_opor[2];
				AA2 = this->B_Luch[im][j]->Yzels_opor[2];
				AA22 = this->B_Luch[imm][j]->Yzels_opor[2];
				AA4 = this->B_Luch[ip][j]->Yzels_opor[2];
				AA44 = this->B_Luch[ipp][j]->Yzels_opor[2];

				if (j == 0)
				{
					int nnn = this->A_Luch[0].size();
					AA1 = this->A_Luch[i][nnn - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][nnn - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}
				else if (j == 1)
				{
					int nnn = this->A_Luch[0].size();
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][nnn - 1]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}
				else if (j == mm - 1)
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->E_Luch[i][0]->Yzels_opor[1];
					AA33 = this->E_Luch[i][1]->Yzels_opor[1];
				}
				else if (j == mm - 2)
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->E_Luch[i][0]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}

				macros1();
			}
		}
	
		// Параболическое сглаживание на E лучах
		nn = this->E_Luch.size();
		mm = this->E_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->E_Luch[i][j]->Yzels_opor[1];
				AA2 = this->E_Luch[im][j]->Yzels_opor[1];
				AA22 = this->E_Luch[imm][j]->Yzels_opor[1];
				AA4 = this->E_Luch[ip][j]->Yzels_opor[1];
				AA44 = this->E_Luch[ipp][j]->Yzels_opor[1];

				if (j == 0)
				{
					int nnn = this->B_Luch[0].size();
					AA1 = this->B_Luch[i][nnn - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][nnn - 2]->Yzels_opor[2];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == 1)
				{
					int nnn = this->B_Luch[0].size();
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->B_Luch[i][nnn - 1]->Yzels_opor[2];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == mm - 1)
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][0]->Yzels_opor[1];
					AA33 = this->D_Luch[i][1]->Yzels_opor[1];
				}
				else if (j == mm - 2)
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][0]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				macros1();
			}
		}

		// Параболическое сглаживание на D лучах
		nn = this->D_Luch.size();
		mm = this->geo->N4 - 1; // this->D_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->D_Luch[i][j]->Yzels_opor[1];
				AA2 = this->D_Luch[im][j]->Yzels_opor[1];
				AA22 = this->D_Luch[imm][j]->Yzels_opor[1];
				AA4 = this->D_Luch[ip][j]->Yzels_opor[1];
				AA44 = this->D_Luch[ipp][j]->Yzels_opor[1];

				if (j == 0)
				{
					int nnn = this->E_Luch[0].size();
					AA1 = this->E_Luch[i][nnn - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][nnn - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == 1)
				{
					int nnn = this->E_Luch[0].size();
					AA1 = this->D_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][nnn - 1]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->D_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->D_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				macros1();
			}
		}

	}


	if (this->phys_param->sglag_BS == true)
	{
		Eigen::Vector3d A, B, V;

		// Сглаживание по серединам граней
		if (false)
		{
			// Лаплас в декартовых
#pragma omp parallel for private(A, B, V)
			for (int i_step = 0; i_step < this->Gran_BS.size(); i_step++)
			{
				auto gr = this->Gran_BS[i_step];
				A << gr->center[now][0], gr->center[now][1], gr->center[now][2];

				B = A;

				if (gr->grans_surf.size() < 4) continue;    // TODO!

				for (auto& j : gr->grans_surf)
				{
					B[0] += j->center[now][0];
					B[1] += j->center[now][1];
					B[2] += j->center[now][2];
				}
				B /= (gr->grans_surf.size() + 1.0);


				V = this->phys_param->velocity_BS * this->phys_param->sglag_BS_k * (B - A) / time;

				for (auto& yz : gr->yzels)
				{
					yz->mut.lock();
					yz->velocity[0] += V(0);
					yz->velocity[1] += V(1);
					yz->velocity[2] += V(2);
					yz->num_velocity++;
					yz->mut.unlock();
				}
			}
		}

		// Сглаживание по узлам
		if (true)
		{
#pragma omp parallel for private(A, B, V) schedule(dynamic)
			for (size_t ili = 0; ili < this->All_Yzel.size(); ili++)
			{
				auto yz = this->All_Yzel[ili];
				if (yz->type != Type_yzel::BS) continue;

				A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];

				B = A;

				if (yz->Yzel_sosed_sglag2.size() <= 5) continue;    // TODO!

				for (auto& j : yz->Yzel_sosed_sglag2)
				{
					B[0] += j->coord[now][0];
					B[1] += j->coord[now][1];
					B[2] += j->coord[now][2];
				}
				B /= (yz->Yzel_sosed_sglag2.size() + 1.0);


				V = this->phys_param->velocity_BS * this->phys_param->sglag_BS_k * (B - A) / time;

				yz->mut.lock();
				yz->velocity[0] += V(0);
				yz->velocity[1] += V(1);
				yz->velocity[2] += V(2);
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}

	// Вычисляем новые координаты узлов на поверхности
	for (auto& yz : this->All_Yzel)
	{
		if (yz->num_velocity > 0)
		{
			yz->velocity[0] /= yz->num_velocity;
			yz->velocity[1] /= yz->num_velocity;
			yz->velocity[2] /= yz->num_velocity;
		}

		// Вычисляем новую координату узла
		if (yz->type == Type_yzel::TS && this->phys_param->move_TS == true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			B = A;
			B.normalize();
			A = A + V.dot(B) * time * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
		else if (yz->type == Type_yzel::HP && this->phys_param->move_HP == true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			if (A(0) >= 0.0)
			{
				B = A;
			}
			else
			{
				B << 0.0, A(1), A(2);
			}
			B.normalize();
			A = A + V.dot(B) * time * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
		else if (yz->type == Type_yzel::BS && this->phys_param->move_BS == true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			B = A;
			B.normalize();
			A = A + V.dot(B) * time * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
	}

	// Остальные узлы на HP (невыделяемой части) надо подвинуть
	if (this->phys_param->move_HP == true)
	{
		short int NN = this->D_Luch[0].size() - 1;
		for (auto& L : this->D_Luch)
		{

			double h1 = norm2(0.0, L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][1],
				L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][2]);

			// Высота контакта в хвосте (остаётся постоянной)
			//double h2 = norm2(0.0, L[NN]->Yzels_opor[1]->coord[now2][1],
			//	L[NN]->Yzels_opor[1]->coord[now2][2]);

			for (short int i = this->geo->N4 - 3; i <= NN; i++)
			{
				//double h = h1 + (i - this->geo->N4 + 2) * (h2 - h1) / (NN - this->geo->N4 + 2);
				double h = h1;
				auto yz = L[i]->Yzels_opor[1];
				double hh = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
				if (hh < 0.0000001 || h < 0.0000001 || std::isnan(hh) || std::isnan(h) ||
					std::fpclassify(h) == FP_SUBNORMAL || std::fpclassify(hh) == FP_SUBNORMAL)
					// || fabs(1.0 - h / hh) > 0.0001)
				{
					cout << "Error 8563529613" << endl;
					whach(h);
					whach(hh);
					whach(h1);
					//whach(h2);
				}


				yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
				yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
			}
		}
	}
	
	// Надо подвинуть узлы, которые продолжают BS
	if (this->phys_param->move_BS == true)
	{
		short int NN = this->A_Luch[0].size() - 1;
		for (short int i = 0; i < this->A_Luch.size(); i++)
		{
			auto yz = this->A_Luch[i][NN]->Yzels_opor[3];
			double H = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
			for (auto& L : this->B_Luch[i])
			{
				auto yyz = L->Yzels_opor[3];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
			for (auto& L : this->E_Luch[i])
			{
				auto yyz = L->Yzels_opor[2];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
			for (auto& L : this->D_Luch[i])
			{
				auto yyz = L->Yzels_opor[2];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
		}


	}
}

void Setka::Smooth_head_TS(void)
{
	cout << "Start Smooth_head_TS" << endl;

	std::vector<double> a; // коэффициенты a1, a2, a3
	std::vector<double> b; // коэффициенты b1, b2, b3
	std::vector<double> c; // коэффициенты b1, b2, b3
	std::vector<double> d; // коэффициенты b1, b2, b3
	std::vector<double> e; // коэффициенты b1, b2, b3
	std::vector<double> f; // коэффициенты b1, b2, b3
	std::vector<double> w; // веса
	Eigen::Vector3d A, B, V;

	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			a.push_back(kv(A[0]));
			b.push_back(kv(A[1]));
			c.push_back(kv(A[2]));
			d.push_back(A[0] * A[1]);
			e.push_back(A[0] * A[2]);
			f.push_back(A[2] * A[1]);
			w.push_back(0.05 + phi);
		}
	}

	unsigned int n = a.size();
	Eigen::MatrixXd AA(n, 6);
	for (int i = 0; i < n; ++i) {
		AA(i, 0) = a[i] * w[i];  // коэффициент при x
		AA(i, 1) = b[i] * w[i];  // коэффициент при x
		AA(i, 2) = c[i] * w[i];  // коэффициент при x
		AA(i, 3) = d[i] * w[i];  // коэффициент при x
		AA(i, 4) = e[i] * w[i];  // коэффициент при x
		AA(i, 5) = f[i] * w[i];  // коэффициент при x
	}

	Eigen::VectorXd b_vec(n);
	for (int i = 0; i < n; ++i) {
		b_vec(i) = 1.0 * w[i];
	}

	// Решаем систему методом наименьших квадратов
	Eigen::VectorXd solution = AA.colPivHouseholderQr().solve(b_vec);

	// Извлекаем решение
	double x = solution(0);
	double y = solution(1);
	double z = solution(2);
	double zz = solution(3);
	double zzz = solution(4);
	double zzzz = solution(5);

	double ak = sqrt(1.0 / x);
	double bk = sqrt(1.0 / y);
	double ck = sqrt(1.0 / z);
	double dk = zz;
	double ek = zzz;
	double fk = zzzz;
	
	cout << "a = " << ak << endl;
	cout << "b = " << bk << endl;
	cout << "c = " << ck << endl;
	cout << "d = " << dk << endl;
	cout << "e = " << ek << endl;
	cout << "f = " << fk << endl;

	// Теперь сглаживаем поверхность

	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		for (auto& yz : gr->yzels)
		{
			A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			//double rr = A.norm();

			double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
			if (phi < this->geo->tetta0 + 0.17)
			{
				double rr = sqrt(kv(A[0]) / kv(ak) + kv(A[1]) / kv(bk)
					+ kv(A[2]) / kv(ck) + A[0] * A[1] * dk + 
					A[0] * A[2] * ek + A[2] * A[1] * fk);
				double aaa = 1.0/rr;
				yz->coord[0][0] *= aaa;
				yz->coord[0][1] *= aaa;
				yz->coord[0][2] *= aaa;
				yz->coord[1][0] = yz->coord[0][0];
				yz->coord[1][1] = yz->coord[0][1];
				yz->coord[1][2] = yz->coord[0][2];
			}
		}
	}


	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_TS" << endl;
}

void Setka::Smooth_head_TS2(void)
{
	cout << "Start Smooth_head_TS2" << endl;
	
	vector<unsigned int> nomera;
	vector<double> RRR;

	this->Renumerate();

	Eigen::Vector3d A, B, V;

	// Для больших r
	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz < r_yyz) goto b1;
					}
				}
				// Если дошли до сюда, значит этот узел действительно максимальный
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new/ i_yz_new);
			b1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr/r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}

	nomera.clear();
	RRR.clear();

	// Для малых R
	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz > r_yyz) goto bb1;
					}
				}
				// Если дошли до сюда, значит этот узел действительно максимальный
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new / i_yz_new);
			bb1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr / r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}



	// ---------------

	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_TS2" << endl;
}

void Setka::Smooth_head_HP2(void)
{
	cout << "Start Smooth_head_HP2" << endl;
	// Функция ручного сглаживания HP
	int now = 0;
	int now2 = 1;
	int M = this->A_Luch.size() - 1;
	int M1 = 0;
	int M2 = M - 1;

	for (auto& yz : this->All_Yzel)
	{
		yz->num_velocity = 0;
		yz->velocity[0] = 0.0;
		yz->velocity[1] = 0.0;
		yz->velocity[2] = 0.0;
	}

	int N = this->A_Luch[0].size() - 1;
	for (int i_step = 0; i_step < this->A_Luch[0].size(); i_step++)
	{
		auto A = this->A_Luch[M][i_step]->Yzels_opor[2];
		auto B = this->A_Luch[M1][i_step]->Yzels_opor[2];
		auto C = this->A_Luch[M2][i_step]->Yzels_opor[2];

		A->velocity[0] = - A->coord[0][0] + (B->coord[0][0] + C->coord[0][0]) / 2.0;
		A->velocity[1] = -A->coord[0][1] + (B->coord[0][1] + C->coord[0][1]) / 2.0;
		A->velocity[2] = -A->coord[0][2] + (B->coord[0][2] + C->coord[0][2]) / 2.0;
		A->num_velocity++;
	}

	N = this->B_Luch[0].size() - 1;
	for (int i_step = 0; i_step < this->B_Luch[0].size(); i_step++)
	{
		auto A = this->B_Luch[M][i_step]->Yzels_opor[2];
		auto B = this->B_Luch[M1][i_step]->Yzels_opor[2];
		auto C = this->B_Luch[M2][i_step]->Yzels_opor[2];

		A->velocity[0] = -A->coord[0][0] + (B->coord[0][0] + C->coord[0][0]) / 2.0;
		A->velocity[1] = -A->coord[0][1] + (B->coord[0][1] + C->coord[0][1]) / 2.0;
		A->velocity[2] = -A->coord[0][2] + (B->coord[0][2] + C->coord[0][2]) / 2.0;
		A->num_velocity++;
	}

	N = this->E_Luch[0].size() - 1;
	for (int i_step = 0; i_step < this->E_Luch[0].size(); i_step++)
	{
		auto A = this->E_Luch[M][i_step]->Yzels_opor[1];
		auto B = this->E_Luch[M1][i_step]->Yzels_opor[1];
		auto C = this->E_Luch[M2][i_step]->Yzels_opor[1];

		A->velocity[0] = -A->coord[0][0] + (B->coord[0][0] + C->coord[0][0]) / 2.0;
		A->velocity[1] = -A->coord[0][1] + (B->coord[0][1] + C->coord[0][1]) / 2.0;
		A->velocity[2] = -A->coord[0][2] + (B->coord[0][2] + C->coord[0][2]) / 2.0;
		A->num_velocity++;
	}

	N = this->D_Luch[0].size() - 1;
	for (int i_step = 0; i_step < this->geo->N4 - 3; i_step++)
	{
		auto A = this->D_Luch[M][i_step]->Yzels_opor[1];
		auto B = this->D_Luch[M1][i_step]->Yzels_opor[1];
		auto C = this->D_Luch[M2][i_step]->Yzels_opor[1];

		A->velocity[0] = -A->coord[0][0] + (B->coord[0][0] + C->coord[0][0]) / 2.0;
		A->velocity[1] = -A->coord[0][1] + (B->coord[0][1] + C->coord[0][1]) / 2.0;
		A->velocity[2] = -A->coord[0][2] + (B->coord[0][2] + C->coord[0][2]) / 2.0;
		A->num_velocity++;
	}


	// Вычисляем новые координаты узлов на поверхности
	for (auto& yz : this->All_Yzel)
	{
		if (yz->num_velocity > 0)
		{
			yz->velocity[0] /= yz->num_velocity;
			yz->velocity[1] /= yz->num_velocity;
			yz->velocity[2] /= yz->num_velocity;
		}
		else
		{
			continue;
		}

		if (true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			if (A(0) >= 0.0)
			{
				B = A;
			}
			else
			{
				B << 0.0, A(1), A(2);
			}
			B.normalize();
			A = A + V.dot(B) * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
	}

	// Остальные узлы на HP (невыделяемой части) надо подвинуть
	if (true)
	{
		short int NN = this->D_Luch[0].size() - 1;
		for (auto& L : this->D_Luch)
		{

			double h1 = norm2(0.0, L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][1],
				L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][2]);

			// Высота контакта в хвосте (остаётся постоянной)
			//double h2 = norm2(0.0, L[NN]->Yzels_opor[1]->coord[now2][1],
			//	L[NN]->Yzels_opor[1]->coord[now2][2]);

			for (short int i = this->geo->N4 - 3; i <= NN; i++)
			{
				//double h = h1 + (i - this->geo->N4 + 2) * (h2 - h1) / (NN - this->geo->N4 + 2);
				double h = h1;
				auto yz = L[i]->Yzels_opor[1];
				double hh = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
				if (hh < 0.0000001 || h < 0.0000001 || std::isnan(hh) || std::isnan(h) ||
					std::fpclassify(h) == FP_SUBNORMAL || std::fpclassify(hh) == FP_SUBNORMAL)
					// || fabs(1.0 - h / hh) > 0.0001)
				{
					cout << "Error 8563529613" << endl;
					whach(h);
					whach(hh);
					whach(h1);
					//whach(h2);
				}


				yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
				yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
			}
		}
	}


	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(1);
	}

	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[0][j] = i->coord[1][j];
		}
	}

	for (auto& yz : this->All_Yzel)
	{
		yz->num_velocity = 0;
		yz->velocity[0] = 0.0;
		yz->velocity[1] = 0.0;
		yz->velocity[2] = 0.0;
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_HP2" << endl;
}

void Setka::Smooth_head_HP(void)
{
	cout << "Start Smooth_head_HP" << endl;

	std::vector<double> a; // коэффициенты a1, a2, a3
	std::vector<double> b; // коэффициенты b1, b2, b3
	std::vector<double> c; // коэффициенты b1, b2, b3
	std::vector<double> d; // коэффициенты b1, b2, b3
	std::vector<double> e; // коэффициенты b1, b2, b3
	std::vector<double> f; // коэффициенты b1, b2, b3
	std::vector<double> w; // веса
	Eigen::Vector3d A, B, V;

	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			a.push_back(kv(A[0]));
			b.push_back(kv(A[1]));
			c.push_back(kv(A[2]));
			d.push_back(A[0] * A[1]);
			e.push_back(A[0] * A[2]);
			f.push_back(A[2] * A[1]);
			w.push_back(0.05 + phi);
		}
	}

	unsigned int n = a.size();
	Eigen::MatrixXd AA(n, 6);
	for (int i = 0; i < n; ++i) {
		AA(i, 0) = a[i] * w[i];  // коэффициент при x
		AA(i, 1) = b[i] * w[i];  // коэффициент при x
		AA(i, 2) = c[i] * w[i];  // коэффициент при x
		AA(i, 3) = d[i] * w[i];  // коэффициент при x
		AA(i, 4) = e[i] * w[i];  // коэффициент при x
		AA(i, 5) = f[i] * w[i];  // коэффициент при x
	}

	Eigen::VectorXd b_vec(n);
	for (int i = 0; i < n; ++i) {
		b_vec(i) = 1.0 * w[i];
	}

	// Решаем систему методом наименьших квадратов
	Eigen::VectorXd solution = AA.colPivHouseholderQr().solve(b_vec);

	// Извлекаем решение
	double x = solution(0);
	double y = solution(1);
	double z = solution(2);
	double zz = solution(3);
	double zzz = solution(4);
	double zzzz = solution(5);

	double ak = sqrt(1.0 / x);
	double bk = sqrt(1.0 / y);
	double ck = sqrt(1.0 / z);
	double dk = zz;
	double ek = zzz;
	double fk = zzzz;

	cout << "a = " << ak << endl;
	cout << "b = " << bk << endl;
	cout << "c = " << ck << endl;
	cout << "d = " << dk << endl;
	cout << "e = " << ek << endl;
	cout << "f = " << fk << endl;

	// Теперь сглаживаем поверхность

	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		for (auto& yz : gr->yzels)
		{
			A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			//double rr = A.norm();

			double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
			if (phi < this->geo->tetta0 + 0.17)
			{
				double rr = sqrt(kv(A[0]) / kv(ak) + kv(A[1]) / kv(bk)
					+ kv(A[2]) / kv(ck) + A[0] * A[1] * dk +
					A[0] * A[2] * ek + A[2] * A[1] * fk);
				double aaa = 1.0 / rr;
				yz->coord[0][0] *= aaa;
				yz->coord[0][1] *= aaa;
				yz->coord[0][2] *= aaa;
				yz->coord[1][0] = yz->coord[0][0];
				yz->coord[1][1] = yz->coord[0][1];
				yz->coord[1][2] = yz->coord[0][2];
			}
		}
	}


	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_HP" << endl;
}

void Setka::Smooth_head_HP3(void)
{
	int max_iterations = 8;
	double tolerance = 1e-6;
	double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
		norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
			this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

	cout << "Start Smooth_head_HP3" << endl;

	for (auto& yz : this->All_Yzel)
	{
		if (yz->type != Type_yzel::HP) continue;

		double phi = polar_angle(yz->coord[0][0], 
			norm2(0.0, yz->coord[0][1], yz->coord[0][2]));

		// Включаем радиально-сферическое сглаживание только в головной зоне
		if (phi > phi_a + const_pi / 180.0) continue;

		std::unordered_set<Yzel*> Yzels;
		Yzels.insert(yz);
		for (auto& gr : yz->grans)
		{
			if (gr->type2 != Type_Gran_surf::HP) continue;
			for (auto& yz2 : gr->yzels)
			{
				Yzels.insert(yz2);
			}
		}

		int n = Yzels.size();

		if (n < 3)
		{
			cout << "Error  " << yz->coord[0][0] << " " << yz->coord[0][1] << " "
				<< yz->coord[0][2] << endl;
			continue;
		}

		yz->Yzel_sosed_sglag2.clear();
		for (auto& i : Yzels)
		{
			if (i == yz) continue;
			yz->Yzel_sosed_sglag2.insert(i);
		}

		std::vector<Eigen::Vector3d> points;
		for (auto& i : Yzels)
		{
			Eigen::Vector3d AA;
			AA << i->coord[0][0],
				i->coord[0][1], i->coord[0][2];
			points.push_back(AA);
		}

		Eigen::MatrixXd A(n, 3);
		Eigen::VectorXd b(n);
		for (int i = 0; i < n; ++i) {
			A.row(i) = 2.0 * points[i];
			b(i) = points[i].squaredNorm();
		}

		// 2. Вычитаем уравнения для уменьшения влияния численных ошибок
		Eigen::Vector3d mean_point = Eigen::Vector3d::Zero();
		for (const auto& p : points) mean_point += p;
		mean_point /= n;

		Eigen::MatrixXd A_centered = A.rowwise() - A.colwise().mean();
		Eigen::VectorXd b_centered = b.array() - b.mean();

		Eigen::Vector3d center = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
		double radius = 0.0;

		// 2. Итеративное уточнение
		Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
		double prev_residual = std::numeric_limits<double>::max();

		for (int iter = 0; iter < max_iterations; ++iter) {
			// 2.1. Вычисляем текущий радиус как среднее расстояние до центра
			radius = 0.0;
			for (const auto& p : points) {
				radius += (p - center).norm();
			}
			radius /= n;

			// 2.2. Вычисляем невязки и веса (исключаем выбросы)
			Eigen::VectorXd residuals(n);
			double total_residual = 0.0;

			for (size_t i = 0; i < n; ++i) {
				double dist = (points[i] - center).norm();
				residuals[i] = std::abs(dist - radius);
				total_residual += residuals[i] * residuals[i];
			}

			double rms_residual = std::sqrt(total_residual / n);

			// 2.3. Проверка критерия остановки
			if (std::abs(prev_residual - rms_residual) < tolerance) {
				break;
			}
			prev_residual = rms_residual;

			// 2.4. Обновляем веса (используем стандартное отклонение вместо MAD)
			double mean = residuals.mean();
			double stddev = std::sqrt((residuals.array() - mean).square().sum() / n);
			double cutoff = 2.5 * stddev;  // 2.5 сигма

			for (size_t i = 0; i < n; ++i) 
			{
				if (residuals[i] < cutoff) 
				{
					double t = residuals[i] / cutoff;
					weights[i] = (1.0 - t * t) * (1.0 - t * t);
				}
				else 
				{
					weights[i] = 0.0;
				}
			}

			// 2.5. Решаем взвешенную систему
			Eigen::MatrixXd WA = weights.asDiagonal() * A;
			Eigen::VectorXd Wb = weights.asDiagonal() * b;
			center = WA.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Wb);
		}

		// 3. Вычисляем финальный радиус и невязку
		radius = 0.0;
		double final_residual = 0.0;
		for (const auto& p : points) 
		{
			double dist = (p - center).norm();
			radius += dist;
		}
		radius /= n;

		for (const auto& p : points)
		{
			double dist = (p - center).norm();
			final_residual += std::abs(dist - radius);
		}
		final_residual /= n;


		yz->parameters["xc"] = center[0];
		yz->parameters["yc"] = center[1];
		yz->parameters["zc"] = center[2];

		//cout << "center = " << center[0] << " " << center[1] <<
		//	" " << center[2] << endl;
		//cout << "final_residual = " << final_residual << endl;
	}

	//exit(-1);
	cout << "End Smooth_head_HP3" << endl;
}

void Setka::Smooth_head_TS3(void)
{
	int max_iterations = 8;
	double tolerance = 1e-6;
	double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
		norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
			this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

	cout << "Start Smooth_head_TS3" << endl;

	for (auto& yz : this->All_Yzel)
	{
		if (yz->type != Type_yzel::TS) continue;

		double phi = polar_angle(yz->coord[0][0],
			norm2(0.0, yz->coord[0][1], yz->coord[0][2]));

		// Включаем радиально-сферическое сглаживание только в головной зоне
		//if (phi > phi_a + const_pi / 180.0) continue;

		std::unordered_set<Yzel*> Yzels;
		Yzels.insert(yz);
		for (auto& gr : yz->grans)
		{
			if (gr->type2 != Type_Gran_surf::TS) continue;
			for (auto& yz2 : gr->yzels)
			{
				Yzels.insert(yz2);
			}
		}

		int n = Yzels.size();

		if (n < 3)
		{
			cout << "Error  " << yz->coord[0][0] << " " << yz->coord[0][1] << " "
				<< yz->coord[0][2] << endl;
			continue;
		}

		yz->Yzel_sosed_sglag2.clear();
		for (auto& i : Yzels)
		{
			if (i == yz) continue;
			yz->Yzel_sosed_sglag2.insert(i);
		}

		std::vector<Eigen::Vector3d> points;
		for (auto& i : Yzels)
		{
			Eigen::Vector3d AA;
			AA << i->coord[0][0],
				i->coord[0][1], i->coord[0][2];
			points.push_back(AA);
		}

		Eigen::MatrixXd A(n, 3);
		Eigen::VectorXd b(n);
		for (int i = 0; i < n; ++i) {
			A.row(i) = 2.0 * points[i];
			b(i) = points[i].squaredNorm();
		}

		// 2. Вычитаем уравнения для уменьшения влияния численных ошибок
		Eigen::Vector3d mean_point = Eigen::Vector3d::Zero();
		for (const auto& p : points) mean_point += p;
		mean_point /= n;

		Eigen::MatrixXd A_centered = A.rowwise() - A.colwise().mean();
		Eigen::VectorXd b_centered = b.array() - b.mean();

		Eigen::Vector3d center = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
		double radius = 0.0;

		// 2. Итеративное уточнение
		Eigen::VectorXd weights = Eigen::VectorXd::Ones(n);
		double prev_residual = std::numeric_limits<double>::max();

		for (int iter = 0; iter < max_iterations; ++iter) {
			// 2.1. Вычисляем текущий радиус как среднее расстояние до центра
			radius = 0.0;
			for (const auto& p : points) {
				radius += (p - center).norm();
			}
			radius /= n;

			// 2.2. Вычисляем невязки и веса (исключаем выбросы)
			Eigen::VectorXd residuals(n);
			double total_residual = 0.0;

			for (size_t i = 0; i < n; ++i) {
				double dist = (points[i] - center).norm();
				residuals[i] = std::abs(dist - radius);
				total_residual += residuals[i] * residuals[i];
			}

			double rms_residual = std::sqrt(total_residual / n);

			// 2.3. Проверка критерия остановки
			if (std::abs(prev_residual - rms_residual) < tolerance) {
				break;
			}
			prev_residual = rms_residual;

			// 2.4. Обновляем веса (используем стандартное отклонение вместо MAD)
			double mean = residuals.mean();
			double stddev = std::sqrt((residuals.array() - mean).square().sum() / n);
			double cutoff = 2.5 * stddev;  // 2.5 сигма

			for (size_t i = 0; i < n; ++i)
			{
				if (residuals[i] < cutoff)
				{
					double t = residuals[i] / cutoff;
					weights[i] = (1.0 - t * t) * (1.0 - t * t);
				}
				else
				{
					weights[i] = 0.0;
				}
			}

			// 2.5. Решаем взвешенную систему
			Eigen::MatrixXd WA = weights.asDiagonal() * A;
			Eigen::VectorXd Wb = weights.asDiagonal() * b;
			center = WA.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Wb);
		}

		// 3. Вычисляем финальный радиус и невязку
		radius = 0.0;
		double final_residual = 0.0;
		for (const auto& p : points)
		{
			double dist = (p - center).norm();
			radius += dist;
		}
		radius /= n;

		for (const auto& p : points)
		{
			double dist = (p - center).norm();
			final_residual += std::abs(dist - radius);
		}
		final_residual /= n;


		yz->parameters["xc"] = center[0];
		yz->parameters["yc"] = center[1];
		yz->parameters["zc"] = center[2];

		//cout << "center = " << center[0] << " " << center[1] <<
		//	" " << center[2] << endl;
		//cout << "final_residual = " << final_residual << endl;
	}

	//exit(-1);
	cout << "End Smooth_head_TS3" << endl;
}


void Setka::Find_Yzel_Sosed_for_BS(void)
{
	for (auto& yz : this->All_Yzel)
	{
		if (yz->type != Type_yzel::BS) continue;

		std::unordered_set<Yzel*> Yzels;
		Yzels.insert(yz);
		for (auto& gr : yz->grans)
		{
			if (gr->type2 != Type_Gran_surf::BS) continue;
			for (auto& yz2 : gr->yzels)
			{
				Yzels.insert(yz2);
			}
		}

		int n = Yzels.size();

		if (n < 3)
		{
			continue;
		}

		yz->Yzel_sosed_sglag2.clear();
		for (auto& i : Yzels)
		{
			if (i == yz) continue;
			yz->Yzel_sosed_sglag2.insert(i);
		}
	}
}


void Setka::Find_Yzel_Sosed_for_sglag(void)
{
	cout << "Start: Find_Yzel_Sosed_for_sglag" << endl;
	this->Renumerate();
	short int now = 0;
	Eigen::Vector3d A;
	Eigen::Vector3d B, V1, V2, VV;
	double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
		norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
			this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

	// Сначала находим для HP
	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[now][0], gr->center[now][1], gr->center[now][2];

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

		// Включаем радиально-сферическое сглаживание только в головной зоне
		if (phi > phi_a + const_pi / 180.0) continue;

		//  Для следующих узлов нужно вводить сглаживание
		for (auto& yz : gr->yzels)
		{
			if (yz->Yzel_sosed_sglag.size() > 0)
			{
				continue;
			}
			Yzel* AA3, * AA1;
			Yzel* AA2, * AA4;
			AA1 = nullptr;
			AA2 = nullptr;
			AA3 = nullptr;
			AA4 = nullptr;
			// AA33  AA3  (AA)  AA1  AA11 - вдоль HP  направление в апвинд
			// A22  AA2  (AA)  AA4  AA44  - вращение HP
			double d1 = -1.0;
			double d2 = 1.0;
			double d3 = -1.0;
			double d4 = 1.0;

			B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			dekard_skorost(B[0], B[1], B[2],
				0.0, 0.0, 1.0,
				V1[0], V1[1], V1[2]);
			dekard_skorost(B[0], B[1], B[2],
				0.0, 1.0, 0.0,
				V2[0], V2[1], V2[2]);
			V1.normalize();
			V2.normalize();
			for (auto& grr : gr->grans_surf)
			{
				int num = yz->number;
				if (std::find_if(
					grr->yzels.begin(),
					grr->yzels.end(),
					[num](Yzel* yzz) { return yzz->number == num; })
					== grr->yzels.end())
				{
					continue;
				}

				for (auto& yzz : grr->yzels)
				{
					if (yzz->number == yz->number) continue;
					VV << yzz->coord[0][0] - B[0], yzz->coord[0][1] - B[1],
						yzz->coord[0][2] - B[2];
					VV.normalize();
					double sk1 = VV.dot(V1);
					if (sk1 > d1)
					{
						d1 = sk1;
						AA1 = yzz;
					}
					if (sk1 < d2)
					{
						d2 = sk1;
						AA2 = yzz;
					}

					double sk2 = VV.dot(V2);
					if (sk2 > d3)
					{
						d3 = sk2;
						AA3 = yzz;
					}
					if (sk2 < d4)
					{
						d4 = sk2;
						AA4 = yzz;
					}
				}
			}

			if (AA1 != nullptr && AA2 != nullptr && AA3 != nullptr && AA4 != nullptr)
			{
				yz->Yzel_sosed_sglag["AA1"] = AA1;
				yz->Yzel_sosed_sglag["AA2"] = AA2;
				yz->Yzel_sosed_sglag["AA3"] = AA3;
				yz->Yzel_sosed_sglag["AA4"] = AA4;
				this->Yzels_HP_sglag.push_back(yz);
			}
		}
	}
	
	for (auto& yz : this->Yzels_HP_sglag)
	{
		auto s1 = yz->Yzel_sosed_sglag["AA1"];
		auto s2 = yz->Yzel_sosed_sglag["AA2"];
		auto s3 = yz->Yzel_sosed_sglag["AA3"];
		auto s4 = yz->Yzel_sosed_sglag["AA4"];


		if (s1->Yzel_sosed_sglag.find("AA1") != s1->Yzel_sosed_sglag.end() &&
			s2->Yzel_sosed_sglag.find("AA2") != s2->Yzel_sosed_sglag.end() &&
			s3->Yzel_sosed_sglag.find("AA3") != s3->Yzel_sosed_sglag.end() &&
			s4->Yzel_sosed_sglag.find("AA4") != s4->Yzel_sosed_sglag.end())
		{
			yz->Yzel_sosed_sglag["AA11"] = s1->Yzel_sosed_sglag["AA1"];
			yz->Yzel_sosed_sglag["AA22"] = s2->Yzel_sosed_sglag["AA2"];
			yz->Yzel_sosed_sglag["AA33"] = s3->Yzel_sosed_sglag["AA3"];
			yz->Yzel_sosed_sglag["AA44"] = s4->Yzel_sosed_sglag["AA4"];
		}
		else
		{
			yz->Yzel_sosed_sglag.clear();
		}

	}
	
	// Удаляем узлы, у которых не нашлось нужных соседей
	this->Yzels_HP_sglag.erase(
		std::remove_if(
			this->Yzels_HP_sglag.begin(),
			this->Yzels_HP_sglag.end(),
			[](Yzel* yz) { return yz->Yzel_sosed_sglag.size() == 0; }
		),
		this->Yzels_HP_sglag.end()
	);


	cout << "End: Find_Yzel_Sosed_for_sglag" << endl;
}

Eigen::Vector2d fitCircleRobust(const std::vector<Eigen::Vector2d>& points, double& radius_, 
	double& error_, int max_iter = 50, double tol = 1e-8)
{
	const int n = points.size();
	if (n < 3) throw std::runtime_error("Need at least 3 points");

	// 1. Initial guess using geometric median
	Eigen::Vector2d center = Eigen::Vector2d::Zero();
	for (const auto& p : points) center += p;
	center /= n;

	// 2. Levenberg-Marquardt optimization
	double radius = 0.0;
	for (const auto& p : points) radius += (p - center).norm();
	radius /= n;

	double lambda = 0.001;
	double prev_error = std::numeric_limits<double>::max();
	Eigen::Vector3d params;
	params << center.x(), center.y(), radius;

	for (int iter = 0; iter < max_iter; ++iter) {
		Eigen::MatrixXd J(n, 3);
		Eigen::VectorXd residuals(n);
		double error = 0.0;

		// Compute Jacobian and residuals
		for (int i = 0; i < n; ++i) {
			double dx = points[i].x() - params(0);
			double dy = points[i].y() - params(1);
			double dist = std::sqrt(dx * dx + dy * dy);

			// Handle near-zero distance case
			if (dist < 1e-10) {
				residuals(i) = -params(2);
				J.row(i) << 0, 0, -1;
			}
			else {
				residuals(i) = dist - params(2);
				J(i, 0) = -dx / dist;
				J(i, 1) = -dy / dist;
				J(i, 2) = -1.0;
			}
			error += residuals(i) * residuals(i);
		}

		error = std::sqrt(error / n);
		if (std::abs(prev_error - error) < tol) break;
		prev_error = error;

		// LM update
		Eigen::Matrix3d H = J.transpose() * J;
		H.diagonal() *= (1.0 + lambda);
		Eigen::Vector3d delta = H.ldlt().solve(J.transpose() * (-residuals));

		// Safeguard against invalid updates
		if (!delta.allFinite()) {
			lambda *= 10.0;
			continue;
		}

		params += delta;
		lambda = (error < prev_error) ? lambda / 2.0 : lambda * 10.0;
	}

	// Final parameters
	center << params(0), params(1);
	radius = params(2);

	// Calculate final error
	double final_error = 0.0;
	for (const auto& p : points) {
		double dist = (p - center).norm();
		final_error += std::pow(dist - radius, 2);
	}
	final_error = std::sqrt(final_error / n);

	error_ = final_error;
	radius_ = radius;

	return center;
}

Eigen::Vector2d fitCircleCenter(const std::vector<Eigen::Vector2d>& points) 
{
	const int n = points.size();
	if (n < 3) throw std::runtime_error("Need at least 3 points to fit a circle");

	// 1. Центрируем точки для улучшения численной устойчивости
	Eigen::Vector2d centroid = Eigen::Vector2d::Zero();
	for (const auto& p : points) centroid += p;
	centroid /= n;

	// 2. Строим систему уравнений для метода наименьших квадратов
	Eigen::MatrixXd A(n, 2);
	Eigen::VectorXd b(n);

	for (int i = 0; i < n; ++i) {
		Eigen::Vector2d p = points[i] - centroid;
		A.row(i) = 2.0 * p;
		b(i) = p.squaredNorm();
	}

	// 3. Решаем систему
	Eigen::Vector2d center_offset = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	// 4. Возвращаем центр в исходной системе координат
	return center_offset + centroid;
}

// Функция для вычисления радиуса (после нахождения центра)
double computeRadius(const std::vector<Eigen::Vector2d>& points, const Eigen::Vector2d& center) 
{
	double radius = 0.0;
	for (const auto& p : points) {
		radius += (p - center).norm();
	}
	return radius / points.size();
}

void Setka::Smooth_angle_HP(void)
{
	cout << "Start Smooth_angle_HP " << endl;
	for (auto& yz : this->All_Yzel)
	{
		yz->num_velocity = 0;
		yz->velocity[0] = 0.0;
		yz->velocity[1] = 0.0;
		yz->velocity[2] = 0.0;
	}

	vector<vector <vector<Luch*>>> VVV;
	vector<int> VVV_int;

	VVV.push_back(this->D_Luch);
	VVV.push_back(this->E_Luch);
	VVV.push_back(this->B_Luch);
	VVV_int.push_back(1);
	VVV_int.push_back(1);
	VVV_int.push_back(2);

	for (size_t ii = 0; ii < VVV.size(); ii++)
	{
		int nn = VVV[ii].size();
		int mm = VVV[ii][0].size();

		#pragma omp parallel for schedule(dynamic)
		for (size_t i = 0; i < nn; i++)
		{
			//if (i == 1) continue;
			if (i > 3 && i < nn - 3) continue;

			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				std::vector<Eigen::Vector2d> points;
				auto AA = VVV[ii][i][j]->Yzels_opor[VVV_int[ii]];
				points.push_back(Eigen::Vector2d(AA->coord[0][1], AA->coord[0][2]));
				Yzel* AA2 = AA;

				AA = VVV[ii][im][j]->Yzels_opor[VVV_int[ii]];
				points.push_back(Eigen::Vector2d(AA->coord[0][1], AA->coord[0][2]));
				AA2->Yzel_sosed_sglag2.insert(AA);

				AA = VVV[ii][imm][j]->Yzels_opor[VVV_int[ii]];
				points.push_back(Eigen::Vector2d(AA->coord[0][1], AA->coord[0][2]));
				AA2->Yzel_sosed_sglag2.insert(AA);

				AA = VVV[ii][ip][j]->Yzels_opor[VVV_int[ii]];
				points.push_back(Eigen::Vector2d(AA->coord[0][1], AA->coord[0][2]));
				AA2->Yzel_sosed_sglag2.insert(AA);

				AA = VVV[ii][ipp][j]->Yzels_opor[VVV_int[ii]];
				points.push_back(Eigen::Vector2d(AA->coord[0][1], AA->coord[0][2]));
				AA2->Yzel_sosed_sglag2.insert(AA);

				//Eigen::Vector2d center = fitCircleCenter(points);
				//double radius = computeRadius(points, center);

				double radius, error;
				Eigen::Vector2d center = fitCircleRobust(points,
					radius, error, 8);

				AA2->parameters["xcc"] = center[0];
				AA2->parameters["ycc"] = center[1];
				

				AA = VVV[ii][i][j]->Yzels_opor[VVV_int[ii]];

				Eigen::Vector2d A, B;
				A << AA->coord[0][1], AA->coord[0][2];
				B = A;

				continue;

				// Дальше алгоритм не используется, мы просто нашли центр для 
				// сглаживания HP в программе

				unsigned int kl = 0;
				while ((B - center).norm() > radius)
				{
					kl++;
					B *= 0.999;
					if (kl > 100000)
					{
						cout << "Error 12ertewrfwr3 " << endl;
						cout << B[0] << " " << B[1] << endl;
						cout << center[0] << " " << center[1] << endl;
						cout << radius << endl;
						cout << "Points: " << endl;
						for (auto& ik : points)
						{
							cout << ik[0] << " " << ik[1] << endl;
						}
						exit(-1);
					}
				}

				kl = 0;
				while ((B - center).norm() < radius)
				{
					B *= 1.001;
					if (kl > 100000)
					{
						cout << "Error gegy45432ty5hre" << endl;
						cout << B[0] << " " << B[1] << endl;
						cout << center[0] << " " << center[1] << endl;
						cout << radius << endl;
						exit(-1);
					}
				}

				AA->velocity[0] = 0.0;
				AA->velocity[1] = (B - A)[0];
				AA->velocity[2] = (B - A)[1];
			}
		}
	}

	// Вычисляем новые координаты узлов на HP
	for (auto& yz : this->All_Yzel)
	{

		// Вычисляем новую координату узла
		if (yz->type == Type_yzel::HP && this->phys_param->move_HP == true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			if (A(0) >= 0.0)
			{
				B = A;
			}
			else
			{
				B << 0.0, A(1), A(2);
			}
			B.normalize();
			A = A + V.dot(B) * B;
			yz->coord[0][0] = A[0];
			yz->coord[0][1] = A[1];
			yz->coord[0][2] = A[2];
		}
		
	}

	// Остальные узлы на HP (невыделяемой части) надо подвинуть
	if (false)
	{
		short int NN = this->D_Luch[0].size() - 1;
		for (auto& L : this->D_Luch)
		{

			double h1 = norm2(0.0, L[this->geo->N4 - 4]->Yzels_opor[1]->coord[0][1],
				L[this->geo->N4 - 4]->Yzels_opor[1]->coord[0][2]);


			for (short int i = this->geo->N4 - 3; i <= NN; i++)
			{
				//double h = h1 + (i - this->geo->N4 + 2) * (h2 - h1) / (NN - this->geo->N4 + 2);
				double h = h1;
				auto yz = L[i]->Yzels_opor[1];
				double hh = norm2(0.0, yz->coord[0][1], yz->coord[0][2]);


				yz->coord[0][1] = yz->coord[0][1] * h / hh;
				yz->coord[0][2] = yz->coord[0][2] * h / hh;
			}
		}
	}


	// Перестраиваем сетку
	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	
	cout << "End Smooth_angle_HP " << endl;
}