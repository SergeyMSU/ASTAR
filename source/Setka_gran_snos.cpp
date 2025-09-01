#include "Setka.h"

void Setka::Snos_on_Gran(Gran* gr, unordered_map<string, double>& par_left,
	unordered_map<string, double>& par_right, short int now)
{
	// Для граничных граней не надо делать ТВД
	if (gr->type != Type_Gran::Us)
	{
		//cout << "D1 " << endl;
		auto C = gr->cells[0];

		int8_t type_gran;
		if (gr->type == Type_Gran::Inner_Hard) type_gran = 0;
		if (gr->type == Type_Gran::Outer_Hard) type_gran = 1;
		if (gr->type == Type_Gran::Outer_Soft) type_gran = 2;

		int8_t plasma_cond = this->phys_param->plasma_condition(type_gran, 0);

		if (plasma_cond == 1)
		{
			for (auto& nam : this->phys_param->plasma_name)
			{
				par_left[nam] = C->parameters[now][nam];
			}

			for (auto& nam : this->phys_param->plasma_name)
			{
				par_right[nam] = par_left[nam];
			}
		}
		else if (plasma_cond == 2)
		{
			for (auto& nam : this->phys_param->plasma_name)
			{
				par_left[nam] = C->parameters[now][nam];
				par_right[nam] = gr->parameters[nam];
			}
		}
		else
		{
			cout << "Error 8965655332" << endl;
			exit(-1);
		}

	}
	else
	{
		//cout << "C1 " << endl;
		if (gr->cells.size() < 2)
		{
			cout << "Error i7ktyrh4567456h45" << endl;
			exit(-1);
		}

		auto A = gr->cells[0];
		auto B = gr->cells[1];

		// Делаем ли ТВД?
		if (this->phys_param->TVD == true)
		{
			// AA  -  A  -|-  B  -  BB
			Eigen::Vector3d Ac, Bc, AAc, BBc, G, vec;
			if (gr->cells_TVD.size() < 2)
			{
				cout << "Error 4t43t4gwt3w" << endl;
				exit(-1);
			}
			auto AA = gr->cells_TVD[0];
			auto BB = gr->cells_TVD[1];
			if (AA == nullptr || BB == nullptr) goto a1; // Для таких не надо делать TVD
		
			G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
			Ac << A->center[now][0], A->center[now][1], A->center[now][2];
			Bc << B->center[now][0], B->center[now][1], B->center[now][2];
			AAc << AA->center[now][0], AA->center[now][1], AA->center[now][2];
			BBc << BB->center[now][0], BB->center[now][1], BB->center[now][2];
			double d1, d2, dd1, dd2;
			vec = Ac - G;
			d1 = vec.norm();
			vec = Bc - G;
			d2 = vec.norm();
			vec = Ac - AAc;
			dd1 = vec.norm();
			vec = Bc - BBc;
			dd2 = vec.norm();

			// AA  -  A  -|-  B  -  BB
			// r3,    r1, rr, r2,   r4
			double r3, r1, rr, r2, r4;

			if (true)
			{
				Eigen::Vector3d VAA, VA, VB, VBB, Vleft, Vright;
				r3 = AAc.norm();
				r1 = Ac.norm();
				rr = G.norm();
				r2 = Bc.norm();
				r4 = BBc.norm();

				// Преобразуем нужные векторные величины
				if (true)
				{
					std::array<string, 2> V1;
					std::array<string, 2> V2;
					std::array<string, 2> V3;

					V1[0] = "Vx";
					V2[0] = "Vy";
					V3[0] = "Vz";

				    V1[1] = "Bx";
					V2[1] = "By";
					V3[1] = "Bz";

					double the1 = acos(G[2] / rr);

					if(the1 > const_pi/9 && the1 < 8 * const_pi / 9)
					{ 
						for (short int ik = 0; ik < 2; ik++)
						{

							// Переводим скорости в сферическую с.к.
							spherical_skorost(AAc[2], AAc[0], AAc[1],
								AA->parameters[now][V3[ik]], AA->parameters[now][V1[ik]],
								AA->parameters[now][V2[ik]], VAA[0], VAA[1], VAA[2]);

							spherical_skorost(Ac[2], Ac[0], Ac[1],
								A->parameters[now][V3[ik]], A->parameters[now][V1[ik]],
								A->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

							spherical_skorost(Bc[2], Bc[0], Bc[1],
								B->parameters[now][V3[ik]], B->parameters[now][V1[ik]],
								B->parameters[now][V2[ik]], VB[0], VB[1], VB[2]);

							spherical_skorost(BBc[2], BBc[0], BBc[1],
								BB->parameters[now][V3[ik]], BB->parameters[now][V1[ik]],
								BB->parameters[now][V2[ik]], VBB[0], VBB[1], VBB[2]);


							// Интерполируем скорости (в сферической с.к.)
							for (short int i = 0; i < 3; i++)
							{
								if (V1[ik] == "Bx")
								{
									if (i == 0)
									{
										Vleft[i] = linear(-dd1 - d1, VAA[i] * kv(r3), -d1, VA[i] * kv(r1),
											d2, VB[i] * kv(r2), 0.0) / kv(rr);
										Vright[i] = linear(-d1, VA[i] * kv(r1), d2, VB[i] * kv(r2),
											d2 + dd2, VBB[i] * kv(r4), 0.0) / kv(rr);
									}
									else
									{
										Vleft[i] = linear(-dd1 - d1, VAA[i] * (r3), -d1, VA[i] * (r1),
											d2, VB[i] * (r2), 0.0) / (rr);
										Vright[i] = linear(-d1, VA[i] * (r1), d2, VB[i] * (r2),
											d2 + dd2, VBB[i] * (r4), 0.0) / (rr);
									}
								}
								else
								{
									Vleft[i] = linear(-dd1 - d1, VAA[i], -d1, VA[i], d2, VB[i], 0.0);
									Vright[i] = linear(-d1, VA[i], d2, VB[i], d2 + dd2, VBB[i], 0.0);
								}
							}

							dekard_skorost(G[2], G[0], G[1],
								Vleft[0], Vleft[1], Vleft[2],
								par_left[V3[ik]], par_left[V1[ik]], par_left[V2[ik]]);

							dekard_skorost(G[2], G[0], G[1],
								Vright[0], Vright[1], Vright[2],
								par_right[V3[ik]], par_right[V1[ik]], par_right[V2[ik]]);
						}
					}
					else
					{
						for (short int ik = 0; ik < 2; ik++)
						{

							// Переводим скорости в сферическую с.к.
							spherical_skorost(AAc[0], AAc[1], AAc[2],
								AA->parameters[now][V1[ik]], AA->parameters[now][V2[ik]],
								AA->parameters[now][V3[ik]], VAA[0], VAA[1], VAA[2]);

							spherical_skorost(Ac[0], Ac[1], Ac[2],
								A->parameters[now][V1[ik]], A->parameters[now][V2[ik]],
								A->parameters[now][V3[ik]], VA[0], VA[1], VA[2]);

							spherical_skorost(Bc[0], Bc[1], Bc[2],
								B->parameters[now][V1[ik]], B->parameters[now][V2[ik]],
								B->parameters[now][V3[ik]], VB[0], VB[1], VB[2]);

							spherical_skorost(BBc[0], BBc[1], BBc[2],
								BB->parameters[now][V1[ik]], BB->parameters[now][V2[ik]],
								BB->parameters[now][V3[ik]], VBB[0], VBB[1], VBB[2]);


							// Интерполируем скорости (в сферической с.к.)
							for (short int i = 0; i < 3; i++)
							{
								if (V1[ik] == "Bx")
								{
									if (i == 0)
									{
										Vleft[i] = linear(-dd1 - d1, VAA[i] * kv(r3), -d1, VA[i] * kv(r1),
											d2, VB[i] * kv(r2), 0.0) / kv(rr);
										Vright[i] = linear(-d1, VA[i] * kv(r1), d2, VB[i] * kv(r2),
											d2 + dd2, VBB[i] * kv(r4), 0.0) / kv(rr);
									}
									else
									{
										Vleft[i] = linear(-dd1 - d1, VAA[i] * (r3), -d1, VA[i] * (r1),
											d2, VB[i] * (r2), 0.0) / (rr);
										Vright[i] = linear(-d1, VA[i] * (r1), d2, VB[i] * (r2),
											d2 + dd2, VBB[i] * (r4), 0.0) / (rr);
									}
								}
								else
								{
									Vleft[i] = linear(-dd1 - d1, VAA[i], -d1, VA[i], d2, VB[i], 0.0);
									Vright[i] = linear(-d1, VA[i], d2, VB[i], d2 + dd2, VBB[i], 0.0);
								}
							}

							dekard_skorost(G[0], G[1], G[2],
								Vleft[0], Vleft[1], Vleft[2],
								par_left[V1[ik]], par_left[V2[ik]], par_left[V3[ik]]);

							dekard_skorost(G[0], G[1], G[2],
								Vright[0], Vright[1], Vright[2],
								par_right[V1[ik]], par_right[V2[ik]], par_right[V3[ik]]);
						}
					}
				}


			
				for (auto& nam : this->phys_param->param_names)
				{
					if (this->phys_param->r2_snos_names.find(nam) != this->phys_param->r2_snos_names.end())
					{
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam] * kv(r3),
							-d1, A->parameters[now][nam] * kv(r1),
							d2, B->parameters[now][nam] * kv(r2), 0.0)/kv(rr);
						par_right[nam] = linear(-d1, A->parameters[now][nam] * kv(r1),
							d2, B->parameters[now][nam] * kv(r2),
							d2 + dd2, BB->parameters[now][nam] * kv(r4), 0.0) / kv(rr);
					}
					else if (this->phys_param->r2g_snos_names.find(nam) != this->phys_param->r2g_snos_names.end())
					{
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam] * kvg(r3),
							-d1, A->parameters[now][nam] * kvg(r1),
							d2, B->parameters[now][nam] * kvg(r2), 0.0) / kvg(rr);
						par_right[nam] = linear(-d1, A->parameters[now][nam] * kvg(r1),
							d2, B->parameters[now][nam] * kvg(r2),
							d2 + dd2, BB->parameters[now][nam] * kvg(r4), 0.0) / kvg(rr);
					}
					else if (nam == "Vx" || nam == "Vy" || nam == "Vz" 
						|| nam == "Bx" || nam == "By" || nam == "Bz")
					{
						// PASS
					}
					else
					{
						cout << "SS = " << nam << endl;
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam],
							-d1, A->parameters[now][nam],
							d2, B->parameters[now][nam], 0.0);
						par_right[nam] = linear(-d1, A->parameters[now][nam],
							d2, B->parameters[now][nam],
							d2 + dd2, BB->parameters[now][nam], 0.0);
					}
				}

			}
			else
			{
				vector<string>* parameters;
				if (this->phys_param->culc_plasma == true)
				{
					parameters = &this->phys_param->param_names;
				}
				else
				{
					parameters = &this->phys_param->H_param_names;
				}

				for (const auto& nam : *parameters)
				{
					par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam],
						-d1, A->parameters[now][nam],
						d2, B->parameters[now][nam], 0.0);
					par_right[nam] = linear(-d1, A->parameters[now][nam],
						d2, B->parameters[now][nam],
						d2 + dd2, BB->parameters[now][nam], 0.0);
				}
			}

			if (this->phys_param->culc_plasma == true)
			{
				if (par_left["rho"] <1e-12) par_left["rho"] = A->parameters[now]["rho"];
				if (par_left["p"] < 1e-12) par_left["p"] = A->parameters[now]["p"];
				if (par_right["rho"] < 1e-12) par_right["rho"] = B->parameters[now]["rho"];
				if (par_right["p"] < 1e-12) par_right["p"] = B->parameters[now]["p"];
			}
			
		}
		else
		{
			a1:

			Eigen::Vector3d VAA, VA, VB, VBB, Vleft, Vright;

			if (false)//(A->type == Type_cell::Zone_1)
			{
				Eigen::Vector3d Ac, G;
				Ac << A->center[now][0], A->center[now][1], A->center[now][2];
				G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
				double r1 = Ac.norm();
				double rr = G.norm();

				// Преобразуем нужные векторные величины
				if (true)
				{
					std::array<string, 2> V1;
					std::array<string, 2> V2;
					std::array<string, 2> V3;

					V1[0] = "Vx";
					V2[0] = "Vy";
					V3[0] = "Vz";

					V1[1] = "Bx";
					V2[1] = "By";
					V3[1] = "Bz";

					for (short int ik = 0; ik < 2; ik++)
					{
						// Переводим скорости в сферическую с.к.
						VA << 0.0, 0.0, 0.0;

						spherical_skorost(Ac[2], Ac[0], Ac[1],
							A->parameters[now][V3[ik]], A->parameters[now][V1[ik]],
							A->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

						if(V1[ik] == "Bx")
						{
							//cout << "B = " << VA[0] << " " << VA[1] << " " << VA[2] << endl;
							VA[0] = VA[0] * (kv(r1 / rr));
							VA[1] = VA[1] * ((r1) / (rr));
							VA[2] = VA[2] * ((r1) / (rr));
							//cout << VA[0] << " " << VA[1] << " " << VA[2] << endl << endl;
						}

						dekard_skorost(G[2], G[0], G[1],
							VA[0], VA[1], VA[2],
							par_left[V3[ik]], par_left[V1[ik]], par_left[V2[ik]]);
					}
				}

				if (this->phys_param->culc_plasma == true)
				{
					par_left["rho"] = A->parameters[now]["rho"] * kv(r1) / kv(rr);
					par_left["p"] = A->parameters[now]["p"] * kvg(r1) / kvg(rr);
				

					for (auto& nam : this->phys_param->pui_name)
					{
						par_left["rho" + nam] = A->parameters[now]["rho" + nam] * kv(r1) / kv(rr);
						par_left["p" + nam] = A->parameters[now]["p" + nam] * kvg(r1) / kvg(rr);
					}
				}

				for (auto& nam : this->phys_param->H_name)
				{
					par_left["rho" + nam] = A->parameters[now]["rho" + nam];
					par_left["p" + nam] = A->parameters[now]["p" + nam];
					if (nam == "_H1") continue;
					par_left["Vx" + nam] = A->parameters[now]["Vx" + nam];
					par_left["Vy" + nam] = A->parameters[now]["Vy" + nam];
					par_left["Vz" + nam] = A->parameters[now]["Vz" + nam];
				}
			}
			else
			{
				if (this->phys_param->culc_plasma == true)
				{
					par_left["rho"] = A->parameters[now]["rho"];
					par_left["p"] = A->parameters[now]["p"];
					par_left["Vx"] = A->parameters[now]["Vx"];
					par_left["Vy"] = A->parameters[now]["Vy"];
					par_left["Vz"] = A->parameters[now]["Vz"];
					par_left["Bx"] = A->parameters[now]["Bx"];
					par_left["By"] = A->parameters[now]["By"];
					par_left["Bz"] = A->parameters[now]["Bz"];


				}

			}

			if (false)//(B->type == Type_cell::Zone_1)
			{
				Eigen::Vector3d Ac, G;
				Ac << B->center[now][0], B->center[now][1], B->center[now][2];
				G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
				double r1 = Ac.norm();
				double rr = G.norm();

				// Преобразуем нужные векторные величины
				if (true)
				{
					std::array<string, 2> V1;
					std::array<string, 2> V2;
					std::array<string, 2> V3;

					V1[0] = "Vx";
					V2[0] = "Vy";
					V3[0] = "Vz";

					V1[1] = "Bx";
					V2[1] = "By";
					V3[1] = "Bz";


					for (short int ik = 0; ik < 2; ik++)
					{
						// Переводим скорости в сферическую с.к.
						spherical_skorost(Ac[2], Ac[0], Ac[1],
							B->parameters[now][V3[ik]], B->parameters[now][V1[ik]],
							B->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

						if (V1[ik] == "Bx")
						{
							//cout << "A = " << VA[0] << " " << VA[1] << " " << VA[2] << endl;
							VA[0] = VA[0] * (kv(r1 / rr));
							VA[1] = VA[1] * ((r1) / (rr));
							VA[2] = VA[2] * ((r1) / (rr));
							//cout << VA[0] << " " << VA[1] << " " << VA[2] << endl << endl;;
						}

						dekard_skorost(G[2], G[0], G[1],
							VA[0], VA[1], VA[2],
							par_right[V3[ik]], par_right[V1[ik]], par_right[V2[ik]]);
					}
				}

				if (this->phys_param->culc_plasma == true)
				{
					par_right["rho"] = B->parameters[now]["rho"] * kv(r1) / kv(rr);
					par_right["Q"] = B->parameters[now]["Q"] * kv(r1) / kv(rr);
					par_right["rho_He"] = B->parameters[now]["rho_He"] * kv(r1) / kv(rr);
					par_right["p"] = B->parameters[now]["p"] * kvg(r1) / kvg(rr);

					for (auto& nam : this->phys_param->pui_name)
					{
						par_right["rho" + nam] = B->parameters[now]["rho" + nam] * kv(r1) / kv(rr);
						par_right["p" + nam] = B->parameters[now]["p" + nam] * kvg(r1) / kvg(rr);
					}
				}
					
				if (this->phys_param->culc_atoms == true)
				{
					for (auto& nam : this->phys_param->H_name)
					{
						par_right["rho" + nam] = B->parameters[now]["rho" + nam];
						par_right["p" + nam] = B->parameters[now]["p" + nam];
						if (nam == "_H1") continue;
						par_right["Vx" + nam] = B->parameters[now]["Vx" + nam];
						par_right["Vy" + nam] = B->parameters[now]["Vy" + nam];
						par_right["Vz" + nam] = B->parameters[now]["Vz" + nam];
					}
				}
			}
			else
			{
				if (this->phys_param->culc_plasma == true)
				{
					par_right["rho"] = B->parameters[now]["rho"];
					par_right["p"] = B->parameters[now]["p"];
					par_right["Vx"] = B->parameters[now]["Vx"];
					par_right["Vy"] = B->parameters[now]["Vy"];
					par_right["Vz"] = B->parameters[now]["Vz"];
					par_right["Bx"] = B->parameters[now]["Bx"];
					par_right["By"] = B->parameters[now]["By"];
					par_right["Bz"] = B->parameters[now]["Bz"];

				}

			}
			
		}

		//cout << "C2 " << endl;
	}
}