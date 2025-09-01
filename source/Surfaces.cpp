#include "Surfaces.h"

Surfaces::~Surfaces()
{
	this->HP_cilindr.resize(boost::extents[0][0]);
	this->x_cilindr.resize(boost::extents[0][0]);
	this->TS_radial2.resize(boost::extents[0][0]);
	this->the_angle2.resize(boost::extents[0][0]);
	this->BS_radial.resize(boost::extents[0][0]);
	this->HP_radial.resize(boost::extents[0][0]);
	this->TS_radial.resize(boost::extents[0][0]);
	this->the_angle.resize(boost::extents[0][0]);
	this->phi_angle.clear();
	this->parameters.clear();
}

void Surfaces::Read_old(string nam)
{
	this->name = nam;

	ifstream fin(this->name, ios::binary | ios::in);
	if (!fin)
	{
		cout << "ERROR  8767654541  Net takogo fajla" << endl;
		exit(-1);
	}

	int i = 0;
	double a = 0.0;
	fin.read((char*)&i, sizeof i);
	parameters["phi"] = i;
	fin.read((char*)&i, sizeof i);
	parameters["mA"] = i;
	fin.read((char*)&i, sizeof i);
	parameters["mBC"] = i; 
	fin.read((char*)&i, sizeof i);
	parameters["mO"] = i;
	fin.read((char*)&i, sizeof i);
	parameters["mK"] = i;

	fin.read((char*)&a, sizeof a);
	fin.read((char*)&a, sizeof a);

	phi_angle.resize(parameters["phi"]);

	TS_radial.resize(boost::extents[parameters["phi"]][parameters["mA"]]);
	HP_radial.resize(boost::extents[parameters["phi"]][parameters["mA"]]);
	BS_radial.resize(boost::extents[parameters["phi"]][parameters["mA"]]);
	the_angle.resize(boost::extents[parameters["phi"]][parameters["mA"]]);

	the_angle2.resize(boost::extents[parameters["phi"]][parameters["mBC"] + parameters["mK"]]);
	TS_radial2.resize(boost::extents[parameters["phi"]][parameters["mBC"] + parameters["mK"]]);
	
	
	x_cilindr.resize(boost::extents[parameters["phi"]][parameters["mBC"] + parameters["mO"]]);
	HP_cilindr.resize(boost::extents[parameters["phi"]][parameters["mBC"] + parameters["mO"]]);


	for(int i = 0; i < parameters["phi"]; i++)
	{
		fin.read((char*)&a, sizeof a);
		phi_angle[i] = a;

		for (int j = 0; j < parameters["mA"]; j++)
		{
			fin.read((char*)&a, sizeof a);
			the_angle[i][j] = a;
			fin.read((char*)&a, sizeof a);
			TS_radial[i][j] = a;
			fin.read((char*)&a, sizeof a);
			HP_radial[i][j] = a;
			fin.read((char*)&a, sizeof a);
			BS_radial[i][j] = a;
		}

		for (int j = 0; j < parameters["mBC"]; j++)
		{
			fin.read((char*)&a, sizeof a);
			the_angle2[i][j] = a;
			fin.read((char*)&a, sizeof a);
			TS_radial2[i][j] = a;
			fin.read((char*)&a, sizeof a);
			x_cilindr[i][j] = a;
			fin.read((char*)&a, sizeof a);
			HP_cilindr[i][j] = a;
		}

		for (int j = 0; j < parameters["mO"]; j++)
		{
			fin.read((char*)&a, sizeof a);
			x_cilindr[i][parameters["mBC"] + j] = a;
			fin.read((char*)&a, sizeof a);
			HP_cilindr[i][parameters["mBC"] + j] = a;
		}

		for (int j = 0; j < parameters["mK"]; j++)
		{
			fin.read((char*)&a, sizeof a);
			the_angle2[i][parameters["mBC"] + j] = a;
			fin.read((char*)&a, sizeof a);
			TS_radial2[i][parameters["mBC"] + j] = a;
		}

	}


	fin.close();
}

double Surfaces::Get_TS(const double& phi, const double& the)
{
	int i1 = 0, i2 = 0; // по углу phi
	int j1 = 0, j2 = 0; // по углу the

	i1 = this->phi_angle.size() - 1;

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
 		if (this->phi_angle[i] > phi)
		{
			i2 = i;
			i1 = i2 - 1;
			if (i1 < 0) i1 = this->phi_angle.size() - 1;
			break;
		}
	}

	

	if (the < const_pi / 2.0)
	{
		// ¬ головной области
		for (int i = 0; i < this->the_angle.shape()[1]; i++)
		{
			if (this->the_angle[i1][i] > the)
			{
				j2 = i;
				j1 = j2 - 1;
				if (j1 < 0) 
				{
					j1 = 0;
					j2 = 1;
				}
				break;
			}
		}

		double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
		double y = (the - this->the_angle[i1][j1]) / (this->the_angle[i1][j2] - this->the_angle[i1][j1]);

		if (i1 == this->phi_angle.size() - 1)
		{
			x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
		}

		double b =   this->TS_radial[i1][j1] * (1 - x) * (1 - y) +
			this->TS_radial[i2][j1] * x * (1 - y) +
			this->TS_radial[i1][j2] * (1 - x) * y +
			this->TS_radial[i2][j2] * x * y;

		return b;

	}
	else if (the < this->the_angle2[0][0])
	{
		double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
		double y = (the - const_pi / 2.0) / (this->the_angle2[0][0] - const_pi / 2.0);

		if (i1 == this->phi_angle.size() - 1)
		{
			x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
		}


		return this->TS_radial[i1][this->the_angle.shape()[1] - 1] * (1 - x) * (1 - y) +
			this->TS_radial[i2][this->the_angle.shape()[1] - 1] * x * (1 - y) +
			this->TS_radial2[i1][0] * (1 - x) * y +
			this->TS_radial2[i2][0] * x * y;
	}
	else
	{

		//cout << this->phi_angle[i1] << " " << phi << " " << this->phi_angle[i2] << endl;


		// ¬ хвостовой области
		j2 = this->the_angle2.shape()[1] - 1;
		j1 = j2 - 1;

		for (int i = 0; i < this->the_angle2.shape()[1]; i++)
		{
			if (this->the_angle2[i1][i] > the)
			{
				j2 = i;
				j1 = j2 - 1;
				if (j1 < 0)
				{
					j1 = 0;
					j2 = 1;
				}
				break;
			}
		}

		double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
		double y = (the - this->the_angle2[i1][j1]) / (this->the_angle2[i1][j2] - this->the_angle2[i1][j1]);

		if (i1 == this->phi_angle.size() - 1)
		{
			x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
		}


		double b =  this->TS_radial2[i1][j1] * (1 - x) * (1 - y) +
			this->TS_radial2[i2][j1] * x * (1 - y) +
			this->TS_radial2[i1][j2] * (1 - x) * y +
			this->TS_radial2[i2][j2] * x * y;

		return b;
	}

	

	return 0.0;
}

double Surfaces::Get_HP(const double& phi, const double& the, int met)
{
	// met = 0  -  дл€ радиального положени€ HP (x > 0)  phi  the 
	// met = 1  -  дл€ x < 0                             phi  x

	if (met == 0)
	{
		int i1 = 0, i2 = 0; // по углу phi
		int j1 = 0, j2 = 0; // по углу the

		j2 = this->the_angle.shape()[1] - 1;
		j1 = j2 - 1;

		i2 = this->phi_angle.size() - 1;
		i1 = i2 - 1;

		for (int i = 0; i < this->phi_angle.size(); i++)
		{
			if (this->phi_angle[i] > phi)
			{
				i2 = i;
				i1 = i2 - 1;
				if (i1 < 0) i1 = this->phi_angle.size() - 1;
				break;
			}
		}

		if (the > const_pi / 2.0 + 0.0000001) 
		{
			cout << "ERROR 0912312603" << endl;
		}
		
		// ¬ головной области
		for (int i = 0; i < this->the_angle.shape()[1]; i++)
		{
			if (this->the_angle[i1][i] > the)
			{
				j2 = i;
				j1 = j2 - 1;
				if (j1 < 0)
				{
					j1 = 0;
					j2 = 1;
				}
				break;
			}
		}

		double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
		double y = (the - this->the_angle[i1][j1]) / (this->the_angle[i1][j2] - this->the_angle[i1][j1]);

		if (i1 == this->phi_angle.size() - 1)
		{
			x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
		}

		return this->HP_radial[i1][j1] * (1 - x) * (1 - y) +
			this->HP_radial[i2][j1] * x * (1 - y) +
			this->HP_radial[i1][j2] * (1 - x) * y +
			this->HP_radial[i2][j2] * x * y;
	}
	else if(met == 1)
	{
		int i1 = 0, i2 = 0; // по углу phi
		int j1 = 0, j2 = 0; // по координате x

		i2 = this->phi_angle.size() - 1;
		i1 = i2 - 1;

		j2 = this->x_cilindr.shape()[1] - 1;
		j1 = j2 - 1;

		for (int i = 0; i < this->phi_angle.size(); i++)
		{
			if (this->phi_angle[i] > phi)
			{
				i2 = i;
				i1 = i2 - 1;
				if (i1 < 0) i1 = this->phi_angle.size() - 1;
				break;
			}
		}

		if (the > 0.0) cout << "ERROR 1964960286" << endl;

		if (the < this->x_cilindr[i1][this->x_cilindr.shape()[1] - 1] || 
			the < this->x_cilindr[i2][this->x_cilindr.shape()[1] - 1])
		{
			double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);

			if (i1 == this->phi_angle.size() - 1)
			{
				x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
			}

			return this->HP_cilindr[i1][this->HP_cilindr.shape()[1] - 1] * (1 - x) + 
				this->HP_cilindr[i2][this->HP_cilindr.shape()[1] - 1] * x;
		}
		else if (the < this->x_cilindr[i1][0])
		{

			for (int i = 0; i < this->x_cilindr.shape()[1]; i++)
			{
				if (this->x_cilindr[i1][i] < the)
				{
					j2 = i;
					j1 = j2 - 1;
					if (j1 < 0)
					{
						j2 = 1;
						j1 = 0;
					}
					break;
				}
			}

			double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
			double y = (the - this->x_cilindr[i1][j1]) / (this->x_cilindr[i1][j2] - this->x_cilindr[i1][j1]);

			if (i1 == this->phi_angle.size() - 1)
			{
				x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
			}


			return this->HP_cilindr[i1][j1] * (1 - x) * (1 - y) +
				this->HP_cilindr[i2][j1] * x * (1 - y) +
				this->HP_cilindr[i1][j2] * (1 - x) * y +
				this->HP_cilindr[i2][j2] * x * y;
		}
		else
		{
			double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
			double y = (the - 0.0) / (this->x_cilindr[i1][0] - 0.0);

			if (i1 == this->phi_angle.size() - 1)
			{
				x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
			}

			return this->HP_radial[i1][this->HP_radial.shape()[1] - 1] * (1 - x) * (1 - y) +
				this->HP_radial[i2][this->HP_radial.shape()[1] - 1] * x * (1 - y) +
				this->HP_cilindr[i1][0] * (1 - x) * y +
				this->HP_cilindr[i2][0] * x * y;
		}
	}
	else
	{
		cout << "ERROR  9897095454" << endl;
		exit(-1);
	}

}

double Surfaces::Get_BS(const double& phi, const double& the)
{
	int i1 = 0, i2 = 0; // по углу phi
	int j1 = 0, j2 = 0; // по углу the

	i1 = this->phi_angle.size() - 1;

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		if (this->phi_angle[i] > phi)
		{
			i2 = i;
			i1 = i2 - 1;
			if (i1 < 0) i1 = this->phi_angle.size() - 1;
			break;
		}
	}

	/*if (the == const_pi / 2.0 && phi == 0.0)
	{
		cout << "edrge" << endl;
	}*/


	if (the < const_pi / 2.0 + 0.000001)
	{
		// ¬ головной области
		j2 = this->the_angle.shape()[1] - 1;
		j1 = j2 - 1;

		for (int i = 0; i < this->the_angle.shape()[1]; i++)
		{
			if (this->the_angle[i1][i] > the)
			{
				j2 = i;
				j1 = j2 - 1;
				if (j1 < 0)
				{
					j1 = 0;
					j2 = 1;
				}
				break;
			}
		}

		double x = (phi - this->phi_angle[i1]) / (this->phi_angle[i2] - this->phi_angle[i1]);
		double y = (the - this->the_angle[i1][j1]) / (this->the_angle[i1][j2] - this->the_angle[i1][j1]);

		if (i1 == this->phi_angle.size() - 1)
		{
			x = (phi - this->phi_angle[i1]) / (2.0 * const_pi + this->phi_angle[i2] - this->phi_angle[i1]);
		}

		double b = this->BS_radial[i1][j1] * (1 - x) * (1 - y) +
			this->BS_radial[i2][j1] * x * (1 - y) +
			this->BS_radial[i1][j2] * (1 - x) * y +
			this->BS_radial[i2][j2] * x * y;

		return b;

	}
	else
	{
		cout << "error 9890901432" << endl;
	}


	return 0.0;
}

void Surfaces::Print_TS()
{
	ofstream fout;
	fout.open("Surface_Print_TS.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		double phi = this->phi_angle[i];
		for (int j = 0; j < this->the_angle.shape()[1]; j++)
		{
			double the = this->the_angle[i][j];
			double r = this->TS_radial[i][j];

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		double phi = this->phi_angle[i];
		for (int j = 0; j < this->the_angle2.shape()[1]; j++)
		{
			double the = this->the_angle2[i][j];
			double r = this->TS_radial2[i][j];

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	fout.close();


	fout.open("Surface_Print_TS_my_point.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int N1 = 100;
	int N2 = 100;

	for (int i = 0; i < N1; i++)
	{
		double phi = i * 2 * const_pi / N1;
		for (int j = 0; j < N2; j++)
		{
			double the = j * const_pi / N2;
			double r = this->Get_TS(phi, the);

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	fout.close();
}

void Surfaces::Print_BS()
{
	ofstream fout;
	fout.open("Surface_Print_BS.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		double phi = this->phi_angle[i];
		for (int j = 0; j < this->the_angle.shape()[1]; j++)
		{
			double the = this->the_angle[i][j];
			double r = this->BS_radial[i][j];

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	fout.close();


	fout.open("Surface_Print_BS_my_point.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int N1 = 100;
	int N2 = 100;

	for (int i = 0; i < N1; i++)
	{
		double phi = i * 2 * const_pi / N1;
		for (int j = 0; j < N2; j++)
		{
			double the = j * const_pi/2 / N2;
			double r = this->Get_BS(phi, the);

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	fout.close();
}

void Surfaces::Print_HP()
{
	ofstream fout;
	fout.open("Surface_Print_HP.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		double phi = this->phi_angle[i];
		for (int j = 0; j < this->the_angle.shape()[1]; j++)
		{
			double the = this->the_angle[i][j];
			double r = this->HP_radial[i][j];

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	for (int i = 0; i < this->phi_angle.size(); i++)
	{
		double phi = this->phi_angle[i];
		for (int j = 0; j < this->x_cilindr.shape()[1]; j++)
		{
			double x = this->x_cilindr[i][j];
			double r = this->HP_cilindr[i][j];

			fout << x << " " << r * cos(phi) << " " << r * sin(phi) << endl;
		}
	}

	fout.close();


	fout.open("Surface_Print_HP_my_point.txt");

	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int N1 = 100;
	int N2 = 50;

	for (int i = 0; i < N1; i++)
	{
		double phi = i * 2 * const_pi / N1;
		for (int j = 0; j < N2; j++)
		{
			double x = -(j + 1) * (50.0)/N2;
			double r = this->Get_HP(phi, x, 1);

			fout << x << " " << r * cos(phi) << " " << r * sin(phi) << endl;
		}
	}

	for (int i = 0; i < N1; i++)
	{
		double phi = i * 2 * const_pi / N1;
		for (int j = 0; j < N2; j++)
		{
			double the = j * const_pi/2.0 / N2;
			double r = this->Get_HP(phi, the, 0);

			fout << r * cos(the) << " " << r * sin(the) * cos(phi) << " " << r * sin(the) * sin(phi) << endl;
		}
	}

	fout.close();
}
