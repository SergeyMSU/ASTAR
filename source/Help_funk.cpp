#include "Help_funk.h"
#include "Yzel.h"


Eigen::Matrix3d rotationMatrixZ(double angle_rad) 
{
	Eigen::Matrix3d R;
	double cos_angle = std::cos(angle_rad);
	double sin_angle = std::sin(angle_rad);

	R << cos_angle, -sin_angle, 0,
		sin_angle, cos_angle, 0,
		0, 0, 1;

	return R;
}

Eigen::Matrix3d rotationMatrixX(double angle_rad) 
{
	Eigen::Matrix3d R;
	double cos_angle = std::cos(angle_rad);
	double sin_angle = std::sin(angle_rad);

	R << 1, 0, 0,
		0, cos_angle, -sin_angle,
		0, sin_angle, cos_angle;

	return R;
}

void pause_seconds(unsigned int n) 
{
	cout << "Pause - " << n << endl;
	std::this_thread::sleep_for(std::chrono::seconds(n));
	cout << "End  Pause - " << n << endl;
}


bool file_exists(const std::string& filename) {
	std::ifstream file(filename);
	return file.good();  // или просто return file.is_open();
}

double maxwell(const double& n_H, const double& cp, const double& u1,
	const double& u2, const double& u3,
	const double& x, const double& y, const double& z)
{
	return n_H * pow((sqrt(const_pi) * cp), -3) * exp((-kv(x - u1)
		- kv(y - u2) - kv(z - u3)) / kv(cp));
}

short int signum(const double& x)
{
	if (x > 0.00000001)
	{
		return 1;
	}
	else if (x < -0.00000001)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}


double minmod(const double& x, const double& y)
{
	if (signum(x) + signum(y) == 0)
	{
		return 0.0;
	}
	else
	{
		return  ((signum(x) + signum(y)) / 2.0) * min(fabs(x), fabs(y));
	}
}

double linear(const double& x1, const double& t1, const double& x2, 
	const double& t2, const double& x3, const double& t3, 
	const double& y)
{
	//Главное значение с параметрами 2
	//Строим линии между 1 и 2, 2 и 3, потом находим минмодом значение в y
	double d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3));
		return (d * (y - x2) + t2);
}

double linear2(const double& x1, const double& t1, const double& x2,
	const double& t2, const double& y)
{
	double d = (t1 - t2) / (x1 - x2);
	return (d * (y - x2) + t2);
}


double rbfKernel(double r, double epsilon) 
{
	return std::sqrt(1 + (epsilon * r) * (epsilon * r));
}

double polar_angle(const double& x, const double& y)
{
	if (fabs(x) + fabs(y) < 0.000001)
	{
		return 0.0;
	}

	if (x < 0)
	{
		return atan(y / x) + 1.0 * const_pi;
	}
	else if (x > 0 && y >= 0)
	{
		return atan(y / x);
	}
	else if (x > 0 && y < 0)
	{
		return atan(y / x) + 2.0 * const_pi;
	}
	else if (y > 0 && x >= 0 && x <= 0)
	{
		return const_pi / 2.0;
	}
	else if (y < 0 && x >= 0 && x <= 0)
	{
		return  3.0 * const_pi / 2.0;
	}
	return 0.0;
}

void dekard_skorost(const double& z, const double& x, const double& y, 
	const double& Vr, const double& Vphi, const double& Vtheta, 
	double& Vz, double& Vx, double& Vy)
{
	double r_2, the_2, phi_2;
	r_2 = sqrt(x * x + y * y + z * z);
	the_2 = acos(z / r_2);
	phi_2 = polar_angle(x, y);

	Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
	Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
	Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
}

void spherical_skorost(const double& z, const double& x, const double& y, 
	const double& Vz, const double& Vx, const double& Vy,
	double& Vr, double& Vphi, double& Vtheta)
{
	double r_1, the_1, phi_1;

	r_1 = sqrt(x * x + y * y + z * z);
	the_1 = acos(z / r_1);
	phi_1 = polar_angle(x, y);

	Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
	Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
	Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
}

double Yzel_distance(Yzel* A, Yzel* B, int time)
{
	return sqrt(kv(A->coord[time][0] - B->coord[time][0]) +
		kv(A->coord[time][1] - B->coord[time][1]) + kv(A->coord[time][2] - B->coord[time][2]));
}

double Yzel_distance_x(Yzel* A, Yzel* B, int time)
{
	return sqrt(kv(A->coord[time][1] - B->coord[time][1]) + kv(A->coord[time][2] - B->coord[time][2]));
}

bool areCellsEqual(const Gran& cell1, const Gran& cell2) {
	if (cell1.yzels.size() != cell2.yzels.size()) {
		return false;
	}

	// Создаем мультимножество для каждого вектора
	std::unordered_multiset<Yzel*> set1(cell1.yzels.begin(), cell1.yzels.end());
	std::unordered_multiset<Yzel*> set2(cell2.yzels.begin(), cell2.yzels.end());

	return set1 == set2;
}

bool areCellsEqual(const Gran* cell1, const Gran* cell2) 
{
	if (cell1->yzels.size() != cell2->yzels.size()) 
	{
		return false;
	}

	// Создаем мультимножество для каждого вектора
	std::unordered_multiset<Yzel*> set1(cell1->yzels.begin(), cell1->yzels.end());
	std::unordered_multiset<Yzel*> set2(cell2->yzels.begin(), cell2->yzels.end());

	return set1 == set2;
}

bool areCellsEqual_my(const Gran* cell1, const Gran* cell2)
{
	if (cell1->yzels.size() != cell2->yzels.size())
	{
		return false;
	}

	for (int i = 0; i < cell1->yzels.size(); i++)
	{
		if (std::find(cell2->yzels.begin(), cell2->yzels.end(), cell1->yzels[i]) == cell2->yzels.end())
		{
			return false;
		}
	}

	return true;
}


double triangleArea3D(
	const double & x1, const double& y1, const double& z1,
	const double& x2, const double& y2, const double& z2,
	const double& x3, const double& y3, const double& z3)
{
	// Вычисляем векторы сторон
	double v1x = x2 - x1;
	double v1y = y2 - y1;
	double v1z = z2 - z1;

	double v2x = x3 - x1;
	double v2y = y3 - y1;
	double v2z = z3 - z1;

	// Вычисляем векторное произведение
	double cross_x = v1y * v2z - v1z * v2y;
	double cross_y = v1z * v2x - v1x * v2z;
	double cross_z = v1x * v2y - v1y * v2x;

	// Вычисляем длину вектора произведения
	double cross_magnitude = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);

	// Площадь равна половине длины векторного произведения
	return 0.5 * cross_magnitude;
}

void get_bazis(const Eigen::Vector3d& n, Eigen::Vector3d& t, Eigen::Vector3d& m)
{
	// 1. Выбираем произвольный вектор 'a', не коллинеарный с 'n'
	Eigen::Vector3d a(1.0, 0.0, 0.0);
	if (std::abs(n.dot(a)) > 0.9) {  // Если n почти параллелен (1,0,0), берём другой a
		a = Eigen::Vector3d(0.0, 1.0, 0.0);
	}

	// 2. Находим t = (a × n) / |a × n|
	t = a.cross(n).normalized();

	// 3. Находим m = n × t (уже нормирован, так как n и t ортогональны и единичные)
	m = n.cross(t);
}

void crossProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	double& rx, double& ry, double& rz)
{
	rx = ay * bz - az * by;
	ry = az * bx - ax * bz;
	rz = ax * by - ay * bx;
}

double scalarProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz)
{
	return ax * bx + ay * by + az * bz;
}

double tetrahedronVolume(
	const double& x1, const double& y1, const double& z1,  // Вершина A
	const double& x2, const double& y2, const double& z2,  // Вершина B
	const double& x3, const double& y3, const double& z3,  // Вершина C
	const double& x4, const double& y4, const double& z4)  // Вершина D
{
	// Вычисляем векторы рёбер из вершины A
	double abx = x2 - x1, aby = y2 - y1, abz = z2 - z1;  // Вектор AB
	double acx = x3 - x1, acy = y3 - y1, acz = z3 - z1;  // Вектор AC
	double adx = x4 - x1, ady = y4 - y1, adz = z4 - z1;  // Вектор AD

	// Вычисляем векторное произведение AB × AC
	double cross_x = aby * acz - abz * acy;
	double cross_y = abz * acx - abx * acz;
	double cross_z = abx * acy - aby * acx;

	// Вычисляем скалярное произведение (AB × AC) · AD
	double dot_product = cross_x * adx + cross_y * ady + cross_z * adz;

	// Объём равен 1/6 модуля этого произведения
	return std::abs(dot_product) / 6.0;
}


Cell* Get_Sosed(Cell* C, Gran* gr)
{
	if (gr->cells.size() == 1) return nullptr;

	if (gr->cells[0] == C)
	{
		return gr->cells[1];
	}
	else
	{
		return gr->cells[0];
	}
	return nullptr;
}

void solveQuadraticEquation(const double& x1, const double& y1,
	const double& x2, const double& y2,
	const double& x3, const double& y3, double& a, double& b, double& c)
{
	double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	if (denom == 0) {
		cout << "Points are collinear or have duplicate x-coordinates." << endl;
		cout << "ERROR  7543851320" << endl;
		exit(-1);
	}

	a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
	c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	return;
}


void Sootnosheniya(const double& rho, const double& p, const double& rho_He, 
	const double& rho_Pui, const double& T_Pui, 
	const short int& zone, 
	double& rho_Th, double& rho_E, double& p_Th, double& p_Pui, 
	double& T_Th, double& T_E)
{
	// Функция, определяющая температуры и концентрации гелия, пикапов и т.д.
	// al - это заряд гелия
	// если al = 1 то вне гелиопаузы
	// если al = 2, то внутри гелиопаузы
	// Th - термальные протоны, He - гелий, Pui - пикапы, E - электроны
	// без параметров, это общие(те, что считаются в МГД)
	short int al;

	if (zone <= 2)
	{
		al = 2;
	}
	else
	{
		al = 1;
	}
	
	rho_Th = -(MF_meDmp * al * rho_He + 4.0 * (-rho + rho_He)) / (4.0 * (1.0 + MF_meDmp)) - rho_Pui;

	rho_E = MF_meDmp * (4.0 * rho + (-4.0 + al) * rho_He) / (4.0 * (1.0 + MF_meDmp));
			

	p_Th = (p - rho_Pui * T_Pui) * (-4.0 * rho + (4.0 + al * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui) / 
		(-8.0 * rho + (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui);
					
	p_Pui = T_Pui * rho_Pui;
						
	T_Th = -4.0 * (1.0 + MF_meDmp) * (p - rho_Pui * T_Pui) /
		(-8.0 * rho + (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He + 4.0 * (1.0 - MF_meDmp) * rho_Pui);
	
	T_E = T_Th;

	return;
}

bool findIntersection(const std::array<double, 3>& P1, const std::array<double, 3>& P2,
	const double& a, const double& b, const double& c, const double& d,
	std::array<double, 3>& outIntersection)
{
	double D = a * P1[0] + b * P1[1] + c * P1[2] + d;
	double N = a * (P2[0] - P1[0]) + b * (P2[1] - P1[1]) + c * (P2[2] - P1[2]);

	// Отрезок параллелен плоскости
	if (std::abs(N) < 1e-10)
	{
		if (std::abs(D) < 1e-10)
		{
			// Отрезок лежит в плоскости (бесконечно много пересечений)
			// Можно вернуть, например, P1 или P2
			outIntersection = P1;
			return true;
		}
		else
		{
			// Нет пересечения
			return false;
		}
	}

	double t = -D / N;

	// Проверяем, что t ∈ [0, 1]
	if (t >= 0.0 && t <= 1.0)
	{
		outIntersection[0] = P1[0] + t * (P2[0] - P1[0]);
		outIntersection[1] = P1[1] + t * (P2[1] - P1[1]);
		outIntersection[2] = P1[2] + t * (P2[2] - P1[2]);
		return true;
	}
	else
	{
		// Пересечение за пределами отрезка
		return false;
	}
}