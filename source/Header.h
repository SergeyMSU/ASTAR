#pragma once


// Описание всех классов
class Geo_param; 
class AMR_f;
class AMR_cell;
class Phys_param;
class Setka;
class Yzel;
class Luch;
class Cell;
class Gran;
class Edge;
class Surfaces;
class Interpol;
class Int_point;
class MK_particle;


#define kv(x) ((x) * (x))
#define pow3(x) ((x) * (x) * (x))
#define kvg(x) (pow(x, 2.0 * this->phys_param->gamma))
#define kyb(x) ((x) * (x) * (x))
#define kvv(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))
#define norm2(x, y, z) (sqrt(kvv(x, y, z)))
#define whach(x) cout << #x <<": " << (x) << endl
#define max3(a, b, c) (max(max( (a) , (b) ), (c) ))
#define min3(a, b, c) (min(min( (a) , (b) ), (c) ))


#include <filesystem>
#include <string>
#include <utility>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <iterator>
#include <cstdlib>
#include <Eigen/Dense>
#include <optional>
#include <omp.h>
#include <chrono>
#include <random>

#pragma warning(disable : 4146)  // Отключаем предупреждение об унарном минусе
#include <gmp.h>  // Теперь это будет работать!
#pragma warning(default : 4146)  // Восстанавливаем настройки

#define a_2 (0.1307345665) 

inline double sigma(double x)
{
    return kv(1.0 - a_2 * log(x));
}

inline double sigma2(double x, double y)
{
    return kv(1.0 - a_2 * log((x) * (y)));
}


// Boost библиотека (надо подключать к компилятору отдельно)

#include <boost/parameter.hpp>
#include "boost/multi_array.hpp"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/barycenter.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel KKexact;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, KKexact> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, KKexact> Vb2;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Triangulation_data_structure_2<Vb2> Tds2;
typedef CGAL::Delaunay_triangulation_3<KKexact, Tds> Delaunay;
typedef CGAL::Delaunay_triangulation_2<KKexact, Tds2> Delaunay2;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay2::Vertex_handle Vertex_handle2;
typedef Delaunay2::Face_handle Face_handle;
typedef KKexact::FT FT;


#define sqrtpi_ 1.77245385
const double const_pi = 3.14159265358979323846;
const double cpi4 = 4.0 * const_pi;
const double cpi8 = 8.0 * const_pi;
const double spi4 = sqrt(cpi4);
const double eps = 1E-12;
const double epsb = 1E-4;
const double eps_p = 1E-6;
const double eps_d = 1E-3;
const double MF_meDmp = (1.0 / 1836.15);  // Отношения массы электрона к массе протона



using namespace std;

#include "sensor.h"
#include "AMR_f.h"
#include "AMR_cell.h"
#include "Phys_param.h"
#include "Geo_param.h"
#include "Setka.h"
#include "Help_funk.h"
#include "Luch.h"
#include "Gran.h"
#include "Edge.h"
#include "Cell.h"
#include "Yzel.h"
#include "Surfaces.h"
#include "Interpol.h"
#include "Int_point.h"
#include "MK_particle.h"






// Описание всех типов комметариев

// &INIT&  -  переменная инициализации (возможно нужно будет изменить)
