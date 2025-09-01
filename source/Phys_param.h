#pragma once
#include "Header.h"


struct PrintOptions {
    std::optional<double> x;
    std::optional<double> y;
    std::optional<double> z;
    std::optional<string> fluid;
};

class Phys_param
{
public:
    double T_all = 0;      // Всё время счёта программы
    double Omega = 2.0 * const_pi/6.47143;    // Угловая скорость вращения солнца   2 pi/T
    double lambda = 0.52359878;    // Угол магнитного диполя   30 градусов
    double theta0 = 0.26179939;    // Разделение низких и высоких широт    15 градусов
    double B_0 = 1.33992;             // Магнитное поле на внутренней границе
    double p_0 = 0.0913058;             // давление на внутренней границе
    double R_0 = 0.163;             // внутренней границе
    double GM = 0.00439795;             // Безразмерная гравитационная постоянная на массу М G = 0.00439795

    vector<string> param_names;
    // Все имена параметров, наличие которых мы ожидаем в ячейках!
    vector<string> plasma_name; // Имена плазменных значений
    vector<string> plasma_pui_name; // Имена плазменных значений
    vector<string> H_param_names;


    vector<string> pui_name;   // имена дополнительных жидкостей пикапов
    vector<string> MK_param;   // дополнительные параметры для MK

    vector<string> p_pui_name; 
    vector<string> H_name;   // имена дополнительных жидкостей водорода

    std::unordered_set<std::string> r2_snos_names;
    std::unordered_set<std::string> r2g_snos_names;

    unordered_map<string, double> perevod_razmer;


    bool is_PUI;       // Считаем ли пикапы?

    uint8_t num_pui;           // Сколько сортов пикапов в ячейках

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>  pui_in_zone;   
    // Матрица показывает в какой зоне какие сорта пикапов присутствуют


    // Граничные условия заданы для каждого типа поверхности (их три)
    // Inner_Hard, Outer_Hard, Outer_Soft
    // 0 - гран условия не нужны, так как там нет этой компоненты
    // 1 - сносовые   cells[0] | cells[0]
    // 2 - жёсткие  cells[0] | gran
    // 3 - центральная ячейка cells[0] | Center 

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  pui_condition;
    // Матрица граничных условий для пикапов

    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>  plasma_condition;
    // Матрица граничных условий для плазмы

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_condition;
    // Матрица граничных условий для водорода

    // Матрицы взаимодействия сортов по областям
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_1;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_1;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_2;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_2;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_3;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_3;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_4;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_4;

    // Проблема таблиц proton_arise в том, что в них нет "пропущенных столбцов"
    // То есть если в области 1 нет pui_2, то в таблице нет такого столбца
    // это сделано просто, так как pui_2 - это последний столбец
    // Но если в области не будет pui, который не в последнем столбце, то надо переделывать
    // алгоритм вычисления источников.

    std::unordered_map<std::string, bool> Culc_hidrogen;   
    // Какой сорт водорода считается, а какой нет
    // Для увеличения шага по времени иногда полезно какой-то занулить


    bool is_div_V_in_cell;       // Считаем ли дивергенцию скорости в ячейках 
    // на каждом шаге? Нужно, например, для учёта пикапов
    // может и без пикапов понадобиться, поэтому тут отдельный маячок



    std::function<void (const short int&,
        unordered_map<string, double>&,
        unordered_map<string, double>&)> Plasma_components;



    double Velosity_inf = -2.54327;   // Значение скорости смеси на бесконечности
    double B_inf = 8.91006;           // Значение МОДУЛЯ магнитного поля на бесконечности
    double alphaB_inf = 0.6981317;     // угол в радианах - Направление магнитного поля на бесконечности
    double par_a_2 = 0.102046;         // Параметр в сечении перезарядки
    double par_n_H_LISM = 3.0;        // б/р Концентрация водорода на бесконечности 
    double par_Kn = 40.0906;
    double mrho_He_0 = 0.035194;        // Множитель концентрации гелия на 1 АЕ
    double mrho_He_inf = 0.15;        // Множитель концентрации гелия на inf
    double mep = 0.000544617;      // Отношение массы электрона к массе протона
    uint8_t num_H = 4;           // Сколько сортов водорода для МК

    // параметры
    double gamma = (5.0 / 3.0);       // показатель адиабаты
    double g1 = (5.0 / 3.0 - 1.0);       // показатель адиабаты - 1

    // Характерные параметры
    double char_rho = 0.06;             // Характерная концентрация в СГС  см^-3
    double char_v = 10.3804;         // Характерная скорость в км/с


    // Настройки расчёта Плазмы
    double KFL;                   // критерий Куранта
    bool TVD;                   // Делаем ли ТВД для плазмы?
    bool TVD_atom = true;                   // Делаем ли ТВД для атомов?

    bool Snos_on_HP = false;       // Делать ли линейный снос на HP с двух сторон (иначе сносить первым порядком)

    bool culc_plasma;     // НЕ АКТИВЕН      // Считаем ли плазму? Можно заморозить плазму для расчёта водорода
    bool culc_atoms;             // Вычисляем ли атомы или оставляем их вмороженными
    bool move_setka;
    
    
    bool move_TS;               // Двигаем ли HP
    bool move_HP;               // Двигаем ли HP
    bool move_BS;

    bool sglag_TS;              // Делаем ли сглаживание TS
    double velocity_TS;
    double sglag_TS_k;         // Сглаживание на высоких широтах
    double sglag_TS_k_sphere;         // Сглаживание на высоких широтах
    //double sglag_TS_k_sphere = 0.001; // 0.3;   // Сглаживание в головной и хвостовой части
    double sglag_TS_k_sphere_head; // 0.08;   // Сглаживание в головной части
    double sglag_TS_k_sphere_tail; // 0.03;   // Сглаживание в хвостовой части
    // Не может быть больше 1

    bool sglag_HP;
    double velocity_HP;
    double sglag_HP_k_sphere;  //0.005 0.002    // Cглаживание в головной части
    double sglag_HP_k_angle;  //0.005 0.002    // Cглаживание в хвостовой части к локальным окружностям
    double sglag_HP_k; // 0.001         // Сглаживание не в головной области
    double sglag_HP_angle;    // 1.2 коэффициент усиления сглаживания по углу
    double sglag_HP_along;    // коэффициент усиления сглаживания вдоль х
    double sglag_HP_sphere;   // коэффициент усиления сглаживания в головной области - НЕ АКТИВНО


    bool sglag_BS;
    double sglag_BS_k;
    double velocity_BS;


    bool null_bn_on_HP;   // Для ячеек рядом с HP обнуляем нормальную компоненту магнитного поля
    bool bn_in_p_on_HP;   // Для ячеек рядом с HP записываем магнитное поле в давление
    bool contact_hard;    // Надо ли на контакте приравнивать скорость грани к скорости контакта
    bool TS_hard;    // Надо ли на контакте приравнивать скорость грани к скорости контакта

     // Настройки расчёта МК  ********************************************************
    
    bool save_AMR;        // Нужно ли сохранять посчитанные функции распределения?
    bool culc_AMR;        // Нужно ли считать функции распределения?
    bool refine_AMR;      // Нужно ли мельчить посчитанные функции распределения?
    unsigned int N_per_gran;  // Сколько в среднем частиц вылетает с каждой грани
    bool culc_cell_moments;    // Нужно ли считать моменты в ячейках?
    bool de_refine_AMR;        // Нужно ли огрублять AMR сетку, если требуется?
    string MK_file;            // Фаил для сохранения моментов в ячейках

    bool MK_source_S = false;         // Нужно ли считать источники S+ и S- в МК
    int pui_nW;               // Число разбиений массивов S+ и S- в ячейках
    double pui_wR;            // Максимальная относительная скорость в массивах S+ и S-
     




    Eigen::Matrix3d Matr;              // Матрица перехода 1
    Eigen::Matrix3d Matr2;             // Матрица перехода 2


    // Данные для считывания файла параметров на начальной сфере
    std::vector<double> heliolat_deg;
    std::vector<double> n_p_cm3;
    std::vector<double> V_kms;
    std::vector<double> T_K;


    Phys_param();
    void set_parameters(void);

    double Get_rho_0(const double& the);
    double Get_v_0(const double& the);
    double Get_T_0(const double& the);

    double Get_razmer(string par);  // Возвращает коэффициент для перевода 
    // параметра в размерные значения

    // Разделение компонент плазмы
    void Plasma_components_1(const short int& zone,
        unordered_map<string, double>& param_in_cell,
        unordered_map<string, double>& param);

    void Plasma_components_2(const short int& zone,
        unordered_map<string, double>& param_in_cell,
        unordered_map<string, double>& param);


    void Get_Potok(const double& rho, const double& p, const double& u, const double& v,
        const double& w, const double& bx, const double& by, const double& bz,
        const double& n1, const double& n2, const double& n3, 
        std::vector<double>& Potok);  // Potok - 8-мерный вектор

    double energy(const double& rho, const double& p, const double& vv, const double& bb);

    void chlld(unsigned short int n_state, // метод
        const double& al, const double& be, const double& ge, // нормаль
        const double& w, // скорость грани
        const std::vector<double>& qqq1, const std::vector<double>& qqq2, // основные параметры с двух сторон 
        std::vector<double>& qqq, // Выходной поток
        bool null_bn,  // если  true, то нужно обнулить поток магнитного поля через поверхность
        unsigned short int n_disc,    // формула для определения скорости
        const std::vector<double>& konvect_left,  // Дополнительные переменные конвективного переноса слева
        const std::vector<double>& konvect_right, // Дополнительные переменные конвективного переноса справа
        std::vector<double>& konvect, // Дополнительные переменные конвективного переноса ПОТОКИ
        double& dsr, double& dsc, double& dsl,
        PrintOptions& opts, bool left_ydar = false, bool contact = false);

    // Блок Годунова ********
    void lev(const double& enI, const double& pI, const double& rI, const double& enII,//
        const double& pII, const double& rII, double& uuu, double& fee);
    void devtwo(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
        const double& w, double& p);
    void newton(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
        const double& w, double& p);
    void Godunov_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
        const std::vector<double>& n, std::vector<double>& qqq,//
        double& dsl, double& dsp, double& dsc, double w, bool contact = false);

    // ***********

    void raspad_testing(void);

private:
    using VarRef = std::variant<
        double*,
        int*,
        uint8_t*,
        unsigned int*,
        bool*,
        std::string*
    >;

    std::unordered_map<std::string, VarRef> varMap;

    void initVarMap();
    void parseAndAssign(VarRef varRef, const std::string& valueStr);
};

