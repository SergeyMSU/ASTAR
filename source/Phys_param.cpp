#include "Phys_param.h"

#define ga (this->gamma)
#define ggg (ga)
#define gg1 (ga - 1.0)
#define g1 (ga - 1.0)
#define g2 (ga + 1.0)
#define gg2 (ga + 1.0)
#define gp ((g2/ga)/2.0)
#define gm ((g1/ga)/2.0)
#define gga ga
#define ER_S std::cout << "\n---------------------\nStandart error in file: Solvers.cpp\n" << endl
#define watch(x) cout << (#x) << " is " << (x) << endl
#define M(x) cout << (#x)  << endl

#define zone_info false

Phys_param::Phys_param()
{
    this->initVarMap();

    this->set_parameters();

    this->num_H = 0;
    this->num_pui = 0;

    if (this->is_PUI == true)
    {
        this->num_H = 9;
        this->is_div_V_in_cell = true;
        this->num_pui = 2;
        this->pui_in_zone.resize(4, this->num_pui);
        this->pui_in_zone << true, false,
                             true, true,
                             true, false,
                             true, false;
        // Пикапы будут задаваться везде (для простоты), просто где их нет будут определены нулями
        // Предупреждение, может быть проблема с тем, что некоторые
        // поверхности выделяются не полностью

        this->pui_condition.resize(3, this->num_pui);
        this->pui_condition << 2, 0,
                               2, 0,
                               1, 1;


        this->hydrogen_condition.resize(3, this->num_H);
        this->hydrogen_condition << 2, 3, 3, 3, 2, 3, 3, 3, 3,
                                    1, 1, 1, 2, 1, 1, 1, 1, 2,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1;

        this->hydrogen_arise_1.resize(this->num_H, this->num_pui + 1 - 1);
        this->hydrogen_arise_1 << 1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5,
                                  1, 5;

        this->hydrogen_arise_2.resize(this->num_H, this->num_pui + 1);
        this->hydrogen_arise_2 <<   2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7,
                                    2, 6, 7;

        this->hydrogen_arise_3.resize(this->num_H, this->num_pui + 1 - 1);
        this->hydrogen_arise_3 <<   3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8,
                                    3, 8;

        this->hydrogen_arise_4.resize(this->num_H, this->num_pui + 1 - 1);
        this->hydrogen_arise_4 <<   4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9,
                                    4, 9;

        proton_arise_1.resize(this->num_H, this->num_pui + 1 - 1);
        this->proton_arise_1 << 0, 0,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1;

        proton_arise_2.resize(this->num_H, this->num_pui + 1);
        this->proton_arise_2 << 2, 2, 2,
                                0, 0, 0,
                                2, 2, 2,
                                2, 2, 2,
                                2, 2, 2,
                                1, 1, 1,
                                2, 2, 2,
                                2, 2, 2,
                                2, 2, 2;

        proton_arise_3.resize(this->num_H, this->num_pui + 1 - 1);
        this->proton_arise_3 << 1, 1,
                                1, 1,
                                0, 0,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1;

        proton_arise_4.resize(this->num_H, this->num_pui + 1 - 1);
        this->proton_arise_4 << 1, 1,
                                1, 1,
                                1, 1,
                                0, 0,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1,
                                1, 1;


    }
    else if(false)
    {
        this->num_H = 0;
        this->num_pui = 0;

        this->hydrogen_condition.resize(3, this->num_H);
        this->hydrogen_condition << 2, 3, 3, 3,
                                    1, 1, 1, 2,
                                    1, 1, 1, 1;

        this->hydrogen_arise_1.resize(this->num_H, this->num_pui + 1);
        this->hydrogen_arise_1 << 1,
                                  1,
                                  1,
                                  1;

        this->hydrogen_arise_2.resize(this->num_H, this->num_pui + 1);
        this->hydrogen_arise_2 << 2,
                                  2,
                                  2,
                                  2;

        this->hydrogen_arise_3.resize(this->num_H, this->num_pui + 1);
        this->hydrogen_arise_3 << 3,
                                  3,
                                  3,
                                  3;

        this->hydrogen_arise_4.resize(this->num_H, this->num_pui + 1);
        this->hydrogen_arise_4 << 4,
                                  4,
                                  4,
                                  4;

        proton_arise_1.resize(this->num_H, this->num_pui + 1);
        this->proton_arise_1 << 0,
                                0,
                                0,
                                0;

        proton_arise_2.resize(this->num_H, this->num_pui + 1);
        this->proton_arise_2 << 0,
                                0,
                                0,
                                0;

        proton_arise_3.resize(this->num_H, this->num_pui + 1);
        this->proton_arise_3 << 0,
                                0,
                                0,
                                0;

        proton_arise_4.resize(this->num_H, this->num_pui + 1);
        this->proton_arise_4 << 0,
                                0,
                                0,
                                0;
        
    }

    // 1 - 4 старые обычные сорта
    // 5 - 8 сорта водорода, рождённые от пикапов сорта 1
    // 9 - сорт водорода от пикапов сорта 2 (такой сорт есть только во внутреннем ударном слое).
    
    this->plasma_condition.resize(3, 1);
    this->plasma_condition << 2,
                              2,
                              1;

    if (this->is_PUI == true)
    {
        this->Plasma_components = [this](const short int& zone,
            unordered_map<string, double>& param_in_cell,
            unordered_map<string, double>& param) {
                this->Plasma_components_2(zone, param_in_cell, param); };
    }
    else
    {
        this->Plasma_components = [this](const short int& zone,
            unordered_map<string, double>& param_in_cell,
            unordered_map<string, double>& param) {
                this->Plasma_components_1(zone, param_in_cell, param); };
    }


    this->param_names.push_back("rho");  this->plasma_name.push_back("rho"); this->plasma_pui_name.push_back("rho");
    this->param_names.push_back("p"); this->plasma_name.push_back("p"); this->plasma_pui_name.push_back("p");
    this->param_names.push_back("Vx"); this->plasma_name.push_back("Vx"); this->plasma_pui_name.push_back("Vx");
    this->param_names.push_back("Vy"); this->plasma_name.push_back("Vy"); this->plasma_pui_name.push_back("Vy");
    this->param_names.push_back("Vz"); this->plasma_name.push_back("Vz"); this->plasma_pui_name.push_back("Vz");
    this->param_names.push_back("Bx"); this->plasma_name.push_back("Bx"); this->plasma_pui_name.push_back("Bx");
    this->param_names.push_back("By"); this->plasma_name.push_back("By"); this->plasma_pui_name.push_back("By");
    this->param_names.push_back("Bz"); this->plasma_name.push_back("Bz"); this->plasma_pui_name.push_back("Bz");
    //this->param_names.push_back("Q"); this->plasma_name.push_back("Q"); this->plasma_pui_name.push_back("Q");
    //this->param_names.push_back("rho_He"); this->plasma_name.push_back("rho_He"); this->plasma_pui_name.push_back("rho_He");




    // Задаём имена дополнительныхъ жидкостей пикапов
    if (this->is_PUI == true)
    {
        this->pui_name.push_back("_Pui_1");
        this->pui_name.push_back("_Pui_2");
    }

    this->p_pui_name.push_back("_p");
    for (auto& nam : this->pui_name)
    {
        this->p_pui_name.push_back(nam);
    }

    // Именя для особого сноса в ТВД
    this->r2_snos_names.insert("rho");

    this->r2g_snos_names.insert("p");



    // Задаём имена дополнительных жидкостей водорода
    for (size_t ii = 1; ii <= this->num_H; ii++)
    {
        string nii = "_H" + to_string(ii);
        this->H_name.push_back(nii);
        this->Culc_hidrogen[nii] = true;
    }

    // TODO!
    //this->Culc_hidrogen["_H6"] = false;
    //this->Culc_hidrogen["_H5"] = false;



    // Перевод в размерные единицы   СГС
    this->perevod_razmer["r"] = 1.0;
    this->perevod_razmer["rho"] = 1.0;
    this->perevod_razmer["V"] = 1.0;
    this->perevod_razmer["p"] = 1.0;
    this->perevod_razmer["B"] = 1.0;    // Магнитное поле в микрогауссах
    this->perevod_razmer["T"] = 1.0;


	this->Matr << -0.9958639688067077,  0.01776569097515556,  0.08910295088675518,
		           0.07561695085992419, 0.7057402284561812,   0.7044237408557894,
		          -0.05036896241933166, 0.7082479157489926,  -0.7041646522383864;

	this->Matr2 << -0.9958639688067080, 0.0756169508599243, -0.0503689624193315,
	            	0.0177656909751554, 0.7057402284561816,  0.7082479157489927,
		            0.0891029508867553, 0.7044237408557898, -0.7041646522383865;

    // Открываем файл для чтения
    std::ifstream file("nVT1au_new_2000-2022av.dat");
    if (!file.is_open()) {
        std::cerr << "Error 6564397608" << std::endl;
        exit(-1);
    }

    std::string line;
    bool is_first_line = true;

    // Читаем файл построчно
    while (std::getline(file, line)) {
        // Пропускаем первую строку (заголовок)
        if (is_first_line) {
            is_first_line = false;
            continue;
        }

        // Разбиваем строку на отдельные значения
        std::istringstream iss(line);
        double hlat, np, v, t;

        // Если строка содержит данные, записываем их в векторы
        if (iss >> hlat >> np >> v >> t) {
            this->heliolat_deg.push_back(hlat);
            this->n_p_cm3.push_back(np);
            this->V_kms.push_back(v);
            this->T_K.push_back(t);

            cout << hlat << " " << np << " " << v << " " << t << endl;
        }
    }

    file.close();
}

void Phys_param::initVarMap() 
{
    this->varMap = 
    {
    {"is_div_V_in_cell", VarRef(&this->is_div_V_in_cell)},
    {"is_PUI", VarRef(&this->is_PUI)},
    {"num_pui", VarRef(&this->num_pui)},
    {"KFL", VarRef(&this->KFL)},
    {"TVD", VarRef(&this->TVD)},
    {"culc_plasma", VarRef(&this->culc_plasma)},
    {"culc_atoms", VarRef(&this->culc_atoms)},
    {"move_setka", VarRef(&this->move_setka)},
    {"move_TS", VarRef(&this->move_TS)},
    {"move_HP", VarRef(&this->move_HP)},
    {"move_BS", VarRef(&this->move_BS)},
    {"sglag_TS", VarRef(&this->sglag_TS)},
    {"velocity_TS", VarRef(&this->velocity_TS)},
    {"sglag_TS_k_sphere", VarRef(&this->sglag_TS_k_sphere)},
    {"sglag_TS_k", VarRef(&this->sglag_TS_k)},
    {"sglag_TS_k_sphere_head", VarRef(&this->sglag_TS_k_sphere_head)},
    {"sglag_TS_k_sphere_tail", VarRef(&this->sglag_TS_k_sphere_tail)},
    {"sglag_HP", VarRef(&this->sglag_HP)},
    {"velocity_HP", VarRef(&this->velocity_HP)},
    {"sglag_HP_k_sphere", VarRef(&this->sglag_HP_k_sphere)},
    {"sglag_HP_k", VarRef(&this->sglag_HP_k)},
    {"sglag_HP_angle", VarRef(&this->sglag_HP_angle)},
    {"sglag_HP_along", VarRef(&this->sglag_HP_along)},
    {"sglag_HP_sphere", VarRef(&this->sglag_HP_sphere)},
    {"sglag_BS", VarRef(&this->sglag_BS)},
    {"sglag_BS_k", VarRef(&this->sglag_BS_k)},
    {"null_bn_on_HP", VarRef(&this->null_bn_on_HP)},
    {"bn_in_p_on_HP", VarRef(&this->bn_in_p_on_HP)},
    {"contact_hard", VarRef(&this->contact_hard)},
    {"TS_hard", VarRef(&this->TS_hard)},
    {"save_AMR", VarRef(&this->save_AMR)},
    {"culc_AMR", VarRef(&this->culc_AMR)},
    {"refine_AMR", VarRef(&this->refine_AMR)},
    {"N_per_gran", VarRef(&this->N_per_gran)},
    {"culc_cell_moments", VarRef(&this->culc_cell_moments)},
    {"de_refine_AMR", VarRef(&this->de_refine_AMR)},
    {"MK_file", VarRef(&this->MK_file)},
    {"velocity_BS", VarRef(&this->velocity_BS)},
    {"sglag_HP_k_angle", VarRef(&this->sglag_HP_k_angle)},
    {"pui_nW", VarRef(&this->pui_nW)},
    {"pui_wR", VarRef(&this->pui_wR)},
    {"TVD_atom", VarRef(&this->TVD_atom)},
    {"Snos_on_HP", VarRef(&this->Snos_on_HP)},
    {"MK_source_S", VarRef(&this->MK_source_S)}
    };
}

// Функция для парсинга значения из строки и записи в переменную
void Phys_param::parseAndAssign(VarRef varRef, const std::string& valueStr)
{
    try {
        if (auto p = std::get_if<double*>(&varRef)) {
            **p = std::stod(valueStr);
        }
        else if (auto p = std::get_if<int*>(&varRef)) {
            **p = std::stoi(valueStr);
        }
        else if (auto p = std::get_if<unsigned int*>(&varRef)) {
            unsigned long val = std::stoul(valueStr);
            if (val > std::numeric_limits<unsigned int>::max()) {
                throw std::out_of_range("Value exceeds unsigned int range");
            }
            **p = static_cast<unsigned int>(val);
        }
        else if (auto p = std::get_if<uint8_t*>(&varRef)) {
            int temp = std::stoi(valueStr);
            if (temp < 0 || temp > 255) {
                throw std::out_of_range("uint8_t value out of range");
            }
            **p = static_cast<uint8_t>(temp);
        }
        else if (auto p = std::get_if<bool*>(&varRef)) {
            if (valueStr == "true") **p = true;
            else if (valueStr == "false") **p = false;
            else throw std::invalid_argument("Not a boolean");
        }
        else if (auto p = std::get_if<std::string*>(&varRef)) {
            **p = valueStr;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error parsing value '" << valueStr << "': " << e.what() << "\n";
    }
}

void Phys_param::set_parameters(void)
{
    // Старый вариант задания параметров (есть не все, он не актуален)
    if (false) // задание параметров в ручную
    {
        this->is_div_V_in_cell = false;

        this->is_PUI = true;       // Считаем ли пикапы?
        this->num_pui = 2;         // Сколько сортов пикапов в ячейках

        // Настройки расчёта Плазмы
        this->KFL = 0.8;                   // критерий Куранта
        this->TVD = false;                  // Делаем ли ТВД?

        this->culc_plasma = true;          // Считаем ли плазму? Можно заморозить плазму для расчёта водорода
        this->culc_atoms = true;           // Вычисляем ли атомы или оставляем их вмороженными
        this->move_setka = true;


        this->move_TS = true;               // Двигаем ли TS
        this->move_HP = true;               // Двигаем ли HP
        this->move_BS = true;

        this->sglag_TS = true;              // Делаем ли сглаживание TS
        this->velocity_TS = 0.1;            // 0.05
        this->sglag_TS_k_sphere = 0.01;  //0.01 0.002    // Cглаживание в головной части
        this->sglag_TS_k = 0.001;            // Сглаживание на высоких широтах
        this->sglag_TS_k_sphere_head = 0.05; // 0.08;   // Сглаживание в головной части
        this->sglag_TS_k_sphere_tail = 0.005; // 0.03;   // Сглаживание в хвостовой части



        this->sglag_HP = true;
        this->velocity_HP = 0.1;  // 0.1
        this->sglag_HP_k_sphere = 0.01;  //0.005     // Cглаживание в головной части
        this->sglag_HP_k = 0.03; // 0.01         // Сглаживание не в головной области
        this->sglag_HP_k_angle = 0.001; // 0.01         // Сглаживание не в головной области
        this->sglag_HP_angle = 1.8;    // 1.2 коэффициент усилинея сглаживания по углу
        this->sglag_HP_along = 1.0;    // коэффициент усилинея сглаживания вдоль х
        this->sglag_HP_sphere = 5.0;   // коэффициент усиления сглаживания в головной области - НЕ АКТИВНО


        this->sglag_BS = false;
        this->sglag_BS_k = 0.05;
        this->velocity_BS = 0.3;  // 0.1


        this->null_bn_on_HP = false;   // Для ячеек рядом с HP обнуляем нормальную компоненту магнитного поля
        this->bn_in_p_on_HP = false;   // Для ячеек рядом с HP записываем магнитное поле в давление и решаем Годунова
        this->contact_hard = true;
        this->TS_hard = false;

        // Парметры настройки MK
        this->save_AMR = true;        // Нужно ли сохранять посчитанные функции распределения?
        this->culc_AMR = true;        // Нужно ли считать функции распределения?
        this->refine_AMR = true;      // Нужно ли мельчить посчитанные функции распределения?
        this->N_per_gran = 40000;  // Сколько в среднем частиц вылетает с каждой грани
        this->culc_cell_moments = true;    // Нужно ли считать моменты в ячейках?
        this->de_refine_AMR = true;        // Нужно ли огрублять AMR сетку, если требуется?
        this->MK_file = "parameters_MK_0002.bin";
    }
    else
    {
        // Считываем параметры с файла
        std::ifstream file("set_parameters.txt");

        if (!file.is_open()) {
            cerr << "Error: Could not open file!   set_parameters.txt" << endl;
            exit(-1);
        }

        std::string line;
        while (std::getline(file, line)) 
        {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            // Пропускаем пустые строки
            if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) {
                continue;
            }

            std::istringstream iss(line);
            std::string varName, valueStr;

            // Извлекаем имя и значение (с поддержкой кавычек)
            if ((iss >> varName >> std::ws) && std::getline(iss >> std::ws, valueStr)) {
                // Удаляем пробелы в начале/конце значения
                valueStr.erase(0, valueStr.find_first_not_of(" \t"));
                valueStr.erase(valueStr.find_last_not_of(" \t") + 1);

                if (varMap.count(varName)) {
                    parseAndAssign(varMap[varName], valueStr);
                }
                else {
                    std::cerr << "Unknown variable: " << varName << "\n";
                }
            }
            else {
                std::cerr << "Invalid line format: " << line << "\n";
            }
        }

        file.close();

        /*whach(this->is_div_V_in_cell);
        whach(this->MK_file);
        whach(this->refine_AMR);
        whach(this->sglag_HP_angle);*/
    }

}


void Phys_param::Plasma_components_1(const short int& zone, 
    unordered_map<string, double>& param_in_cell,
    unordered_map<string, double>& param)
{
    // Без пикапов, только протоны и гелий
    // Te == Tth
    // Функция, определяющая температуры и концентрации гелия, 
    // al - это заряд гелия
    // если al = 1 то вне гелиопаузы
    // если al = 2, то внутри гелиопаузы
    short int al;

    double rho = param_in_cell["rho"];
    double p = param_in_cell["p"];
    double rho_He = param_in_cell["rho_He"];

    if (zone <= 2)
    {
        al = 2;
    }
    else
    {
        al = 1;
    }

    param["rho_Th"] = -(MF_meDmp * al * rho_He + 4.0 * (-rho + rho_He))
        / (4.0 * (1.0 + MF_meDmp));


    param["p_Th"] = p * (4.0 * rho - (4.0 + al * MF_meDmp) * rho_He) /
        (8.0 * rho + (-7.0 + al + (1.0 - al) * MF_meDmp) * rho_He);


    param["T_Th"] = 8.0 * (1.0 + MF_meDmp) * p  /
        (8.0 * rho - (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He);

    return;
}

void Phys_param::Plasma_components_2(const short int& zone,
    unordered_map<string, double>& param_in_cell,
    unordered_map<string, double>& param)
{
    // Протоны, гелий, пикапы (два сорта)
    // Te == Tth

    // al - это заряд гелия
    // если al = 1 то вне гелиопаузы
    // если al = 2, то внутри гелиопаузы
    short int al;

    double rho = param_in_cell["rho"];
    double p = param_in_cell["p"];
    double rho_He = param_in_cell["rho_He"];
    double rho_Pui_1, rho_Pui_2;
    double p_Pui_1, p_Pui_2;

    if (zone <= 2)
    {
        al = 2;
    }
    else
    {
        al = 1;
    }


    if (zone == 1 || zone == 3 || zone == 4)
    {
        // Здесь один сорт пикапов
        rho_Pui_1 = param_in_cell["rho_Pui_1"];
        p_Pui_1 = param_in_cell["p_Pui_1"];

        param["rho_Th"] = -(-4.0 * rho + 4.0 * rho_He + al * MF_meDmp * rho_He + 
            (4.0  + 4.0 * MF_meDmp) * rho_Pui_1)
            / (4.0 * (1.0 + MF_meDmp));

        param["T_Th"] = 8.0 * (1.0 + MF_meDmp) * (p - p_Pui_1) /
            (8.0 * rho - (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He - 
                4.0 * (1.0 + MF_meDmp) * rho_Pui_1);
    }
    else
    {
        // Здесь два сорта пикапов
        rho_Pui_1 = param_in_cell["rho_Pui_1"];
        p_Pui_1 = param_in_cell["p_Pui_1"];

        rho_Pui_2 = param_in_cell["rho_Pui_2"];
        p_Pui_2 = param_in_cell["p_Pui_2"];

        param["rho_Th"] = -(-4.0 * rho + 4.0 * rho_He + al * MF_meDmp * rho_He +
            (4.0 + 4.0 * MF_meDmp) * (rho_Pui_1 + rho_Pui_2))
            / (4.0 * (1.0 + MF_meDmp));

        param["T_Th"] = 8.0 * (1.0 + MF_meDmp) * (p - p_Pui_1 - p_Pui_2) /
            (8.0 * rho - (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He -
                4.0 * (1.0 + MF_meDmp) * (rho_Pui_1 + rho_Pui_2));
    }

    return;
}

double Phys_param::Get_rho_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->n_p_cm3;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1));
}

double Phys_param::Get_T_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->T_K;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1));
}

double Phys_param::Get_v_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->V_kms;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) 
    {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1));

    return 0.0;
}

void Phys_param::Get_Potok(const double& rho, const double& p, const double& u, const double& v, 
    const double& w, const double& bx, const double& by, const double& bz, 
    const double& n1, const double& n2, const double& n3, std::vector<double>& Potok)
{
    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d m;

    n << n1, n2, n3;

    if (Potok.size() < 8 || std::abs(n.norm() - 1.0) > 1e-6)
    {
        cout << "Error 8765086751   " << n.norm() << "  " << Potok.size() << endl;
        exit(-1);
    }

    get_bazis(n, t, m);


    double vn = scalarProductFast(u, v, w, n1, n2, n3);
    double vt = scalarProductFast(u, v, w, t(0), t(1), t(2));
    double vm = scalarProductFast(u, v, w, m(0), m(1), m(2));

    double Bn = scalarProductFast(bx, by, bz, n1, n2, n3) / spi4;
    double Bt = scalarProductFast(bx, by, bz, t(0), t(1), t(2)) / spi4;
    double Bm = scalarProductFast(bx, by, bz, m(0), m(1), m(2)) / spi4;

    double bb = kvv(bx, by, bz) / spi4 / spi4;
    double vv = kvv(u, v, w);

    double e = this->energy(rho, p, vv, bb);
    double pT = p + bb / 2.0;

    Potok[0] = rho * vn; 
    Potok[7] = (e + pT) * vn - Bn * scalarProductFast(bx, by, bz, u, v, w) / spi4;
    double Pvn = rho * kv(vn) + pT - kv(Bn);
    double Pvt = rho * vn * vt - Bn * Bt;
    double Pvm = rho * vn * vm - Bn * Bm;
    double Pbt = Bt * vn - Bn * vt;
    double Pbm = Bm * vn - Bn * vm;
    Pbt *= spi4;
    Pbm *= spi4;
    Potok[1] = Pvn * n(0) + Pvt * t(0) + Pvm * m(0);
    Potok[2] = Pvn * n(1) + Pvt * t(1) + Pvm * m(1);
    Potok[3] = Pvn * n(2) + Pvt * t(2) + Pvm * m(2);
    Potok[4] = Pbt * t(0) + Pbm * m(0);
    Potok[5] = Pbt * t(1) + Pbm * m(1);
    Potok[6] = Pbt * t(2) + Pbm * m(2);
}

double Phys_param::energy(const double& rho, const double& p, const double& vv, const double& bb)
{
    return p / (this->gamma - 1) + rho * vv / 2.0 + bb / 2.0;
}


void Phys_param::chlld(unsigned short int n_state, // метод
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
    PrintOptions& opts, bool left_ydar, bool contact)  // Дополнительные опциональные параметры
    // n_state = 0 Лакс, // 1 HLL, // 2 HLLC,  3 HLLD
    // Конвективные переменные добавил только для HLL
{

    double FL[8];
    double FR[8];
    double UL[8];
    double UR[8];
    double UZ[8];

    vector<double> FLk(konvect_left.size());
    vector<double> FRk(konvect_left.size());
    vector<double> ULk(konvect_left.size());
    vector<double> URk(konvect_left.size());
    vector<double> UZk(konvect_left.size());

    double wv = w;

    double r1 = qqq1[0];
    double u1 = qqq1[1];
    double v1 = qqq1[2];
    double w1 = qqq1[3];
    double p1 = qqq1[4];
    double bx1 = qqq1[5] / spi4;
    double by1 = qqq1[6] / spi4;
    double bz1 = qqq1[7] / spi4;

    double r2 = qqq2[0];
    double u2 = qqq2[1];
    double v2 = qqq2[2];
    double w2 = qqq2[3];
    double p2 = qqq2[4];
    double bx2 = qqq2[5] / spi4;
    double by2 = qqq2[6] / spi4;
    double bz2 = qqq2[7] / spi4;

    double ro = (r2 + r1) / 2.0;
    double au = (u2 + u1) / 2.0;
    double av = (v2 + v1) / 2.0;
    double aw = (w2 + w1) / 2.0;
    double ap = (p2 + p1) / 2.0;
    double abx = (bx2 + bx1) / 2.0;
    double aby = (by2 + by1) / 2.0;
    double abz = (bz2 + bz1) / 2.0;

    double bk = abx * al + aby * be + abz * ge;
    double b2 = kv(abx) + kv(aby) + kv(abz);

    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d m;

    n << al, be, ge;

    get_bazis(n, t, m);

    //t << 0.0, 1.0, 0.0;
    //m << 0.0, 0.0, 1.0;


    Eigen::Vector3d vL;
    Eigen::Vector3d vR;
    Eigen::Vector3d bL;
    Eigen::Vector3d bR;


    vL(0) = scalarProductFast(n(0), n(1), n(2), u1, v1, w1);
    vL(1) = scalarProductFast(t(0), t(1), t(2), u1, v1, w1);
    vL(2) = scalarProductFast(m(0), m(1), m(2), u1, v1, w1);

    vR(0) = scalarProductFast(n(0), n(1), n(2), u2, v2, w2);
    vR(1) = scalarProductFast(t(0), t(1), t(2), u2, v2, w2);
    vR(2) = scalarProductFast(m(0), m(1), m(2), u2, v2, w2);

    bL(0) = scalarProductFast(n(0), n(1), n(2), bx1, by1, bz1);
    bL(1) = scalarProductFast(t(0), t(1), t(2), bx1, by1, bz1);
    bL(2) = scalarProductFast(m(0), m(1), m(2), bx1, by1, bz1);


    bR(0) = scalarProductFast(n(0), n(1), n(2), bx2, by2, bz2);
    bR(1) = scalarProductFast(t(0), t(1), t(2), bx2, by2, bz2);
    bR(2) = scalarProductFast(m(0), m(1), m(2), bx2, by2, bz2);


    double aaL = bL(0) / sqrt(r1);
    double b2L = bL.dot(bL);
    double b21 = b2L / r1;
    double cL = sqrt(this->gamma * p1 / r1);
    double qp = sqrt(b21 + cL * (cL + 2.0 * aaL));
    double qm = sqrt(b21 + cL * (cL - 2.0 * aaL));
    double cfL = (qp + qm) / 2.0;
    double ptL = p1 + b2L / 2.0;

    double aaR = bR(0) / sqrt(r2);
    double b2R = bR.dot(bR);
    double b22 = b2R / r2;
    double cR = sqrt(this->gamma * p2 / r2);
    qp = sqrt(b22 + cR * (cR + 2.0 * aaR));
    qm = sqrt(b22 + cR * (cR - 2.0 * aaR));
    double cfR = (qp + qm) / 2.0;
    double ptR = p2 + b2R / 2.0;

    double aC = (aaL + aaR) / 2.0;
    double b2o = (b22 + b21) / 2.0;
    double cC = sqrt(this->gamma * ap / ro);
    qp = sqrt(b2o + cC * (cC + 2.0 * aC));
    qm = sqrt(b2o + cC * (cC - 2.0 * aC));
    double cfC = (qp + qm) / 2.0;
    double vC1 = (vL(0) + vR(0)) / 2.0;

    double SL = 0.0, SR = 0.0;

    if (left_ydar == true)  n_disc = 2;

    if (n_disc == 0)
    {
        SL = min((vL(0) - cfL), (vR(0) - cfR));
        SR = max((vL(0) + cfL), (vR(0) + cfR));
    }
    else if (n_disc == 1)
    {
        SL = min((vL(0) - cfL), (vC1 - cfC));
        SR = max((vR(0) + cfR), (vC1 + cfC));
    }
    else if (n_disc == 2)
    {
        double SL_1 = min((vL(0) - cfL), (vC1 - cfC));
        double SR_1 = max((vR(0) + cfR), (vC1 + cfC));
        double SL_2 = min((vL(0) - cfL), (vR(0) - cfR));
        double SR_2 = max((vL(0) + cfL), (vR(0) + cfR));
        double oo = 0.6;
        double oo1 = 1.0 - oo;
        SL = oo * SL_1 + oo1 * SL_2;
        SR = oo * SR_1 + oo1 * SR_2;
    }
    else if (n_disc == 3)
    {
        SL = min(vL[0], vR[0]) - max(cfL, cfR);
        SR = max(vL[0], vR[0]) + max(cfL, cfR);
    }
    else
    {
        cout << "Error 8654354789" << endl;
    }

    dsr = SR;
    dsl = SL;

    double suR = SR - vR(0);
    double suL = SL - vL(0);
    double SM = (suR * r2 * vR(0) - ptR + ptL - suL * r1 * vL(0))
        / (suR * r2 - suL * r1);
    dsc = SM;

    if (contact == true)
    {
        wv = SM;
    }

    if (n_state == 0)
    {
        double TR0 = fabs(vL(0) + vR(0)) / 2.0 + cfC;
        double TL0 = -TR0;
        SR = TR0;
        SL = TL0;
    }


    double upt1 = (kv(u1) + kv(v1) + kv(w1)) / 2.0;
    double sbv1 = u1 * bx1 + v1 * by1 + w1 * bz1;

    double upt2 = (kv(u2) + kv(v2) + kv(w2)) / 2.0;
    double sbv2 = u2 * bx2 + v2 * by2 + w2 * bz2;

    double e1 = p1 / g1 + r1 * upt1 + b2L / 2.0;
    double e2 = p2 / g1 + r2 * upt2 + b2R / 2.0;

    FL[0] = r1 * vL[0];
    FL[1] = r1 * vL[0] * vL[0] + ptL - kv(bL[0]);
    FL[2] = r1 * vL[0] * vL[1] - bL[0] * bL[1];
    FL[3] = r1 * vL[0] * vL[2] - bL[0] * bL[2];
    FL[4] = (e1 + ptL) * vL[0] - bL[0] * sbv1;
    FL[5] = 0.0;
    FL[6] = vL[0] * bL[1] - vL[1] * bL[0];
    FL[7] = vL[0] * bL[2] - vL[2] * bL[0];
    for (short int i = 0; i < konvect_left.size(); i++)
    {
        FLk[i] = konvect_left[i] * vL[0];
    }

    FR[0] = r2 * vR[0];
    FR[1] = r2 * vR[0] * vR[0] + ptR - kv(bR[0]);
    FR[2] = r2 * vR[0] * vR[1] - bR[0] * bR[1];
    FR[3] = r2 * vR[0] * vR[2] - bR[0] * bR[2];
    FR[4] = (e2 + ptR) * vR[0] - bR[0] * sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0] * bR[1] - vR[1] * bR[0];
    FR[7] = vR[0] * bR[2] - vR[2] * bR[0];   // ***
    for (short int i = 0; i < konvect_right.size(); i++)
    {
        FRk[i] = konvect_right[i] * vR[0];
    }

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;

    for (short int i = 0; i < konvect_right.size(); i++)
    {
        ULk[i] = konvect_left[i];
        URk[i] = konvect_right[i];
    }



    for (short unsigned int i = 0; i < 3; i++)
    {
        UL[i + 1] = r1 * vL[i];
        UL[i + 5] = bL[i];
        UR[i + 1] = r2 * vR[i];
        UR[i + 5] = bR[i];
    }


    for (short unsigned int ik = 0; ik < 8; ik++)
    {
        UZ[ik] = (SR * UR[ik] - SL * UL[ik] + FL[ik] - FR[ik]) / (SR - SL);
    }

    for (short int ik = 0; ik < konvect_right.size(); ik++)
    {
        UZk[ik] = (SR * URk[ik] - SL * ULk[ik] + FLk[ik] - FRk[ik]) / (SR - SL);
    }


    if (null_bn == true) UZ[5] = 0.0;

   
    if (n_state <= 1)  // Лакс или HLL
    {
        double dq[8];
        double FW[8];

        vector<double> dqk(konvect_left.size());
        vector<double> FWk(konvect_left.size());

        for (short unsigned int ik = 0; ik < 8; ik++)
        {
            dq[ik] = UR[ik] - UL[ik];
        }
        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            dqk[ik] = URk[ik] - ULk[ik];
        }

        if (left_ydar == true)
        {
            wv = SL;
        }

        double TL = SL;
        double TR = SR;
        if (SL > wv || left_ydar == true)
        {
            TL = 0.0;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UL[ik];
            }

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * ULk[ik];
            }
        }
        else if (SL <= wv && wv <= SR)
        {
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UZ[ik];
            }

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * UZk[ik];
            }
        }
        else if (SR < wv)
        {
            TR = 0.0;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UR[ik];
            }
            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * URk[ik];
            }
        }
        else
        {
            cout << "error  8767854390  " << endl;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                cout << "ik: " << qqq1[ik] << "   " << qqq2[ik] << ";" << endl;
            }
            if (opts.x) whach(*opts.x);
            if (opts.y) whach(*opts.y);
            if (opts.z) whach(*opts.z);
            if (opts.fluid) whach(*opts.fluid);
            exit(-1);
        }

        double a = TR * TL;
        double b = TR - TL;

        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        qqq[0] = (TR * FL[0] - TL * FR[0] + a * dq[0]) / b - FW[0];
        qqq[4] = (TR * FL[4] - TL * FR[4] + a * dq[4]) / b - FW[4];
        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            konvect[ik] = (TR * FLk[ik] - TL * FRk[ik] + a * dqk[ik]) / b - FWk[ik];
        }


        for (short unsigned int ik = 1; ik < 4; ik++)
        {
            qv(ik - 1) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }
        for (short unsigned int ik = 5; ik < 8; ik++)
        {
            qb(ik - 5) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }

        double SN = max(fabs(SL), fabs(SR));
        double wbn = 0.0;
        if (wv >= SR)
        {
            wbn = wv * bR[0];
        }
        else if (wv <= SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / 2.0;
        }

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;
        if (null_bn == true) qb[0] = 0.0;

        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
            qqq[ik + 5] = spi4 * qqq[ik + 5];
        }

        return;
    }
    if (n_state == 2) // HLLC
    {
        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = r2 * suRm;
        double rzL = r1 * suLm;

        Eigen::Vector3d vzR, vzL, bzR, bzL;

        vzR[0] = SM;
        vzL[0] = SM;

        double ptz = (suR * r2 * ptL - suL * r1 * ptR +
            r1 * r2 * suR * suL * (vR[0] - vL[0]))
            / (suR * r2 - suL * r1);

        bzR[0] = UZ[5];   // bR[0]
        bzL[0] = UZ[5];   // bL[0]

        vzR[1] = UZ[2] / UZ[0];
        vzR[2] = UZ[3] / UZ[0];
        vzL[1] = vzR[1];
        vzL[2] = vzR[2];

        vzR[1] = vR[1] + UZ[5] * (bR[1] - UZ[6]) / suR / r2;
        vzR[2] = vR[2] + UZ[5] * (bR[2] - UZ[7]) / suR / r2;
        vzL[1] = vL[1] + UZ[5] * (bL[1] - UZ[6]) / suL / r1;
        vzL[2] = vL[2] + UZ[5] * (bL[2] - UZ[7]) / suL / r1;

        bzR[1] = UZ[6];
        bzR[2] = UZ[7];
        bzL[1] = bzR[1];
        bzL[2] = bzR[2];

        double sbvz = (UZ[5] * UZ[1] + UZ[6] * UZ[2] + UZ[7] * UZ[3]) / UZ[0];


        double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + UZ[5] * (sbv2 - sbvz)) / (SR - SM);
        double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + UZ[5] * (sbv1 - sbvz)) / (SL - SM);


        if (fabs(UZ[5]) <= epsb)
        {
            vzR[1] = vR[1];
            vzR[2] = vR[2];
            vzL[1] = vL[1];
            vzL[2] = vL[2];
            bzR[1] = bR[1] * suRm;
            bzR[2] = bR[2] * suRm;
            bzL[1] = bL[1] * suLm;
            bzL[2] = bL[2] * suLm;
        }

        // Сгладим Скорость     Я добавил, не было изначально !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // vzR[1] = UZ[2] / UZ[0]
        // vzR[2] = UZ[3] / UZ[0]
        // vzL[1] = vzR[1]
        // vzL[2] = vzR[2]

        double UZL[8];
        double UZR[8];
        vector<double> UZLk(konvect_left.size());
        vector<double> UZRk(konvect_left.size());

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            UZLk[ik] = konvect_left[ik] * suLm;
            UZRk[ik] = konvect_right[ik] * suRm;
        }

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;

        for (size_t ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik] * rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik] * rzR;
            UZR[ik + 5] = bzR[ik];
        }

        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        vector<double> FWk(konvect_left.size());

        if (left_ydar == true)
        {
            wv = SL;
        }

        if (SL >= wv || left_ydar == true)
        {
            qqq[0] = FL[0] - wv * UL[0];
            qqq[4] = FL[4] - wv * UL[4];

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FLk[ik] - wv * ULk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv * UL[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv * UL[ik];
            }
        }
        else if (SL <= wv && SM >= wv)
        {
            qqq[0] = FL[0] + SL * (rzL - r1) - wv * UZL[0];
            qqq[4] = FL[4] + SL * (ezL - e1) - wv * UZL[4];

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FLk[ik] + SL * (UZLk[ik] - konvect_left[ik]) - wv * UZLk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
        }
        else if (SM <= wv && SR > wv)
        {
            qqq[0] = FR[0] + SR * (rzR - r2) - wv * UZR[0];
            qqq[4] = FR[4] + SR * (ezR - e2) - wv * UZR[4];

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FRk[ik] + SR * (UZRk[ik] - konvect_right[ik]) - wv * UZRk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
        }
        else if (SR <= wv)
        {
            qqq[0] = FR[0] - wv * UR[0];
            qqq[4] = FR[4] - wv * UR[4];

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FRk[ik] - wv * URk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv * UR[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv * UR[ik];
            }
        }
        else
        {
            cout << "error  3423167452  " << endl;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                cout << "ik: " << qqq1[ik] << "   " << qqq2[ik] << ";" << endl;
            }
            whach(SR);
            whach(SL);
            whach(SM);
            whach(wv);
            exit(-1);
        }

        // Bn
        double SN = max(fabs(SL), fabs(SR));

        double wbn = 0.0;
        if (wv > SR)
        {
            wbn = wv * bR[0];
        }
        else if(wv < SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / 2.0;
        }

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;

        if (null_bn == true) qb[0] = 0.0;


        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
            qqq[ik + 5] = spi4 * qqq[ik + 5];
        }

        return;
    }
    if (n_state == 3) // HLLD
    {
        Eigen::Vector3d vzR, vzL, vzzL, vzzR, bzR, bzL;

        double ptz = (suR * r2 * ptL - suL * r1 * ptR +
            r1 * r2 * suR * suL * (vR[0] - vL[0]))
            / (suR * r2 - suL * r1);

        vzL[0] = SM;
        vzR[0] = SM;
        vzzL[0] = SM;
        vzzR[0] = SM;
        double ptzL = ptz;
        double ptzR = ptz;
        double ptzzL = ptz;
        double ptzzR = ptz;

        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = r2 * suRm;
        double rzL = r1 * suLm;

        vector<double> qzL(konvect_left.size());
        vector<double> qzR(konvect_left.size());

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            qzL[ik] = konvect_left[ik] * suLm;
            qzR[ik] = konvect_right[ik] * suRm;
        }

        Eigen::Vector3d bzzL, bzzR;

        double bn = UZ[5];
        double bn2 = bn * bn;
        bzL[0] = bn;
        bzR[0] = bn;
        bzzL[0] = bn;
        bzzR[0] = bn;

        double tvR, tbR;
        double ttR = r2 * suR * (SR - SM) - bn2;
        if (fabs(ttR) < 0.000000001)
        {
            tvR = 0.0;
            tbR = 0.0;
        }
        else
        {
            tvR = (SM - vR[0]) / ttR;
            tbR = (r2 * suR * suR - bn2) / ttR;
        }

        double tvL, tbL;
        double ttL = r1 * suL * (SL - SM) - bn2;
        if (fabs(ttL) < 0.000000001)
        {
            tvL = 0.0;
            tbL = 0.0;
        }
        else
        {
            tvL = (SM - vL[0]) / ttL;
            tbL = (r1 * suL * suL - bn2) / ttL;
        }

        vzL[1] = vL[1] - bn * bL[1] * tvL;
        vzL[2] = vL[2] - bn * bL[2] * tvL;
        vzR[1] = vR[1] - bn * bR[1] * tvR;
        vzR[2] = vR[2] - bn * bR[2] * tvR;

        bzL[1] = bL[1] * tbL;
        bzL[2] = bL[2] * tbL;
        bzR[1] = bR[1] * tbR;
        bzR[2] = bR[2] * tbR;

        double sbvL = bzL[0] * vzL[0] + bzL[1] * vzL[1] + bzL[2] * vzL[2];
        double sbvR = bzR[0] * vzR[0] + bzR[1] * vzR[1] + bzR[2] * vzR[2];

        double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + bn * (sbv2 - sbvR)) / (SR - SM);
        double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + bn * (sbv1 - sbvL)) / (SL - SM);


        vector<double> qzzL(konvect_left.size());
        vector<double> qzzR(konvect_left.size());

        vector<double> qzLs(konvect_left.size());
        vector<double> qzRs(konvect_left.size());

        double rzzR = rzR;
        double rzzL = rzL;
        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            qzzR[ik] = qzR[ik];
            qzzL[ik] = qzL[ik];
        }

        double rzRs = sqrt(rzR);
        double rzLs = sqrt(rzL);

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            qzRs[ik] = sqrt(qzR[ik]);
            qzLs[ik] = sqrt(qzL[ik]);
        }

        double rzss = rzRs + rzLs;
        double rzps = rzRs * rzLs;

        vector<double> qzss(konvect_left.size());
        vector<double> qzps(konvect_left.size());

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            qzss[ik] = qzRs[ik] + qzLs[ik];
            qzps[ik] = qzRs[ik] * qzLs[ik];
        }

        double SZL = SM - fabs(bn) / rzLs;
        double SZR = SM + fabs(bn) / rzRs;


        double sbn;
        short int ibn = 0;
        if (fabs(bn) > epsb)
        {
            sbn = fabs(bn) / bn;
            ibn = 1;
        }
        else
        {
            sbn = 0.0;
            ibn = 0;
            SZL = SM;
            SZR = SM;
        }

        vzzL[1] = (rzLs * vzL[1] + rzRs * vzR[1]
            + sbn * (bzR[1] - bzL[1])) / rzss;
        vzzL[2] = (rzLs * vzL[2] + rzRs * vzR[2]
            + sbn * (bzR[2] - bzL[2])) / rzss;
        vzzR[1] = vzzL[1];
        vzzR[2] = vzzL[2];

        bzzL[1] = (rzLs * bzR[1] + rzRs * bzL[1]
            + sbn * rzps * (vzR[1] - vzL[1])) / rzss;
        bzzL[2] = (rzLs * bzR[2] + rzRs * bzL[2]
            + sbn * rzps * (vzR[2] - vzL[2])) / rzss;
        bzzR[1] = bzzL[1];
        bzzR[2] = bzzL[2];

        double sbzz = bzzL[0] * vzzL[0] + bzzL[1] * vzzL[1] + bzzL[2] * vzzL[2];

        double ezzR = ezR + rzRs * sbn * (sbvR - sbzz);
        double ezzL = ezL - rzLs * sbn * (sbvL - sbzz);

        double UZL[8];
        double UZR[8];
        vector<double> UZLk(konvect_left.size());
        vector<double> UZRk(konvect_left.size());

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            UZLk[ik] = qzL[ik];
            UZRk[ik] = qzR[ik];
        }

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;
        for (size_t ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik] * rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik] * rzR;
            UZR[ik + 5] = bzR[ik];
        }

        double UZZL[8];
        double UZZR[8];
        vector<double> UZZLk(konvect_left.size());
        vector<double> UZZRk(konvect_left.size());

        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            UZZLk[ik] = qzzL[ik];
            UZZRk[ik] = qzzR[ik];
        }

        UZZL[0] = rzzL;
        UZZL[4] = ezzL;
        UZZR[0] = rzzR;
        UZZR[4] = ezzR;
        for (size_t ik = 0; ik < 3; ik++)
        {
            UZZL[ik + 1] = vzzL[ik] * rzzL;
            UZZL[ik + 5] = bzzL[ik];
            UZZR[ik + 1] = vzzR[ik] * rzzR;
            UZZR[ik + 5] = bzzR[ik];
        }

        short int j_ccs = -1;
        short int ik;
        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        if (left_ydar == true)
        {
            wv = SL;
        }

        if (SL >= wv || left_ydar == true)
        {
            if(zone_info) cout << "Zone 1" << endl;
            qqq[0] = FL[0] - wv * UL[0];
            qqq[4] = FL[4] - wv * UL[4];
            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv * UL[ik];
            }
            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv * UL[ik];
            }
            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FLk[ik] - wv * ULk[ik];
            }
            j_ccs = 1;
        }
        else if (SL <= wv && SZL >= wv)
        {
            if (zone_info) cout << "Zone 2" << endl;
            ik = 0;
            qqq[ik] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            ik = 4;
            qqq[ik] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];


            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FLk [ik] + SL * (UZLk[ik] - ULk[ik]) - wv * UZLk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
            j_ccs = 2;
        }

        if (ibn == 1)
        {
            if (SZL <= wv && SM >= wv)
            {
                if (zone_info) cout << "Zone 3" << endl;
                ik = 0;
                qqq[ik] = FL[ik] + SZL * (UZZL[ik] - UZL[ik])
                    + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                ik = 4;
                qqq[ik] = FL[ik] + SZL * (UZZL[ik] - UZL[ik])
                    + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];

                for (short int ik = 0; ik < konvect_right.size(); ik++)
                {
                    konvect[ik] = FLk[ik] + SZL * (UZZLk[ik] - UZLk[ik])
                        + SL * (UZLk[ik] - ULk[ik]) - wv * UZZLk[ik];
                }

                for (size_t ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FL[ik] + SZL * (UZZL[ik] - UZL[ik])
                        + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                }
                for (size_t ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FL[ik] + SZL * (UZZL[ik] - UZL[ik])
                        + SL * (UZL[ik] - UL[ik]) - wv * UZZL[ik];
                }
                j_ccs = 3;
            }
            else if (SM <= wv && SZR >= wv)
            {
                if (zone_info) cout << "Zone 4" << endl;
                ik = 0;
                qqq[ik] = FR[ik] + SZR * (UZZR[ik] - UZR[ik])
                    + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                ik = 4;
                qqq[ik] = FR[ik] + SZR * (UZZR[ik] - UZR[ik])
                    + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                for (short int ik = 0; ik < konvect_right.size(); ik++)
                {
                    konvect[ik] = FRk[ik] + SZR * (UZZRk[ik] - UZRk[ik])
                        + SR * (UZRk[ik] - URk[ik]) - wv * UZZRk[ik];
                }
                for (size_t ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FR[ik] + SZR * (UZZR[ik] - UZR[ik])
                        + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                }
                for (size_t ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FR[ik] + SZR * (UZZR[ik] - UZR[ik])
                        + SR * (UZR[ik] - UR[ik]) - wv * UZZR[ik];
                }
                j_ccs = 4;
            }
        }

        if (SZR <= wv && SR >= wv)
        {
            if (zone_info) cout << "Zone 5" << endl;
            ik = 0;
            qqq[ik] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            ik = 4;
            qqq[ik] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FRk[ik] + SR * (UZRk[ik] - URk[ik]) - wv * UZRk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
            j_ccs = 5;
        }
        else if (SR <= wv)
        {
            if (zone_info) cout << "Zone 6" << endl;
            qqq[0] = FR[0] - wv * UR[0];
            qqq[4] = FR[4] - wv * UR[4];
            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                konvect[ik] = FRk[ik] - wv * URk[ik];
            }

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv * UR[ik];
            }
            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv * UR[ik];
            }
            j_ccs = 6;
        }

        if (j_ccs == -1)
        {
            cout << "ERROR chlld_Q  HLLD solver, nstate=3, wrong choise!!!!  " << j_ccs << endl;
            cout << "w SL SZL SM SZR SR" << endl;
            exit(-1);
        }

        double SN = max(fabs(SL), fabs(SR));

        double wbn = 0.0;
        if (wv > SR)
        {
            wbn = wv * bR[0];
        }
        else if(wv < SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / 2.0;
        }

        // wbn = 0.0; //   !Korolkov

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;
        if (null_bn == true) qb[0] = 0.0;

        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
            qqq[ik + 5] = spi4 * qqq[ik + 5];
        }

        return;
    }
    else
    {
        cout << "Error  8562085512" << endl;
        exit(-1);
    }

    return;
}

void Phys_param::raspad_testing(void)
{
    std::vector<double> qqq, qqq1, qqq2;
    qqq.resize(8);
    qqq1.resize(8);
    qqq2.resize(8);

    std::vector<double> konvect_left, konvect_right, konvect;
    PrintOptions Option = PrintOptions{};

    double dsr, dsc, dsl;

    qqq1[0] = 1.0;

    qqq1[1] = 1.0;
    qqq1[2] = 1.0;
    qqq1[3] = 1.0;

    qqq1[4] = 1.0;

    qqq1[5] = 1.1;
    qqq1[6] = 1.2;
    qqq1[7] = 1.3;


    qqq2[0] = 2.0;

    qqq2[1] = -1.0;
    qqq2[2] = 1.0;
    qqq2[3] = 1.0;

    qqq2[4] = 6.0;

    qqq2[5] = 3.1;
    qqq2[6] = 1.2;
    qqq2[7] = 1.3;

    double w = -4.0;

    konvect_left.push_back(10.0);
    konvect_right.push_back(20.0);
    konvect.push_back(0.0);

    this->chlld(2, 1.0, 0.0, 0.0,
        w, qqq1, qqq2, qqq, false, 3,
        konvect_left, konvect_right, konvect, dsr, dsc, dsl,
        Option);

    for (size_t ff = 0; ff < 8; ff++)
    {
        whach(ff);
        whach(qqq[ff]);
    }

    whach(konvect[0]);
}


void Phys_param::Godunov_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
    const std::vector<double>& n, std::vector<double>& qqq,//
    double& dsl, double& dsp, double& dsc, double w, bool contact)
{
    double al = n[0];
    double be = n[1];
    double ge = n[2];
    double time = 0.0;

    Eigen::Vector3d nn;
    Eigen::Vector3d t;
    Eigen::Vector3d m;

    nn << al, be, ge;

    get_bazis(nn, t, m);

    double al2 = t[0];
    double be2 = t[1];
    double ge2 = t[2];
    double al3 = m[0];
    double be3 = m[1];
    double ge3 = m[2];

    double enI = al * qqq1[1] + be * qqq1[2] + ge * qqq1[3];
    double teI2 = al2 * qqq1[1] + be2 * qqq1[2] + ge2 * qqq1[3];
    double teI3 = al3 * qqq1[1] + be3 * qqq1[2] + ge3 * qqq1[3];
    double enII = al * qqq2[1] + be * qqq2[2] + ge * qqq2[3];
    double teII2 = al2 * qqq2[1] + be2 * qqq2[2] + ge2 * qqq2[3];
    double teII3 = al3 * qqq2[1] + be3 * qqq2[2] + ge3 * qqq2[3];

    double pI = qqq1[4];
    double pII = qqq2[4];
    double rI = qqq1[0];
    double rII = qqq2[0];

    int ipiz = 0;
    if (pI > pII)   // Смена местами величин
    {
        double eno2 = enII;;
        double teo22 = teII2;
        double teo23 = teII3;
        double p2 = pII;
        double r2 = rII;

        double eno1 = enI;
        double teo12 = teI2;
        double teo13 = teI3;
        double p1 = pI;
        double r1 = rI;

        enI = -eno2;
        teI2 = teo22;
        teI3 = teo23;
        pI = p2;
        rI = r2;

        enII = -eno1;
        teII2 = teo12;
        teII3 = teo13;
        pII = p1;
        rII = r1;
        w = -w;
        ipiz = 1;                                                                // ???? Он точно здесь должен быть?
    }

    double cI = 0.0;
    double cII = 0.0;
    if (pI != 0.0)
    {
        cI = sqrt(this->gamma * pI / rI);
    }
    if (pII != 0.0)
    {
        cII = sqrt(this->gamma * pII / rII);
    }

    double a = sqrt(rI * (g2 * pII + g1 * pI) / 2.0);
    double Uud = (pII - pI) / a;
    double Urz = -2.0 * cII / g1 * (1.0 - pow((pI / pII), gm));
    double Uvk = -2.0 * (cII + cI) / g1;
    double Udf = enI - enII;


    int il, ip;
    double p, r, te2, te3, en;

    if (Udf < Uvk)
    {
        il = -1;
        ip = -1;
    }
    else if ((Udf >= Uvk) && (Udf <= Urz))
    {
        p = pI * pow(((Udf - Uvk) / (Urz - Uvk)), (1.0 / gm));
        il = 0;
        ip = 0;
    }
    else if ((Udf > Urz) && (Udf <= Uud))
    {
        devtwo(enI, pI, rI, enII, pII, rII, w, p);
        il = 1;
        ip = 0;
    }
    else if (Udf > Uud)
    {
        newton(enI, pI, rI, enII, pII, rII, w, p);
        il = 1;
        ip = 1;
    }

    //*********TWO SHOCKS**********************************************
    if ((il == 1) && (ip == 1))
    {
        //cout << "TWO SHOCKS" << endl;
        double aI = sqrt(rI * (g2 / 2.0 * p + g1 / 2.0 * pI));
        double aII = sqrt(rII * (g2 / 2.0 * p + g1 / 2.0 * pII));

        double u = (aI * enI + aII * enII + pI - pII) / (aI + aII);
        double dI = enI - aI / rI;
        double dII = enII + aII / rII;
        dsl = dI;
        dsp = dII;
        dsc = u;

        if (contact == true)
        {
            w = u;
        }

        double UU = max(fabs(dsl), fabs(dsp));


        if (w <= dI)
        {
            en = enI;
            p = pI;
            r = rI;
            te2 = teI2;
            te3 = teI3;
        }
        else if ((w > dI) && (w <= u))
        {
            en = u;
            p = p;
            r = rI * aI / (aI - rI * (enI - u));
            te2 = teI2;
            te3 = teI3;
        }
        else if ((w > u) && (w < dII))
        {
            en = u;
            p = p;
            r = rII * aII / (aII + rII * (enII - u));
            te2 = teII2;
            te3 = teII3;
        }
        else if (w >= dII)
        {
            en = enII;
            p = pII;
            r = rII;
            te2 = teII2;
            te3 = teII3;
        }
    }


    //*********LEFT - SHOCK, RIGHT - EXPANSION FAN*******************
    if ((il == 1) && (ip == 0))
    {
        //cout << "LEFT - SHOCK, RIGHT - EXPANSION FAN" << endl;
        double aI = sqrt(rI * (g2 / 2.0 * p + g1 / 2.0 * pI));
        double aII;
        if (fabs(p - pII) < eps)
        {
            aII = rII * cII;
        }
        else
        {
            aII = gm * rII * cII * (1.0 - p / pII) / (1.0 - pow((p / pII), gm));
        }

        double u = (aI * enI + aII * enII + pI - pII) / (aI + aII);
        double dI = enI - aI / rI;
        double dII = enII + cII;
        double ddII = u + cII - g1 * (enII - u) / 2.0;
        dsl = dI;
        dsp = dII;
        dsc = u;

        if (contact == true)
        {
            w = u;
        }

        double UU = max(fabs(dsl), fabs(dsp));
        UU = max(UU, fabs(ddII));

        if (w <= dI)
        {
            en = enI;
            p = pI;
            r = rI;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w > dI) && (w <= u))
        {
            en = u;
            p = p;
            r = rI * aI / (aI - rI * (enI - u));
            te2 = teI2;
            te3 = teI3;
        }
        if ((w > u) && (w <= ddII))
        {
            double ce = cII - g1 / 2.0 * (enII - u);
            en = u;
            p = p;
            r = ga * p / ce / ce;
            te2 = teII2;
            te3 = teII3;
        }
        if ((w > ddII) && (w < dII))
        {
            double ce = -g1 / g2 * (enII - w) + 2.0 / g2 * cII;
            en = w - ce;
            p = pII * pow((ce / cII), (1.0 / gm));
            r = ga * p / ce / ce;
            te2 = teII2;
            te3 = teII3;
        }
        if (w >= dII)
        {
            en = enII;
            p = pII;
            r = rII;
            te2 = teII2;
            te3 = teII3;
        }
    }
    //*********TWO EXPANSION FANS**************************************
    if ((il == 0) && (ip == 0))
    {
        //cout << "TWO EXPANSION FANS" << endl;
        //printf("p = %lf\n", p);
        double aI;
        if (fabs(p - pI) < eps)
        {
            aI = rI * cI;
        }
        else
        {
            aI = gm * rI * cI * (1.0 - p / pI) / (1.0 - pow((p / pI), gm));
        }

        //printf("aI = %lf\n", aI);
        double aII;
        if (fabs(p - pII) < eps)
        {
            aII = rII * cII;
        }
        else
        {
            aII = gm * rII * cII * (1.0 - p / pII) / (1.0 - pow((p / pII), gm));
        }

        //printf("aII = %lf\n", aI);
        double u = (aI * enI + aII * enII + pI - pII) / (aI + aII);
        double dI = enI - cI;
        double ddI = u - cI - g1 * (enI - u) / 2.0;
        double dII = enII + cII;
        double ddII = u + cII - g1 * (enII - u) / 2.0;
        dsl = dI;
        dsp = dII;
        dsc = u;

        //whach(dI);
        //whach(dII);
        //whach(u);

        if (contact == true)
        {
            w = u;
        }

        double UU = max(fabs(dsl), fabs(dsp));
        UU = max(UU, fabs(ddII));
        UU = max(UU, fabs(ddI));


        if (w <= dI)
        {
            en = enI;
            p = pI;
            r = rI;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w > dI) && (w < ddI))
        {
            double ce = g1 / g2 * (enI - w) + 2.0 / g2 * cI;
            en = w + ce;
            p = pI * pow((ce / cI), (1.0 / gm));
            r = ga * p / ce / ce;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w >= ddI) && (w <= u))
        {
            double ce = cI + g1 / 2.0 * (enI - u);
            en = u;
            p = p;
            r = ga * p / ce / ce;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w > u) && (w <= ddII))
        {
            double ce = cII - g1 / 2.0 * (enII - u);
            en = u;
            p = p;
            r = ga * p / ce / ce;
            te2 = teII2;
            te3 = teII3;
        }
        if ((w > ddII) && (w < dII))
        {
            double ce = -g1 / g2 * (enII - w) + 2.0 / g2 * cII;
            en = w - ce;
            p = pII * pow((ce / cII), (1.0 / gm));
            r = ga * p / ce / ce;
            te2 = teII2;
            te3 = teII3;
        }
        if (w >= dII)
        {
            en = enII;
            p = pII;
            r = rII;
            te2 = teII2;
            te3 = teII3;
        }
    }

    //*********VAKUUM ************************************************
    if ((il == -1) && (ip == -1))
    {
        //cout << "VAKUUM" << endl;
        double dI = enI - cI;
        double ddI = enI + 2.0 / gg1 * cI;
        double dII = enII + cII;
        double ddII = enII - 2.0 / gg1 * cII;

        dsl = dI;
        dsp = dII;
        dsc = (dI + dII) / 2.0;

        if (contact == true)
        {
            w = dsc;
        }

        double UU = max(fabs(dsl), fabs(dsp));
        UU = max(UU, fabs(ddII));
        UU = max(UU, fabs(ddI));


        if (w <= dI)
        {
            en = enI;
            p = pI;
            r = rI;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w > dI) && (w < ddI))
        {
            double ce = gg1 / gg2 * (enI - w) + 2.0 / gg2 * cI;
            en = w + ce;
            p = pI * pow((ce / cI), (1.0 / gm));
            r = gga * p / ce / ce;
            te2 = teI2;
            te3 = teI3;
        }
        if ((w >= ddI) && (w <= ddII))
        {
            en = w;
            p = 0.0;
            r = 0.0;
            te2 = 0.0;
            te3 = 0.0;
        }
        if ((w > ddII) && (w < dII))
        {
            double ce = -gg1 / gg2 * (enII - w) + 2.0 / gg2 * cII;
            en = w - ce;
            p = pII * pow((ce / cII), (1.0 / gm));
            r = gga * p / ce / ce;
            te2 = teII2;
            te3 = teII3;
        }
        if (w >= dII)
        {
            en = enII;
            p = pII;
            r = rII;
            te2 = teII2;
            te3 = teII3;
        }
    }


    if (ipiz == 1)
    {
        en = -en;
        double dsl1 = dsl;
        double dsp1 = dsp;
        dsl = -dsp1;
        dsp = -dsl1;
        dsc = -dsc;
        w = -w;
    }

    double uo = al * en + al2 * te2 + al3 * te3;
    double vo = be * en + be2 * te2 + be3 * te3;
    double wo = ge * en + ge2 * te2 + ge3 * te3;


    double eo = p / g1 + 0.5 * r * (uo * uo + vo * vo + wo * wo);
    en = al * uo + be * vo + ge * wo;

    if (contact == true)
    {
        w = en;
    }

    qqq[0] = (r * (en - w));
    qqq[1] = (r * (en - w) * uo + al * p);
    qqq[2] = (r * (en - w) * vo + be * p);
    qqq[3] = (r * (en - w) * wo + ge * p);
    qqq[4] = ((en - w) * eo + en * p);


    return;

}

void Phys_param::newton(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
    const double& w, double& p)
{
    double fI, fIs, fII, fIIs;
    double cI = sqrt(ga * pI / rI);
    double cII = sqrt(ga * pII / rII);
    double pn = pI * rII * cII + pII * rI * cI + (enI - enII) * rI * cI * rII * cII;
    pn = pn / (rI * cI + rII * cII);

    double pee = pn;

    int kiter = 0;
a1:
    p = pn;
    if (p <= 0.0)
    {
        ER_S;
        cout << "negative pressure, newton" << endl;
        exit(-1);
    }

    kiter = kiter + 1;

    fI = (p - pI) / (rI * cI * sqrt(gp * p / pI + gm));
    fIs = (ga + 1.0) * p / pI + (3.0 * ga - 1.0);
    fIs = fIs / (4.0 * ga * rI * cI * pow((gp * p / pI + gm), (3.0 / 2.0)));

    fII = (p - pII) / (rII * cII * sqrt(gp * p / pII + gm));
    fIIs = (ga + 1.0) * p / pII + (3.0 * ga - 1.0);
    fIIs = fIIs / (4.0 * ga * rII * cII * pow((gp * p / pII + gm), (3.0 / 2.0)));


    if (kiter == 1100)
    {
        ER_S;
        cout << "zaciklilsya v raspade,i,j,k,KOBL,kdir:" << endl;
        watch(enI);
        watch(pI);
        watch(rI);
        watch(enII);
        watch(pII);
        watch(rII);
        exit(-1);
    }

    pn = p - (fI + fII - (enI - enII)) / (fIs + fIIs);

    if (fabs(pn / pee - p / pee) >= eps)
    {
        goto a1;
    }

    p = pn;

    return;
}

void Phys_param::devtwo(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
    const double& w, double& p)
{
    const double epsil = 10e-10;
    double kl, kp, kc, ksi, ksir, um, ksit;
    int kpizd;

    kl = pI;
    kp = pII;


    this->lev(enI, pI, rI, enII, pII, rII, kl, ksi);
    this->lev(enI, pI, rI, enII, pII, rII, kp, ksir);

    if (fabs(ksi) <= epsil)
    {
        um = kl;
        goto a1;
    }

    if (fabs(ksir) <= epsil)
    {
        um = kp;
        goto a1;
    }

    kpizd = 0;

a2:
    kpizd = kpizd + 1;

    if (kpizd == 1100)
    {
        ER_S;
        cout << "zaciklilsya, devtwo.f i,j,k,KOBL,kdir:" << endl;
        cout << "Error 0978453645" << endl;
        watch(enI);
        watch(pI);
        watch(rI);
        watch(enII);
        watch(pII);
        watch(rII);
        exit(-1);
    }


    kc = (kl + kp) / 2.0;

    this->lev(enI, pI, rI, enII, pII, rII, kc, ksit);

    if (fabs(ksit) <= epsil)
    {
        goto a3;
    }

    if ((ksi * ksit) <= 0.0)
    {
        kp = kc;
        ksir = ksit;
    }
    else
    {
        kl = kc;
        ksi = ksit;
    }

    goto a2;

a3:
    um = kc;
a1:

    p = um;

    return;
}

void Phys_param::lev(const double& enI, const double& pI, const double& rI, const double& enII,//
    const double& pII, const double& rII, double& uuu, double& fee)
{
    double cI = sqrt(ga * pI / rI);
    double cII = sqrt(ga * pII / rII);

    double fI = (uuu - pI) / (rI * cI * sqrt(gp * uuu / pI + gm));

    double fII = 2.0 / g1 * cII * (pow((uuu / pII), gm) - 1.0);

    double f1 = fI + fII;
    double f2 = enI - enII;
    fee = f1 - f2;
    return;
}


double Phys_param::Get_razmer(string par)
{
    if (par == "r") return this->perevod_razmer["r"];

    if (par == "rho" || par == "rho_He") return this->perevod_razmer["rho"];

    if (par == "V" || par == "Vx" || par == "Vy"
        || par == "Vz") return this->perevod_razmer["V"];

    if (par == "p") return this->perevod_razmer["p"];

    if (par == "B" || par == "Bx" || par == "By" || 
        par == "Bz") return this->perevod_razmer["B"];

    if (par == "T") return this->perevod_razmer["T"];

    return 1.0;
}
