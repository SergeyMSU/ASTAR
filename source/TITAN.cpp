// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

Eigen::Vector3d shortestVectorToAxis(const Eigen::Vector3d& point,
    const Eigen::Vector3d& axisDirection) {
    // Нормируем направляющий вектор оси
    Eigen::Vector3d unitAxis = axisDirection.normalized();

    // Находим проекцию точки на ось
    double t = point.dot(unitAxis);
    Eigen::Vector3d projection = t * unitAxis;

    // Возвращаем вектор от проекции к точке (перпендикулярный оси)
    return point - projection;
}

int main()
{
    /*Eigen::Vector3d X, V, Omega;
    Eigen::Vector3d eX, eV, eOmega;
    double alpha = const_pi/2;
    double lambda = 0.0;
    X << 1.0, 0.0, 0.0;
    Omega << 0.0, 0.0, 0.9;
    V = X * 3.0;

    eX = rotationMatrixZ(alpha) * rotationMatrixX(lambda) * X;
    cout << "eX = " << eX(0) << " " << eX(1) << " " << eX(2) << endl;
    eOmega = rotationMatrixZ(alpha) * rotationMatrixX(lambda) * Omega;
    eV = rotationMatrixZ(alpha) * rotationMatrixX(lambda) * V - eOmega.cross(shortestVectorToAxis(eX, eOmega));
    cout << "eOmega = " << eOmega(0) << " " << eOmega(1) << " " << eOmega(2) << endl;
    cout << "eV = " << eV(0) << " " << eV(1) << " " << eV(2) << endl;
    cout << "V = " << V(0) << " " << V(1) << " " << V(2) << endl;

    X(2) = 0.0;
    V = rotationMatrixX(-lambda) * rotationMatrixZ(-alpha) * eV + Omega.cross(X);
    cout << "V = " << V(0) << " " << V(1) << " " << V(2) << endl;

    return 0;*/


    cout << "Start Programm" << endl;
    if (false)
    {
        Interpol SS = Interpol("For_intertpolate_11.bin");

        Setka S1 = Setka(false);
        S1.phys_param->T_all = 50.0012;

        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0),
            -20.0, 20.0, -20.0, 20.0, "_XY_", false);
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            -20.0, 20.0, -20.0, 20.0, "_YZ_", false);

        return 0;
    }

    Setka S1 = Setka();
    S1.Calculating_measure(0);
    S1.Calculating_measure(1);

    //S1.Tecplot_print_krug_yzel_in_3D(1);
    S1.Tecplot_print_all_yzel_with_condition();
    //S1.Tecplot_print_all_cell_in_3D();

    S1.Init_boundary_grans(); 
    S1.Init_TVD();

    //S1.Download_cell_parameters("paramet_0001.bin");   // 23
    S1.phys_param->T_all = 0;
    S1.Init_physics();

    //S1.Write_file_for_FCMHD();
    S1.Read_file_for_FCMHD();

    cout << "S1.phys_param->T_all  " << S1.phys_param->T_all << endl;
    //return 0;

    while (S1.phys_param->T_all < 0.0)
    {
        auto start = std::chrono::high_resolution_clock::now();
        S1.Go(false, 2000, 1); // 400   1
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Execution time: " << duration.count() / 1000.0 / 60.0 << " minutes" << std::endl;
        cout << "ALL T  =  " << S1.phys_param->T_all << endl;

        S1.Tecplot_print_2D_sphere(10, "move_Sphere_13", false, true);
        S1.Tecplot_print_2D_sphere(10, "Sphere_13", false, false);

        S1.Tecplot_print_2D_sphere(121, "move_Sphere_75", false, true);
        S1.Tecplot_print_2D_sphere(121, "Sphere_75", false, false);

        S1.Save_cell_parameters("paramet_0002-.bin");
    }

    cout << "TIME = " << S1.phys_param->T_all << endl;

    //S1.Save_cell_parameters("paramet_0002.bin");

    // Для экономии памяти почистим некоторые массивы

    S1.Save_for_interpolate("For_intertpolate_11.bin", true);

    if (true)
    {
        for (auto& i : S1.All_Yzel)
        {
            delete i;
        }
        S1.All_Yzel.clear();
    }

    Interpol SS = Interpol("For_intertpolate_11.bin");

    if (false)
    {
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0),
            -10.0, 10.0, -10.0, 10.0, "_XY_", false);
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            -10.0, 10.0, -10.0, 10.0, "_YZ_", false);

        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0),
            -10.0, 10.0, -10.0, 10.0, "_XY_move_", true);
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            -10.0, 10.0, -10.0, 10.0, "_YZ_move_", true);

        S1.Tecplot_print_spherik(&SS, 1.0, 200, 100, "move_Sphere_1ae", true);
        S1.Tecplot_print_spherik(&SS, 1.0, 200, 100, "Sphere_1ae", false);

        S1.Tecplot_print_spherik(&SS, 5.0, 200, 100, "move_Sphere_5ae", true);
        S1.Tecplot_print_spherik(&SS, 5.0, 200, 100, "Sphere_5ae", false);

        S1.Tecplot_print_spherik(&SS, 9.5, 200, 100, "move_Sphere_9.5ae", true);
        S1.Tecplot_print_spherik(&SS, 9.5, 200, 100, "Sphere_9.5ae", false);
    }
    else
    {
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0),
            -5.0, 5.0, -5.0, 5.0, "_XY_", false);
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            -5.0, 5.0, -5.0, 5.0, "_YZ_", false);

        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0),
            -5.0, 5.0, -5.0, 5.0, "_XY_move_", true);
        S1.Tecplot_print_2D_dekard(&SS, Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            -5.0, 5.0, -5.0, 5.0, "_YZ_move_", true);

        S1.Tecplot_print_spherik(&SS, 1.0, 200, 100, "move_Sphere_1ae", true);
        S1.Tecplot_print_spherik(&SS, 1.0, 200, 100, "Sphere_1ae", false);

        S1.Tecplot_print_spherik(&SS, 2.0, 200, 100, "move_Sphere_2ae", true);
        S1.Tecplot_print_spherik(&SS, 2.0, 200, 100, "Sphere_2ae", false);

        S1.Tecplot_print_spherik(&SS, 4.5, 200, 100, "move_Sphere_4.5ae", true);
        S1.Tecplot_print_spherik(&SS, 4.5, 200, 100, "Sphere_4.5ae", false);
    }

    return 0;

    


   

    //S1.Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "XYZ_2d_(0, 0, 1, 0)_", false, false);
    //S1.Tecplot_print_2D(&SS, 0.0, 1.0, 0.0, -0.00001, "XYZ_2d_(0, 1, 0, 0)_", false, false);

    S1.Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");
    S1.Tecplot_print_2D(&SS, 0.0, 1.0, 0.0, -0.00001, "_2d_(0, 1, 0, 0)_");

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0), "_(1, 0, 0)_", 25.0);

    


    //S1.Download_cell_parameters("parameters_promeg_116.bin");

    

    return 0;

    S1.Read_old_surface("ASurf_Save00591.bin");
    S1.Move_to_surf(S1.Surf1);

    S1.auto_set_luch_geo_parameter(0);
    cout << "A " << endl;
    S1.Calculating_measure(0);
    cout << "B " << endl;
    S1.Calculating_measure(1);
    cout << "B2 " << endl;

    S1.Init_boundary_grans();
    cout << "C " << endl;

    S1.Download_cell_parameters("parameters_promeg_116.bin");   // 23
    // 19 стартовая точка от которой две параллели с пикапами и без
    // 32 с пикапами
    // 62 включи TVD

    S1.geo->R0 = 0.237455;

    // 23 полностью установленное решение без Пикапов (у контакта есть артефакт нужно сглаживание по
    // углу увеличить)

    cout << "C2 " << endl;

    S1.auto_set_luch_geo_parameter(0);

    //S1.Smooth_head_HP2(); // Ручное сглаживыание

    S1.Init_TVD();
    cout << "D2 " << endl;

    S1.Init_physics();

    cout << "E " << endl;

    //S1.Smooth_head_TS2();
    //S1.Smooth_head_HP2();

    S1.Tecplot_print_cell_plane_parameters();
    S1.Tecplot_print_all_lush_in_2D();
    S1.Tecplot_print_all_cell_in_3D();

    //S1.Algoritm(2);
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");

    

    S1.Find_Yzel_Sosed_for_BS();

    S1.Smooth_angle_HP();
    S1.Smooth_head_HP3();
    S1.Smooth_head_TS3();


    for (int i = 1; i <= 4 * 3; i++) // 6 * 2
    {
        auto start = std::chrono::high_resolution_clock::now();
        cout << "IIIII = " << i << endl;
        S1.Go(false, 400, 1); // 400   1
        S1.Smooth_head_HP3();
        S1.Smooth_head_TS3();

        S1.Tecplot_print_cell_plane_parameters();
        S1.Tecplot_print_all_lush_in_2D();
        S1.Tecplot_print_all_gran_in_surface("TS");
        S1.Tecplot_print_all_gran_in_surface("HP");
        S1.Tecplot_print_all_gran_in_surface("BS");
        //S1.Go(true, 100, 1);
        //S1.Tecplot_print_cell_plane_parameters();

        //S1.Init_physics();

        if (i % 6 == 0)
        {
            string namn = "parameters_promeg_11" + to_string(i) + ".bin";
            S1.Save_cell_parameters(namn);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Execution time: " << duration.count()/1000.0/60.0 << " minutes" << std::endl;
    }


    if (false)
    {
        S1.~Setka();
        std::cout << "Setka delete\n";
        return 0;
    }

    S1.Save_cell_parameters("parameters_0077.bin");
    //S1.Save_cell_pui_parameters("parameters_0026.bin");

    //S1.Edges_create();
    //S1.Culc_divergence_in_cell();
    //S1.Culc_rotors_in_cell();

    if (false) // Проверка интерполятора
    {
        std::unordered_map<string, double> parameters;
        Cell_handle next_cell;
        Cell_handle prev_cell = Cell_handle();
        // Открываем файл для записи
        std::ofstream outfile("angles.txt");
        if (!outfile.is_open()) {
            return 1;
        }
        const int N = 500;         // Количество шагов
        const double step = const_pi / N;  // Размер шага

        if (false)
        {
            // Основной цикл
            for (double angle = 0.0; angle <= const_pi / 2.0 + 1e-6; angle += step)
            {
                for (double r = 10.0; r < 300.0; r += 0.01)
                {
                    double x = r * cos(angle);
                    double y = r * sin(angle);
                    double z = 0.0;
                    short int zoon = 0;
                    // Записываем в файл
                    SS.Get_param(x, y, z, parameters, prev_cell, next_cell, zoon);
                    if (zoon == 3)
                    {
                        outfile << x << " " << y << std::endl;
                        break;
                    }
                }
            }

            // Закрываем файл
            outfile.close();
        }

        // Открываем файл для записи
        outfile = std::ofstream("angles2.txt");
        if (!outfile.is_open()) {
            return 1;
        }

        // Основной цикл
        cout << "AA1" << endl;
        for (double angle = 0.0; angle <= const_pi/2; angle += step)
        {
            //cout << "angle = " << angle << endl;
            double x = cos(angle);
            double y = sin(angle);
            double z = 0.0;
            short int zoon = 0;
            // Записываем в файл
            SS.Get_HP(x, y, z, parameters);
            double r = parameters["r"];
            outfile << r * cos(angle) << " " << r * sin(angle) << std::endl;
        }
        cout << "AA3" << endl;

        // Закрываем файл
        outfile.close();

        std::unordered_map<string, double> param;
        SS.Get_TS(13.94, 0.0, 0.0, param);
        for (const auto& [key, value] : param) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "A " << endl;


        outfile = std::ofstream("2D.txt");

        for (double x = -200.0; x < 400.0; x = x + 0.4)
        {
            for (double y = -400.0; y < 400.0; y = y + 0.3)
            {
                bool vv = SS.Get_param(x, y, 0.0, parameters, prev_cell, next_cell);
                if (vv == false) continue;
                outfile << x << " " << y << " " << parameters["rho"] << std::endl;
            }
        }

        outfile.close();

        SS.Get_param(0.00001, -350, 0.0, parameters, prev_cell, next_cell);
        for (const auto& [key, value] : parameters) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "B1 " << endl;

        SS.Get_param(-0.00001, -350, 0.0, parameters, prev_cell, next_cell);
        for (const auto& [key, value] : parameters) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "B2 " << endl;
        exit(-1);
    }


    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0), "_(1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(-1.0, 0.0, 0.0), "_(-1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0), "_(0, 1, 0)_", 500.0);

    S1.Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");


    cout << "F " << endl;

    S1.Tecplot_print_cell_plane_parameters();


    cout << "YSPEX" << endl;


    //S1.Tecplot_print_all_yzel_in_3D("SDK1");
    
    

    //S1.Tecplot_print_krug_yzel_in_3D(1);
    //S1.Tecplot_print_krug_yzel_in_3D(2);

    //S1.Tecplot_print_all_lush_in_2D();
    //S1.Tecplot_print_All_surfase_in_2D();
    //S1.Tecplot_print_plane_lush(0);
    //S1.Tecplot_print_plane_surfase(0);
    //S1.Tecplot_print_all_gran_in_cell();
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");
    //S1.Tecplot_print_all_yzel_with_condition();


    std::cout << "Hello World!\n";

    S1.~Setka();
    std::cout << "Setka delete\n";


    SS.~Interpol();
    std::cout << "Interpol delete\n";


}

