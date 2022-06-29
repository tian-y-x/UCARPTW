//
// Created by 段心童 on 2021/12/21.
//

#include "read_ucarp_tw.h"
#include <fstream>

extern int INSTANCE_NUM;
std::string& trim_tw(std::string& s) {
    if (s.empty()) {
        return s;
    }
    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
}


Instance* data_process_macOS_tw(string dir_path) {
    std::string str;
    std::ifstream in;

    std::string str_const;
    int instance_num = 30;
    Instance* instance_list = new Instance[instance_num];

    for (int i = 0; i < 30; ++i) {
        std::string path = dir_path + std::to_string(i + 1) + ".dat";
        in.open(path, std::ios::in);
        if (!in) {
            std::cout << "open failed!" << std::endl;
            exit(1);
        }

        double** cost;
        vector<Task> task;
        vector<Task> task_rev;
        double** demand;
        int points_number;//节点的数量
        int cost_mul;
        int demand_mul;
        double Q = 0;
        int len = 0; //|UT|
        int task_cnt = 1;

        bool read_non_required=false;

        while (!in.eof()) {
            getline(in, str);
            str = trim_tw(str);

            if (str==""){
                break;
            }

            str += "\r";

            int index = 0;
            std::string line;
            while (str[index] != ' ') {
                line += str[index];
                index++;
            }

            if (line=="COST_MUL"){
                index+=3;
                std::string str_cost;
                while (str[index] != '\r') {
                    str_cost += str[index];
                    index++;
                }
                cost_mul = std::stoi(str_cost);
            }

            if (line=="DEMAND_MUL"){
                index+=3;
                std::string str_demand;
                while (str[index] != '\r') {
                    str_demand += str[index];
                    index++;
                }
                demand_mul = std::stoi(str_demand);
            }

            str_const = "VERTICES";
            if (line == str_const) {
                index += 3;

                std::string temp_vertices;
                while (str[index] != '\r') {
                    temp_vertices += str[index];
                    index++;
                }
                points_number = std::stoi(temp_vertices);
                //定义cost数组
                cost = new double* [points_number];
                for (int j = 0; j < points_number; ++j) {
                    cost[j] = new double[points_number];
                    for (int k = 0; k < points_number; ++k) {
                        if (j != k) cost[j][k] = -1;
                        else cost[j][k] = 0;
                    }
                }
                //定义demand数组
                demand = new double* [points_number];
                for (int j = 0; j < points_number; ++j) {
                    demand[j] = new double[points_number];
                    for (int k = 0; k < points_number; ++k) {
                        demand[j][k] = 0;
                    }
                }
            }

            str_const = "ARISTAS_REQ";
            if (line == str_const) {
                index += 3;
                std::string temp_task_count;
                while (str[index] != '\r') {
                    temp_task_count += str[index];
                    index++;
                }
                len = std::stoi(temp_task_count);
            }

            str_const = "CAPACIDAD";
            if (line == str_const) {
                index += 3;
                std::string temp_capacity_count;
                while (str[index] != '\r') {
                    temp_capacity_count += str[index];
                    index++;
                }
                Q = std::stod(temp_capacity_count);
            }


            str_const = "LISTA_ARISTAS_REQ";
            if (line == str_const) {
                read_non_required=true;
                int x, y, st, et;
                double c, d;
                std::string temp_x, temp_y, temp_c, temp_d, temp_st, temp_et;
                getline(in, str);
                str = trim_tw(str);
                while (str[0] == '(') {
                    str += "\r";
                    int line_index = 2;
                    temp_x = "";
                    temp_y = "";
                    temp_c = "";
                    temp_d = "";
                    temp_st = "";
                    temp_et = "";

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_x += str[line_index];
                        line_index++;
                    }
                    x = std::stoi(temp_x);

                    line_index += 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_y += str[line_index];
                        line_index++;
                    }
                    y = std::stoi(temp_y);

                    line_index += 11;
                    while (str[line_index] != ' ') {
                        temp_c += str[line_index];
                        line_index++;
                    }
                    c = std::stod(temp_c);

                    line_index += 10;

                    while (str[line_index] != ' ') {
                        line_index++;
                    }
                    line_index += 11;

                    while (str[line_index] != ' ') {
                        temp_d += str[line_index];
                        line_index++;
                    }
                    d = std::stod(temp_d);

                    line_index += 9;
                    while (str[line_index] != ' ') {
                        temp_st += str[line_index];
                        line_index++;
                    }
                    st = std::stod(temp_st);


                    line_index += 7;
                    while (str[line_index] != '\r') {
                        temp_et += str[line_index];
                        line_index++;
                    }
                    et = std::stod(temp_et);


                    Task* task1 = new Task(task_cnt, c, d, x, y,st, et);
                    task.push_back(*task1);
                    delete task1;
//                    task_cnt++;
                    Task* task2 = new Task(task_cnt, c, d, y, x, st, et);
                    task_rev.push_back(*task2);
                    delete task2;
                    task_cnt++;


                    cost[x - 1][y - 1] = c;
                    cost[y - 1][x - 1] = c;
                    demand[x - 1][y - 1] = d;
                    demand[y - 1][x - 1] = d;

                    getline(in, str);
                    str = trim_tw(str);
                }
                read_non_required=true;
            }

            str_const = "LISTA_ARISTAS_NOREQ : ";
            if (read_non_required) {
                read_non_required=false;
                getline(in, str);
                str = trim_tw(str);
                int x, y;
                double c;
                std::string temp_x, temp_y, temp_c;
                while (str[0] == '(') {
                    str += "\r";
                    temp_x = "";
                    temp_y = "";
                    temp_c = "";
                    int line_index = 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_x += str[line_index];
                        line_index++;
                    }
                    x = std::stoi(temp_x);

                    line_index += 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_y += str[line_index];
                        line_index++;
                    }
                    y = std::stoi(temp_y);

                    line_index += 10;
                    while (str[line_index] != ' ') {
                        temp_c += str[line_index];
                        line_index++;
                    }
                    c = std::stod(temp_c);


                    cost[x - 1][y - 1] = c;
                    cost[y - 1][x - 1] = c;

                    getline(in, str);
                    str = trim_tw(str);
                }
            }


            //NOMBRE 所属UCARP的文件夹
            //COMENTARIO 不知道
            //VERTICES 节点个数1-77
            //ARISTAS_REQ: 一个CARP中Task的数量
            //ARISTAS_NOREQ : 一个CARP中不是Task的数量（但是联通不包含-1）
            //VEHICULOS : 车的数量
            //CAPACIDAD : 一辆车的容量
            //TIPO_COSTES_ARISTAS : 不知
            //COSTE_TOTAL_REQ : 1395
            //LISTA_ARISTAS_REQ :不知
            //...
            //LISTA_ARISTAS_NOREQ : 不是task的边
            //...
            //DEPOSITO 回收站
        }
        for (auto & j : task_rev) {
            j.task_number += (task_cnt-1);
        }
        for (const auto & j : task_rev) {
            task.push_back(j);
        }
        Instance* instance = new Instance(task, cost, demand, points_number, task_cnt-1, Q,cost_mul,demand_mul);
        instance_list[i] = *instance;

//                cout<<i<<endl;
//        cout<<"VERTICES : "<<instance_list[i].point_number<<endl;
//        cout<<"CAPACIDAD : "<<instance_list[i].capacity<<endl;
//        cout<<"UT : "<<instance_list[i].UT<<endl;
//        Instance cur=instance_list[i];
//        for (int kk = 0; kk < cur.point_number; ++kk) {
//            for (int j = 0; j < cur.point_number; ++j) {
//                cout<<cur.cost[kk][j]<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<endl;
//        for (int kk = 0; kk < cur.UT; ++kk) {
//            cout<<cur.tasklist[kk].demand<<" ";
//        }
//        cout<<endl;

        in.close();
    }

    return instance_list;
}

Instance* data_process_win_tw(string dir_path) {
    std::string str;
    std::ifstream in;

    std::string str_const;
    int instance_num = 30;
    Instance* instance_list = new Instance[instance_num];

    for (int i = 0; i < 30; ++i) {
        std::string path = dir_path + std::to_string(i + 1) + ".dat";
        in.open(path, std::ios::in);
        if (!in) {
            std::cout << "open failed!" << std::endl;
            exit(1);
        }

        double** cost;
        vector<Task> task;
        vector<Task> task_rev;
        double** demand;
        int points_number;//节点的数量
        int cost_mul;
        int demand_mul;
        double Q = 0;
        int len = 0; //|UT|
        int task_cnt = 1;


        while (!in.eof()) {
            getline(in, str);
            str = trim_tw(str);

            if (str==""){
                break;
            }

            str += "\r";

            int index = 0;
            std::string line;
            while (str[index] != ' ') {
                line += str[index];
                index++;
            }

            if (line=="COST_MUL"){
                index+=3;
                std::string str_cost;
                while (str[index] != '\r') {
                    str_cost += str[index];
                    index++;
                }
                cost_mul = std::stoi(str_cost);
            }

            if (line=="DEMAND_MUL"){
                index+=3;
                std::string str_demand;
                while (str[index] != '\r') {
                    str_demand += str[index];
                    index++;
                }
                demand_mul = std::stoi(str_demand);
            }

            str_const = "VERTICES";
            if (line == str_const) {
                index += 3;
                std::string temp_vertices;
                while (str[index] != '\r') {
                    temp_vertices += str[index];
                    index++;
                }
                points_number = std::stoi(temp_vertices);
                //定义cost数组
                cost = new double* [points_number];
                for (int j = 0; j < points_number; ++j) {
                    cost[j] = new double[points_number];
                    for (int k = 0; k < points_number; ++k) {
                        if (j != k) cost[j][k] = -1;
                        else cost[j][k] = 0;
                    }
                }
                //定义demand数组
                demand = new double* [points_number];
                for (int j = 0; j < points_number; ++j) {
                    demand[j] = new double[points_number];
                    for (int k = 0; k < points_number; ++k) {
                        demand[j][k] = 0;
                    }
                }
            }

            str_const = "ARISTAS_REQ";
            if (line == str_const) {
                index += 3;
                std::string temp_task_count;
                while (str[index] != '\r') {
                    temp_task_count += str[index];
                    index++;
                }
                len = std::stoi(temp_task_count);
            }

            str_const = "CAPACIDAD";
            if (line == str_const) {
                index += 3;
                std::string temp_capacity_count;
                while (str[index] != '\r') {
                    temp_capacity_count += str[index];
                    index++;
                }
                Q = std::stod(temp_capacity_count);
            }

            str_const = "LISTA_ARISTAS_REQ";
            if (line == str_const) {
                int x, y, st, et;
                double c, d;
                std::string temp_x, temp_y, temp_c, temp_d, temp_st, temp_et;
                getline(in, str);
                str = trim_tw(str);
                while (str[0] == '(') {
                    str += "\r";
                    int line_index = 2;
                    temp_x = "";
                    temp_y = "";
                    temp_c = "";
                    temp_d = "";
                    temp_st = "";
                    temp_et = "";

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_x += str[line_index];
                        line_index++;
                    }
                    x = std::stoi(temp_x);

                    line_index += 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_y += str[line_index];
                        line_index++;
                    }
                    y = std::stoi(temp_y);

                    line_index += 11;
                    while (str[line_index] != ' ') {
                        temp_c += str[line_index];
                        line_index++;
                    }
                    c = std::stod(temp_c);

                    line_index += 10;

                    while (str[line_index] != ' ') {
                        line_index++;
                    }
                    line_index += 11;

                    while (str[line_index] != ' ') {
                        temp_d += str[line_index];
                        line_index++;
                    }
                    d = std::stod(temp_d);

                    line_index += 9;
                    while (str[line_index] != ' ') {
                        temp_st += str[line_index];
                        line_index++;
                    }
                    st = std::stod(temp_st);


                    line_index += 7;
                    while (str[line_index] != '\r') {
                        temp_et += str[line_index];
                        line_index++;
                    }
                    et = std::stod(temp_et);


                    Task* task1 = new Task(task_cnt, c, d, x, y,st, et);
                    task.push_back(*task1);
                    delete task1;
//                    task_cnt++;
                    Task* task2 = new Task(task_cnt, c, d, y, x, st, et);
                    task_rev.push_back(*task2);
                    delete task2;
                    task_cnt++;


                    cost[x - 1][y - 1] = c;
                    cost[y - 1][x - 1] = c;
                    demand[x - 1][y - 1] = d;
                    demand[y - 1][x - 1] = d;

                    getline(in, str);
                    str = trim_tw(str);
                }
            }

            str_const = "LISTA_ARISTAS_NOREQ :";
            if (str == str_const) {
                getline(in, str);
                str = trim_tw(str);
                int x, y;
                double c;
                std::string temp_x, temp_y, temp_c;
                while (str[0] == '(') {
                    str += "\r";
                    temp_x = "";
                    temp_y = "";
                    temp_c = "";
                    int line_index = 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_x += str[line_index];
                        line_index++;
                    }
                    x = std::stoi(temp_x);

                    line_index += 2;

                    while (str[line_index] <= '9' && str[line_index] >= '0') {
                        temp_y += str[line_index];
                        line_index++;
                    }
                    y = std::stoi(temp_y);

                    line_index += 11;
                    while (str[line_index] != ' ') {
                        temp_c += str[line_index];
                        line_index++;
                    }
                    c = std::stod(temp_c);


                    cost[x - 1][y - 1] = c;
                    cost[y - 1][x - 1] = c;

                    getline(in, str);
                    str = trim_tw(str);
                }
            }


            //NOMBRE 所属UCARP的文件夹
            //COMENTARIO 不知道
            //VERTICES 节点个数1-77
            //ARISTAS_REQ: 一个CARP中Task的数量
            //ARISTAS_NOREQ : 一个CARP中不是Task的数量（但是联通不包含-1）
            //VEHICULOS : 车的数量
            //CAPACIDAD : 一辆车的容量
            //TIPO_COSTES_ARISTAS : 不知
            //COSTE_TOTAL_REQ : 1395
            //LISTA_ARISTAS_REQ :不知
            //...
            //LISTA_ARISTAS_NOREQ : 不是task的边
            //...
            //DEPOSITO 回收站

        }
        for (auto & j : task_rev) {
            j.task_number += (task_cnt-1);
        }
        for (const auto & j : task_rev) {
            task.push_back(j);
        }
        Instance* instance = new Instance(task, cost, demand, points_number, task_cnt-1, Q,cost_mul,demand_mul);
        instance_list[i] = *instance;

        in.close();
    }

    return instance_list;
}