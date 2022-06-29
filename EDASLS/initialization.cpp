#include <iostream>
#include <random>
#include <vector>
#include "Task.h"
#include <chrono>
#include <ctime>
#include "initialization.h"
#include <algorithm>
#include <climits>

using namespace std;

//extern int deadheading_multiplier[32];
//extern int serving_multiplier[32];

void PS(Instance &instance,
        vector<int> &chromosome, int kind) {

    int len = instance.UT;
    double left_capacity = instance.capacity;
    int current_node = 0;
    while (chromosome.size() != len) {
        int next_task;
        bool isfull = true;
        if (kind == 0) {
            double criteria = INT_MAX;
            for (int i = 0; i < 2 * len; ++i) {
                if (instance.tasklist[i].enable && left_capacity >= instance.tasklist[i].demand) {
                    if (instance.tasklist[i].demand <= 0.000000001 && instance.tasklist[i].demand >= -0.000000001) {
                        next_task = i;
                        isfull = false;
                        break;
                    }
                    double ratio = instance.tasklist[i].cost / instance.tasklist[i].demand;
                    if (ratio < criteria) {
                        criteria = ratio;
                        next_task = i;
                        isfull = false;
                    }
                }
            }
        } else if (kind == 1) {
            double biggest = 0;
            for (int i = 0; i < 2 * len; ++i) {
                if (instance.tasklist[i].enable && left_capacity >= instance.tasklist[i].demand) {
                    if (instance.tasklist[i].demand <= 0.000000001 && instance.tasklist[i].demand >= -0.000000001) {
                        next_task = i;
                        isfull = false;
                        break;
                    }
                    double ratio = instance.tasklist[i].cost / instance.tasklist[i].demand;
                    if (ratio > biggest) {
                        biggest = ratio;
                        next_task = i;
                        isfull = false;
                    }
                }
            }
        } else if (kind == 2 || (kind == 4 && left_capacity <= instance.capacity / 2)) {
            double smallest = INT_MAX;
            for (int i = 0; i < 2 * len; ++i) {
                if (instance.tasklist[i].enable && left_capacity >= instance.tasklist[i].demand) {
                    if (instance.tasklist[i].demand <= 0.000000001 && instance.tasklist[i].demand >= -0.000000001) {
                        next_task = i;
                        isfull = false;
                        break;
                    }
                    double ratio = instance.map[0][instance.tasklist[i].end_point - 1];
                    if (ratio < smallest) {
                        smallest = ratio;
                        next_task = i;
                        isfull = false;
                    }
                }
            }
        } else if (kind == 3 || (kind == 4 && left_capacity > instance.capacity / 2)) {
            double biggest = 0;
            for (int i = 0; i < 2 * len; ++i) {
                if (instance.tasklist[i].enable && left_capacity >= instance.tasklist[i].demand) {
                    if (instance.tasklist[i].demand <= 0.000000001 && instance.tasklist[i].demand >= -0.000000001) {
                        next_task = i;
                        isfull = false;
                        break;
                    }
                    double ratio = instance.map[0][instance.tasklist[i].end_point - 1];
                    if (ratio > biggest) {
                        biggest = ratio;
                        next_task = i;
                        isfull = false;
                    }
                }
            }
        } else {
            double smallest = INT_MAX;
            for (int i = 0; i < 2 * len; ++i) {
                if (instance.tasklist[i].enable && left_capacity >= instance.tasklist[i].demand) {
                    if (instance.tasklist[i].demand <= 0.000000001 && instance.tasklist[i].demand >= -0.000000001) {
                        next_task = i;
                        isfull = false;
                        break;
                    }
                    double ratio = instance.map[current_node][instance.tasklist[i].start_point - 1];
                    if (ratio < smallest) {
                        smallest = ratio;
                        next_task = i;
                        isfull = false;
                    }
                }
            }
        }
        if (isfull) {
            current_node = 0;
            left_capacity = instance.capacity;
        } else {
            chromosome.push_back(instance.tasklist[next_task].task_number);
            instance.tasklist[next_task].enable = false;
            instance.tasklist[(next_task + len) % (2 * len)].enable = false;
            current_node = instance.tasklist[next_task].end_point - 1;
            left_capacity -= instance.tasklist[next_task].demand;
        }
    }
}

int **PS_BIH(Instance *inslist, int ins_size, int pop_size, int fileInd) {
    int cost_mul = inslist->cost_mul;
    int demand_mul = inslist->demand_mul;
    for (int i = 0; i < ins_size; i++) {
        inslist[i].map = new double *[inslist[i].point_number];
        for (int j = 0; j < inslist[i].point_number; ++j) {
            inslist[i].map[j] = new double[inslist[i].point_number];
        }
//        if (i==20){
//            for (int j = 0; j < inslist[i].point_number; ++j) {
//                for (int k = 0; k < inslist[i].point_number; ++k) {
//                    cout<<j<<" : "<<k<<" : "<<inslist[i].cost[j][k]<<endl;
//                }
//            }
//        }
    }

    for (int i = 0; i < ins_size; i++) {
        for (int j = 0; j < inslist[i].point_number; ++j) {
            for (int k = 0; k < inslist[i].point_number; ++k) {
//                if (i==20){
//                    cout<<j<<" : "<<k<<" : "<<shortest_distance(j, k, inslist[i].cost, inslist[i].point_number)<<endl;
//                }

                inslist[i].map[j][k] = shortest_distance(j, k, inslist[i].cost, inslist[i].point_number);
            }
        }
    }

    for (int ins_cur = 0; ins_cur < ins_size; ++ins_cur) {
        int n = inslist[ins_cur].tasklist.size();
        inslist[ins_cur].e_map = new double *[n];
        for (int i = 0; i < n; ++i) {
            inslist[ins_cur].e_map[i] = new double[n];
        }
    }

    for (int ins_cur = 0; ins_cur < ins_size; ++ins_cur) {
        int n = inslist[ins_cur].tasklist.size();
        for (int i = 0; i < n; ++i) {
            Task ti = inslist[ins_cur].tasklist[i];
            for (int j = 0; j < n; ++j) {
                Task tj = inslist[ins_cur].tasklist[j];
                if (i == j) {
                    inslist[ins_cur].e_map[i][j] = inslist[ins_cur].map[0][tj.start_point - 1]
                                                   + tj.cost
                                                   + inslist[ins_cur].map[tj.end_point - 1][0];
                } else {
                    inslist[ins_cur].e_map[i][j] = inslist[ins_cur].map[ti.end_point - 1][tj.start_point - 1]
                                                   + tj.cost
                                                   + inslist[ins_cur].map[tj.end_point - 1][0]
                                                   - inslist[ins_cur].map[ti.end_point - 1][0];
                }
            }
        }
    }

    for (int ins_cur = 0; ins_cur < ins_size; ++ins_cur) {
        int n = inslist[ins_cur].tasklist.size();
        inslist[ins_cur].t_map = new double *[n];
        for (int i = 0; i < n; ++i) {
            inslist[ins_cur].t_map[i] = new double[n];
        }
    }

    for (int ins_cur = 0; ins_cur < ins_size; ++ins_cur) {
        int n = inslist[ins_cur].tasklist.size();
        for (int i = 0; i < n; ++i) {
            Task ti = inslist[ins_cur].tasklist[i];
            for (int j = 0; j < n; ++j) {
                Task tj = inslist[ins_cur].tasklist[j];
                if (i == j) {
                    inslist[ins_cur].t_map[i][j]
                            = inslist[ins_cur].map[0][tj.start_point - 1]
                                    * cost_mul + tj.cost * demand_mul;
                } else {
                    inslist[ins_cur].t_map[i][j]
                    = inslist[ins_cur].map[ti.end_point - 1][tj.start_point - 1]
                            * cost_mul + tj.cost * demand_mul;
                }
            }
        }
    }

    int *a = new int[ins_size];
    for (int i = 0; i < ins_size; ++i) {
        a[i] = i;
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(a, a + ins_size, std::default_random_engine(seed));

    vector<int> chromosome_set[ins_size];
    for (int i = 0; i < ins_size; ++i) {
        PS(inslist[a[i]], chromosome_set[i], i);
    }

    int **res = new int *[pop_size];
    int len = inslist[0].UT;
    for (int i = 0; i < pop_size; ++i) {
        res[i] = new int[len];
    }
    for (int m = 0; m < ins_size; ++m) {
        for (int i = 0; i < len; ++i) {
            res[m][i] = chromosome_set[m][i];
        }
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<double> dis(0, 1);

    for (int i = ins_size; i < pop_size; ++i) {
        for (int j = 0; j < len; ++j) {
//            int rand_num = rand() % (999 + 1) / (float)(999 + 1);
            double rand_num = dis(generator);
            res[i][j] = j + 1;
            if (rand_num < 0.5) {
                res[i][j] += len;
            }
        }
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle(res[i], res[i] + len, std::default_random_engine(seed));
        bool flag_break;
        for (int j = 0; j < i; ++j) {
            flag_break = true;
            for (int k = 0; k < len; ++k) {
                if (res[j][k] != res[i][k]) {
                    flag_break = false;
                    break;
                }
            }
            if (flag_break) {
                i--;
                break;
            }
        }
    }
    delete[] a;
    return res;
}

double shortest_distance(int point_1, int point_2, double **cost, int point_number) {
    double *dist = new double[point_number];
    for (int i = 0; i < point_number; ++i) {
        if (i == point_1) {
            dist[i] = 0;
        } else {
            dist[i] = INT_MAX;
        }
    }
    bool *checking = new bool[point_number];
    for (int i = 0; i < point_number; ++i) {
        checking[i] = true;
    }
    while (checking[point_2]) {
        double smallest_dist = INT_MAX;
        int smallpoint;
        for (int i = 0; i < point_number; ++i) {
            if (checking[i] && dist[i] < smallest_dist) {
                smallpoint = i;
                smallest_dist = dist[i];
            }
        }
        checking[smallpoint] = false;
        for (int i = 0; i < point_number; ++i) {
            if (cost[i][smallpoint] != -1 && cost[i][smallpoint] + dist[smallpoint] < dist[i]) {
                dist[i] = cost[i][smallpoint] + dist[smallpoint];
            }
        }
    }

    double kk = dist[point_2];
    delete[] checking;
    delete[] dist;
    return kk;
}
