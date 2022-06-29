#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <ctime>
#include <math.h>
#include "instance.h"
#include "initialization.h"
#include "ehbsa.h"
#include "sls.h"
#include "split_v1.h"
#include "read_ucarp_tw.h"


double *health;
double Bratio = 0.01;
int fileInd;
const int size_of_UCARP_TW_instance = 32;
const int INSTANCE_NUM = 30;
const int FEm = 300000;
const int executions = 20;
int pop_size = 120;
double Pls = 0.1;
extern double best_rs;
double curve[FEm];
double time_curve[FEm];
double start_time=0;

extern int routes[INSTANCE_NUM][500];
extern int evaluate_times;

bool compare(const int a, const int b) {
    return health[a] < health[b];
}

int main() {

    double scores[executions];
    double running_time[executions];

    string dataset[size_of_UCARP_TW_instance] = {
            "TW-A10A",
            "TW-A13A",
            "TW-A13B",
            "TW-A13C",
            "TW-A20B",
            "TW-A40C",
            "TW-A40D",
            "TW-A60A",
            "TW-B10A",
            "TW-B13A",
            "TW-B13B",
            "TW-B13C",
            "TW-B20B",
            "TW-B40C",
            "TW-B40D",
            "TW-B60A",
            "TW-C10A",
            "TW-C13A",
            "TW-C13B",
            "TW-C13C",
            "TW-C20B",
            "TW-C40C",
            "TW-C40D",
            "TW-C60A",
            "TW-D10A",
            "TW-D13A",
            "TW-D13B",
            "TW-D13C",
            "TW-D20B",
            "TW-D40C",
            "TW-D40D",
            "TW-D60A"
    };

    ofstream ofs;
    ofstream t_ofs;
    for (fileInd = 0; fileInd < size_of_UCARP_TW_instance; ++fileInd) {
        cout << dataset[fileInd] << endl;

        string path = "../Dataset/UncertainInstances_from_tw/"+ dataset[fileInd]+"/";
//        string path = "..\\Dataset\\UncertainInstances_from_tw\\" + dataset[fileInd] + "\\";

        ofs.open(dataset[fileInd] + "_curve.csv", ios::out);
        ofs << "FEm,best" << endl;

        t_ofs.open(dataset[fileInd] + "time_curve.csv", ios::out);
        t_ofs << "time,best" << endl;

        for (double & i : curve) {
            i=0;
        }
        for (double & i : time_curve) {
            i=0;
        }

        for (int times = 0; times < executions; ++times) {
//            cout << times << endl;
            best_rs=9999999;
            Instance *instance_list = data_process_macOS_tw(path);
//            Instance *instance_list = data_process_win_tw(path);

            start_time = clock();
            std::vector<Task> task = instance_list[0].tasklist;
            int UT = instance_list[0].UT; //|UT|
            int **pop = PS_BIH(instance_list, INSTANCE_NUM, pop_size, fileInd);
            health = new double[pop_size];

            for (int i = 0; i < pop_size; ++i) {
                for (int j = 0; j < instance_list[0].UT; ++j) {
                    pop[i][j] -= 1;
                }
            }
            int *index = new int[pop_size];


            //test if every initialization is similar
            for (int i = 0; i < pop_size; ++i) {
                health[i] = CEvaluate(instance_list, vector<int>(pop[i], pop[i] + UT));
                index[i] = i;
            }
            sort(index, index + pop_size, compare);

            int cur = 0;
            int **pop_temp = new int *[pop_size / 2];
            for (int i = 0; i < pop_size / 2; ++i) {
                pop_temp[i] = new int[UT];
            }
            evaluate_times=0;//
            while (evaluate_times < FEm) {
                for (int i = 0; i < pop_size / 2; ++i) {
                    for (int j = 0; j < UT; ++j) {
                        pop_temp[i][j] = pop[index[i]][j];
                    }
                }
                double epsilon = ((pop_size / 2) * Bratio) / (UT - 1);
                double **ehm = build_ehm(pop_temp, pop_size / 2, UT, epsilon);
                int rand_num = rand() % pop_size;
                int *child_chormosome = EHBSA(pop[rand_num], ehm, 2, UT);

                double rand_prob = (rand() % 100) / (float) (100);
                vector<int> child_chromosome_vec = vector<int>(child_chormosome, child_chormosome + UT);
                if (rand_prob < Pls && evaluate_times < FEm - 1) {
                    StochasticSLS(child_chromosome_vec, ehm, instance_list);
                    bool flag_break = true;
                    for (int i = 0; i < pop_size; ++i) {
                        flag_break = true;
                        for (int k = 0; k < UT; ++k) {
                            if (pop[i][k] != child_chromosome_vec[k]) {
                                flag_break = false;
                                break;
                            }
                        }
                        if (flag_break) {
                            break;
                        }
                    }
                    if (!flag_break) {
                        for (int i = 0; i < UT; ++i) {
                            child_chormosome[i] = child_chromosome_vec[i];
                        }
                    }
                }
                for (int s = 0; s < 2 * UT; s++)
                    delete[] ehm[s];
                delete[] ehm;


                double child_value = CEvaluate(instance_list, vector<int>(child_chormosome, child_chormosome + UT));

                if (child_value < health[rand_num]) {
                    bool flag_break = true;
                    for (int i = 0; i < pop_size; ++i) {
                        flag_break = true;
                        for (int k = 0; k < UT; ++k) {
                            if (pop[i][k] != child_chormosome[k]) {
                                flag_break = false;
                                break;
                            }
                        }
                        if (flag_break) {
                            break;
                        }
                    }

                    if (!flag_break) {
                        for (int i = 0; i < UT; ++i) {
                            pop[rand_num][i] = child_chormosome[i];
                        }
                        health[rand_num] = child_value;
                    }
                }
                sort(index, index + pop_size, compare);

                delete[] child_chormosome;
                cur++;
            }

            for (int s = 0; s < pop_size / 2; s++)
                delete[] pop_temp[s];
            delete[] pop_temp;

            double end_time = clock();
            cout << "best Rs: " << best_rs << endl;
//            scores[times] = health[index[0]];
            scores[times] = best_rs;
            running_time[times] = (end_time - start_time) / CLOCKS_PER_SEC;
//            cout << "time: " << running_time[times] << endl;
            for (int s = 0; s < pop_size; ++s) {
                delete[] pop[s];
            }
            delete[] pop;
            delete[] health;
            delete[] index;
            delete[] instance_list;

            extern int robust_solution[];
            int len=robust_solution[0];
//            for (int i = 1; i <= len; ++i) {
//                cout<<robust_solution[i]<<" ";
//            }
//            cout<<endl;

        }
        for (int i = 0; i < FEm; i++) {
            ofs << i+1 << "," << curve[i] / executions << std::endl;
        }
        for (int i = 0; i < FEm; i++) {
            t_ofs << i+1 << "," << time_curve[i] / executions << std::endl;
        }
        double sum_score = 0;
        double sum_time = 0;
        double var = 0;
        for (int i = 0; i < executions; ++i) {
            sum_score += scores[i];
            sum_time += running_time[i];
        }
        double average_score = sum_score / executions;
        double average_time = sum_time / executions;

        for (double score : scores) {
            var += pow((score - average_score), 2);
        }
        double std = pow((var / executions), 0.5);
        cout<<"avg_cost: "<<average_score<<endl;
        cout<<"avg_time: "<<average_time<<endl;
        cout<<"std: "<<std<<endl;
        t_ofs.close();
        ofs.close();
    }
}
