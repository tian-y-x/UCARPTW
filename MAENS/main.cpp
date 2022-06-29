#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include "functions.h"
#include "heutistic.h"
#include <string>
#include "MAENS/MAENS.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "ins30.h"
#define tes
using namespace std;

int evaluate_time;
int req_arc_num = 0;
int req_edge_num;
int nonreq_arc_num = 0;
int nonreq_edge_num;
int vertex_num;
int vehicle_num;
int capacity;
int task_num;
int total_arc_num;
int DEPOT;
int cost_mul;
int demand_mul;
int iter_time;
double start_time;
double trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
double serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
double min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
double shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
ofstream outfile;
int* res[33];
int ind;
double eva_cost[31][FEM];
double curve[FEM];
double time_curve[FEM];
int FEm = 300000;
int execution = 20;
void ShowMatrix(int(*Min)[MAX_NODE_TAG_LENGTH], int dim);
void init_flush();
const string dir_path = "../UCARP_MAENS_HS/UncertainInstances_from_tw/";
const string rpath = "../UCARP_MAENS_HS/UncertainInstances_from_tw/";
Individual best_rob;
string filename[32] = {
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
ins30* instances;
int main()
{
    best_rob.TotalCost = INF;
    best_rob.TotalCost = INF;
    Task inst_tasks[MAX_TASKS_TAG_LENGTH];
    Arc inst_arcs[MAX_ARCS_TAG_LENGTH];

    for (ind = 0; ind < 32; ind++)
    {
        instances = new ins30[31];
        double result_fitness[30];
        memset(best_rob.Sequence, -1, sizeof(int) * 250);
        ofstream ofs;
        ofs.open(filename[ind] + "_MAENS_curve.csv", ios::out);
        ofs << "FEm,best" << endl;
        for (int ins = 1; ins <= 30; ins++)
        {
            init_flush();
            double res[31];
            string path = rpath + filename[ind] + "/" + to_string(ins) + ".dat";
            //            cout << path << endl;
            data_process(instances[ins].inst_tasks, instances[ins].inst_arcs, path);
            instances[ins].init_ins();
            instances[ins].mod_dijkstra();

        }
        //        cout << "dataset " << ind << endl;
        res[ind] = new int[100];
        init_flush();

        //"S:\\MyCppSpace\\UCARPTW_MAENS\\UncertainInstances_tw\\TW-B40C\\"; 
        string path = dir_path + filename[ind] + "/" + to_string(ind) + ".dat";

        path = "../UCARP_MAENS_HS/TWdata/" + filename[ind] + ".dat";

        //        cout << path << endl;
        //        cout << filename[ind] << endl;
        readFile(inst_tasks, inst_arcs, path);
        for (double& i : curve) {
            i = 0;
        }
        for (double& i : time_curve) {
            i = 0;
        }


        iter_time = 1;
        double sum_time = 0;
        while (iter_time <= execution)
        {
            start_time = clock();
            best_rob.TotalCost = INF;
            best_rob.TotalCost = INF;
            memset(best_rob.Sequence, -1, sizeof(int) * 250);
            evaluate_time = 0;
            //data_process_win_tw(inst_tasks, inst_arcs, path);
            for (int i = 1; i <= vertex_num; i++)
            {
                for (int j = 1; j <= vertex_num; j++)
                {
                    trav_cost[i][j] = INF;
                    serve_cost[i][j] = 0;
                }
            }

            trav_cost[1][0] = 0;
            trav_cost[0][1] = 0;

            for (int i = 1; i <= task_num; i++)
            {
                serve_cost[inst_tasks[i].head_node][inst_tasks[i].tail_node] = inst_tasks[i].serv_cost;
            }

            for (int i = 1; i <= total_arc_num; i++)
            {
                trav_cost[inst_arcs[i].head_node][inst_arcs[i].tail_node] = inst_arcs[i].trav_cost;
            }

            mod_dijkstra();
#ifdef test
            //            for (int e = 0; e < 50; e++)
            //            {
            //                cout << instances[1].inst_arcs[e].trav_cost << endl;
            //            }


            //            for (int d1 = 0; d1 < 30; d1++)
            //            {
            //                cout << endl;
            //                for (int d2 = 0; d2 < 30; d2++)
            //                {
            //                    cout << min_cost[d1][d2] << ",";
            //                }
            //
            //            }
            //            for (int d = 1; d < 31; d++)
            //            {
            //                cout << d << endl;
            //                for (int d1 = 0; d1 < 30; d1++)
            //                {
            //                    cout << endl;
            //                    for (int d2 = 0; d2 < 30; d2++)
            //                    {
            //                        cout << instances[d].min_cost[d1][d2] << ",";
            //                    }
            //
            //                }
            //            }
#endif
            //int test_seq[50] = { 15,0,6,11,10,9,0,5,7,19,0,4,14,1,2,0

            //};

            //double test_cost = get_total_cost(test_seq, inst_tasks);
            ////cout << test_cost << endl;

            MAENS(inst_tasks);

            outfile.open(filename[ind] + "_res.txt");

            //            for (int i = 0; i <= res[ind][0]; i++)
            //            {
            //                //cout << best_fsb_solution.Sequence[i] << " " << inst_tasks[best_fsb_solution.Sequence[i]].demand << endl;
            //                cout << res[ind][i] << "," << endl;
            //
            //                if (i == res[ind][0])
            //                {
            //                    outfile << res[ind][i] << endl;
            //                    get_total_cost(res[ind], inst_tasks);
            //                }
            //                else
            //                    outfile << res[ind][i] << ",";
            //            }
            sum_time += (clock() - start_time) / CLOCKS_PER_SEC;
            for (int i = 0; i <= best_rob.Sequence[0]; i++)
            {
                if (i == best_rob.Sequence[0])
                {
                    outfile << best_rob.Sequence[i] << endl;
                    cout << best_rob.Sequence[i] << endl;
                    //get_total_cost(res[ind], inst_tasks);
                }
                else
                    cout << best_rob.Sequence[i] << ",";
                outfile << best_rob.Sequence[i] << ",";
            }
            result_fitness[iter_time - 1] = best_rob.TotalCost;
            iter_time++;

        }
        for (int i = 0; i < FEm; i++) {
            ofs << i << "," << curve[i] / execution << ","<< time_curve[i]  << std::endl;
        }
        double variance_exc = variance(result_fitness, execution);
        //double variance_pop = variance(curve, 210);
        ofs << "variance_exc_std" << "," << variance_exc << endl;
        //ofs << "variance_pop_std" << "," << variance_pop << endl;
        ofs << "times_results" << "," << sum_time / execution << endl;
        for (int exce_times = 0; exce_times < execution; exce_times++)
        {
            ofs << exce_times << "std" << "," << result_fitness[exce_times] << endl;
        }

        ofs.close();

        delete[] instances;

    }
    printf("Hello, World!\n");
    return 0;
}


void init_flush()
{
    req_arc_num = 0; //NRA
    req_edge_num = 0; //NRE
    nonreq_arc_num = 0;
    nonreq_edge_num = 0;
    vertex_num = 0;
    vehicle_num = 0;
    capacity = 0;
    task_num = 0;
    total_arc_num = 0;
    DEPOT = 1;

    memset(&trav_cost[0][0], 0, sizeof(double) * MAX_NODE_TAG_LENGTH * MAX_NODE_TAG_LENGTH);
    memset(&serve_cost[0][0], 0, sizeof(double) * MAX_NODE_TAG_LENGTH * MAX_NODE_TAG_LENGTH);
    memset(&min_cost[0][0], 0, sizeof(double) * MAX_NODE_TAG_LENGTH * MAX_NODE_TAG_LENGTH);
    memset(&shortest_path[0][0][0], 0, sizeof(double) * MAX_NODE_TAG_LENGTH * MAX_NODE_TAG_LENGTH * MAX_NODE_TAG_LENGTH);

    costlb = INF;
    costub = 0;
    demandlb = INF;
    demandub = 0;
}

void ShowMatrix(int(*Min)[MAX_NODE_TAG_LENGTH], int dim)
{
    for (int i = 0; i <= dim; i++)
    {
        for (int j = 0; j <= dim; j++)
        {
            printf("%d \t", Min[i][j]);
        }
        printf("\n");
    }
}

double variance(double* x, int len)
{
    double sum = 0;
    for (int i = 0; i < len; i++)
        sum += x[i];
    double average = sum / len;
    double num = 0;
    for (int i = 0; i < len; i++)
        num += pow(x[i] - average, 2);
    return sqrt(num / len);
}

