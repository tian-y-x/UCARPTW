
#ifndef CARP_FUNCTIONS_H
#define CARP_FUNCTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#define MAX_TASKS_TAG_LENGTH 500
#define MAX_ARCS_TAG_LENGTH 1001
#define MAX_NODE_TAG_LENGTH 100
#define INF 1000000000

#define MAX_TASK_SEG_LENGTH 550
#define MAX_ROUTE_TAG_LENGTH 50

#define MAX_TASK_TAG_LENGTH 381
#define MAX_TASK_SEQ_LENGTH 250
#define MIN_VALUE 1e-8
#define IS_DOUBLE_ZERO(d)  (abs(d) < MIN_VALUE)

extern std::string filename[];
extern int cost_mul;
extern int demand_mul;
extern int evaluate_time;
extern int req_arc_num;
extern int nonreq_arc_num;
extern int vertex_num;
extern int req_edge_num;
extern int nonreq_edge_num;
extern int vehicle_num;
extern int capacity;
extern int task_num;
extern int total_arc_num;
extern int DEPOT;
extern int ind;
extern int* res[];
extern int iter_time;
extern int FEm;
extern double start_time;
extern double curve[];
extern double trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern double serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

extern double min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern double shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

extern double costlb;
extern double costub;
extern double demandlb;
extern double demandub;
extern std::ofstream outfile;
const int FEM = 300000;
extern double eva_cost[31][FEM];
extern double time_curve[FEM];

typedef struct task
{
    int head_node;
    int tail_node;
    double dead_cost;
    double serv_cost;
    double demand;
    int inverse;
    double an;
    double bn;
} Task;

// format: head ----> tail
typedef struct arc
{
    int head_node;
    int tail_node;
    double trav_cost;
    unsigned int change;
    int link;
} Arc;



typedef struct individual
{
    int Sequence[250]; // task id with depot inside
    int Assignment[250]; // when the representation is chromosome,this is used as the chromosome.
    double TotalCost;
    double Loads[50];
    double TotalVioLoad;
    double Fitness;
//    int start[50]; // for DCARP
} Individual;

extern Individual best_rob;
void mod_dijkstra();


class ins30 {
public:
    Task inst_tasks[MAX_TASKS_TAG_LENGTH];
    Arc inst_arcs[MAX_ARCS_TAG_LENGTH];
    double trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
    double serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
    double min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
    double shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
    void mod_dijkstra();

    void init_ins();
};
extern ins30* instances;

void readUncertainFile(Task* inst_tasks, Arc* inst_arcs, std::string path);

double eva_fitness(const int* task_seq, ins30* instances, const Task* inst_tasks);

void data_process(Task* inst_tasks, Arc* inst_arcs, std::string char_dir_path);

void readFile(Task* inst_tasks, Arc* inst_arcs, std::string char_dir_path);

double get_total_cost(int* task_seq, const Task* inst_tasks);

int rand_choose(int num);

void rand_selection(int *id1, int *id2, int popsize);

int Max(int *Array);

void find_ele_positions(int *positions, int *a, int e);

void delete_element(int *a, int k);

void delete_element(double* a, int k);

void add_element(int *a, int e, int k);

void add_element(double* a, int e, int k);

double get_task_seq_total_cost(int *task_seq, const Task *inst_tasks);

double get_total_vio_load(double *route_seg_load);

int FindTask(int a, int b, const Task *ARPTask, int NO_Task);

void AssignArray(int *Array1, int *Array2);

void AssignSubArray(int *Array1, int k1, int k2, int *Array2);

void JoinArray(int *JointArray, int *Array);

int RandChoose(int n);

void ReverseDirection(int *Array, int k1, int k2);

double variance(double* x, int len);

#endif 
