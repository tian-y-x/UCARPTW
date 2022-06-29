#include <iostream>
#include <vector>
#include "split_v1.h"
#include "instance.h"
#include "initialization.h"
#include "Task.h"
#include <ctime>

using namespace std;

const int INSTANCE_NUM = 30;
const int INF = 0x7fffffff;
const int DEPOT = 1;
const int MAX_TASK_TAG_LENGTH = 500;
double w=1;
double v=30;
double w2v=w/v;
int evaluate_times=0;
double V[MAX_TASK_TAG_LENGTH];
int P[MAX_TASK_TAG_LENGTH];
int routes[INSTANCE_NUM][MAX_TASK_TAG_LENGTH];
int robust_solution[MAX_TASK_TAG_LENGTH];
double best_rs=INF;
extern int fileInd;
extern int FEm;
extern double curve[];
extern double time_curve[];
extern double start_time;

//int deadheading_multiplier[32] = {
//    1,2,2,2,1,2,2,1,
//    1,2,2,2,1,2,2,1,
//    1,2,2,2,1,2,2,1,
//    1,2,2,2,1,2,2,1
//};
//
//int serving_multiplier[32] = {
//    2,5,4,3,3,3,5,2,
//    2,5,4,3,3,3,5,2,
//    2,5,4,3,3,3,5,2,
//    2,5,4,3,3,3,5,2
//};

inline bool is_zero(double tar){
    return tar>-0.000001&&tar<0.000001;
}

double execute(const Instance& instance, int rt){
//    int cost_mul = instance.cost_mul;
//    int demand_mul = instance.demand_mul;
    int length = routes[rt][0];
    double current_capacity = 0;
    double current_time = 0;
    double total_cost = 0;
    int last_pos = -2;
    for (int i = 2; i <= length; ++i) {
        if (routes[rt][i] == -1)
        {
            current_capacity = 0;
            current_time = 0;
            if (last_pos>=0)
                total_cost += instance.map[0][instance.tasklist[last_pos].end_point - 1];
            last_pos = -2;
        } else {
            Task task = instance.tasklist[routes[rt][i]];
            if (is_zero(task.demand)){
                continue;
            }
            if (current_capacity + task.demand > instance.capacity){
                current_capacity = task.demand;
                current_time = 0;
                if (last_pos>=0)
                    total_cost += instance.map[0][instance.tasklist[last_pos].end_point - 1];
                last_pos = -2;
            }else{
                current_capacity += task.demand;
            }
            if (last_pos==-2){
//                current_time += instance.map[0][task.start_point - 1] * cost_mul + task.cost * demand_mul;
                current_time +=instance.t_map[routes[rt][i]][routes[rt][i]];
                total_cost += instance.map[0][task.start_point - 1]
                        + task.cost
                        + penalize(task, current_time);
//                total_cost += instance.map[0][task.start_point - 1] + task.cost;
            }else{
//                current_time += instance.map[instance.tasklist[last_pos].end_point - 1][task.start_point - 1] * cost_mul + task.cost * demand_mul;
                current_time += instance.t_map[last_pos][routes[rt][i]];
                total_cost += instance.map[instance.tasklist[last_pos].end_point - 1][task.start_point - 1]
                              + task.cost
                              + penalize(task, current_time);
//                total_cost += instance.map[instance.tasklist[last_pos].end_point - 1][task.start_point - 1]
//                              + task.cost;
            }

            last_pos = routes[rt][i];
        }
    }
    return total_cost;
}

double CEvaluate(const Instance *instances, const vector<int> &chrom)
{

    //chech_chromosome(chrom,instances);
    double Fc = 0;
//    double start=clock();
    for (int i = 0; i < INSTANCE_NUM; i++)
    {
        double tmp_Fc= split(chrom, instances[i], routes[i]);
        if (Fc < tmp_Fc){
            Fc = tmp_Fc;
        }
    }
//    cout<<"split_time: "<<clock()-start<<endl;

    for (int i = 0; i < INSTANCE_NUM; i++)
    {
        double rs=0;
        for (int j = 0; j < INSTANCE_NUM; j++)
        {
            double res_ij=execute(instances[j],i);
            if (res_ij>rs)
                rs=res_ij;
        }
        if (rs<best_rs){
            best_rs=rs;
            int len=routes[i][0];
            for (int k = 0; k <= len; ++k) {
                robust_solution[k]=routes[i][k];
            }
        }
//        evaluate_times+=1;
    }
    if (evaluate_times<FEm){
        curve[evaluate_times] += best_rs;
        time_curve[evaluate_times] += (clock() - start_time) / CLOCKS_PER_SEC;
    }
    evaluate_times+=1;

//    if (smallest_rs<best_rs){
//        best_rs=smallest_rs;
//        int len=routes[rs_s][0];
//        for (int i = 0; i <= len; ++i) {
//            robust_solution[i]=routes[rs_s][i];
//        }
//    }
//    cout<<"total_time: "<<clock()-start<<endl;
    return Fc;
}

double split(const vector<int> &chrom, const Instance &instance, int *split_task_seq)
{
    int task_len = chrom.size();
    V[0] = 0;
    P[0] = 0;
    for (int i = 1; i <= task_len; i++){
        V[i] = INF;
    }
    for (int i = 0; i < task_len; i++)
    {
        double load = 0, cost = 0, time=0;
        int cur=i,last =i;
        while (cur < task_len)
        {
            Task tj = instance.tasklist[chrom[cur]];
            Task pre_tj = instance.tasklist[chrom[last]];
            load += tj.demand;
            time += instance.t_map[pre_tj.task_number-1][tj.task_number-1];
            cost += instance.e_map[pre_tj.task_number-1][tj.task_number-1]+penalize(tj, time);
//            cost += instance.e_map[pre_tj.task_number-1][tj.task_number-1];

            if (load <= instance.capacity)
            {
                double V_new = V[i] + cost;
                if (V_new < V[cur + 1])
                {
                    V[cur + 1] = V_new;
                    P[cur + 1] = i;
                }
                last=cur;
                cur++;
            }else{
                break;
            }
        }
    }
    split_task_seq[0] = 1;
    split_task_seq[1] = -1;
    int ptr = chrom.size();
    do {
        for (int k = P[ptr]; k < ptr; k++)
        {
            split_task_seq[++split_task_seq[0]] = chrom[k];
        }
        split_task_seq[++split_task_seq[0]] = -1;
        ptr = P[ptr];
    }while (ptr > 0);
    return V[task_len];
}

double penalize(Task tj, double time) {
    if (time < tj.an) {
        return w2v * (tj.an - time) * (tj.an - time) /  (tj.bn - tj.an);
    }
    else if (time > tj.bn) {
        return w2v * (time - tj.bn) * (time - tj.bn) / (tj.bn - tj.an);
    }
    else {
        return 0;
    }
}
