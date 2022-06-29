#include "functions.h"
#include <iostream>
#include <string>
#include <ctime>
using namespace std;

int rand_choose(int num)
{
    if (num == 0)
        exit(-1);
    int k = rand()%num;

    k++;

    return k;
}

void rand_selection(int *id1, int *id2, int popsize)
/* pop has been sorted increasingly already */
{
    int k1, k2;
    int* candi = (int*)malloc(sizeof(int)*(popsize + 1));
    candi[0] = popsize;
    for (int i = 1; i <= popsize; i++)
    {
        candi[i] = i-1;
    }

    k1 = rand_choose(candi[0]);
    *id1 = candi[k1];
    delete_element(candi, k1);
    k2 = rand_choose(candi[0]);
    *id2 = candi[k2];
}

int Max(int *Array)
{
    int Length = Array[0];
    int i, maximum;

    maximum = 0;
    for (i = 1; i <= Length; i++)
    {
        if (Array[i] > maximum)
        {
            maximum = Array[i];
        }
    }

    return maximum;
}

void ReverseDirection(int *Array, int k1, int k2)
// reverse the subarray of Array[1:n] from the position k1 to position k2
{
    int i, k, tmp;
    double m = (k2-k1+1)/2;
    k = (int)m;

    for (i = k1; i < k1+k; i++)
    {
        tmp = Array[i];
        Array[i] = Array[k1+k2-i];
        Array[k1+k2-i] = tmp;
    }
}


void AssignArray(int *Array1, int *Array2)
// assign array1 to array2
{
    int i;

    for (i = 0; i <= Array1[0]; i++)
    {
        Array2[i] = Array1[i];
    }
}

void AssignSubArray(int *Array1, int k1, int k2, int *Array2)
// assign array1[k1:k2] to array2
{
    int i;

    Array2[0] = k2-k1+1;

    for (i = k1; i <= k2; i++)
    {
        Array2[i-k1+1] = Array1[i];
    }
}

void JoinArray(int *JointArray, int *Array)
{
    int i;
    for (i = 1; i <= Array[0]; i++)
    {
        JointArray[0] ++;
        JointArray[JointArray[0]] = Array[i];
    }
}

int RandChoose(int n)
// randomly choose a number between 1 and n
{
    int random;

    random = rand();

    int k = random%n;

    k++;

    return k;
}

void find_ele_positions(int *positions, int *a, int e)
{
    positions[0] = 0;
    for (int i = 1; i <= a[0]; i++)
    {
        if (a[i] == e)
        {
            positions[0] ++;
            positions[positions[0]] = i;
        }
    }
}

void delete_element(int *a, int k)
{
    if (k < 1 || k > a[0])
    {
        printf("%d %d \n", a[0], k);
        printf("the deleting position is wrong!\n");
        return;
        exit(-1);
    }

    for (int i = k; i < a[0]; i++)
    {
        a[i] = a[i+1];
    }
    a[0] --;
}
void delete_element(double* a, int k)
{
    if (k < 1 || k > a[0])
    {
        printf("%lf %d \n", a[0], k);
        printf("the deleting position is wrong!\n");
        exit(-1);
    }

    for (int i = k; i < a[0]; i++)
    {
        a[i] = a[i + 1];
    }
    a[0] --;
}

void add_element(int *a, int e, int k)
// add element e in k position of a
{
    if (k < 1 || k > a[0]+1)
    {
        printf("the inserting position is wrong!\n");
        exit(-1);
        return;
    }

    a[0] ++;
    for (int i = a[0]; i > k; i--)
    {
        a[i] = a[i-1];
    }
    a[k] = e;
}

double get_task_seq_total_cost(int *task_seq, const Task *inst_tasks)
{
    double total_cost =  0;
    for (int i = 1; i < task_seq[0]; i++)
    {
        total_cost += min_cost[inst_tasks[task_seq[i]].tail_node][inst_tasks[task_seq[i+1]].head_node]+inst_tasks[task_seq[i]].serv_cost;
    }

    return total_cost;
}

double penalize(Task tj, double time) {
    double w = 1;
    double v = 30;
    double w2v = w / v;
    if (IS_DOUBLE_ZERO(tj.an) == 1&& IS_DOUBLE_ZERO(tj.bn) == 1)
        return 0;
    if (time < tj.an) {
        return w2v * (tj.an - time) * (tj.an - time) / (tj.bn - tj.an);
    }
    else if (time > tj.bn) {
        return w2v * (time - tj.bn) * (time - tj.bn) / (tj.bn - tj.an);
    }
    else {
        return 0;
    }
}
int is_valid(const int* task_seq, const Task* inst_tasks)
{
    double load = 0;
    for (int i = 1; i <= task_seq[0]; i++)
    {
        load += inst_tasks[i].demand;
        if (task_seq[i] == 0)
        {
            if (load <= capacity)
            {
                return 0;
            }
            load = 0;
        }
    }
    return 1;
}
double eva_fitness(const int* task_seq, ins30* instances, const Task* inst_tasks)
{
//    double start=clock();
    //if (!is_valid(task_seq, inst_tasks))
    //{
    //    eva_cost[iter_time][evaluate_time] = best_rob.TotalCost;
    //    return -1;
    //}
    
    double max_cost = 0;
    for (int ins = 1; ins <= 30; ins++)
    {

        int task_copy[250];
        memcpy(task_copy, task_seq, sizeof(int) * 250);
        for (int i = 1; i < task_copy[0]; i++)
        {
            if (IS_DOUBLE_ZERO(instances[ins].inst_tasks[task_copy[i]].demand) == 1 && task_copy[i] != 0)
            {
                delete_element(task_copy, i);
            }
        }
        double total_cost = 0;
        double distance = 0;
        double time = 0;
        double time_cost = 0;
        int last_point = -1;
        if (!is_valid(task_copy, instances[ins].inst_tasks))
        {
            double load = 0;
            for (int i = 1; i < task_copy[0]; i++)
            {
                if (task_copy[i] == 0)
                {
                    time = 0;
                    load = 0;
                }
                load += instances[ins].inst_tasks[task_copy[i]].demand;
                if (load > capacity)
                {
                    add_element(task_copy, 0, i);
                    load = 0;
                }
            }
        }

        for (int i = 1; i < task_copy[0]; i++)
        {
            /*if (IS_DOUBLE_ZERO(instances[ins].inst_tasks[task_seq[i]].demand) == 1 && task_seq[i] != 0)
            {
                continue;
            }*/
            if (task_copy[i] == 0)
                time = 0;

            distance += instances[ins].min_cost[instances[ins].inst_tasks[task_copy[i]].tail_node]
                [instances[ins].inst_tasks[task_copy[i + 1]].head_node]
            + instances[ins].inst_tasks[task_copy[i]].serv_cost;

            //time += instances[ins].min_cost[instances[ins].inst_tasks[task_copy[i]].tail_node]
            //    [instances[ins].inst_tasks[task_copy[i + 1]].head_node] * cost_mul
            //    + instances[ins].inst_tasks[task_copy[i]].serv_cost * demand_mul;
            time += instances[ins].inst_tasks[task_copy[i]].serv_cost * demand_mul;

            time_cost += penalize(instances[ins].inst_tasks[task_copy[i]], time);


            time += instances[ins].min_cost[instances[ins].inst_tasks[task_copy[i]].tail_node]
                [instances[ins].inst_tasks[task_copy[i + 1]].head_node] * cost_mul;



            //cout << "distance" << distance << endl;
            //std::cout << inst_tasks[task_copy[i]].head_node << " to " << inst_tasks[task_copy[i]].tail_node << ", " <<
            //    inst_tasks[task_copy[i]].tail_node << " to " <<
            //    inst_tasks[task_copy[i + 1]].head_node << std::endl;

        }
        total_cost = distance + time_cost;
        /*cout << ins << ":";
        cout << total_cost << endl;*/
        if (total_cost > max_cost)
        {
            max_cost = total_cost;
            //cout << ins << endl;
        }
    }
    
    //cout << max_cost << endl;
    if (max_cost < best_rob.TotalCost)
    {
        best_rob.TotalCost = max_cost;
        memcpy(best_rob.Sequence, task_seq, sizeof(task_seq)*250);
    }
    eva_cost[iter_time][evaluate_time] = best_rob.TotalCost;
    //cout << evaluate_time << "," << best_rob.TotalCost << endl;
//    cout<<"rs_time: "<<clock()-start<<endl;
    return max_cost;
}
double get_total_cost(int* task_seq, const Task* inst_tasks)
{
//    int task_seq2[25] = { 15,0,8,12,3,15,0,7,11,10,9,13,0,5,17,0 };
//    eva_fitness(task_seq2, instances, inst_tasks);
    
    double fitness=eva_fitness(task_seq, instances, inst_tasks);
    if (evaluate_time < FEm) {
        curve[evaluate_time] += best_rob.TotalCost;
        time_curve[evaluate_time] += (clock() - start_time) / CLOCKS_PER_SEC;
    }
    
    evaluate_time++;
    return fitness;
//    double total_cost = 0;
//    double distance = 0;
//    double time = 0;
//    double time_cost = 0;
//    for (int i = 1; i < task_seq[0]; i++)
//    {
//        if (task_seq[i] == 0)
//            time = 0;
//
//        distance += min_cost[inst_tasks[task_seq[i]].tail_node][inst_tasks[task_seq[i + 1]].head_node]
//            + inst_tasks[task_seq[i]].serv_cost;
//        time += min_cost[inst_tasks[task_seq[i]].tail_node][inst_tasks[task_seq[i + 1]].head_node] * cost_mul
//            + inst_tasks[task_seq[i]].serv_cost * demand_mul;
//        time_cost += penalize(inst_tasks[task_seq[i]], time);
//    }
//    total_cost = distance + time_cost;
//    return total_cost;
}


double get_total_vio_load(double *route_seg_load)
{
    double total_vio_load = 0;
    for (int i = 1; i <= (int)route_seg_load[0]; i++)
    {
        if (route_seg_load[i] > capacity)
            total_vio_load += route_seg_load[i]-capacity;
    }

    return total_vio_load;
}

int FindTask(int a, int b, const Task *inst_tasks, int NO_Task)
{
    int i;

    for (i = 1; i <= NO_Task; i++)
    {
        if (inst_tasks[i].head_node == a && inst_tasks[i].tail_node == b)
            return i;
    }

    return 0;
}