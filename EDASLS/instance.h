

//
// Created by 段心童 on 2021/9/23.
//

#ifndef DATA_INSTANCE_H
#define DATA_INSTANCE_H
#include "Task.h"
#include <vector>
#include <iostream>
using namespace std;
class Instance
{
public:
    vector<Task> tasklist;
    double** cost;
    double** demand;
    int point_number;
    double** map;
    double** e_map;
    double** t_map;
    int UT;
    int cost_mul;
    int demand_mul;
    double capacity;
    Instance();
    Instance(vector<Task> tasklist, double** cost, double** demand, int point_number, int UT, double capacity,int cost_mul,int demand_mul);
    ~Instance();
    Instance& operator =(Instance& a);
};

#endif //DATA_INSTANCE_H