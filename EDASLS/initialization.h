//
// Created by lenovo on 2021/9/7.
//

#ifndef RESEARCH_INITIALIZATION_H
#define RESEARCH_INITIALIZATION_H

#include "Task.h"
#include "instance.h"
double shortest_distance(int point_1, int point_2, double **cost, int point_number);
int** PS_BIH(Instance* inslist, int ins_size, int pop_size, int fileInd);
#endif //RESEARCH_INITIALIZATION_H