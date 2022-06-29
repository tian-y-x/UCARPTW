//
// Created by lenovo on 2021/9/14.
//

#ifndef RESEARCH_SLS_H
#define RESEARCH_SLS_H
#include <iostream>
#include <vector>
#include "instance.h"

using namespace std;

void s_insert(vector<int> &a, int ith, int target_pos);

void d_insert(vector<int> &a, int ith, int jth);

void exchange(vector<int> &a, int ith, int jth);

void reverse(vector<int> &a, int i);

void d_reverse(vector<int> &a, int i, int j);

vector<int> SI(vector<int> &c_ori, double **ehm, double standard, Instance *inslist);

vector<int> DI(vector<int> &c_ori, double **ehm, double standard, Instance *inslist);

vector<int> Swap(vector<int> &c_ori, double **ehm, double standard, Instance *inslist);

vector<int> Two_opt(vector<int> &c_ori, double **ehm, double standard, Instance *inslist);

double StochasticSLS(vector<int> &c, double **ehm, Instance* inslist);

#endif //RESEARCH_SLS_H
