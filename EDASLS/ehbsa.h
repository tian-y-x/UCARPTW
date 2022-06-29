//
// Created by lenovo on 2021/9/11.
//

#ifndef RESEARCH_EHBSA_H
#define RESEARCH_EHBSA_H

#include <vector>

using namespace std;

double **build_ehm(int **pop, int pop_size, int len, double epsilon);
int *EHBSA(const int* T, double **ehm, int n,int len);

#endif //RESEARCH_EHBSA_H
