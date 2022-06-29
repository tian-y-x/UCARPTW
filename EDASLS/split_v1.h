#pragma once
#include <vector>
#include "instance.h"
using namespace std;

double CEvaluate(const Instance *instances, const vector<int> &chrom);
double split(const vector<int> &chrom, const Instance &instance, int *split_task_seq);
double penalize(Task tj, double time);
