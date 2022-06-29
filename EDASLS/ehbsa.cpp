#include <iostream>
#include <ctime>

#include "ehbsa.h"

//len is the number of task, the size of ehm is (2*len) * (2*len)

using namespace std;

int times=0;

double **build_ehm(int **pop, int pop_size, int len, double epsilon) {//epsilon should be proportional
    double** ehm = new double* [2 * len];
    for (int i = 0; i < 2 * len; ++i) {
        ehm[i]=new double [2 * len];
        for (int j = 0; j < 2 * len; ++j) {
            ehm[i][j] = epsilon;
        }
    }
    for (int i = 0; i < pop_size; ++i) {
        ehm[pop[i][len - 1]][pop[i][0]]++;
        ehm[(pop[i][0] + len) % (2 * len)][(pop[i][len - 1] + len) % (2 * len)]++;
        for (int j = 1; j < len; ++j) {
            ehm[pop[i][j - 1]][pop[i][j]]++;
            ehm[(pop[i][j] + len) % (2 * len)][(pop[i][j - 1] + len) % (2 * len)]++;
        }
    }
    return ehm;
}

void ring(int len, int slot, double *rw, bool *is_waiting, double **ehm, int* res){
    double sum = 0;
    for (int j = 0; j < 2 * len; ++j) {
        rw[j] = is_waiting[j % len] ? ehm[res[(slot + len - 1) % len]][j] : 0;
        sum += rw[j];
    }
    rw[0] /= sum;
    for (int j = 1; j < 2 * len; ++j) {
        rw[j] /= sum;
        rw[j] += rw[j - 1];
    }
    double random = (rand() % 9999 + 1) / (double) 10000;
    for (int j = 0; j < 2 * len; ++j) {
        if (random <= rw[j]) {
            res[slot] = j;
            is_waiting[j % len] = false;
            break;
        }
    }
}

int *EHBSA(const int* T, double **ehm, int n, int len) {// n is the amount of divided parts
    //chech_chromosome(T);
    if (n < 1 || n > len) {
        cerr << "error in range" << endl;
        exit(1);
    }
    int start_pos;
    int end_pos;
    if (n==2){
        start_pos=rand() % len;
        end_pos = (start_pos + rand() % (len-1) + 1) % len;
    }else{
        int* divide_plan=new int[n];// divide T into n segments
        for (int i = 0; i < n; ++i) {
            divide_plan[i] = 1;
        }
        int left_point = len - n;
        while (left_point > 0) {
            divide_plan[rand() % n]++;
            left_point--;
        }
        int pick = rand() % n;// pick one segment
        int bias = rand() % len;
        start_pos = bias; // where to start in T
        end_pos = (bias + divide_plan[pick] - 1) % len;// where to end in T
        delete[] divide_plan;
    }
    times++;

    double* rw = new double[2 * len];
    bool is_waiting[len];
    int* res = new int[len];
    for (int i = 0; i < len; ++i) {
        res[i] = T[i];
    }
    if (start_pos<=end_pos){
        for (int i = 0; i < len; ++i) {
            is_waiting[T[i] % len] = i >= start_pos && i <= end_pos;
        }
        for (int i = start_pos; i <= end_pos; ++i) {
            ring(len, i, rw, is_waiting, ehm, res);
        }
    }else{
        for (int i = 0; i < len; ++i) {
            is_waiting[T[i] % len] = i >= start_pos || i <= end_pos;
        }
        for (int i = start_pos; i < len; ++i) {
            ring(len, i, rw, is_waiting, ehm, res);
        }
        for (int i = 0; i <= end_pos; ++i) {
            ring(len, i, rw, is_waiting, ehm, res);
        }
    }
    delete[] rw;
    return res;
}