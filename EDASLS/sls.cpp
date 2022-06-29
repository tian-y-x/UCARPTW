//
// Created by lenovo on 2021/8/19.
//
#include <iostream>
#include <vector>
#include "sls.h"
#include "split_v1.h"
#include "instance.h"

const int INSTANCE_NUM = 30;
const int MAX_TASK_TAG_LENGTH = 500;
extern int routes[INSTANCE_NUM][MAX_TASK_TAG_LENGTH];
extern int evaluate_times;
//extern int sls_times;
int FEm = 300000;

using namespace std;

double StochasticSLS(vector<int> &c, double **ehm, Instance *inslist) {
//    sls_times++;
    double standard = CEvaluate(inslist, c);
    double base;
    do {
        base = standard;
        if (evaluate_times>=FEm-4){
            return standard;
        }
        vector<int> si(SI(c, ehm, standard, inslist));
        double fitness_si = CEvaluate(inslist, si);
        if (fitness_si > standard) {
            cout << "bug in si" << endl;
        }

        vector<int> di(DI(c, ehm, standard, inslist));
        double fitness_di = CEvaluate(inslist, di);
        if (fitness_di > standard) {
            cout << "bug in di" << endl;
        }

        vector<int> swap(Swap(c, ehm, standard, inslist));
        double fitness_swap = CEvaluate(inslist, swap);
        if (fitness_swap > standard) {
            cout << "bug in swap" << endl;
        }

        vector<int> topt(Two_opt(c, ehm, standard, inslist));
        double fitness_topt = CEvaluate(inslist, topt);
        if (fitness_topt > standard) {
            cout << "bug in two opt" << endl;
        }

        if (fitness_di < standard) {
            standard = fitness_di;
            c = di;
        }
        if (fitness_swap < standard) {
            standard = fitness_swap;
            c = swap;
        }
        if (fitness_topt < standard) {
            standard = fitness_topt;
            c = topt;
        }
        if (fitness_si < standard) {
            standard = fitness_si;
            c = si;
        }
    } while (standard < base);
    return standard;
}

vector<int> Two_opt(vector<int> &c_ori, double **ehm, double standard, Instance *inslist) {
    vector<int> c(c_ori);
    int len = c.size();
    for (int i = 0; i < len - 1; ++i) {
        for (int j = i + 1; j < len; ++j) {
            int inv_ti = (c[i] + len) % (2 * len);
            int inv_tj = (c[j] + len) % (2 * len);
            double old_ehm = 0;
            double new_ehm = 0;
            if (i != 0) {
                old_ehm += ehm[c[i - 1]][c[i]];
                new_ehm += ehm[c[i - 1]][inv_tj];
            }
            if (j != len - 1) {
                old_ehm += ehm[c[j]][c[j + 1]];
                new_ehm += ehm[inv_ti][c[j + 1]];
            }
            if (new_ehm > old_ehm) {
                int a = i, b = j;
                while (a < b) {
                    exchange(c, a, b);
                    reverse(c, a);
                    reverse(c, b);
                    a++;
                    b--;
                }
                if (a == b)
                    reverse(c, a);
                if (CEvaluate(inslist, c) < standard) {
                    return c;
                }
                a = i, b = j;
                while (a < b) {
                    exchange(c, a, b);
                    reverse(c, a);
                    reverse(c, b);
                    a++;
                    b--;
                }
                if (a == b)
                    reverse(c, a);
            }
            if (evaluate_times>=FEm){
                return c;
            }
        }
    }
    return c;
}

vector<int> Swap(vector<int> &c_ori, double **ehm, double standard, Instance *inslist) {
    vector<int> c(c_ori);
    int len = c.size();
    for (int i = 0; i < len - 1; ++i) {//swap i and j
        for (int j = i + 1; j < len; ++j) {
            int inv_ti = (c[i] + len) % (2 * len);
            int inv_tj = (c[j] + len) % (2 * len);
            double new_ehm, old_ehm, inv_i_ehm, inv_j_ehm, inv_both_ehm;
            if (i == 0 && j == len - 1) {
                new_ehm = ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][c[i]];
                old_ehm = +ehm[c[i]][c[i + 1]] + ehm[c[j - 1]][c[j]];
                inv_i_ehm = +ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][inv_ti];
                inv_j_ehm = ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][c[i]];
                inv_both_ehm = ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][inv_ti];
            } else if (i == 0) {
                new_ehm = ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][c[i]]
                          + ehm[c[i]][c[j + 1]];
                old_ehm = ehm[c[i]][c[i + 1]] + ehm[c[j - 1]][c[j]]
                          + ehm[c[j]][c[j + 1]];
                inv_i_ehm = ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][inv_ti]
                            + ehm[inv_ti][c[j + 1]];
                inv_j_ehm = ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][c[i]]
                            + ehm[c[i]][c[j + 1]];
                inv_both_ehm = ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][inv_ti]
                               + ehm[inv_ti][c[j + 1]];
            } else if (j == len - 1) {
                new_ehm = ehm[c[i - 1]][c[j]] + ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][c[i]];

                old_ehm = ehm[c[i - 1]][c[i]] + ehm[c[i]][c[i + 1]] + ehm[c[j - 1]][c[j]];

                inv_i_ehm = ehm[c[i - 1]][c[j]] + ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][inv_ti];

                inv_j_ehm = ehm[c[i - 1]][inv_tj] + ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][c[i]];

                inv_both_ehm = ehm[c[i - 1]][inv_tj] + ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][inv_ti];

            } else {
                new_ehm = ehm[c[i - 1]][c[j]] + ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][c[i]]
                          + ehm[c[i]][c[j + 1]];
                old_ehm = ehm[c[i - 1]][c[i]] + ehm[c[i]][c[i + 1]] + ehm[c[j - 1]][c[j]]
                          + ehm[c[j]][c[j + 1]];
                inv_i_ehm = ehm[c[i - 1]][c[j]] + ehm[c[j]][c[i + 1]] + ehm[c[j - 1]][inv_ti]
                            + ehm[inv_ti][c[j + 1]];
                inv_j_ehm = ehm[c[i - 1]][inv_tj] + ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][c[i]]
                            + ehm[c[i]][c[j + 1]];
                inv_both_ehm = ehm[c[i - 1]][inv_tj] + ehm[inv_tj][c[i + 1]] + ehm[c[j - 1]][inv_ti]
                               + ehm[inv_ti][c[j + 1]];
            }

            if (new_ehm > old_ehm) {
                exchange(c, i, j);
                double temp = CEvaluate(inslist, c);// jth ith
                if (temp < standard) {
                    return c;
                } else {
                    exchange(c, i, j);
                }
            }
            if (evaluate_times>=FEm){
                return c;
            }
            if (inv_i_ehm > old_ehm) {
                reverse(c, i);
                exchange(c, i, j);
                double temp = CEvaluate(inslist, c);// jth ith
                if (temp < standard) {
                    return c;
                } else {
                    exchange(c, i, j);
                    reverse(c, i);
                }
            }
            if (evaluate_times>=FEm){
                return c;
            }
            if (inv_j_ehm > old_ehm) {
                reverse(c, j);
                exchange(c, i, j);
                double temp = CEvaluate(inslist, c);// jth ith
                if (temp < standard) {
                    return c;
                } else {
                    exchange(c, i, j);
                    reverse(c, j);
                }
            }
            if (evaluate_times>=FEm){
                return c;
            }
            if (inv_both_ehm > old_ehm) {
                reverse(c, i);
                reverse(c, j);
                exchange(c, i, j);
                double temp = CEvaluate(inslist, c);// jth ith
                if (temp < standard) {
                    return c;
                } else {
                    exchange(c, i, j);
                    reverse(c, j);
                    reverse(c, i);
                }
            }
            if (evaluate_times>=FEm){
                return c;
            }
        }
    }
    return c;
}

vector<int> DI(vector<int> &c_ori, double **ehm, double standard, Instance *inslist) {
    vector<int> c(c_ori);
    int len = c.size();
    for (int i = 0; i < len - 1; ++i) {//let the ith task be the task after the jth task
        for (int j = 0; j < len; ++j) {
            if (i - j >= 2 || j - i >= 2) {// They can not be adjacent.
                int inv_ti = (c[i] + len) % (2 * len);
                int inv_af_i = (c[i + 1] + len) % (2 * len);
                double new_ehm_value = ehm[c[j]][c[i]], old_ehm_value = 0, inv_ehm_value = ehm[c[j]][inv_af_i];
                if (i != 0) {
                    old_ehm_value += ehm[c[i - 1]][c[i]];
                    if (i != len - 2) {
                        new_ehm_value += ehm[c[i - 1]][c[i + 2]];
                        inv_ehm_value += ehm[c[i - 1]][c[i + 2]];
                    }
                }
                if (i != len - 2) {
                    old_ehm_value += ehm[c[i + 1]][c[i + 2]];
                }
                if (j != len - 1) {
                    new_ehm_value += ehm[c[i + 1]][c[j + 1]];
                    old_ehm_value += ehm[c[j]][c[j + 1]];
                    inv_ehm_value += ehm[inv_ti][c[j + 1]];
                }

//                double new_ehm_value = ehm[c[bf_i]][c[af_i_2]] + ehm[c[j]][c[i]] + ehm[c[i+1]][c[af_j]];
//                double old_ehm_value = ehm[c[bf_i]][c[i]] + ehm[c[i+1]][c[af_i_2]] + ehm[c[j]][c[af_j]];
//                double inv_ehm_value = ehm[c[bf_i]][c[af_i_2]] + ehm[c[j]][inv_af_i] + ehm[inv_ti][c[af_j]];


                if (new_ehm_value > old_ehm_value) {
                    if (i < j) {
                        s_insert(c, i, j);
                        s_insert(c, i, j);
                    } else {
                        s_insert(c, i + 1, j + 1);
                        s_insert(c, i + 1, j + 1);
                    }
                    double temp = CEvaluate(inslist, c);
                    if (temp < standard) {
                        return c;
                    }
                    if (i < j) {
                        s_insert(c, j, i);
                        s_insert(c, j, i);
                    } else {
                        s_insert(c, j + 1, i + 1);
                        s_insert(c, j + 1, i + 1);
                    }
                }
                if (evaluate_times>=FEm){
                    return c;
                }
                if (inv_ehm_value > old_ehm_value) {
                    d_reverse(c, i, i + 1);
                    if (i < j) {
                        s_insert(c, i, j);
                        s_insert(c, i, j);
                    } else {
                        s_insert(c, i + 1, j + 1);
                        s_insert(c, i + 1, j + 1);
                    }

                    double temp = CEvaluate(inslist, c);
                    if (temp < standard) {
                        return c;
                    }
                    if (i < j) {
                        s_insert(c, j, i);
                        s_insert(c, j, i);
                    } else {
                        s_insert(c, j + 1, i + 1);
                        s_insert(c, j + 1, i + 1);
                    }
                    d_reverse(c, i, i + 1);
                }
                if (evaluate_times>=FEm){
                    return c;
                }
            }
        }
    }
    return c;
}

vector<int> SI(vector<int> &c_ori, double **ehm, double standard, Instance *inslist) {
    double threshold = standard;
    vector<int> c(c_ori);
    int len = c.size();
    if (evaluate_times>=FEm){
        return c;
    }
    for (int i = 0; i < len; ++i) {//let the ith task be the task after the jth task
        for (int j = 0; j < len; ++j) {
            //If 2 tasks are adjacent, SI is meaningless,
            // cause Swap would task the responsibility
            if (i != j && i != j - 1 && i != j + 1) {
                int inv_ti = (c[i] + len) % (2 * len);
                double new_ehm_value = ehm[c[j]][c[i]], old_ehm_value = 0, inv_new_ehm_value = ehm[c[j]][inv_ti];
                if (i != 0) {
                    old_ehm_value += ehm[c[i - 1]][c[i]];
                    if (i != len - 1) {
                        new_ehm_value += ehm[c[i - 1]][c[i + 1]];
                        inv_new_ehm_value += ehm[c[i - 1]][c[i + 1]];
                    }
                }
                if (i != len - 1) {
                    old_ehm_value += ehm[c[i]][c[i + 1]];
                }
                if (j != len - 1) {
                    new_ehm_value += ehm[c[i]][c[j + 1]];
                    old_ehm_value += ehm[c[j]][c[j + 1]];
                    inv_new_ehm_value += ehm[inv_ti][c[j + 1]];
                }

//                double new_ehm_value = ehm[c[bf_i]][c[af_i]] + ehm[c[j]][c[i]] + ehm[c[i]][c[af_j]];
//                double old_ehm_value = ehm[c[bf_i]][c[i]] + ehm[c[i]][c[af_i]] + ehm[c[j]][c[af_j]];
//                double inv_new_ehm_value = ehm[c[bf_i]][c[af_i]] + ehm[c[j]][inv_ti]
//                                           + ehm[inv_ti][c[af_j]];
                if (new_ehm_value > old_ehm_value) {
                    if (i > j) {// j----i
                        s_insert(c, i, j + 1);// ji----
                    } else {// i----j
                        s_insert(c, i, j);// ----ji
                    }
                    double temp = CEvaluate(inslist, c);
                    if (temp < threshold) {
                        return c;
                    } else if (i > j) {
                        s_insert(c, j + 1, i);
                    } else {
                        s_insert(c, j, i);
                    }
                    if (evaluate_times>=FEm){
                        return c;
                    }
                }
                if (inv_new_ehm_value > old_ehm_value) {
                    reverse(c, i);// inv(i)----j
                    if (i > j) {
                        s_insert(c, i, j + 1);
                    } else {
                        s_insert(c, i, j);
                    }
                    double temp = CEvaluate(inslist, c);
                    if (temp < threshold) {
                        return c;
                    } else {
                        if (i > j) {// ji----i'
                            s_insert(c, j + 1, i);
                        } else { // i'---ji
                            s_insert(c, j, i);
                        }
                        reverse(c, i);
                    }
                    if (evaluate_times>=FEm){
                        return c;
                    }
                }
            }
        }
    }
    return c;
}

inline void s_insert(vector<int> &a, int ith, int target_pos) {
    int ith_task = a[ith];
    a.erase(a.cbegin() + ith);
    a.insert(a.cbegin() + target_pos, ith_task);
}

inline void d_insert(vector<int> &a, int ith, int jth) {
    if (ith < jth) {
        s_insert(a, ith, jth);
        s_insert(a, ith, jth);
    } else {
        s_insert(a, ith + 1, jth + 1);
        s_insert(a, ith + 1, jth + 1);
    }
}

inline void exchange(vector<int> &a, int ith, int jth) {
    int temp = a[ith];
    a[ith] = a[jth];
    a[jth] = temp;
}

inline void reverse(vector<int> &a, int i) {
    a[i] += a.size();
    a[i] %= 2 * a.size();
}

inline void d_reverse(vector<int> &a, int i, int j) {
    reverse(a, i);
    reverse(a, j);
    exchange(a, i, j);
}