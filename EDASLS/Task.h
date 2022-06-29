//
// Created by lenovo on 2021/9/3.
//

#ifndef RESEARCH_TASK_H
#define RESEARCH_TASK_H


class Task {
public:
    double cost;
    double demand;
    bool enable= true;
    int task_number = -1;
    int start_point = -1;
    int end_point = -1;
    int an=0;
    int bn=100000000;

    Task();
    Task(int task_number,double cost,double demand,int start_point,int end_point);
    Task(int task_number,double cost,double demand,int start_point,int end_point,int an,int bn);
};


#endif //RESEARCH_TASK_H
