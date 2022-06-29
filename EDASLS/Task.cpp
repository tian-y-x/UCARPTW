//
// Created by lenovo on 2021/9/3.
//

#include "Task.h"

Task::Task(int task_number,double cost,double demand,int start_point,int end_point) {
    this->task_number=task_number;
    this->cost=cost;
    this->demand=demand;
    this->start_point=start_point;
    this->end_point=end_point;
}

Task::Task(int task_number, double cost, double demand, int start_point, int end_point, int an, int bn) {
    this->task_number=task_number;
    this->cost=cost;
    this->demand=demand;
    this->start_point=start_point;
    this->end_point=end_point;
    this->an=an;
    this->bn=bn;
}

Task::Task() = default;
