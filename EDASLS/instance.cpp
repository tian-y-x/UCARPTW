#include "instance.h"
Instance::Instance() = default;

Instance::Instance(vector<Task> tasklist, double** cost, double** demand, int point_number, int UT, double capacity,int cost_mul,int demand_mul)
{
    this->tasklist = tasklist;
    this->cost = cost;
    this->demand = demand;
    this->point_number = point_number;
    this->UT = UT;
    this->capacity = capacity;
    this->map = NULL;
    this->e_map = NULL;
    this->t_map = NULL;
    this->demand_mul=demand_mul;
    this->cost_mul=cost_mul;
}
Instance& Instance:: operator =(Instance& a)
{
    this->tasklist = a.tasklist;
    this->cost = a.cost;
    this->demand = a.demand;
    this->point_number = a.point_number;
    this->UT = a.UT;
    this->capacity = a.capacity;
    this->map = a.map;
    this->e_map=a.e_map;
    this->t_map=a.t_map;
    this->cost_mul=a.cost_mul;
    this->demand_mul=a.demand_mul;
    return *this;
}
Instance::~Instance()
{
    for (int i = 0; i < point_number; ++i)
    {
        delete[] cost[i];
        delete[] demand[i];
    }
    if (map != NULL)
    {
        for (int i = 0; i < point_number; ++i)
            delete[] map[i];
        delete[] map;
    }

    if (e_map != NULL)
    {
        for (int i = 0; i < tasklist.size(); ++i)
            delete[] e_map[i];
        delete[] e_map;
    }

    if (t_map != NULL)
    {
        for (int i = 0; i < tasklist.size(); ++i)
            delete[] t_map[i];
        delete[] t_map;
    }
    delete[] cost;
    delete[] demand;

}
