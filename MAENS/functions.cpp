#include <stdio.h>
#include "functions.h"

void mod_dijkstra()
{
    int i, j, k, m;
    double minimum;
    printf("Dijastra\n");
    for (i = 1; i <= vertex_num; i++)
    {
        for (j = 1; j <= vertex_num; j++)
        {
            if (j == i)
                continue;

            shortest_path[i][j][0] = 1;
            shortest_path[i][j][1] = i;
            min_cost[i][j] = INF;
        }
    }

    int mark[MAX_NODE_TAG_LENGTH], nearest_neighbor[MAX_NODE_TAG_LENGTH];
    double dist[MAX_NODE_TAG_LENGTH], dist1[MAX_NODE_TAG_LENGTH];

    for (i = 1; i <= vertex_num; i++)
    {
        mark[i] = 1;

        for (j = 1; j <= vertex_num; j++)
        {
            if (j == i)
                continue;

            mark[j] = 0;
            dist[j] = trav_cost[i][j];
            dist1[j] = dist[j];
        }

        for (k = 1; k < vertex_num; k++)
        {
            minimum = INF;
            nearest_neighbor[0] = 0;

            for (j = 1; j <= vertex_num; j++)
            {
                if (mark[j])
                    continue;

                if (IS_DOUBLE_ZERO(dist1[j]-INF) == 1)
                    continue;

                if (dist1[j] < minimum)
                    minimum = dist1[j];
            }

            if (IS_DOUBLE_ZERO(minimum-INF) == 1)
                continue;

            for (j = 1; j <= vertex_num; j++)
            {
                if (mark[j])
                    continue;

                if (dist1[j] == minimum)
                {
                    nearest_neighbor[0] ++;
                    nearest_neighbor[nearest_neighbor[0]] = j;
                }
            }

            int v = nearest_neighbor[1];
            dist1[v] = INF;
            mark[v] = 1;

            if (shortest_path[i][v][0] == 0 || (shortest_path[i][v][0] > 0 && shortest_path[i][v][(int)shortest_path[i][v][0]] != v))
            {
                shortest_path[i][v][0] ++;
                shortest_path[i][v][(int)shortest_path[i][v][0]] = v;
            }

            for (j = 1; j <= vertex_num; j++)
            {
                if (mark[j])
                    continue;

                if (minimum+trav_cost[v][j] < dist[j])
                {
                    dist[j] = minimum+trav_cost[v][j];
                    dist1[j] = minimum+trav_cost[v][j];
                    for (m = 0; m <= shortest_path[i][v][0]; m++)
                    {
                        shortest_path[i][j][m] = shortest_path[i][v][m];
                    }
                }
            }

            for (j = 1; j <= vertex_num; j++)
            {
                if (j == i)
                    continue;

                min_cost[i][j] = dist[j];
            }
        }
    }

    for (i = 1; i <= vertex_num; i++)
    {
        for (j = 1; j <= vertex_num; j++)
        {
            if (shortest_path[i][j][0] == 1)
                shortest_path[i][j][0] = 0;
        }
    }

    for (i = 1; i <= vertex_num; i++)
    {
        shortest_path[i][i][0] = 1;
        shortest_path[i][i][1] = i;
        min_cost[i][i] = 0;
    }
}


void path_scanning(Individual* ps_indi, const Task* inst_tasks, int* serve_mark)
{

    // min_cost, NRE, NRA, NVeh, capacity, is the extern variables.
    int i, j, k;
    int serve_task_num = 0;
    for (i = req_edge_num + 1; i <= task_num; i++)
    {
        if (serve_mark[i])
        {
             serve_task_num++;
        }
    }
    int trial;
    double load, mindist;
    int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
    int nearest_isol_task[MAX_TASK_TAG_LENGTH], nearest_inci_task[MAX_TASK_TAG_LENGTH], sel_task[MAX_TASK_TAG_LENGTH];
    int current_task, next_task;

    int positions[MAX_TASK_SEQ_LENGTH];

    ps_indi->TotalCost = INF;
    ps_indi->TotalCost = INF;
    Individual tmp_indi1, tmp_indi2, tmp_indi3, tmp_indi4, tmp_indi5;

    double dep_dist[MAX_TASK_TAG_LENGTH], max_dep_dist, min_dep_dist;
    double yield[MAX_TASK_TAG_LENGTH], max_yield, min_yield;

    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;

        dep_dist[i] = min_cost[inst_tasks[i].tail_node][DEPOT];
        yield[i] = 1.0 * inst_tasks[i].demand / inst_tasks[i].serv_cost;
    }

    /*
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * Use Rule 1 to obtain a solution
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * */
    tmp_indi1.Sequence[0] = 1;
    tmp_indi1.Sequence[1] = 0;
    tmp_indi1.Loads[0] = 0;

    unserved_task[0] = 0;
    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;
        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = tmp_indi1.Sequence[tmp_indi1.Sequence[0]]; // get the current task's id
        candi_task[0] = 0;

        for (i = 1; i <= unserved_task[0]; i++)
        {
            if (inst_tasks[unserved_task[i]].demand <= capacity - load)
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            tmp_indi1.Sequence[0] ++;
            tmp_indi1.Sequence[tmp_indi1.Sequence[0]] = 0;
            tmp_indi1.Loads[0] ++;
            tmp_indi1.Loads[(int)tmp_indi1.Loads[0]] = load;
            load = 0;
            continue;
        }

        mindist = INF;
        nearest_task[0] = 0;
        for (i = 1; i <= candi_task[0]; i++)
        {
            if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
            {
                mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
                nearest_task[0] = 1;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
            else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        nearest_inci_task[0] = 0;
        nearest_isol_task[0] = 0;
        for (i = 1; i <= nearest_task[0]; i++)
        {
            if (inst_tasks[nearest_task[i]].tail_node == 1)
            {
                nearest_inci_task[0] ++;
                nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
            }
            else
            {
                nearest_isol_task[0] ++;
                nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
            }
        }
        if (nearest_isol_task[0] == 0)
        {
            memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0] + 1) * sizeof(int));
        }

        // for 5 five phase, the above part is the same
        max_dep_dist = -1;
        sel_task[0] = 0;
        for (i = 1; i <= nearest_isol_task[0]; i++)
        {
            if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
            {
                max_dep_dist = dep_dist[nearest_isol_task[i]];
                sel_task[0] = 1;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
            else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
            {
                sel_task[0] ++;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
        }
        k = 1;
        next_task = sel_task[k];

        trial++;
        tmp_indi1.Sequence[0]++;
        tmp_indi1.Sequence[tmp_indi1.Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;

        // delete the served task in unserved_task array
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    tmp_indi1.Sequence[0] ++;
    tmp_indi1.Sequence[tmp_indi1.Sequence[0]] = 0;
    tmp_indi1.Loads[0] ++;
    tmp_indi1.Loads[(int)tmp_indi1.Loads[0]] = load;

    tmp_indi1.TotalCost = get_task_seq_total_cost(tmp_indi1.Sequence, inst_tasks);
    //    tmp_indi1.TotalVioLoad = get_total_vio_load(tmp_indi1.Loads);

    if (tmp_indi1.TotalCost < ps_indi->TotalCost)
    {
        memcpy(ps_indi->Sequence, tmp_indi1.Sequence, (tmp_indi1.Sequence[0] + 1) * sizeof(int));
        memcpy(ps_indi->Loads, tmp_indi1.Loads, ((int)tmp_indi1.Loads[0] + 1) * sizeof(double));
        ps_indi->TotalCost = tmp_indi1.TotalCost;
    }

    /*
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * Use Rule 2 to obtain a solution
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * */
    tmp_indi2.Sequence[0] = 1;
    tmp_indi2.Sequence[1] = 0;
    tmp_indi2.Loads[0] = 0;

    unserved_task[0] = 0;
    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;
        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = tmp_indi2.Sequence[tmp_indi2.Sequence[0]]; // get the current task's id
        candi_task[0] = 0;

        for (i = 1; i <= unserved_task[0]; i++)
        {
            if (inst_tasks[unserved_task[i]].demand <= capacity - load)
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            tmp_indi2.Sequence[0] ++;
            tmp_indi2.Sequence[tmp_indi2.Sequence[0]] = 0;
            tmp_indi2.Loads[0] ++;
            tmp_indi2.Loads[(int)tmp_indi2.Loads[0]] = load;
            load = 0;
            continue;
        }

        mindist = INF;
        nearest_task[0] = 0;
        for (i = 1; i <= candi_task[0]; i++)
        {
            if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
            {
                mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
                nearest_task[0] = 1;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
            else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        nearest_inci_task[0] = 0;
        nearest_isol_task[0] = 0;
        for (i = 1; i <= nearest_task[0]; i++)
        {
            if (inst_tasks[nearest_task[i]].tail_node == 1)
            {
                nearest_inci_task[0] ++;
                nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
            }
            else
            {
                nearest_isol_task[0] ++;
                nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
            }
        }
        if (nearest_isol_task[0] == 0)
        {
            memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0] + 1) * sizeof(int));
        }

        // for 5 five phase, the above part is the same
        min_dep_dist = INF;
        sel_task[0] = 0;
        for (i = 1; i <= nearest_isol_task[0]; i++)
        {
            if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
            {
                min_dep_dist = dep_dist[nearest_isol_task[i]];
                sel_task[0] = 1;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
            else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
            {
                sel_task[0] ++;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
        }


        k = 1;
        next_task = sel_task[k];

        trial++;
        tmp_indi2.Sequence[0]++;
        tmp_indi2.Sequence[tmp_indi2.Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;

        // delete the served task in unserved_task array
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    tmp_indi2.Sequence[0] ++;
    tmp_indi2.Sequence[tmp_indi2.Sequence[0]] = 0;
    tmp_indi2.Loads[0] ++;
    tmp_indi2.Loads[(int)tmp_indi2.Loads[0]] = load;

    tmp_indi2.TotalCost = get_task_seq_total_cost(tmp_indi2.Sequence, inst_tasks);
    //    tmp_indi2.TotalVioLoad = get_total_vio_load(tmp_indi2.Loads);

    if (tmp_indi2.TotalCost < ps_indi->TotalCost)
    {
        memcpy(ps_indi->Sequence, tmp_indi2.Sequence, (tmp_indi2.Sequence[0] + 1) * sizeof(int));
        memcpy(ps_indi->Loads, tmp_indi2.Loads, ((int)tmp_indi2.Loads[0] + 1) * sizeof(double));
        ps_indi->TotalCost = tmp_indi2.TotalCost;
    }

    /*
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * Use Rule 3 to obtain a solution
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * */
    tmp_indi3.Sequence[0] = 1;
    tmp_indi3.Sequence[1] = 0;
    tmp_indi3.Loads[0] = 0;

    unserved_task[0] = 0;
    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;
        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = tmp_indi3.Sequence[tmp_indi3.Sequence[0]]; // get the current task's id
        candi_task[0] = 0;

        for (i = 1; i <= unserved_task[0]; i++)
        {
            if (inst_tasks[unserved_task[i]].demand <= capacity - load)
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            tmp_indi3.Sequence[0] ++;
            tmp_indi3.Sequence[tmp_indi3.Sequence[0]] = 0;
            tmp_indi3.Loads[0] ++;
            tmp_indi3.Loads[(int)tmp_indi3.Loads[0]] = load;
            load = 0;
            continue;
        }

        mindist = INF;
        nearest_task[0] = 0;
        for (i = 1; i <= candi_task[0]; i++)
        {
            if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
            {
                mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
                nearest_task[0] = 1;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
            else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        nearest_inci_task[0] = 0;
        nearest_isol_task[0] = 0;
        for (i = 1; i <= nearest_task[0]; i++)
        {
            if (inst_tasks[nearest_task[i]].tail_node == 1)
            {
                nearest_inci_task[0] ++;
                nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
            }
            else
            {
                nearest_isol_task[0] ++;
                nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
            }
        }
        if (nearest_isol_task[0] == 0)
        {
            memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0] + 1) * sizeof(int));
        }

        // for 5 five phase, the above part is the same
        max_yield = -1;
        sel_task[0] = 0;
        for (i = 1; i <= nearest_isol_task[0]; i++)
        {
            if (yield[nearest_isol_task[i]] > max_yield)
            {
                max_yield = yield[nearest_isol_task[i]];
                sel_task[0] = 1;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
            else if (yield[nearest_isol_task[i]] == max_yield)
            {
                sel_task[0] ++;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
        }


        k = 1;
        next_task = sel_task[k];

        trial++;
        tmp_indi3.Sequence[0]++;
        tmp_indi3.Sequence[tmp_indi3.Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;

        // delete the served task in unserved_task array
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    tmp_indi3.Sequence[0] ++;
    tmp_indi3.Sequence[tmp_indi3.Sequence[0]] = 0;
    tmp_indi3.Loads[0] ++;
    tmp_indi3.Loads[(int)tmp_indi3.Loads[0]] = load;

    tmp_indi3.TotalCost = get_task_seq_total_cost(tmp_indi3.Sequence, inst_tasks);
    //    tmp_indi3.TotalVioLoad = get_total_vio_load(tmp_indi3.Loads);

    if (tmp_indi3.TotalCost < ps_indi->TotalCost)
    {
        memcpy(ps_indi->Sequence, tmp_indi3.Sequence, (tmp_indi3.Sequence[0] + 1) * sizeof(int));
        memcpy(ps_indi->Loads, tmp_indi3.Loads, ((int)tmp_indi3.Loads[0] + 1) * sizeof(double));
        ps_indi->TotalCost = tmp_indi3.TotalCost;
    }

    /*
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * Use Rule 4 to obtain a solution
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * */
    tmp_indi4.Sequence[0] = 1;
    tmp_indi4.Sequence[1] = 0;
    tmp_indi4.Loads[0] = 0;

    unserved_task[0] = 0;
    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;
        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = tmp_indi4.Sequence[tmp_indi4.Sequence[0]]; // get the current task's id
        candi_task[0] = 0;

        for (i = 1; i <= unserved_task[0]; i++)
        {
            if (inst_tasks[unserved_task[i]].demand <= capacity - load)
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            tmp_indi4.Sequence[0] ++;
            tmp_indi4.Sequence[tmp_indi4.Sequence[0]] = 0;
            tmp_indi4.Loads[0] ++;
            tmp_indi4.Loads[(int)tmp_indi4.Loads[0]] = load;
            load = 0;
            continue;
        }

        mindist = INF;
        nearest_task[0] = 0;
        for (i = 1; i <= candi_task[0]; i++)
        {
            if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
            {
                mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
                nearest_task[0] = 1;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
            else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        nearest_inci_task[0] = 0;
        nearest_isol_task[0] = 0;
        for (i = 1; i <= nearest_task[0]; i++)
        {
            if (inst_tasks[nearest_task[i]].tail_node == 1)
            {
                nearest_inci_task[0] ++;
                nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
            }
            else
            {
                nearest_isol_task[0] ++;
                nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
            }
        }
        if (nearest_isol_task[0] == 0)
        {
            memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0] + 1) * sizeof(int));
        }

        // for 5 five phase, the above part is the same
        min_yield = INF;
        sel_task[0] = 0;
        for (i = 1; i <= nearest_isol_task[0]; i++)
        {
            if (yield[nearest_isol_task[i]] < min_yield)
            {
                min_yield = yield[nearest_isol_task[i]];
                sel_task[0] = 1;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
            else if (yield[nearest_isol_task[i]] == min_yield)
            {
                sel_task[0] ++;
                sel_task[sel_task[0]] = nearest_isol_task[i];
            }
        }


        k = 1;
        next_task = sel_task[k];

        trial++;
        tmp_indi4.Sequence[0]++;
        tmp_indi4.Sequence[tmp_indi4.Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;

        // delete the served task in unserved_task array
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    tmp_indi4.Sequence[0] ++;
    tmp_indi4.Sequence[tmp_indi4.Sequence[0]] = 0;
    tmp_indi4.Loads[0] ++;
    tmp_indi4.Loads[(int)tmp_indi4.Loads[0]] = load;

    tmp_indi4.TotalCost = get_task_seq_total_cost(tmp_indi4.Sequence, inst_tasks);
    //    tmp_indi4.TotalVioLoad = get_total_vio_load(tmp_indi4.Loads);

    if (tmp_indi4.TotalCost < ps_indi->TotalCost)
    {
        memcpy(ps_indi->Sequence, tmp_indi4.Sequence, (tmp_indi4.Sequence[0] + 1) * sizeof(int));
        memcpy(ps_indi->Loads, tmp_indi4.Loads, ((int)tmp_indi4.Loads[0] + 1) * sizeof(double));
        ps_indi->TotalCost = tmp_indi4.TotalCost;
    }

    /*
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * Use Rule 5 to obtain a solution
     * / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
     * */
    tmp_indi5.Sequence[0] = 1;
    tmp_indi5.Sequence[1] = 0;
    tmp_indi5.Loads[0] = 0;

    unserved_task[0] = 0;
    for (i = 1; i <= task_num; i++)
    {
        if (!serve_mark[i])
            continue;
        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = tmp_indi5.Sequence[tmp_indi5.Sequence[0]]; // get the current task's id
        candi_task[0] = 0;

        for (i = 1; i <= unserved_task[0]; i++)
        {
            if (inst_tasks[unserved_task[i]].demand <= capacity - load)
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            tmp_indi5.Sequence[0] ++;
            tmp_indi5.Sequence[tmp_indi5.Sequence[0]] = 0;
            tmp_indi5.Loads[0] ++;
            tmp_indi5.Loads[(int)tmp_indi5.Loads[0]] = load;
            load = 0;
            continue;
        }

        mindist = INF;
        nearest_task[0] = 0;
        for (i = 1; i <= candi_task[0]; i++)
        {
            if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
            {
                mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
                nearest_task[0] = 1;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
            else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        nearest_inci_task[0] = 0;
        nearest_isol_task[0] = 0;
        for (i = 1; i <= nearest_task[0]; i++)
        {
            if (inst_tasks[nearest_task[i]].tail_node == 1)
            {
                nearest_inci_task[0] ++;
                nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
            }
            else
            {
                nearest_isol_task[0] ++;
                nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
            }
        }
        if (nearest_isol_task[0] == 0)
        {
            memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0] + 1) * sizeof(int));
        }

        // for 5 five phase, the above part is the same
        if (load < capacity / 2)
        {
            max_dep_dist = -1;
            sel_task[0] = 0;
            for (i = 1; i <= nearest_isol_task[0]; i++)
            {
                if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
                {
                    max_dep_dist = dep_dist[nearest_isol_task[i]];
                    sel_task[0] = 1;
                    sel_task[sel_task[0]] = nearest_isol_task[i];
                }
                else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
                {
                    sel_task[0] ++;
                    sel_task[sel_task[0]] = nearest_isol_task[i];
                }
            }
        }
        else
        {
            min_dep_dist = INF;
            sel_task[0] = 0;
            for (i = 1; i <= nearest_isol_task[0]; i++)
            {
                if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
                {
                    min_dep_dist = dep_dist[nearest_isol_task[i]];
                    sel_task[0] = 1;
                    sel_task[sel_task[0]] = nearest_isol_task[i];
                }
                else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
                {
                    sel_task[0] ++;
                    sel_task[sel_task[0]] = nearest_isol_task[i];
                }
            }
        }


        k = 1;
        next_task = sel_task[k];

        trial++;
        tmp_indi5.Sequence[0]++;
        tmp_indi5.Sequence[tmp_indi5.Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;

        // delete the served task in unserved_task array
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    tmp_indi5.Sequence[0] ++;
    tmp_indi5.Sequence[tmp_indi5.Sequence[0]] = 0;
    tmp_indi5.Loads[0] ++;
    tmp_indi5.Loads[(int)tmp_indi5.Loads[0]] = load;

    tmp_indi5.TotalCost = get_task_seq_total_cost(tmp_indi5.Sequence, inst_tasks);
    //    tmp_indi5.TotalVioLoad = get_total_vio_load(tmp_indi5.Loads);

    if (tmp_indi5.TotalCost < ps_indi->TotalCost)
    {
        memcpy(ps_indi->Sequence, tmp_indi5.Sequence, (tmp_indi5.Sequence[0] + 1) * sizeof(int));
        memcpy(ps_indi->Loads, tmp_indi5.Loads, ((int)tmp_indi5.Loads[0] + 1) * sizeof(double));
        ps_indi->TotalCost = tmp_indi5.TotalCost;
    }

    ps_indi->TotalVioLoad = 0;
}