#include "functions.h"


void ins30::mod_dijkstra()
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

            this->shortest_path[i][j][0] = 1;
            this->shortest_path[i][j][1] = i;
            this->min_cost[i][j] = INF;
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

                if (IS_DOUBLE_ZERO(dist1[j] - INF) == 1)
                    continue;

                if (dist1[j] < minimum)
                    minimum = dist1[j];
            }

            if (IS_DOUBLE_ZERO(minimum - INF) == 1)
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

            if (this->shortest_path[i][v][0] == 0 || (this->shortest_path[i][v][0] > 0 && this->shortest_path[i][v][(int)this->shortest_path[i][v][0]] != v))
            {
                this->shortest_path[i][v][0] ++;
                this->shortest_path[i][v][(int)this->shortest_path[i][v][0]] = v;
            }

            for (j = 1; j <= vertex_num; j++)
            {
                if (mark[j])
                    continue;

                if (minimum + trav_cost[v][j] < dist[j])
                {
                    dist[j] = minimum + this->trav_cost[v][j];
                    dist1[j] = minimum + this->trav_cost[v][j];
                    for (m = 0; m <= this->shortest_path[i][v][0]; m++)
                    {
                        this->shortest_path[i][j][m] = this->shortest_path[i][v][m];
                    }
                }
            }

            for (j = 1; j <= vertex_num; j++)
            {
                if (j == i)
                    continue;

                this->min_cost[i][j] = dist[j];
            }
        }
    }

    for (i = 1; i <= vertex_num; i++)
    {
        for (j = 1; j <= vertex_num; j++)
        {
            if (this->shortest_path[i][j][0] == 1)
                this->shortest_path[i][j][0] = 0;
        }
    }

    for (i = 1; i <= vertex_num; i++)
    {
        this->shortest_path[i][i][0] = 1;
        this->shortest_path[i][i][1] = i;
        this->min_cost[i][i] = 0;
    }
}
void ins30::init_ins()
{
    for (int i = 1; i <= vertex_num; i++)
    {
        for (int j = 1; j <= vertex_num; j++)
        {
            this->trav_cost[i][j] = INF;
            this->serve_cost[i][j] = 0;
        }
    }
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            this->min_cost[i][j] = -1;
        }
    }
    this->trav_cost[1][0] = 0;
    this->trav_cost[0][1] = 0;

    for (int i = 1; i <= task_num; i++)
    {
        this->serve_cost[this->inst_tasks[i].head_node][this->inst_tasks[i].tail_node] = this->inst_tasks[i].serv_cost;
    }

    for (int i = 1; i <= total_arc_num; i++)
    {
        this->trav_cost[this->inst_arcs[i].head_node][this->inst_arcs[i].tail_node] = this->inst_arcs[i].trav_cost;
    }
}
