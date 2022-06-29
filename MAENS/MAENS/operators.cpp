#include "MAENS.h"

int task_routes[50][MAX_TASK_SEG_LENGTH];

void rand_scanning(Individual *rs_indi, const Task *inst_tasks, const int *serve_mark)
{
    int i, k;
//    for(i =0 ;i<200;i++)
//        printf("head %d, end %d, serve_cost %d, dead_cost %d, demand %d, reverse %d\n",
//               inst_tasks[i].head_node,inst_tasks[i].tail_node,inst_tasks[i].serv_cost,inst_tasks[i].dead_cost,
//               inst_tasks[i].demand,inst_tasks[i].inverse);

    int serve_task_num = 0; // the number of tasks required to be served.
    for (i = req_edge_num+1; i <= task_num; i++)
    {
        if(serve_mark[i])
            serve_task_num++;
    }

    double load, mindist;
    int trial;

    int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
    int current_task, next_task;

    int positions[MAX_TASK_SEG_LENGTH];


    rs_indi->Sequence[0] = 1;
    rs_indi->Sequence[1] = 0;
    rs_indi->Loads[0] = 0;

    unserved_task[0] = 0;

    for (i = 1; i <= task_num; i++)
    {
        if ( !serve_mark[i] )
            continue;

        unserved_task[0] ++;
        unserved_task[unserved_task[0]] = i;
    }

    load = 0;
    trial = 0;
    while (trial < serve_task_num)
    {
        current_task = rs_indi->Sequence[rs_indi->Sequence[0]];

        candi_task[0] = 0;
        for( i = 1; i <= unserved_task[0]; i++)
        {
            if( inst_tasks[unserved_task[i]].demand <= capacity - load )
            {
                candi_task[0] ++;
                candi_task[candi_task[0]] = unserved_task[i];
            }
        }

        if (candi_task[0] == 0)
        {
            rs_indi->Sequence[0] ++;
            rs_indi->Sequence[rs_indi->Sequence[0]] = 0;
            rs_indi->Loads[0] ++;
            rs_indi->Loads[(int)rs_indi->Loads[0]] = load;
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
            }else if(min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
            {
                nearest_task[0] ++;
                nearest_task[nearest_task[0]] = candi_task[i];
            }
        }

        k = rand_choose(nearest_task[0]);
        next_task = nearest_task[k];

        trial ++;
        rs_indi->Sequence[0] ++;
        rs_indi->Sequence[rs_indi->Sequence[0]] = next_task;
        load += inst_tasks[next_task].demand;
        find_ele_positions(positions, unserved_task, next_task);
        delete_element(unserved_task, positions[1]);

        if (inst_tasks[next_task].inverse > 0)
        {
            find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
            delete_element(unserved_task, positions[1]);
        }
    }

    rs_indi->Sequence[0] ++;
    rs_indi->Sequence[rs_indi->Sequence[0]] = 0;
    rs_indi->Loads[0] ++;
    rs_indi->Loads[(int)rs_indi->Loads[0]] = load;

    rs_indi->TotalCost = get_total_cost(rs_indi->Sequence, inst_tasks);
    rs_indi->TotalVioLoad = get_total_vio_load(rs_indi->Loads);

}


void indi_copy(Individual *target, Individual *source)
{
    memcpy(target->Sequence, source->Sequence, (source->Sequence[0]+1)*sizeof(int));
    memcpy(target->Loads, source->Loads, ((int)source->Loads[0]+1)*sizeof(double));
    target->TotalCost = source->TotalCost;
    target->TotalVioLoad = source->TotalVioLoad;
    target->Fitness = source->Fitness;
//    memcpy(target->pred, source->pred, sizeof(source->pred));
//    memcpy(target->succ, source->succ, sizeof(source->succ));
}

void SBX(Individual *xed_child, Individual *p1, Individual *p2, const Task *inst_tasks)
// the SBX crossover operator
{
    /*
     * P1 -> R1, P2 -> R2
     * R1 -> {R11, R12}, R2 -> {R21, R22}
     * NewR = {R11, R22}
     * remove duplicate task in NewR and insert missing task to NewR -> offspring
     */
//    for(int i = 0;i<=p1->Sequence[0];i++)
//        printf("%d ",p1->Sequence[i]);
//    printf("\n");
//    for(int i = 0;i<=p2->Sequence[0];i++)
//        printf("%d ",p2->Sequence[i]);
//    printf("\n");
//    printf("\n");
    int i, j, k;
    int SubPath1[MAX_TASK_SEQ_LENGTH], SubPath2[MAX_TASK_SEQ_LENGTH];
    int Routes1[50][MAX_TASK_SEQ_LENGTH], Routes2[50][MAX_TASK_SEQ_LENGTH];
    double XCLds[50];
    int LeftTasks[MAX_TASK_SEQ_LENGTH], Positions[50];

    struct CandSelection
    {
        int RouteID;
        int Pos;
    };

    // CandSLCTList record which routes and the postion of each task located
    struct CandSelection CandSLCTList1[MAX_TASK_SEQ_LENGTH], CandSLCTList2[MAX_TASK_SEQ_LENGTH];
    int LLength1 = 0, LLength2 = 0;

    find_ele_positions(Positions, p1->Sequence, 0);
    Routes1[0][0] = Positions[0] - 1;
    for (i=1; i < Positions[0]; i++)
    {
        AssignSubArray(p1->Sequence, Positions[i], Positions[i+1], Routes1[i]);
    }

    find_ele_positions(Positions, p2->Sequence, 0);
    Routes2[0][0] = Positions[0] - 1;
    for (i=1; i < Positions[0]; i++)
    {
        AssignSubArray(p2->Sequence, Positions[i], Positions[i+1], Routes2[i]);
    }

    memcpy(XCLds, p1->Loads, sizeof(p1->Loads));

    for (i=1; i <= Routes1[0][0]; i++)
    {
        for (j = 2; j < Routes1[i][0]; j++)
        {
            LLength1 ++;
            CandSLCTList1[LLength1].RouteID = i;
            CandSLCTList1[LLength1].Pos = j;
        }
    }

    for (i=1; i <= Routes2[0][0]; i++)
    {
        for (j = 2; j < Routes2[i][0]; j++)
        {
            LLength2 ++;
            CandSLCTList2[LLength2].RouteID = i;
            CandSLCTList2[LLength2].Pos = j;
        }
    }

    int k1 = rand_choose(LLength1);
    int k2 = rand_choose(LLength2);

    AssignSubArray(Routes1[CandSLCTList1[k1].RouteID], 1, CandSLCTList1[k1].Pos, SubPath1);
    AssignSubArray(Routes2[CandSLCTList2[k2].RouteID], CandSLCTList2[k2].Pos,Routes2[CandSLCTList2[k2].RouteID][0], SubPath2);
    AssignSubArray(Routes1[CandSLCTList1[k1].RouteID], CandSLCTList1[k1].Pos+1,Routes1[CandSLCTList1[k1].RouteID][0]-1, LeftTasks);

    // remove duplicated tasks for Route1
    int Checked[MAX_TASK_SEQ_LENGTH];
    memset(Checked, 0, sizeof(Checked));

    for (i = 1; i < SubPath2[0]; i++)
    {
        if (Checked[i])
            continue;

        for (j = SubPath1[0]; j > 1; j--)
        {
            if (SubPath1[j] == SubPath2[i] || SubPath1[j] == inst_tasks[SubPath2[i]].inverse)
            {
                delete_element(SubPath1, j);
                Checked[i] = 1;
                break;
            }
        }
    }

    for (i = 1; i < SubPath2[0]; i++)
    {
        if (Checked[i])
            continue;

        for (j = LeftTasks[0]; j > 0; j--)
        {
            if (LeftTasks[j] == SubPath2[i] || LeftTasks[j] == inst_tasks[SubPath2[i]].inverse)
            {
                delete_element(LeftTasks, j);
                Checked[i] = 1;
                break;
            }
        }
    }

    for (i = 1; i < SubPath2[0]; i ++)
    {
        if (Checked[i])
            continue;

        for (j=1; j <= Routes1[0][0]; j++)
        {
            if (j == CandSLCTList1[k1].RouteID)
                continue;

            for (k = Routes1[j][0]; k > 1; k --)
            {
                if (Routes1[j][k] == SubPath2[i] || Routes1[j][k] == inst_tasks[SubPath2[i]].inverse)
                {
                    delete_element(Routes1[j], k);
                    XCLds[j] -= inst_tasks[SubPath2[i]].demand;
                    Checked[i] = 1;
                }
            }
            if (Checked[i])
                break;
        }
    }

    JoinArray(SubPath1, SubPath2);
    memcpy(Routes1[CandSLCTList1[k1].RouteID], SubPath1, sizeof(SubPath1));
    XCLds[CandSLCTList1[k1].RouteID] = 0;
    for (i = 2; i < Routes1[CandSLCTList1[k1].RouteID][0]; i++)
    {
        XCLds[CandSLCTList1[k1].RouteID] += inst_tasks[Routes1[CandSLCTList1[k1].RouteID][i]].demand;
    }

    int NO_LeftTasks = LeftTasks[0];

    // insert missing tasks
    struct Insert{
        int InsertedTask;
        int InsertRouteID;
        int InsertPos;
        int InsertCost;
        int InsertVioLoad;
    };

    struct Insert CandInsertions[6000];
    int NO_CandInsertions;
    double IVLoad;


    struct Insert ParetoSetInsertions[6000];
    int ParetoSetSize;
    int Out[6000], Add;

    struct Insert BestInsertion;

    int NO_CurrRoutes = 0;
    for (j = 1; j <= Routes1[0][0]; j++)
    {
        if (Routes1[j][0] == 2)
            continue;

        NO_CurrRoutes ++;
    }

    int n, m, l, w;
    int CurrTask;
    for (n = 1; n <= NO_LeftTasks; n++)
    {
        NO_CandInsertions = 0;
        ParetoSetSize = 0;

        for (j = 1; j <= Routes1[0][0]; j++)
        {
            if (Routes1[j][0] == 2)
                continue;

            CurrTask = LeftTasks[n];
            if (XCLds[j] > capacity)
            {
                IVLoad = inst_tasks[CurrTask].demand;
            } else if (XCLds[j] > (capacity - inst_tasks[CurrTask].demand)) {
                IVLoad = XCLds[j] + inst_tasks[CurrTask].demand - capacity;
            } else {
                IVLoad = 0;
            }

            for (k = 2; k <= Routes1[j][0]; k++)
            {
                NO_CandInsertions ++;
                CandInsertions[NO_CandInsertions].InsertedTask = CurrTask;
                CandInsertions[NO_CandInsertions].InsertRouteID = j;
                CandInsertions[NO_CandInsertions].InsertPos = k;
                CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[CurrTask].head_node] +
                        min_cost[inst_tasks[CurrTask].tail_node][inst_tasks[Routes1[j][k]].head_node] -
                        min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];

                CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;

                Out[0] = 0;
                Add = 1;

                for (m = 1; m <= ParetoSetSize; m++)
                {
                    if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost &&
                            CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
                    {
                        Add = 0;
                        break;
                    } else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost &&
                               CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
                    {
                        Out[0] ++;
                        Out[Out[0]] = m;
                    }
                }

                if (Add)
                {
                    for ( m = Out[0]; m > 0; m--)
                    {
                        for (l = Out[m]; l < ParetoSetSize; l++)
                        {
                            ParetoSetInsertions[l] = ParetoSetInsertions[l + 1];
                        }
                        ParetoSetSize --;
                    }

                    ParetoSetSize ++;
                    ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
                }

                w = inst_tasks[CurrTask].inverse;
                NO_CandInsertions ++;
                CandInsertions[NO_CandInsertions].InsertedTask = w;
                CandInsertions[NO_CandInsertions].InsertRouteID = j;
                CandInsertions[NO_CandInsertions].InsertPos = k;
                CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[CurrTask].head_node] +
                                                               min_cost[inst_tasks[w].tail_node][inst_tasks[Routes1[j][k]].head_node] -
                                                               min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];

                CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;

                Out[0] = 0;
                Add = 1;

                for (m = 1; m <= ParetoSetSize; m++)
                {
                    if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost &&
                        CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
                    {
                        Add = 0;
                        break;
                    } else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost &&
                               CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
                    {
                        Out[0] ++;
                        Out[Out[0]] = m;
                    }
                }

                if (Add)
                {
                    for ( m = Out[0]; m > 0; m--)
                    {
                        for (l = Out[m]; l < ParetoSetSize; l++)
                        {
                            ParetoSetInsertions[l] = ParetoSetInsertions[l + 1];
                        }
                        ParetoSetSize --;
                    }

                    ParetoSetSize ++;
                    ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
                }
            }
        }

        // insert as a new route;

        NO_CandInsertions ++;
        CandInsertions[NO_CandInsertions].InsertedTask = CurrTask;
        CandInsertions[NO_CandInsertions].InsertRouteID = 0;
        CandInsertions[NO_CandInsertions].InsertPos = 2;
        CandInsertions[NO_CandInsertions].InsertCost = min_cost[DEPOT][inst_tasks[CurrTask].head_node] + min_cost[inst_tasks[CurrTask].tail_node][DEPOT];
        CandInsertions[NO_CandInsertions].InsertVioLoad = 0;

        Out[0] = 0;
        Add = 1;

        for (m = 1; m <= ParetoSetSize; m++)
        {
            if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost &&
                CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
            {
                Add = 0;
                break;
            } else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost &&
                       CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
            {
                Out[0] ++;
                Out[Out[0]] = m;
            }
        }

        if (Add)
        {
            for ( m = Out[0]; m > 0; m--)
            {
                for (l = Out[m]; l < ParetoSetSize; l++)
                {
                    ParetoSetInsertions[l] = ParetoSetInsertions[l + 1];
                }
                ParetoSetSize --;
            }

            ParetoSetSize ++;
            ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
        }

        k = rand_choose(ParetoSetSize);
        BestInsertion = ParetoSetInsertions[k];
        if (BestInsertion.InsertRouteID == 0)
        {
            Routes1[0][0] ++;
            Routes1[Routes1[0][0]][0] = 3;
            Routes1[Routes1[0][0]][1] = 0;
            Routes1[Routes1[0][0]][2] = BestInsertion.InsertedTask;
            Routes1[Routes1[0][0]][3] = 0;

            XCLds[0] ++;
            XCLds[(int)XCLds[0]] = inst_tasks[BestInsertion.InsertedTask].demand;
        } else {
            add_element(Routes1[BestInsertion.InsertRouteID], BestInsertion.InsertedTask, BestInsertion.InsertPos);
            XCLds[BestInsertion.InsertRouteID] += inst_tasks[BestInsertion.InsertedTask].demand;
        }
    }


    // transfer Routes1 to sequence
    xed_child->Sequence[0] = 1;
    for (i = 1; i <= Routes1[0][0]; i++)
    {
        if (Routes1[i][0] == 2)
            continue;

        xed_child->Sequence[0] --;
        JoinArray(xed_child->Sequence, Routes1[i]);
    }

    xed_child->TotalCost = get_total_cost(xed_child->Sequence, inst_tasks);
    memcpy(xed_child->Loads, XCLds, sizeof(XCLds));

    for (i = xed_child->Loads[0]; i > 0; i --)
    {
        if (IS_DOUBLE_ZERO(xed_child->Loads[i]) == 1)
            delete_element(xed_child->Loads, i);
    }
    xed_child->TotalVioLoad = get_total_vio_load(xed_child->Loads);

    if (51 + xed_child->Loads[0] + 1 != xed_child->Sequence[0])
    {
        int kk = 0;
        for (int i = 1; i <= Routes1[0][0]; i++)
        {
            for (int j = 1; j <= Routes1[i][0]; j++)
            {
                //printf("%d \t", Routes1[i][j]);
                kk ++;
            }
            //printf("\n");
        }
        //printf("k: %d\n", kk);
        int a = 0;
    }

}

void lns_mut(Individual *c, Individual *p, Individual *best_fsb_solution, Task *inst_tasks)
{

    
    indi_copy(c, p);
    double coef = 1.0*best_fsb_solution->TotalCost/capacity*(1.0*best_fsb_solution->TotalCost/p->TotalCost+1.0*p->TotalVioLoad/capacity+1.0);
    c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
    double oldcost = best_fsb_solution->TotalCost;
    int count = 0;
    int count_fsb = 0, count_infsb = 0;

    int imp = 1; // improved
    while (imp)
    {
        count ++;
        imp = 0;

        if (c->TotalVioLoad == 0)
        {
            count_fsb ++;
        } else {
            count_infsb ++;
        }

        if (count % 5 == 0)
        {
            if (count_fsb == 5)
            {
                coef /= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            } else if (count_infsb == 5) {
                coef *= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            }
            count_fsb = 0;
            count_infsb = 0;
        }

        Individual tmp_indi;
        indi_copy(&tmp_indi, c);

        // trditional move
        lns(c, coef, 1, inst_tasks);

        if (c->Fitness < tmp_indi.Fitness)
            imp = 1;

        if (c->TotalVioLoad == 0 && c->TotalCost < best_fsb_solution->TotalCost)
        {
            indi_copy(best_fsb_solution, c); // check if assign or not
        }
    }

    imp = 1;
    while (imp)
    {
        imp = 0;
        count ++;

        if (c->TotalVioLoad == 0)
        {
            count_fsb ++;
        } else {
            count_infsb ++;
        }

        if (count % 5 == 0)
        {
            if (count_fsb == 5)
            {
                coef /= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            } else if (count_infsb == 5) {
                coef *= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            }
            count_fsb = 0;
            count_infsb = 0;
        }

        Individual tmp_indi;
        indi_copy(&tmp_indi, c);

        // MS
        lns(&tmp_indi, coef, 2, inst_tasks);

        if (c->Fitness < tmp_indi.Fitness)
            imp = 1;

        if (c->TotalVioLoad == 0 && c->TotalCost < best_fsb_solution->TotalCost)
        {
            indi_copy(best_fsb_solution, c); 
        }

    }

    imp = 1;
    while (imp)
    {
        imp = 0;
        count ++;

        if (c->TotalVioLoad == 0)
        {
            count_fsb ++;
        } else {
            count_infsb ++;
        }

        if (count % 5 == 0)
        {
            if (count_fsb == 5)
            {
                coef /= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            } else if (count_infsb == 5) {
                coef *= 5;
                c->Fitness = c->TotalCost + coef * c->TotalVioLoad;
            }
            count_fsb = 0;
            count_infsb = 0;
        }

        Individual tmp_indi;
        indi_copy(&tmp_indi, c);

        // trditional move

        lns(c, coef, 1, inst_tasks);

        if (c->Fitness < tmp_indi.Fitness)
            imp = 1;

        if (c->TotalVioLoad == 0 && c->TotalCost < best_fsb_solution->TotalCost)
        {
            indi_copy(best_fsb_solution, c); // check if assign or not
        }
    }
}

void lns(Individual *indi, double coef, int nsize, const Task *inst_tasks)
{
    indi->Fitness = indi->TotalCost + coef * indi->TotalVioLoad;

    if (nsize == 1) // traditional move
    {
        Move si_move, di_move, swap_move, next_move;
        next_move.fitness = INF;
        next_move.type = -1;
        single_insertion(&si_move, indi, coef, inst_tasks);
        double_insertion(&di_move, indi, coef, inst_tasks);
        swap(&swap_move, indi, coef, inst_tasks);

        if (si_move.fitness < next_move.fitness)
            next_move = si_move;

        if (di_move.fitness < next_move.fitness)
            next_move = di_move;

        if (swap_move.fitness < next_move.fitness)
            next_move = swap_move;

        if (next_move.type == -1)
            return;

        int orig_ptr = 0, targ_ptr = 0;
        int seg_ptr1 = 0, seg_ptr2 = 0;
        for (int i = 1; i < indi->Sequence[0]; i++)
        {
            if (indi->Sequence[i] == 0)
            {
                if (seg_ptr1 < next_move.orig_seg)
                    seg_ptr1 ++;

                if (seg_ptr2 < next_move.targ_seg)
                    seg_ptr2 ++;

                if (seg_ptr1 == next_move.orig_seg && orig_ptr == 0)
                    orig_ptr = i + next_move.orig_pos - 1;

                if (seg_ptr2 == next_move.targ_seg && targ_ptr == 0)
                    targ_ptr = i + next_move.targ_pos - 1;
            }
            if (orig_ptr != 0 && targ_ptr != 0)
                break;
        }
        switch (next_move.type) {
            case SI:
                {
                    delete_element(indi->Sequence, orig_ptr);
                    if (targ_ptr > orig_ptr)
                        targ_ptr --;
                    indi->Loads[next_move.orig_seg] -= inst_tasks[next_move.task1].demand;
                    if (next_move.targ_seg > indi->Loads[0])
                    {
                        indi->Sequence[0] ++;
                        indi->Sequence[indi->Sequence[0]] = next_move.task1;
                        indi->Sequence[0] ++;
                        indi->Sequence[indi->Sequence[0]] = 0;
                        indi->Loads[0] ++;
                        indi->Loads[(int)indi->Loads[0]] = inst_tasks[next_move.task1].demand;
                    } else {
                        add_element(indi->Sequence, next_move.task1, targ_ptr);
                        indi->Loads[next_move.targ_seg] += inst_tasks[next_move.task1].demand;
                    }
                }
                break;
            case DI:
                {
                    delete_element(indi->Sequence, orig_ptr+1);
                    delete_element(indi->Sequence, orig_ptr);
                    if (targ_ptr > orig_ptr)
                        targ_ptr -= 2;

                    indi->Loads[next_move.orig_seg] -= inst_tasks[next_move.task1].demand + inst_tasks[next_move.task2].demand;
                    if (next_move.targ_seg > indi->Loads[0])
                    {
                        indi->Sequence[0] ++;
                        indi->Sequence[indi->Sequence[0]] = next_move.task1;
                        indi->Sequence[0] ++;
                        indi->Sequence[indi->Sequence[0]] = next_move.task2;
                        indi->Sequence[0] ++;
                        indi->Sequence[indi->Sequence[0]] = 0;
                        indi->Loads[0] ++;
                        indi->Loads[(int)indi->Loads[0]] = inst_tasks[next_move.task1].demand + inst_tasks[next_move.task2].demand;
                    } else {
                        add_element(indi->Sequence, next_move.task2, targ_ptr);
                        add_element(indi->Sequence, next_move.task1, targ_ptr);
                        indi->Loads[next_move.targ_seg] += inst_tasks[next_move.task1].demand + inst_tasks[next_move.task2].demand;
                    }
                }
                break;
            case SWAP:
                {
                    indi->Sequence[targ_ptr] = next_move.task1;
                    indi->Sequence[orig_ptr] = next_move.task2;
                    indi->Loads[next_move.orig_seg] -= inst_tasks[next_move.task1].demand - inst_tasks[next_move.task2].demand;
                    indi->Loads[next_move.targ_seg] += inst_tasks[next_move.task1].demand - inst_tasks[next_move.task2].demand;
                }
                break;
        }

        
        indi->TotalVioLoad = next_move.total_vio_load;
        indi->Fitness = next_move.fitness;
        if (IS_DOUBLE_ZERO(indi->Loads[next_move.orig_seg]) == 1)
        {
            if (next_move.type == DI && next_move.orig_seg > next_move.targ_seg)
            {
                delete_element(indi->Sequence, orig_ptr+1);
            } else {
                delete_element(indi->Sequence, orig_ptr);
            }

            delete_element(indi->Loads, next_move.orig_seg);
        }
        indi->TotalCost = get_total_cost(indi->Sequence,inst_tasks);

    }
    else
    {
        int i, j, k;
        // merge and split
        task_routes[0][0] = 1;
        task_routes[1][0] = 1;
        task_routes[1][1] = 0;
        for (i = 2; i <= indi->Sequence[0]; i++)
        {
            task_routes[task_routes[0][0]][0] ++;
            task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->Sequence[i];
            if (indi->Sequence[i] == 0 && i < indi->Sequence[0])
            {
                task_routes[0][0] ++;
                task_routes[task_routes[0][0]][0] = 1;
                task_routes[task_routes[0][0]][1] = 0;
            }
        }

        if (task_routes[0][0] < nsize)
            return;

        int multi = task_routes[0][0];
        long int ub_trial = task_routes[0][0];
        for (i = 1; i < nsize; i++)
        {
            multi--;
            ub_trial *= multi;
        }

        multi = nsize;
        for (i = 1; i < nsize; i++)
        {
            ub_trial /= multi;
            multi--;
        }

        int maxcount = ub_trial;
        if (maxcount > MAX_ENSSIZE)
            maxcount = MAX_ENSSIZE;

        typedef struct lns_comb {
            int ids[MAX_NSIZE + 1];
        } lns_comb;

        lns_comb cand_combs[MAX_ENSSIZE + 1];

        int pointers[MAX_NSIZE + 1];
        for (i = 1; i <= nsize; i++)
        {
            pointers[i] = i;
        }
        int curr_ptr;
        for (i = 1; i <= maxcount; i++)
        {
            cand_combs[i].ids[0] = nsize;
            for (j = 1; j <= nsize; j++)
            {
                cand_combs[i].ids[j] = pointers[j];
            }
            curr_ptr = nsize;

            while (pointers[curr_ptr] == task_routes[0][0] - nsize + curr_ptr)
            {
                curr_ptr--;
            }
            if (curr_ptr == 0)
                break;

            pointers[curr_ptr] ++;
            for (j = curr_ptr + 1; j <= nsize; j++)
            {
                pointers[j] = pointers[j - 1] + 1;
            }
        }

        Individual tmp_indi, next_indi;
        next_indi.Fitness = INF;

        int lns_routes[MAX_NSIZE + 1];
        double sel_total_load;

        for (i = 1; i <= maxcount; i++)
        {
            memcpy(lns_routes, cand_combs[i].ids, sizeof(cand_combs[i].ids));

            sel_total_load = 0;
            for (j = 1; j <= lns_routes[0]; j++)
            {
                sel_total_load += indi->Loads[lns_routes[j]];
            }

            if (sel_total_load > nsize * capacity)
                continue;

            int serve_mark[MAX_TASK_TAG_LENGTH];
            memset(serve_mark, 0, sizeof(serve_mark));

            for (j = 1; j <= lns_routes[0]; j++)
            {
                for (k = 2; k < task_routes[lns_routes[j]][0]; k++)
                {
                    serve_mark[task_routes[lns_routes[j]][k]] = 1;
                    serve_mark[inst_tasks[task_routes[lns_routes[j]][k]].inverse] = 1;
                }
            }

            path_scanning(&tmp_indi, inst_tasks, serve_mark);

            for (j = 1; j <= task_routes[0][0]; j++)
            {
                int lnused = 0;
                for (k = 1; k <= lns_routes[0]; k++)
                {
                    if (j == lns_routes[k])
                    {
                        lnused = 1;
                        break;
                    }
                }
                if (lnused)
                    continue;

                tmp_indi.Sequence[0] --;
                JoinArray(tmp_indi.Sequence, task_routes[j]);
                tmp_indi.Loads[0] ++;
                tmp_indi.Loads[(int)tmp_indi.Loads[0]] = indi->Loads[j];
                if (indi->Loads[j] > capacity)
                    tmp_indi.TotalVioLoad += (indi->Loads[j] - capacity);  // bug in the original program
            }

            tmp_indi.TotalCost = get_task_seq_total_cost(tmp_indi.Sequence, inst_tasks);
            //            tmp_indi.TotalVioLoad = get_total_vio_load(tmp_indi.Loads);
            tmp_indi.Fitness = tmp_indi.TotalCost + coef * tmp_indi.TotalVioLoad;
            if (tmp_indi.Fitness < next_indi.Fitness)
                indi_copy(&next_indi, &tmp_indi);
        }

        if (next_indi.Fitness < indi->Fitness)
            indi_copy(indi, &next_indi);

    }

}
void single_insertion(Move *best_move, Individual *indi, double coef, const Task *inst_tasks)
{
    best_move->type = SI;
    best_move->fitness = INF;

    int i, j, s1, s2;

    task_routes[0][0] = 1;
    task_routes[1][0] = 1;
    task_routes[1][1] = 0;
    for (i = 2; i <= indi->Sequence[0]; i++)
    {
        task_routes[task_routes[0][0]][0] ++;
        task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->Sequence[i];
        if (indi->Sequence[i] == 0 && i < indi->Sequence[0])
        {
            task_routes[0][0] ++;
            task_routes[task_routes[0][0]][0] = 1;
            task_routes[task_routes[0][0]][1] = 0;
        }
    }

    Move tmp_move;
    for (s1 = 1; s1 <= task_routes[0][0]; s1++)
    {
        tmp_move.orig_seg = s1;
        for (i = 2; i < task_routes[s1][0]; i++)
        {
            tmp_move.orig_pos = i;

            for (s2 = 1; s2 <= task_routes[0][0]; s2++)
            {

                if (s2 == s1)
                    continue;

                tmp_move.targ_seg = s2;
                if (s2 > task_routes[0][0])
                {
                    tmp_move.total_vio_load = indi->TotalVioLoad;
                    if (indi->Loads[s1] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s1]-capacity;
                    if (indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;

                    tmp_move.task1 = task_routes[s1][i];

                    tmp_move.total_cost = indi->TotalCost + min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]
                            + min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]
                            + min_cost[inst_tasks[tmp_move.task1].tail_node][DEPOT]
                            - min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]
                            - min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }
                    continue;
                }

                for (j = 2; j <= task_routes[s2][0]; j++)
                {
                    if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node)
                        continue;

                    tmp_move.targ_pos = j;
                    tmp_move.total_vio_load = indi->TotalVioLoad;
                    if (indi->Loads[s1] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s1]-capacity;
                    if (indi->Loads[s2] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s2]-capacity;
                    if (indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;
                    if (indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand-capacity;


                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.total_cost = indi->TotalCost
                            + min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]
                            - min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]
                            - min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]
                            + min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]
                            + min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]
                            - min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;
                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }


                    tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
                    tmp_move.total_cost = indi->TotalCost
                            + min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]
                            - min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]
                            - min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]
                            + min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]
                            + min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]
                            - min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;
                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }
                }
            }
        }
    }
}

void double_insertion(Move *best_move, Individual *indi, double coef, const Task *inst_tasks)
{
    best_move->type = DI;
    best_move->fitness = INF;

    int i, j, s1, s2;

    task_routes[0][0] = 1;
    task_routes[1][0] = 1;
    task_routes[1][1] = 0;
    for (i = 2; i <= indi->Sequence[0]; i++)
    {
        task_routes[task_routes[0][0]][0] ++;
        task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->Sequence[i];
        if (indi->Sequence[i] == 0 && i < indi->Sequence[0])
        {
            task_routes[0][0] ++;
            task_routes[task_routes[0][0]][0] = 1;
            task_routes[task_routes[0][0]][1] = 0;
        }
    }

    Move tmp_move;

    for (s1 = 1; s1 <= task_routes[0][0]; s1++)
    {
        if (task_routes[s1][0] < 4)
            continue;

        tmp_move.orig_seg = s1;
        for (i = 2; i < task_routes[s1][0]-1; i++)
        {
            tmp_move.orig_pos = i;
            for (s2 = 1; s2 <= task_routes[0][0]; s2++)  /* s2 > task_routes[0][0] --> create a new route */
            {
                if (s2 == s1)
                    continue;

                tmp_move.targ_seg = s2;
                if (s2 > task_routes[0][0] /* && task_routes[0][0] < vehicle_num*/)
                {
                    if (task_routes[s1][0] <= 4)
                        continue;

                    tmp_move.total_vio_load = indi->TotalVioLoad;
                    if (indi->Loads[s1] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s1]-capacity;
                    if (indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = task_routes[s1][i+1];

                    tmp_move.total_cost = indi->TotalCost +
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node] -
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }
                    continue;
                }

                for (j = 2; j <= task_routes[s2][0]; j++)
                {
                    if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node)
                        continue;

                    tmp_move.targ_pos = j;
                    tmp_move.total_vio_load = indi->TotalVioLoad;
                    if (indi->Loads[s1] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s1]-capacity;
                    if (indi->Loads[s2] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s2]-capacity;
                    if (indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;
                    if (indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand-capacity;

                    //

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = task_routes[s1][i+1];

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
                    tmp_move.task2 = task_routes[s1][i+1];

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
                    tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }
                }
            }
        }
    }
}

void swap(Move *best_move, Individual *indi, double coef, const Task *inst_tasks)
{
    best_move->type = SWAP;
    best_move->fitness = INF;

    int i, j, s1, s2;

    task_routes[0][0] = 1;
    task_routes[1][0] = 1;
    task_routes[1][1] = 0;
    for (i = 2; i <= indi->Sequence[0]; i++) {
        task_routes[task_routes[0][0]][0]++;
        task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->Sequence[i];
        if (indi->Sequence[i] == 0 && i < indi->Sequence[0]) {
            task_routes[0][0]++;
            task_routes[task_routes[0][0]][0] = 1;
            task_routes[task_routes[0][0]][1] = 0;
        }
    }

    Move tmp_move;
    for (s1 = 1; s1 < task_routes[0][0]; s1++)
    {
        tmp_move.orig_seg = s1;
        for (i = 2; i < task_routes[s1][0]; i++)
        {
            if (inst_tasks[task_routes[s1][i-1]].tail_node == inst_tasks[task_routes[s1][i]].head_node &&
                inst_tasks[task_routes[s1][i]].tail_node == inst_tasks[task_routes[s1][i+1]].head_node)
                continue;

            tmp_move.orig_pos = i;
            for (s2 = s1+1; s2 <= task_routes[0][0]; s2++)
            {
                tmp_move.targ_seg = s2;
                for (j = 2; j < task_routes[s2][0]; j++)
                {
                    if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node &&
                        inst_tasks[task_routes[s2][j]].tail_node == inst_tasks[task_routes[s2][j+1]].head_node)
                        continue;

                    tmp_move.targ_pos = j;

                    tmp_move.total_vio_load = indi->TotalVioLoad;
                    if (indi->Loads[s1] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s1]-capacity;
                    if (indi->Loads[s2] > capacity)
                        tmp_move.total_vio_load -= indi->Loads[s2]-capacity;
                    if (indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand-capacity;
                    if (indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand > capacity)
                        tmp_move.total_vio_load += indi->Loads[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand-capacity;

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = task_routes[s2][j];

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
                    tmp_move.task2 = task_routes[s2][j];

                    tmp_move.total_cost = indi->TotalCost +
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = task_routes[s1][i];
                    tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }

                    //

                    tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
                    tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;

                    tmp_move.total_cost = indi->TotalCost+
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
                                          min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
                                          min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
                                          min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
                                          min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

                    tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

                    if (tmp_move.fitness < best_move->fitness)
                    {
                        best_move->task1 = tmp_move.task1;
                        best_move->task2 = tmp_move.task2;
                        best_move->orig_seg = tmp_move.orig_seg;
                        best_move->targ_seg = tmp_move.targ_seg;
                        best_move->orig_pos = tmp_move.orig_pos;
                        best_move->targ_pos = tmp_move.targ_pos;
                        best_move->total_cost = tmp_move.total_cost;
                        best_move->total_vio_load = tmp_move.total_vio_load;
                        best_move->fitness = tmp_move.fitness;
                    }
                }
            }
        }
    }

}

























