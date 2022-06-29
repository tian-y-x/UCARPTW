#include "MAENS.h"
#include <fstream>
#include <iostream>
#include "../functions.h"
using namespace std;
void ShowIndi(Individual *pop, int popsize);

void MAENS(Task *inst_tasks)
{
    int i, j;
    int popsize = 30;
    
    Individual pop[MAX_TOTALSIZE];
    Individual best_fsb_solution;
    Individual best_rob_solution;
    best_fsb_solution.TotalCost = INF;
    // initilization
    int tmp_popsize = 0;
    int used;
    while (tmp_popsize < popsize)
    {
        int trial = 0;
        Individual init_indi;
        while (trial < M_trial)
        {
            trial ++;
            int serve_mark[MAX_TASK_TAG_LENGTH];
            memset(serve_mark, 0, sizeof(serve_mark));
            for (i = 1; i <= task_num; i++)
            {
                serve_mark[i] = 1;
            }

            rand_scanning(&init_indi, inst_tasks, serve_mark);
            used = 0;
            for (i = 0; i < tmp_popsize; i++)
            {
                if (init_indi.TotalCost == pop[i].TotalCost && init_indi.TotalVioLoad == pop[i].TotalVioLoad)
                {
                    used = 1;
                    break;
                }
            }
            if ( !used )
                break;
        }

        if (trial == M_trial && used == 1)
            break;

//        indi_copy(&pop[tmp_popsize], &init_indi);

        pop[tmp_popsize] = init_indi;
        tmp_popsize ++;
        if(init_indi.TotalVioLoad == 0 && init_indi.TotalCost < best_fsb_solution.TotalCost)
        {
            best_fsb_solution = init_indi;
        }
    }
    popsize = tmp_popsize;
    //for(int i=0;i< popsize;i++)
    //    cout<< pop[i].TotalCost<<" ";
    // main loop

    int ite, wite;
    Individual parent1, parent2, xed_child, mted_child, child;

    int offsize = 6*popsize;
    int totalsize = popsize + offsize;

    ite = 0;
    wite = 0;
    while (evaluate_time < FEm)
    {
        ite ++;
        wite ++;

        int ptr = popsize;
        while (ptr < totalsize)
        {
            child.TotalCost = 0;

            // randomly select two parents
            int par_id1, par_id2;
            rand_selection(&par_id1, &par_id2, popsize);
            parent1 = pop[par_id1];
            parent2 = pop[par_id2];

            // crossover
            SBX(&xed_child, &parent1, &parent2, inst_tasks);
            if (xed_child.TotalVioLoad == 0 && xed_child.TotalCost < best_fsb_solution.TotalCost)
            {
                best_fsb_solution = xed_child;
                wite = 0;
            }
            used = 0;
            for (i = 0; i < ptr; i++)
            {
                if (i == par_id1 || i == par_id2)
                    continue;

                if (xed_child.TotalCost == pop[i].TotalCost && xed_child.TotalVioLoad == pop[i].TotalVioLoad)
                {
                    used = 1;
                    break;
                }
            }

            if (!used)
            {
                child = xed_child;
            }

            // Local Search with Probability
            double random = 1.0 * rand() / RAND_MAX;
            if (random < M_PROB)
            {
                // do the local search.
//                printf("Local search %d \n", inst_tasks[1].inverse);
                lns_mut(&mted_child, &xed_child, &best_fsb_solution, inst_tasks);

                used = 0;
                for (i = 0; i < ptr; i++)
                {
                    if (i == par_id1 || i == par_id2)
                        continue;

                    if (mted_child.TotalCost == pop[i].TotalCost && mted_child.TotalVioLoad == pop[i].TotalVioLoad)
                    {
                        used = 1;
                        break;
                    }
                }

                if (!used)
                {
                    child = mted_child;
                }
            }

            if (child.TotalCost == parent1.TotalCost && child.TotalVioLoad == parent1.TotalVioLoad)
            {
                pop[par_id1] = child;
            } else if (child.TotalCost == parent2.TotalCost && child.TotalVioLoad == parent2.TotalVioLoad)
            {
                pop[par_id2] = child;
            } else if (child.TotalCost > 0)
            {
                pop[ptr] = child;
                ptr ++;
            }

        }

        float Pf = 0.45;
        Individual tmp_indi;


        for (i = 0; i < totalsize; i++)
        {
            for (j = 0; j < i; j++)
            {
                double random = 1.0 * rand() / RAND_MAX;
                if ( (pop[j].TotalVioLoad == 0 && pop[j+1].TotalVioLoad == 0) || random < Pf )
                {
                    if (pop[j].TotalCost > pop[j+1].TotalCost)
                    {
                        tmp_indi = pop[j];
                        pop[j] = pop[j+1];
                        pop[j+1] = tmp_indi;
                    }
                } else {
                    if (pop[j].TotalVioLoad > pop[j+1].TotalVioLoad)
                    {
                        tmp_indi = pop[j];
                        pop[j] = pop[j+1];
                        pop[j+1] = tmp_indi;
                    }
                }
            }
        }
        /*for(int i = 0;i<offsize;i++)
            cout << pop[i].TotalCost << ",";*/
        //for (int i = 0; i < offsize; i++)
        //{
        //    
        //    pop[i].TotalCost = eva_fitness(pop[i].Sequence, instances);
        //    if (pop[i].TotalCost < best_rob.TotalCost)
        //    {
        //        best_rob = pop[i];
        //    }
        //    
        //}
        cout << ite << ": " << best_fsb_solution.TotalCost << " | " << endl;
        if (best_fsb_solution.TotalVioLoad > 0)
        {
            cout << "\n*****" << endl;
        }


        if (ite == M_ite)
        {
            for (int i = 0; i <= best_fsb_solution.Sequence[0]; i++)
            {
                res[ind][i] = best_fsb_solution.Sequence[i];
                //cout << best_fsb_solution.Sequence[i] << " " << inst_tasks[best_fsb_solution.Sequence[i]].demand << endl;
                //cout << res[ind][i] << " demand: "<< inst_tasks[best_fsb_solution.Sequence[i]].demand << endl;

            }
        }
    }


}

void ShowIndi(Individual *pop, int popsize)
{
    for (int i= 0; i < popsize; i++)
    {
        printf("%d, TC: %lf, TVIO: %lf \n", i, pop[i].TotalCost, pop[i].TotalVioLoad);
    }
}
















