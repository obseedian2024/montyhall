/*
Monty Hall Problem Simulator
Copyright (C) 2024  Obseedian

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
https://www.gnu.org/licenses/
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#if RAND_MAX == 0x7FFF 
static const int RAND_MAX32 = (RAND_MAX << 15) | RAND_MAX;

// Return a random number between 0 and m-1
static unsigned getrand(unsigned m)
{
    int r = (rand() << 15) | rand();
    return (unsigned)((double)r / ((double)(RAND_MAX32 + 1)) * m);
}
#else
static unsigned getrand(unsigned m)
{
    return (unsigned)((double)rand() / ((double)RAND_MAX + 1) * m);
}
#endif

static double mean(double* x, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += x[i];
    return sum / (double)n;
}

static double sample_stddev(double* x, double u, int n)
{
    double sum2 = 0.0;
    for (int i = 0; i < n; i++)
        sum2 += (u - x[i]) * (u - x[i]);
    return sqrt(sum2 / (double)(n - 1));
}

const unsigned door_mask = 7;

int main()
{
    char buf[100];

    int
        switch_win,
        switch_lose,
        stay_win,
        stay_lose,
        n,
        runs;

    puts("Monty Hall Simulator v0.0001 Copyright (C) 2024  Obseedian\n"
        "Licensed under GNU GPL v2\n");
    printf("Enter number of runs...: ");
    fgets((char*)&buf, sizeof(buf), stdin);
    runs = atoi((char*)&buf);
    if (runs <= 0)
        return 0;

    printf("Enter N (games per run): ");
    fgets((char*)&buf, sizeof(buf), stdin);
    n = atoi((char*)&buf);
    if (n <= 0)
        return 0;

    printf("\nRuns = %u\n", runs);
    printf("N    = %u\n", n);

    unsigned cnt_stays_winpct = 0, cnt_switches_winpct = 0;
    double* arr_stays_winpct = calloc(runs, sizeof(double));
    double* arr_switches_winpct = calloc(runs, sizeof(double));
    if ((arr_stays_winpct == NULL) || (arr_switches_winpct == NULL))
    {
        puts("calloc failed");
        return -3;
    }

    srand((unsigned)time(NULL));
    for (int i = 0; i < runs; i++)
    {
        switch_win = switch_lose = stay_win = stay_lose = 0;

        for (int j = 0; j < n; j++)
        {
            unsigned car = getrand(3);  //which door the car is behind
            unsigned pick1 = getrand(3); //initial pick of contestant
            unsigned door;

            if (car == pick1)
            {
                //contestant picks the car. Monty can show any other door
                //Randomly pick the next door or second door, wrap around if needed
                door = 1 << ((car + 1 + getrand(2)) % 3);

                //Alternatively Monty could always pick the next door
                //It will not alter the outcome
                //door = 1 << ((car + 1) % 3);
            }
            else
            {
                //Contestant picked goat. Monty has to show the other goat
                //convert door number to bit, not car or pick
                door = door_mask & ~((1 << car) | (1 << pick1));
            }
            //bit 0 = door 0, bit 1 = door 1, bit 2 = door 2
            car = 1 << car;
            pick1 = 1 << pick1;

            //Let's randomly decide if contestant switches or stays
            unsigned stay = getrand(2);

            if (stay == 0)
            {
                // Case when contestant switches. The second pick is not pick1 or door
                unsigned pick2 = door_mask & ~(pick1 | door);
                //check if he won
                if (pick2 == car)
                    switch_win++;
                else
                    switch_lose++;
            }
            else {
                // Case when contestant stays
                if (pick1 == car)
                    stay_win++;
                else
                    stay_lose++;
            }
        }
        unsigned total_switches = switch_win + switch_lose;
        unsigned total_stays = stay_win + stay_lose;
        double switch_win_pct = 0.0, stay_win_pct = 0.0;

        if (total_switches > 0)
        {
            switch_win_pct = ((double)switch_win / (double)total_switches) * 100.0;
            arr_switches_winpct[cnt_switches_winpct++] = switch_win_pct;
        }

        if (total_stays > 0)
        {
            stay_win_pct = ((double)stay_win / (double)total_stays) * 100.0;
            arr_stays_winpct[cnt_stays_winpct++] = stay_win_pct;
        }
        if (runs == 1)
        {
            printf("\nContestant switches\n");
            printf("Number of wins  : %u\n", switch_win);
            printf("Number of losses: %u\n", switch_lose);
            printf("Win percentage  : %3.4f%%\n\n", switch_win_pct);

            printf("Contestant stays\n");
            printf("Number of wins  : %u\n", stay_win);
            printf("Number of losses: %u\n", stay_lose);
            printf("Win percentage  : %3.4f%%\n\n", stay_win_pct);
        }
    }

    if (runs > 1)
    {
        double mean_switch = mean(arr_switches_winpct, cnt_switches_winpct);
        if (cnt_switches_winpct < 2)
        {
            puts("Not enough sample size to calculate mean & standard deviation \n"
                "for switch win percentage. Try increasing runs.\n");
        }
        else
        {
            printf("\nMean of win percentage for switching : %3.4f\n", mean_switch);
            printf("Sigma of win percentage for switching: %3.4f\n",
                sample_stddev(arr_switches_winpct, mean_switch, cnt_switches_winpct));
        }

        double mean_stay = mean(arr_stays_winpct, cnt_stays_winpct);
        if (cnt_stays_winpct < 2)
        {
            puts("Not enough sample size to calculate mean & standard deviation \n"
                "for stay win percentage. Try increasing runs.\n");
        }
        else
        {
            printf("\nMean of win percentage for staying   : %3.4f\n", mean_stay);
            printf("Sigma of win percentage for staying  : %3.4f\n",
                sample_stddev(arr_stays_winpct, mean_stay, cnt_stays_winpct));
        }
    }
    free(arr_stays_winpct);
    free(arr_switches_winpct);
}
