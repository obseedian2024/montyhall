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
#include <stdint.h>
#include <stdbool.h>

#if (defined(__clang__) || defined(__GNUC__)) && defined(__x86_64__)
#define HAS_UINT128 1
#endif

#define LO(x) ((x) & 0xFFFFFFFFU)
#define HI(x) ((x) >> 32)

#if defined(_MSC_VER)
/* Available on all platforms in MSVC */
#define ror64 _rotr64
#define rol64 _rotl64
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <x86intrin.h>
#define ror64 __rorq
#define rol64 __rolq
#else
static uint64_t ror64(const uint64_t value, int shift)
{
    shift &= 63;
    return (value >> shift) | (value << (64 - shift));
}
static uint64_t rol64(const uint64_t value, int shift)
{
    shift &= 63;
    return (value << shift) | (value >> (64 - shift));
}
#endif

static unsigned num_doors, door_mask;
static const unsigned max_doors = sizeof(unsigned) * 8;

/* Based on Chris Doty-Humphrey’s Small Fast Chaotic PRNG */
typedef struct rng_sfc {
    uint64_t a, b, c, counter;
} rng_sfc;

static rng_sfc sfc_state;

static uint64_t sfc64_rand(void)
{   
    enum { BARREL_SHIFT = 24, RSHIFT = 11, LSHIFT = 3 };//good sets include {24,11,3},{25,12,3},{},{} ; older versions used {25,12,3}, which is decent
    rng_sfc *rng = &sfc_state;
    uint64_t tmp = rng->a + rng->b + rng->counter++;
    rng->a = rng->b ^ (rng->b >> RSHIFT);
    rng->b = rng->c + (rng->c << LSHIFT);
    rng->c = rol64(rng->c, BARREL_SHIFT) + tmp;
    return tmp;
}

static void sfc64_seed(uint64_t s)
{
    rng_sfc *rng = &sfc_state;
    rng->a = rng->b = rng->c = s;
    rng->counter = 1;
    for (int i = 0; i < 12; i++) 
        (void)sfc64_rand();
}

/* Slightly modified SFC that replaces counter with PCG */
/* Is it better? */
static uint64_t sfc64_tf_rand(void)
{
    enum { BARREL_SHIFT = 24, RSHIFT = 11, LSHIFT = 3 };//good sets include {24,11,3},{25,12,3},{},{} ; older versions used {25,12,3}, which is decent
    rng_sfc *rng = &sfc_state;
    uint64_t tmp = rng->a + rng->b + ror64(rng->counter , rng->counter >> 58); 
    rng->counter = rng->counter * 3202034522624059733ull + 11; // Multiplier from L'Ecuyer
    rng->a = rng->b ^ (rng->b >> RSHIFT);
    rng->b = rng->c + (rng->c << LSHIFT);
    rng->c = rol64(rng->c, BARREL_SHIFT) + tmp;
    return tmp;
}

static void sfc64_tf_seed(uint64_t s)
{
    rng_sfc *rng = &sfc_state;
    rng->a = rng->b = rng->c = s;
    rng->counter = 1;
    for (int i = 0; i < 12; i++) 
        (void)sfc64_tf_rand();
}

/* From 'Fast Random Integer Generation in an Interval', Lemire
    https://arxiv.org/pdf/1805.10941 
*/
static const uint64_t RNG_A = 15750249268501108917ULL;

/* MCG RNG with 128-bit state */
typedef union rng_uint128
{
    struct {
        uint32_t dw0, dw1, dw2, dw3;
    };
    struct {
        uint64_t qw0, qw1;
    };
#if defined(HAS_UINT128)
    __uint128_t ow0;
#endif
} rng_uint128;

rng_uint128 pcg_state;

static uint64_t pcg64_rand(void)
{
#if defined(HAS_UINT128)
    pcg_state.ow0 = pcg_state.ow0 * RNG_A;
#else 
    uint64_t t;
    uint32_t cy;
    uint32_t r0, r1, r2, r3;
    const uint32_t a0 = LO(RNG_A);
    const uint32_t a1 = HI(RNG_A);

    t = (uint64_t)pcg_state.dw0 * a0;
    cy = HI(t);
    r0 = LO(t);

    t = (uint64_t)pcg_state.dw1 * a0 + cy;
    cy = HI(t);
    r1 = LO(t);

    t = (uint64_t)pcg_state.dw2 * a0 + cy;
    cy = HI(t);
    r2 = LO(t);

    t = (uint64_t)pcg_state.dw3 * a0 + cy;
    r3 = LO(t);

    t = (uint64_t)pcg_state.dw0 * a1 + r1;
    r1 = LO(t);
    cy = HI(t);

    t = (uint64_t)pcg_state.dw1 * a1 + r2 + cy;
    r2 = LO(t);
    cy = HI(t);

    t = (uint64_t)pcg_state.dw2 * a1 + r3 + cy;
    r3 = LO(t);

    pcg_state.dw0 = r0;
    pcg_state.dw1 = r1;
    pcg_state.dw2 = r2;
    pcg_state.dw3 = r3;
#endif
/* Randomly rotate xor-shifted output
   From 'PCG: A Family of Simple Fast Space-Efficient
   Statistically Good Algorithms for Random Number Generation', O'Neill
   https://www.pcg-random.org/pdf/hmc-cs-2014-0905.pdf
*/
    return ror64(pcg_state.qw0 ^ pcg_state.qw1, pcg_state.dw3 >> 26);
}

static void pcg64_seed(uint64_t seed)
{
    pcg_state.qw0 = seed | 1;
    pcg_state.qw1 = seed;
    for (int i = 0; i < 20; ++i) {
        (void)pcg64_rand();
    }
}

static uint64_t (*rand64)(void);
static void     (*rseed)(uint64_t);

static uint64_t (*rng_rand64[3])(void)      = { pcg64_rand, sfc64_rand, sfc64_tf_rand };
static void     (*rng_seed[3])(uint64_t)    = { pcg64_seed, sfc64_seed, sfc64_tf_seed };
static char*    rng_name[3]                 = { "PCG64", "SFC64", "SFC64-TF" };

enum { RNG_PCG64 = 1, RNG_SFC64 = 2, RNG_SFC64_TF = 3 };

static inline uint32_t rand32(void)
{
    return HI(rand64());
}

static inline uint64_t rand52(void)
{
    return rand64() >> 12;
}

static double frand(void)
{
    union {
        uint64_t i;
        double   d;
    } r; /* type punning this way is legal in C (but NOT C++) */

    /* Generate an IEEE754 double in the interval [1,2) as an int64 */
    r.i = rand52() | 0x3FF0000000000000;
    return r.d - 1.0; /* bitcast & return number in interval [0,1) */
}

static uint32_t bounded_rand(uint32_t range)
{
    uint32_t x = rand32();
    uint64_t m = (uint64_t)x * (uint64_t)range;
    uint32_t l = (uint32_t)m;
    if (l < range) {
        uint32_t t = (uint32_t)(-(int32_t)range) % range;
        while (l < t) {
            x = rand32();
            m = (uint64_t)x * (uint64_t)range;
            l = (uint32_t)m;
        }
    }
    return (uint32_t)(m >> 32);
}

/* Return a random number between 0 and m - 1 and cumulative distribution function
     specified as an array in cdf or NULL for even distribution*/
static unsigned getrand_cdf(unsigned m, double* cdf)
{
    if (cdf == NULL) {
        return bounded_rand(m); /* fallback to even probabilities */
    }
    else {
        unsigned i;
        double r = frand();
        for (i = 0; i < m - 1; i++) {
            if (r < cdf[i])
                return i;
        }
        return i;
    }
}

static unsigned getrand(unsigned m)
{   
    return bounded_rand(m);
}

static unsigned getuint(char* prompt, unsigned def_val)
{
    char buf[100] = "";
    printf(prompt, def_val);
    fgets((char*)&buf, sizeof(buf), stdin);
    if (buf[0] == '\n') 
        return def_val;
    int v = atoi((char*)&buf);
    if (v <= 0)
        return 0;
    return v;
}

static double sample_mean(double* x, int n)
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

static unsigned pick_rnd_door(unsigned doors, unsigned df, double* cdf)
{
    unsigned pick = getrand_cdf(df, cdf);
    for (unsigned i = 1; i <= door_mask; i <<= 1)
    {
        if (doors & i)
        {
            if (pick == 0)
                return i;
            pick--;
        }
    }
    return 0;
}

int main()
{
    int
        switch_win,
        switch_lose,
        stay_win,
        stay_lose,
        n,
        runs;
 
    unsigned rng;

    puts("Monty Hall Simulator v0.1 Copyright (C) 2024  Obseedian\n"
        "Licensed under GNU GPL v2\n");

    if (!(num_doors = getuint("Enter number of doors [%u]: ", 3)))
        return 0;
    if (num_doors < 3 || num_doors > max_doors)
    {
        printf("Number of doors must be between 3 and %u.\n", max_doors);
        return -1;
    }
    if (!(runs = getuint("Enter number of runs [%u]: ", 1)))
        return 0;
    if (!(n =    getuint("Enter N (games per run) [%u]: ", 1000)))
        return 0;
    rng = getuint("Select random number generator (1 - PCG64, 2 - SFC64, 3 - SFC64-TF) [%u]: ", 1);
    if (rng < 1 || rng > 3) {
        return 0;
    }
    door_mask = (1 << num_doors) - 1;
    rand64 = rng_rand64[rng - 1];
    rseed = rng_seed[rng - 1];

    rseed(0xC68FA87D83F72455);

    printf("\nRuns = %u\n", runs);
    printf("N    = %u\n", n);
    printf("RNG  = %s\n", rng_name[rng - 1]);

    unsigned cnt_stays_winpct = 0, cnt_switches_winpct = 0;
    double* arr_stays_winpct = calloc(runs, sizeof(double));
    double* arr_switches_winpct = calloc(runs, sizeof(double));
    if ((arr_stays_winpct == NULL) || (arr_switches_winpct == NULL))
    {
        puts("calloc failed");
        free(arr_stays_winpct);
        free(arr_switches_winpct);
        return -3;
    }

    for (int i = 0; i < runs; i++)
    {;
        switch_win = switch_lose = stay_win = stay_lose = 0;

        for (int j = 0; j < n; j++)
        {
            unsigned car = getrand_cdf(num_doors, NULL);  //which door the car is behind
            unsigned pick1 = getrand_cdf(num_doors, NULL); //initial pick of contestant
            unsigned door;

            if (car == pick1)
            {
                //contestant picks the car. Monty can show any other door
                //Randomly pick the next door or second door, wrap around if needed
                door = 1 << ((car + 1 + getrand(num_doors - 1)) % num_doors);

                //Alternatively Monty could always pick the next door
                //It will not alter the outcome
                //door = 1 << ((car + 1) % NUM_DOORS);
            }
            else
            {
                //Contestant picked goat. Monty has to show the other goat
                //convert door number to bit, not car or pick
                door = door_mask & ~((1 << car) | (1 << pick1));
                //printf("car = %u, pick1 = %u, door = %u\n", car, pick1, door);
                if (num_doors > 3)
                    door = pick_rnd_door(door, num_doors - 2, NULL);
                //printf("xcar = %u, pick1 = %u, door = %u\n",1<<car, 1<<pick1, door);
                
            }
            //bit 0 = door 0, bit 1 = door 1, bit 2 = door 2, etc.
            car = 1 << car;
            pick1 = 1 << pick1;

            //Let's randomly decide if contestant switches or stays
            unsigned stay = getrand(2);

            if (stay == 0)
            {
                // Case when contestant switches. The second pick is not pick1 or door
                unsigned pick2 = door_mask & ~(pick1 | door);

                if (num_doors > 3)
                    pick2 = pick_rnd_door(pick2, num_doors - 2, NULL);

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
        double mean_switch = sample_mean(arr_switches_winpct, cnt_switches_winpct);
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

        double mean_stay = sample_mean(arr_stays_winpct, cnt_stays_winpct);
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
