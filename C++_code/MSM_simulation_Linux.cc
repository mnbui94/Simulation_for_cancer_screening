#include <R.h>
#include <Rinternals.h>
#include <random>
#include <iostream>
#include </home/zcahnmb/LE_3D/eigen/Eigen/Eigenvalues>
#include </home/zcahnmb/LE_3D/eigen/Eigen/Dense>
using namespace Eigen;

extern "C"
{
    void transition_simulation_R_seed(double *tmpt, int i, double *t, double lambda, double beta, double random_seed, int distr, int index)
    {
        switch (distr)
        {
        case 0: // Exponential hazard
            // Rprintf("exp case \n");
            tmpt[i] = t[index] + random_seed;
            // Rprintf("random_seed = %f \n", random_seed);
            break;
        case 1: // Weibull distribution
            tmpt[i] = exp(exp(-beta) * log(-exp(-lambda) * log(random_seed) + pow(t[index], exp(beta))));
            break;
        default: // Gompertz distribution
            if (beta > 0)
            {
                tmpt[i] = (log(exp(lambda + beta * t[index]) - beta * log(random_seed)) - lambda) / beta;
            }
            else
            {
                // int j = 0;
                // while (j < 50)
                // {
                //         // Rprintf("i =  %d", i);
                //         if (U > exp(exp(lambda + beta * t[index]) / beta))
                //     {
                //         tmpt[i] = (log(exp(lambda + beta * t[index]) - beta * log(random_seed[j])) - lambda) / beta;
                //         break;
                //     }
                //     j++;
                // }
                // if (j == 50)
                // {
                // Rprintf("i %d", i);
                tmpt[i] = t[index] + (1.0 / 365);
                // }
            }

            break;
        }

        if ((tmpt[i] - t[index]) < (1.0 / 365))
        {
            tmpt[i] = t[index] + (1.0 / 365);
        }
    }
    void one_natural_history_R_seed(int i, int *cum_index, int *many_state, double *many_t, double *lambda, double *beta, double *t0, double *random_seed, int *distr)
    {
        // int max_trans = 4;
        int *state = new int[4];
        double *t = new double[4];
        double *tmpt = new double[2];
        int *index = new int[1];

        state[0] = 1;
        t[0] = t0[0];

        {
            index[0] = 0;
            transition_simulation_R_seed(tmpt, 0, t, lambda[0], beta[0], random_seed[i * 6 + 0], distr[0], index[0]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[i * 6 + 1], distr[1], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 2;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 4;
            }
        }
        if (state[index[0] + 1] == 2)
        {
            index[0] = 1;
            transition_simulation_R_seed(tmpt, 0, t, lambda[2], beta[2], random_seed[i * 6 + 2], distr[2], index[0]);
            // Rprintf("random_seed[i * 6 + 2]= %f \n", random_seed[i * 6 + 2]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[i * 6 + 3], distr[1], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 3;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 4;
            }
        }
        if (state[index[0] + 1] == 3)
        {
            index[0] = 2;
            transition_simulation_R_seed(tmpt, 0, t, lambda[3], beta[3], random_seed[i * 6 + 4], distr[3], index[0]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[4], beta[4], random_seed[i * 6 + 5], distr[4], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 4;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 5;
            }
        }
        cum_index[i + 1] = cum_index[i] + index[0] + 2;
        for (int row = 0; row < (index[0] + 2); row++)
        {
            many_state[cum_index[i] + row] = state[row];
            many_t[cum_index[i] + row] = t[row];
        }
        delete[] state;
        delete[] t;

        // output
        // SEXP ans;
        // PROTECT(ans = allocMatrix(REALSXP, index[0] + 2, 2));
        // double *rans = REAL(ans);

        // UNPROTECT(1);
        // return (ans);
    }

    int check_no_cancer(int *index, int *state, int *year, int Y1)
    {
        int result = 0; // FALSE
        int min_year = year[index[0] + 1];
        if (state[index[0]] == 3)
        {
            min_year = year[index[0]];
        }

        if (min_year > Y1)
        {
            result = 1; // TRUE
        }
        return (result);
    }
    void one_natural_history_conditioned_R_seed(int i, int *cum_index, int *many_state, double *many_t, int *many_years, double *lambda, double *beta, double *t0, double *random_seed, int *distr, int Y0, int Y1, int *M, int max_seeds)
    {
        // int max_trans = 4;
        int *state = new int[4];
        int *year = new int[4];
        double *t = new double[4];
        double *tmpt = new double[2];
        int *index = new int[1];

        state[0] = 1;
        t[0] = t0[0];

        {
            index[0] = 0;
            transition_simulation_R_seed(tmpt, 0, t, lambda[0], beta[0], random_seed[6 * M[0] + 0], distr[0], index[0]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[6 * M[0] + 1], distr[1], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 2;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 4;
            }
        }
        if (state[index[0] + 1] == 2)
        {
            index[0] = 1;
            transition_simulation_R_seed(tmpt, 0, t, lambda[2], beta[2], random_seed[6 * M[0] + 2], distr[2], index[0]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[6 * M[0] + 3], distr[1], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 3;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 4;
            }
        }
        if (state[index[0] + 1] == 3)
        {
            index[0] = 2;
            transition_simulation_R_seed(tmpt, 0, t, lambda[3], beta[3], random_seed[6 * M[0] + 4], distr[3], index[0]);
            transition_simulation_R_seed(tmpt, 1, t, lambda[4], beta[4], random_seed[6 * M[0] + 5], distr[4], index[0]);
            if (tmpt[0] < tmpt[1])
            {
                t[index[0] + 1] = tmpt[0];
                state[index[0] + 1] = 4;
            }
            else
            {
                t[index[0] + 1] = tmpt[1];
                state[index[0] + 1] = 5;
            }
        }
        year[0] = Y0;
        for (int j = 1; j < (index[0] + 2); j++)
        {
            year[j] = year[0] + round(t[j] - t[0]);
        }
        M[0]++;
        while (check_no_cancer(index, state, year, Y1) == 0)
        {
            {
                index[0] = 0;
                transition_simulation_R_seed(tmpt, 0, t, lambda[0], beta[0], random_seed[6 * M[0] + 0], distr[0], index[0]);
                transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[6 * M[0] + 1], distr[1], index[0]);
                if (tmpt[0] < tmpt[1])
                {
                    t[index[0] + 1] = tmpt[0];
                    state[index[0] + 1] = 2;
                }
                else
                {
                    t[index[0] + 1] = tmpt[1];
                    state[index[0] + 1] = 4;
                }
            }
            if (state[index[0] + 1] == 2)
            {
                index[0] = 1;
                transition_simulation_R_seed(tmpt, 0, t, lambda[2], beta[2], random_seed[6 * M[0] + 2], distr[2], index[0]);
                transition_simulation_R_seed(tmpt, 1, t, lambda[1], beta[1], random_seed[6 * M[0] + 3], distr[1], index[0]);
                if (tmpt[0] < tmpt[1])
                {
                    t[index[0] + 1] = tmpt[0];
                    state[index[0] + 1] = 3;
                }
                else
                {
                    t[index[0] + 1] = tmpt[1];
                    state[index[0] + 1] = 4;
                }
            }
            if (state[index[0] + 1] == 3)
            {
                index[0] = 2;
                transition_simulation_R_seed(tmpt, 0, t, lambda[3], beta[3], random_seed[6 * M[0] + 4], distr[3], index[0]);
                transition_simulation_R_seed(tmpt, 1, t, lambda[4], beta[4], random_seed[6 * M[0] + 5], distr[4], index[0]);
                if (tmpt[0] < tmpt[1])
                {
                    t[index[0] + 1] = tmpt[0];
                    state[index[0] + 1] = 4;
                }
                else
                {
                    t[index[0] + 1] = tmpt[1];
                    state[index[0] + 1] = 5;
                }
            }
            year[0] = Y0;
            for (int j = 1; j < (index[0] + 2); j++)
            {
                year[j] = year[0] + round(t[j] - t[0]);
            }
            M[0]++;
            if (M[0] == max_seeds)
            {
                break;
            }
        }

        cum_index[i + 1] = cum_index[i] + index[0] + 2;
        for (int row = 0; row < (index[0] + 2); row++)
        {
            many_state[cum_index[i] + row] = state[row];
            many_t[cum_index[i] + row] = t[row];
            many_years[cum_index[i] + row] = year[row];
        }
        delete[] state;
        delete[] t;
        delete[] year;
    }

    SEXP R_many_natural_history_R_seed(SEXP R_lambda, SEXP R_beta, SEXP left_age, SEXP initial_years, SEXP hazard_distr, SEXP size, SEXP id_seeds, SEXP R_seeds)
    {
        // input
        double *lambda = REAL(R_lambda);
        double *beta = REAL(R_beta);
        double *t0 = REAL(left_age);
        int *Y0 = INTEGER(initial_years);
        int *id = INTEGER(id_seeds);
        double *random_seed = REAL(R_seeds);
        int *distr = INTEGER(hazard_distr);
        int *N = INTEGER(size);
        // N[0] = 1;

        int *many_years = new int[4 * N[0]];
        int *many_state = new int[4 * N[0]];
        double *many_t = new double[4 * N[0]];
        int *many_id = new int[4 * N[0]];
        int *cum_index = new int[N[0] + 1];
        cum_index[0] = 0;
        for (int i = 0; i < N[0]; i++)
        {
            one_natural_history_R_seed(i, cum_index, many_state, many_t, lambda, beta, t0, random_seed, distr);
            many_years[cum_index[i]] = Y0[i];
            many_id[cum_index[i]] = id[i];
            for (int j = (cum_index[i] + 1); j < cum_index[i + 1]; j++)
            {
                many_id[j] = id[i];
                many_years[j] = many_years[cum_index[i]] + round(many_t[j] - many_t[cum_index[i]]);
            }
        }

        // output
        SEXP ans;
        PROTECT(ans = allocMatrix(REALSXP, cum_index[N[0]], 4));
        double *rans = REAL(ans);
        for (int row = 0; row < (cum_index[N[0]]); row++)
        {
            rans[(cum_index[N[0]]) * 0 + row] = many_id[row];
            rans[(cum_index[N[0]]) * 1 + row] = many_state[row];
            rans[(cum_index[N[0]]) * 2 + row] = many_years[row];
            rans[(cum_index[N[0]]) * 3 + row] = many_t[row];
        }
        // delete[] m]state;
        // delete[] t;
        delete[] many_years;
        delete[] many_state;
        delete[] many_t;
        delete[] many_id;
        delete[] cum_index;
        UNPROTECT(1);
        return (ans);
    }

    SEXP R_many_simulation_control_R_seed(SEXP R_lambda, SEXP R_beta, SEXP left_age, SEXP initial_years, SEXP study_years, SEXP hazard_distr, SEXP size, SEXP id_seeds, SEXP R_seeds, SEXP max_no_seeds)
    {
        // input
        double *lambda = REAL(R_lambda);
        double *beta = REAL(R_beta);
        double *t0 = REAL(left_age);
        int *Y0 = INTEGER(initial_years);
        int *id = INTEGER(id_seeds);
        double *random_seed = REAL(R_seeds);
        int *distr = INTEGER(hazard_distr);
        int *N = INTEGER(size);
        int *Y1 = INTEGER(study_years);
        int *max_seeds = INTEGER(max_no_seeds);

        // nat_hist informations
        int *many_years = new int[4 * N[0]];
        int *many_state = new int[4 * N[0]];
        double *many_t = new double[4 * N[0]];
        int *many_id = new int[4 * N[0]];
        int *cum_index = new int[N[0] + 1];
        cum_index[0] = 0;

        // full_record
        int *full_many_years = new int[3 * N[0]];
        int *full_many_state = new int[3 * N[0]];
        double *full_many_t = new double[3 * N[0]];
        int *full_many_id = new int[3 * N[0]];
        int *full_cum_index = new int[N[0] + 1];
        full_cum_index[0] = 0;

        // study_record
        int *record_many_years = new int[3 * N[0]];
        int *record_many_state = new int[3 * N[0]];
        double *record_many_t = new double[3 * N[0]];
        int *record_many_id = new int[3 * N[0]];
        int *record_cum_index = new int[N[0] + 1];
        record_cum_index[0] = 0;

        int s = 0;
        int *M = new int[1]; // count number of seeds were used
        M[0] = 0;
        for (int i = 0; i < N[0]; i++)
        {
            // keep only natural history simulation if no_cancer and alive by study_years (Y1)
            one_natural_history_conditioned_R_seed(i, cum_index, many_state, many_t, many_years, lambda, beta, t0, random_seed, distr, Y0[i], Y1[0], M, max_seeds[0] * N[0]);
            many_id[cum_index[i]] = id[i];
            for (int j = (cum_index[i] + 1); j < cum_index[i + 1]; j++)
            {
                many_id[j] = id[i];
            }
            if (M[0] < (N[0] * max_seeds[0] + 1))
            {
                // extracting full record
                full_many_t[full_cum_index[i]] = many_t[cum_index[i]] + Y1[0] - Y0[i];
                if (full_many_t[full_cum_index[i]] == t0[0])
                {
                    full_many_state[full_cum_index[i]] = 1;
                }
                else
                {
                    full_many_state[full_cum_index[i]] = -1;
                }

                if (many_state[cum_index[i] + 1] == 2)
                {
                    full_cum_index[i + 1] = full_cum_index[i] + (cum_index[i + 1] - cum_index[i] - 1);
                    for (int j = 1; j < (cum_index[i + 1] - cum_index[i] - 1); j++)
                    {
                        full_many_state[full_cum_index[i] + j] = many_state[cum_index[i] + j + 1];
                        full_many_t[full_cum_index[i] + j] = many_t[cum_index[i] + j + 1];
                    }
                }
                else
                {
                    full_cum_index[i + 1] = full_cum_index[i] + (cum_index[i + 1] - cum_index[i]);
                    for (int j = 1; j < (cum_index[i + 1] - cum_index[i]); j++)
                    {
                        full_many_state[full_cum_index[i] + j] = many_state[cum_index[i] + j];
                        full_many_t[full_cum_index[i] + j] = many_t[cum_index[i] + j];
                    }
                }

                full_many_years[full_cum_index[i]] = Y1[0];
                full_many_id[full_cum_index[i]] = id[i];
                for (int j = (full_cum_index[i] + 1); j < full_cum_index[i + 1]; j++)
                {
                    full_many_id[j] = id[i];
                    full_many_years[j] = full_many_years[full_cum_index[i]] + round(full_many_t[j] - full_many_t[full_cum_index[i]]);
                }

                // extracting study_record
                s = 0;
                // Rprintf("year = %d", Y1[1]); %d = integer, %f = double
                while ((s < (full_cum_index[i + 1] - full_cum_index[i])) & (full_many_years[full_cum_index[i] + s] < Y1[1]))
                {
                    record_many_t[record_cum_index[i] + s] = full_many_t[full_cum_index[i] + s];
                    record_many_state[record_cum_index[i] + s] = full_many_state[full_cum_index[i] + s];
                    record_many_years[record_cum_index[i] + s] = full_many_years[full_cum_index[i] + s];
                    s++;
                }
                // Rprintf("s = %d \n", s);
                switch (record_many_state[record_cum_index[i] + s - 1])
                {
                case 5:
                    record_cum_index[i + 1] = record_cum_index[i] + s;
                    break;
                case 4:
                    record_cum_index[i + 1] = record_cum_index[i] + s;
                    break;
                case 3:
                    record_cum_index[i + 1] = record_cum_index[i] + s + 1;
                    record_many_state[record_cum_index[i] + s] = 3;
                    record_many_years[record_cum_index[i] + s] = Y1[1];
                    record_many_t[record_cum_index[i] + s] = many_t[cum_index[i]] + (Y1[1] - many_years[cum_index[i]]);
                    break;

                default:
                    record_cum_index[i + 1] = record_cum_index[i] + s + 1;
                    record_many_state[record_cum_index[i] + s] = -2;
                    record_many_years[record_cum_index[i] + s] = Y1[1];
                    record_many_t[record_cum_index[i] + s] = many_t[cum_index[i]] + (Y1[1] - many_years[cum_index[i]]);
                    break;
                }

                for (int j = (record_cum_index[i]); j < record_cum_index[i + 1]; j++)
                {
                    record_many_id[j] = id[i];
                }
            }
            else
            {
                Rprintf("Warning: maximum numbers of seeds need to be increased !!! \n");
                break;
            }
        }
        // Rprintf("total number of seeds used = %d \n", M[0]);
        // output

        SEXP ans = PROTECT(allocVector(VECSXP, 3));
        SEXP natural_history = PROTECT(allocMatrix(REALSXP, cum_index[N[0]], 4));
        SEXP life_time_record = PROTECT(allocMatrix(REALSXP, full_cum_index[N[0]], 4));
        SEXP study_record = PROTECT(allocMatrix(REALSXP, record_cum_index[N[0]], 4));
        double *nat_hist = REAL(natural_history);
        double *full_record = REAL(life_time_record);
        double *record = REAL(study_record);

        for (int row = 0; row < (cum_index[N[0]]); row++)
        {
            nat_hist[(cum_index[N[0]]) * 0 + row] = many_id[row];
            nat_hist[(cum_index[N[0]]) * 1 + row] = many_state[row];
            nat_hist[(cum_index[N[0]]) * 2 + row] = many_years[row];
            nat_hist[(cum_index[N[0]]) * 3 + row] = many_t[row];
        }
        for (int row = 0; row < (full_cum_index[N[0]]); row++)
        {
            full_record[(full_cum_index[N[0]]) * 0 + row] = full_many_id[row];
            full_record[(full_cum_index[N[0]]) * 1 + row] = full_many_state[row];
            full_record[(full_cum_index[N[0]]) * 2 + row] = full_many_years[row];
            full_record[(full_cum_index[N[0]]) * 3 + row] = full_many_t[row];
        }

        for (int row = 0; row < (record_cum_index[N[0]]); row++)
        {
            record[(record_cum_index[N[0]]) * 0 + row] = record_many_id[row];
            record[(record_cum_index[N[0]]) * 1 + row] = record_many_state[row];
            record[(record_cum_index[N[0]]) * 2 + row] = record_many_years[row];
            record[(record_cum_index[N[0]]) * 3 + row] = record_many_t[row];
        }

        delete[] M;
        delete[] many_years;
        delete[] many_state;
        delete[] many_t;
        delete[] many_id;
        delete[] cum_index;

        delete[] full_many_years;
        delete[] full_many_state;
        delete[] full_many_t;
        delete[] full_many_id;
        delete[] full_cum_index;

        delete[] record_many_years;
        delete[] record_many_state;
        delete[] record_many_t;
        delete[] record_many_id;
        delete[] record_cum_index;

        SET_VECTOR_ELT(ans, 0, natural_history);
        SET_VECTOR_ELT(ans, 1, life_time_record);
        SET_VECTOR_ELT(ans, 2, study_record);
        UNPROTECT(4);
        return (ans);
    }

    void screen_without_state2(int i, int s, double *t0, int *full_cum_index, int *full_many_state, double *full_many_t, int *cum_index, int *many_state, double *many_t, double *when_to_screen, double attend_prob, int no_seed, int no_screen)
    {
        auto generator = std::mt19937{no_seed};
        std::bernoulli_distribution ber_distribution_attend(attend_prob);

        if (full_many_t[full_cum_index[i]] < when_to_screen[0])
        {
            full_cum_index[i + 1] = full_cum_index[i] + 1;
            if (full_many_t[full_cum_index[i]] == t0[0])
            {
                full_many_state[full_cum_index[i]] = 1;
            }
            else
            {
                full_many_state[full_cum_index[i]] = -1;
            }
        }
        else
        {
            if (full_many_t[full_cum_index[i]] != when_to_screen[s - 1])
            {
                if ((full_many_t[full_cum_index[i]] == when_to_screen[s]) & (s < no_screen))
                {
                    full_cum_index[i + 1] = full_cum_index[i];
                }
                else
                {
                    full_cum_index[i + 1] = full_cum_index[i] + 1;
                    full_many_state[full_cum_index[i + 1] - 1] = -1;
                    s++;
                }
            }
            else
            {
                full_cum_index[i + 1] = full_cum_index[i];
            }
        }
        int l = s;
        while (l < no_screen)
        {
            if (ber_distribution_attend(generator) == 1)
            {
                full_cum_index[i + 1]++;
                full_many_t[full_cum_index[i + 1] - 1] = when_to_screen[l];
                if (full_many_t[full_cum_index[i + 1] - 1] < many_t[cum_index[i + 1] - 1])
                {
                    full_many_state[full_cum_index[i + 1] - 1] = 1;
                }
                else
                {
                    full_many_state[full_cum_index[i + 1] - 1] = 4;
                    // Rprintf("case 2\n");
                    break;
                }
            }

            l++;
        }
        if (full_many_state[full_cum_index[i + 1] - 1] == 4)
        {
            full_many_t[full_cum_index[i + 1] - 1] = many_t[cum_index[i + 1] - 1];
        }
        else
        {
            full_cum_index[i + 1]++;
            full_many_t[full_cum_index[i + 1] - 1] = many_t[cum_index[i + 1] - 1];
            full_many_state[full_cum_index[i + 1] - 1] = 4;
        }
    }

    void screen_with_state2_R_seed(int i, int s, double *t0, int *full_cum_index, int *full_many_state, double *full_many_t, int *cum_index, int *many_state, double *many_t, double *when_to_screen, int freq, double *misc, double attend_prob, double *sub_lambda, double *sub_beta, double *random_seed_sub, int *sub_distr, int no_screen, int *id)
    {
        double misc_prob = 0; // prob of getting true underlying state

        if (full_many_t[full_cum_index[i]] < when_to_screen[0])
        {
            full_cum_index[i + 1] = full_cum_index[i] + 1;
            if (full_many_t[full_cum_index[i]] == t0[0])
            {
                full_many_state[full_cum_index[i]] = 1;
            }
            else
            {
                full_many_state[full_cum_index[i]] = -1;
            }
        }
        else
        {
            if ((full_many_t[full_cum_index[i]] != when_to_screen[s - 1]))
            {
                if ((full_many_t[full_cum_index[i]] == when_to_screen[s]) & (s < no_screen))
                {
                    full_cum_index[i + 1] = full_cum_index[i];
                }
                else
                {
                    full_cum_index[i + 1] = full_cum_index[i] + 1;
                    full_many_state[full_cum_index[i + 1] - 1] = -1;
                    s++;
                }
            }
            else
            {
                full_cum_index[i + 1] = full_cum_index[i];
            }
        }

        int l = s;
        auto generator = std::mt19937{id[i]};
        std::bernoulli_distribution ber_distribution_attend(attend_prob);

        // step 1: determine when cancer is screen-detected.
        if (many_t[cum_index[i] + 1] > when_to_screen[no_screen - 1])
        // state 2 occurs after the end of screening study.
        {
            // Rprintf("case1\n");
            while (l < no_screen)
            {
                if (ber_distribution_attend(generator) == 1)
                {
                    full_cum_index[i + 1]++;
                    full_many_t[full_cum_index[i + 1] - 1] = when_to_screen[l];
                    if (full_many_t[full_cum_index[i + 1] - 1] < many_t[cum_index[i + 1] - 1])
                    {
                        full_many_state[full_cum_index[i + 1] - 1] = 1;
                    }
                }

                l++;
            }
            for (int j = cum_index[i] + 2; j < cum_index[i + 1]; j++)
            {
                full_cum_index[i + 1]++;
                full_many_t[full_cum_index[i + 1] - 1] = many_t[j];
                full_many_state[full_cum_index[i + 1] - 1] = many_state[j];
            }
        }
        // state 2 occurs before the end of screening study.
        else
        {
            // Rprintf("case2\n");
            while (l < no_screen)
            {
                if (ber_distribution_attend(generator) == 1)
                {
                    full_cum_index[i + 1]++;
                    full_many_t[full_cum_index[i + 1] - 1] = when_to_screen[l];

                    if (full_many_t[full_cum_index[i + 1] - 1] < many_t[cum_index[i] + 2])
                    {
                        // apply screening
                        if (full_many_t[full_cum_index[i + 1] - 1] < many_t[cum_index[i] + 1])
                        {
                            full_many_state[full_cum_index[i + 1] - 1] = 1;
                        }
                        else
                        {
                            // applying misclassification
                            misc_prob = 1 / (1 + exp(misc[0] + misc[1] * full_many_t[full_cum_index[i + 1] - 1]));
                            std::bernoulli_distribution ber_distribution_misc(misc_prob);
                            if (ber_distribution_misc(generator) == 1)
                            {
                                full_many_state[full_cum_index[i + 1] - 1] = 2;
                                break;
                            }
                            else
                            {
                                full_many_state[full_cum_index[i + 1] - 1] = 1;
                            }
                        }
                    }
                    else
                    {
                        full_cum_index[i + 1] = full_cum_index[i + 1] - 1;
                        for (int j = cum_index[i] + 2; j < cum_index[i + 1]; j++)
                        {
                            full_cum_index[i + 1]++;
                            full_many_t[full_cum_index[i + 1] - 1] = many_t[j];
                            full_many_state[full_cum_index[i + 1] - 1] = many_state[j];
                        }
                        break;
                    }
                }
                l++;
            }

            if (full_many_state[full_cum_index[i + 1] - 1] == 1)
            {
                for (int j = cum_index[i] + 2; j < cum_index[i + 1]; j++)
                {
                    full_cum_index[i + 1]++;
                    full_many_t[full_cum_index[i + 1] - 1] = many_t[j];
                    full_many_state[full_cum_index[i + 1] - 1] = many_state[j];
                }
            }
            else
            {
                if (full_many_state[full_cum_index[i + 1] - 1] == 2)
                {
                    double *tmpt = new double[2];
                    // std::uniform_int_distribution<int> unif_int_distribution(2 * (no_seed - 1), 2 * no_seed);
                    transition_simulation_R_seed(tmpt, 0, full_many_t, sub_lambda[0], sub_beta[0], random_seed_sub[2 * i + 0], sub_distr[0], full_cum_index[i + 1] - 1);
                    transition_simulation_R_seed(tmpt, 1, full_many_t, sub_lambda[1], sub_beta[1], random_seed_sub[2 * i + 1], sub_distr[1], full_cum_index[i + 1] - 1);
                    full_cum_index[i + 1]++;
                    if (tmpt[0] < tmpt[1])
                    {
                        full_many_t[full_cum_index[i + 1] - 1] = tmpt[0];
                        full_many_state[full_cum_index[i + 1] - 1] = 4;
                    }
                    else
                    {
                        full_many_t[full_cum_index[i + 1] - 1] = tmpt[1];
                        full_many_state[full_cum_index[i + 1] - 1] = 5;
                    }
                    delete[] tmpt;
                }
            }
        }
        if (full_cum_index[i + 1] == (full_cum_index[i] + 1))
        {
            full_cum_index[i + 1]++;
            full_many_t[full_cum_index[i + 1] - 1] = many_t[cum_index[i + 1] - 1];
            full_many_state[full_cum_index[i + 1] - 1] = many_state[cum_index[i + 1] - 1];
        }
    }

    SEXP R_many_simulation_screen_freq_R_seed(SEXP R_lambda, SEXP R_beta, SEXP R_sub_lambda, SEXP R_sub_beta, SEXP R_misclassification, SEXP left_age, SEXP initial_years, SEXP study_years, SEXP screening_time_range, SEXP frequency_screening, SEXP max_no_screen, SEXP screen_take_up, SEXP hazard_distr, SEXP hazard_sub_distr, SEXP size, SEXP id_seeds, SEXP R_seeds, SEXP R_seeds_sub, SEXP max_no_seeds)
    {
        // input
        double *lambda = REAL(R_lambda);
        double *beta = REAL(R_beta);
        double *sub_lambda = REAL(R_sub_lambda);
        double *sub_beta = REAL(R_sub_beta);
        double *misc = REAL(R_misclassification);
        double *t0 = REAL(left_age);
        int *Y0 = INTEGER(initial_years);
        int *Y1 = INTEGER(study_years);
        double *screening_t = REAL(screening_time_range);
        int *freq = INTEGER(frequency_screening);
        int *no_screen = INTEGER(max_no_screen);
        // no_screen[0] += 3;
        double *attend_prob = REAL(screen_take_up);
        int *id = INTEGER(id_seeds);
        double *random_seed = REAL(R_seeds);
        double *random_seed_sub = REAL(R_seeds_sub);
        int *distr = INTEGER(hazard_distr);
        int *sub_distr = INTEGER(hazard_sub_distr);
        int *N = INTEGER(size);
        int *max_seeds = INTEGER(max_no_seeds);

        // at which time screening is offered
        double *when_to_screen = new double[no_screen[0]];
        for (int s = 0; s < no_screen[0]; s++)
        {
            when_to_screen[s] = screening_t[0] + s * freq[0];
        }
        // Rprintf("when_to_screen[0] = %f \n", when_to_screen[0]);
        // Rprintf("when_to_screen[no_screen[0] - 1] = %f \n", when_to_screen[no_screen[0] - 1]);

        // nat_hist informations
        int *many_years = new int[4 * N[0]];
        int *many_state = new int[4 * N[0]];
        double *many_t = new double[4 * N[0]];
        int *many_id = new int[4 * N[0]];
        int *cum_index = new int[N[0] + 1];
        cum_index[0] = 0;

        // full_record
        int *full_many_years = new int[(no_screen[0] + 3) * N[0]];
        int *full_many_state = new int[(no_screen[0] + 3) * N[0]];
        double *full_many_t = new double[(no_screen[0] + 3) * N[0]];
        int *full_many_id = new int[(no_screen[0] + 3) * N[0]];
        int *full_cum_index = new int[N[0] + 1];
        full_cum_index[0] = 0;
        int *full_screen = new int[(no_screen[0] + 3) * N[0]];

        // study_record
        int *record_many_years = new int[(no_screen[0] + 3) * N[0]];
        int *record_many_state = new int[(no_screen[0] + 3) * N[0]];
        double *record_many_t = new double[(no_screen[0] + 3) * N[0]];
        int *record_many_id = new int[(no_screen[0] + 3) * N[0]];
        int *record_cum_index = new int[N[0] + 1];
        record_cum_index[0] = 0;
        int *record_screen = new int[(no_screen[0] + 3) * N[0]];

        int s = 0;
        int *M = new int[1]; // count number of seeds were used
        M[0] = 0;
        for (int i = 0; i < N[0]; i++)
        {
            // keep only natural history simulation if no_cancer and alive by study_years (Y1)
            one_natural_history_conditioned_R_seed(i, cum_index, many_state, many_t, many_years, lambda, beta, t0, random_seed, distr, Y0[i], Y1[0], M, max_seeds[0] * N[0]);
            for (int j = (cum_index[i]); j < cum_index[i + 1]; j++)
            {
                many_id[j] = id[i];
            }
            if (M[0] < (N[0] * max_seeds[0] + 1))
            {
                // implement screening between testing age range and follow over life-time
                full_many_t[full_cum_index[i]] = many_t[cum_index[i]] + Y1[0] - Y0[i];
                // if (full_many_t[full_cum_index[i]] == t0[0])
                // {
                //     full_many_state[full_cum_index[i]] = 1;
                // }
                // else
                // {
                //     full_many_state[full_cum_index[i]] = -1;
                // }
                if (full_many_t[full_cum_index[i]] > screening_t[1])
                // after eligible screening age
                {
                    // Rprintf("case1\n");

                    if (many_state[cum_index[i] + 1] == 2)
                    {
                        full_cum_index[i + 1] = full_cum_index[i] + (cum_index[i + 1] - cum_index[i] - 1);
                        for (int j = 1; j < (cum_index[i + 1] - cum_index[i] - 1); j++)
                        {
                            full_many_state[full_cum_index[i] + j] = many_state[cum_index[i] + j + 1];
                            full_many_t[full_cum_index[i] + j] = many_t[cum_index[i] + j + 1];
                        }
                    }
                    else
                    {
                        full_cum_index[i + 1] = full_cum_index[i] + (cum_index[i + 1] - cum_index[i]);
                        for (int j = 1; j < (cum_index[i + 1] - cum_index[i]); j++)
                        {
                            full_many_state[full_cum_index[i] + j] = many_state[cum_index[i] + j];
                            full_many_t[full_cum_index[i] + j] = many_t[cum_index[i] + j];
                        }
                    }
                }
                else // before eligible screening
                {
                    s = 0; // find when to start screening
                    while (s < no_screen[0])
                    {
                        if ((full_many_t[full_cum_index[i]] == when_to_screen[s]) || (full_many_t[full_cum_index[i]] < when_to_screen[0]) || (full_many_t[full_cum_index[i]] < when_to_screen[s + 1]))
                        {
                            break;
                        }
                        else
                        {
                            s++;
                        }
                    }

                    // Rprintf("s= %d\n", s);

                    switch (many_state[cum_index[i] + 1])
                    {
                    case 4:
                        // Rprintf("Case1\n");
                        screen_without_state2(i, s, t0, full_cum_index, full_many_state, full_many_t, cum_index, many_state, many_t, when_to_screen, attend_prob[0], id[i], no_screen[0]);
                        break;

                    default: // state 2
                        // Rprintf("Case2\n");
                        screen_with_state2_R_seed(i, s, t0, full_cum_index, full_many_state, full_many_t, cum_index, many_state, many_t, when_to_screen, freq[0], misc, attend_prob[0], sub_lambda, sub_beta, random_seed_sub, sub_distr, no_screen[0], id);
                        break;
                    }
                    // }
                }

                full_many_id[full_cum_index[i]] = id[i];
                full_many_years[full_cum_index[i]] = Y0[i] + round(full_many_t[full_cum_index[i]] - round(t0[0]));
                for (int j = (full_cum_index[i] + 1); j < full_cum_index[i + 1]; j++)
                {
                    full_many_id[j] = id[i];
                    full_many_years[j] = full_many_years[full_cum_index[i]] + round(full_many_t[j] - round(full_many_t[full_cum_index[i]]));
                }

                // extracting study_record
                s = 0;
                // Rprintf("year = %d", Y1[1]); %d = integer, %f = double
                while ((s < (full_cum_index[i + 1] - full_cum_index[i])) & (full_many_years[full_cum_index[i] + s] < Y1[1]))
                {
                    record_many_t[record_cum_index[i] + s] = full_many_t[full_cum_index[i] + s];
                    record_many_state[record_cum_index[i] + s] = full_many_state[full_cum_index[i] + s];
                    record_many_years[record_cum_index[i] + s] = full_many_years[full_cum_index[i] + s];
                    s++;
                }
                // Rprintf("s = %d \n", s);
                switch (record_many_state[record_cum_index[i] + s - 1])
                {
                case 5:
                    record_cum_index[i + 1] = record_cum_index[i] + s;
                    break;
                case 4:
                    record_cum_index[i + 1] = record_cum_index[i] + s;
                    break;
                case 1:
                    record_cum_index[i + 1] = record_cum_index[i] + s + 1;
                    record_many_state[record_cum_index[i] + s] = -2;
                    record_many_years[record_cum_index[i] + s] = Y1[1];
                    record_many_t[record_cum_index[i] + s] = many_t[cum_index[i]] + (Y1[1] - many_years[cum_index[i]]);
                    break;
                default:
                    record_cum_index[i + 1] = record_cum_index[i] + s + 1;
                    record_many_state[record_cum_index[i] + s] = record_many_state[record_cum_index[i] + s - 1];
                    record_many_years[record_cum_index[i] + s] = Y1[1];
                    record_many_t[record_cum_index[i] + s] = many_t[cum_index[i]] + (Y1[1] - many_years[cum_index[i]]);
                    break;
                }

                for (int j = (record_cum_index[i]); j < record_cum_index[i + 1]; j++)
                {
                    record_many_id[j] = id[i];
                }

                // adding information about whether screening is performed
                switch (full_many_state[full_cum_index[i]])
                {
                case -1:
                    if ((full_many_state[full_cum_index[i] + 1] == 1) || (full_many_state[full_cum_index[i] + 1] == 2))
                    {
                        full_screen[full_cum_index[i]] = 1;
                        record_screen[record_cum_index[i]] = 1;
                    }
                    else
                    {
                        full_screen[full_cum_index[i]] = 0;
                        record_screen[record_cum_index[i]] = 0;
                    }
                    break;
                default:
                    full_screen[full_cum_index[i]] = 1;
                    record_screen[record_cum_index[i]] = 1;
                    break;
                }
                for (int j = full_cum_index[i] + 1; j < full_cum_index[i + 1]; j++)
                {
                    full_screen[j] = full_screen[full_cum_index[i]];
                }
                for (int j = record_cum_index[i] + 1; j < record_cum_index[i + 1]; j++)
                {
                    record_screen[j] = record_screen[record_cum_index[i]];
                }
            }
            else
            {
                Rprintf("Warning: maximum numbers of seeds need to be increased !!! \n");
                break;
            }
        }
        // Rprintf("total number of seeds used = %d \n", M[0]);
        // output

        SEXP ans = PROTECT(allocVector(VECSXP, 3));
        SEXP natural_history = PROTECT(allocMatrix(REALSXP, cum_index[N[0]], 4));
        SEXP life_time_record = PROTECT(allocMatrix(REALSXP, full_cum_index[N[0]], 5));
        SEXP study_record = PROTECT(allocMatrix(REALSXP, record_cum_index[N[0]], 5));
        double *nat_hist = REAL(natural_history);
        double *full_record = REAL(life_time_record);
        double *record = REAL(study_record);

        for (int row = 0; row < (cum_index[N[0]]); row++)
        {
            nat_hist[(cum_index[N[0]]) * 0 + row] = many_id[row];
            nat_hist[(cum_index[N[0]]) * 1 + row] = many_state[row];
            nat_hist[(cum_index[N[0]]) * 2 + row] = many_years[row];
            nat_hist[(cum_index[N[0]]) * 3 + row] = many_t[row];
        }
        for (int row = 0; row < (full_cum_index[N[0]]); row++)
        {
            full_record[(full_cum_index[N[0]]) * 0 + row] = full_many_id[row];
            full_record[(full_cum_index[N[0]]) * 1 + row] = full_many_state[row];
            full_record[(full_cum_index[N[0]]) * 2 + row] = full_many_years[row];
            full_record[(full_cum_index[N[0]]) * 3 + row] = full_many_t[row];
            full_record[(full_cum_index[N[0]]) * 4 + row] = full_screen[row];
        }
        for (int row = 0; row < (record_cum_index[N[0]]); row++)
        {
            record[(record_cum_index[N[0]]) * 0 + row] = record_many_id[row];
            record[(record_cum_index[N[0]]) * 1 + row] = record_many_state[row];
            record[(record_cum_index[N[0]]) * 2 + row] = record_many_years[row];
            record[(record_cum_index[N[0]]) * 3 + row] = record_many_t[row];
            record[(record_cum_index[N[0]]) * 4 + row] = record_screen[row];
        }

        delete[] many_years;
        delete[] many_state;
        delete[] many_t;
        delete[] many_id;
        delete[] cum_index;

        delete[] full_many_years;
        delete[] full_many_state;
        delete[] full_many_t;
        delete[] full_many_id;
        delete[] full_cum_index;
        delete[] full_screen;

        delete[] record_many_years;
        delete[] record_many_state;
        delete[] record_many_t;
        delete[] record_many_id;
        delete[] record_cum_index;
        delete[] record_screen;

        SET_VECTOR_ELT(ans, 0, natural_history);
        SET_VECTOR_ELT(ans, 1, life_time_record);
        SET_VECTOR_ELT(ans, 2, study_record);
        UNPROTECT(4);
        return (ans);
    }
}

// open command promt= ctrl+shift+C
// compile : R CMD SHLIB Biometrical/C++/MSM_simulation_Linux.cc -shared -o Biometrical/C++/MSM_simulation_Linux.so
// in R
// .Call("myFunc",body(f),as.double(x),new.env())}