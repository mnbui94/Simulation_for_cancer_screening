#include <R.h>
#include <Rinternals.h>
#include </home/zcahnmb/LE_3D/eigen/Eigen/Eigenvalues>
#include </home/zcahnmb/LE_3D/eigen/Eigen/Dense>
using namespace Eigen;

extern "C"
{
    void test_false_grid_4States(double *T_matrix, double *a, double *b, int m1, int m2, double *t, int *S, int k, double grid_pw, double t0, double min_del)
    {
        ArrayXXd L(k, k);
        MatrixXd Q = MatrixXd::Zero(k, k);
        MatrixXcd P;
        MatrixXcd V(k, k);
        MatrixXcd D(k, k);
        MatrixXcd tmpt = MatrixXd::Zero(k, k);
        MatrixXcd T_j = MatrixXd::Identity(k, k);
        double *tmpt_t = new double[1];

        P = MatrixXd::Identity(k, k);
        tmpt_t[0] = t[m1];
        while ((t[m1 + 1] - tmpt_t[0]) > (grid_pw + min_del))
        {
            Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
            Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(0, 0) = -(Q(0, 1) + Q(0, 3));

            Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
            Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(1, 1) = -(Q(1, 2) + Q(1, 3));

            EigenSolver<MatrixXd> es(Q);
            V = es.eigenvectors();
            D = es.eigenvalues().asDiagonal();
            for (int i = 0; i < k; i++)
                D(i, i) = exp(D(i, i) * grid_pw);
            P = P * V * D * V.inverse();
            tmpt_t[0] += grid_pw;
        }

        Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
        Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(0, 0) = -(Q(0, 1) + Q(0, 3));

        Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
        Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(1, 1) = -(Q(1, 2) + Q(1, 3));

        EigenSolver<MatrixXd> es(Q);
        V = es.eigenvectors();
        D = es.eigenvalues().asDiagonal();
        for (int i = 0; i < k; i++)
            D(i, i) = exp(D(i, i) * (t[m1 + 1] - tmpt_t[0]));
        // P = P * V * D * V.inverse();
        P = P * V * D * V.inverse();

        tmpt = MatrixXd::Zero(k, k);
        switch (S[m2 - 1])
        {
        case -2:
            for (int r = 0; r < 2; r++)
            {
                for (int s = 0; s < 2; s++)
                    tmpt(r, s) = P(r, s);
            }
            T_j = tmpt;
            // L = real(tmpt.array());
            break;
        case 3:
            tmpt(0, 0) = 0;
            // tmpt(1, 1) = exp(a[2] + b[2] * t[m2 - 1]);
            tmpt(1, 1) = Q(1, 2);
            T_j = P * tmpt;
            break;
        default:
            // tmpt(0, 0) = exp(a[1] + b[1] * t[m2 - 1]);
            // tmpt(1, 1) = exp(a[1] + b[1] * t[m2 - 1]);
            tmpt(0, 0) = Q(0, 3);
            tmpt(1, 1) = Q(1, 3);
            T_j = P * tmpt;
            break;
        }

        tmpt = MatrixXd::Zero(k, k);
        P = MatrixXd::Identity(k, k);
        tmpt_t[0] = t0;
        while ((t[m1] - tmpt_t[0]) > (grid_pw + min_del))
        {
            Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
            Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(0, 0) = -(Q(0, 1) + Q(0, 3));

            Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
            Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(1, 1) = -(Q(1, 2) + Q(1, 3));
            EigenSolver<MatrixXd> es(Q);
            V = es.eigenvectors();
            D = es.eigenvalues().asDiagonal();
            for (int i = 0; i < k; i++)
                D(i, i) = exp(D(i, i) * grid_pw);
            P = P * V * D * V.inverse();
            tmpt_t[0] += grid_pw;
        }

        Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
        Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(0, 0) = -(Q(0, 1) + Q(0, 3));

        Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
        Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(1, 1) = -(Q(1, 2) + Q(1, 3));

        EigenSolver<MatrixXd> n_es(Q);
        V = n_es.eigenvectors();
        D = n_es.eigenvalues().asDiagonal();
        for (int i = 0; i < k; i++)
            D(i, i) = exp(D(i, i) * (t[m1] - tmpt_t[0]));
        // P = P * V * D * V.inverse();
        P = P * V * D * V.inverse();

        tmpt(0, 0) = P(0, 0) / (P(0, 0) + P(0, 1));
        tmpt(1, 1) = P(0, 1) / (P(0, 0) + P(0, 1));
        T_j = tmpt * T_j;

        L = real(T_j.array());
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < k; j++)
            {
                T_matrix[k * i + j] = L(i, j);
            }
        }
        delete[] tmpt_t;
    }
    void test_true_grid_4States(double *T_matrix, double *a, double *b, double *c, int m1, int m2, double *t, int *S, int k, double grid_pw, double t0, double min_del)
    {
        ArrayXXd L(k, k);
        MatrixXcd T_j = MatrixXd::Identity(k, k);
        MatrixXd Q = MatrixXd::Zero(k, k);
        MatrixXcd P(k, k);
        MatrixXcd V(k, k);
        MatrixXcd D(k, k);
        MatrixXcd tmpt = MatrixXd::Zero(k, k);
        MatrixXd E = MatrixXd::Identity(k, k);
        double *tmpt_t = new double[1];

        if ((S[m2 - 1] == 2) || (S[m2 - 1] == 1))
        {
            for (int j = (m1 + 1); j < (m2); j++)
            {
                tmpt_t[0] = t[j - 1];
                E(1, 1) = 1 / (1 + exp(c[0] + c[1] * t[j]));
                E(1, 0) = 1 - E(1, 1);

                P = MatrixXd::Identity(k, k);
                while ((t[j] - tmpt_t[0]) > (grid_pw + min_del))
                {
                    Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
                    Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                    Q(0, 0) = -(Q(0, 1) + Q(0, 3));

                    Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
                    Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                    Q(1, 1) = -(Q(1, 2) + Q(1, 3));

                    EigenSolver<MatrixXd> es(Q);
                    V = es.eigenvectors();
                    D = es.eigenvalues().asDiagonal();
                    for (int i = 0; i < k; i++)
                        D(i, i) = exp(D(i, i) * grid_pw);
                    P = P * V * D * V.inverse();
                    tmpt_t[0] += grid_pw;
                }

                Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
                Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(0, 0) = -(Q(0, 1) + Q(0, 3));

                Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
                Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(1, 1) = -(Q(1, 2) + Q(1, 3));

                EigenSolver<MatrixXd> es(Q);
                V = es.eigenvectors();
                D = es.eigenvalues().asDiagonal();
                for (int i = 0; i < k; i++)
                    D(i, i) = exp(D(i, i) * (t[j] - tmpt_t[0]));
                P = P * V * D * V.inverse();

                for (int i = 0; i < 4; i++)
                {
                    tmpt(i, i) = E(i, S[j] - 1);
                }
                T_j = T_j * P * tmpt;
                // Rprintf("T_j[1,1] = %f", T_j(0,0));
                // Rprintf("T_j[1,2] = %f", T_j(0,1));
            }
            // L = real(L.array());
        }
        else
        {

            for (int j = (m1 + 1); j < (m2 - 1); j++)
            {
                tmpt_t[0] = t[j - 1];
                E(1, 1) = 1 / (1 + exp(c[0] + c[1] * t[j]));
                E(1, 0) = 1 - E(1, 1);
                P = MatrixXd::Identity(k, k);
                while ((t[j] - tmpt_t[0]) > (grid_pw + min_del))
                {
                    Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
                    Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                    Q(0, 0) = -(Q(0, 1) + Q(0, 3));

                    Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
                    Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                    Q(1, 1) = -(Q(1, 2) + Q(1, 3));

                    EigenSolver<MatrixXd> es(Q);
                    V = es.eigenvectors();
                    D = es.eigenvalues().asDiagonal();
                    for (int i = 0; i < k; i++)
                        D(i, i) = exp(D(i, i) * grid_pw);
                    P = P * V * D * V.inverse();
                    tmpt_t[0] += grid_pw;
                }

                Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
                Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(0, 0) = -(Q(0, 1) + Q(0, 3));

                Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
                Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(1, 1) = -(Q(1, 2) + Q(1, 3));

                EigenSolver<MatrixXd> es(Q);
                V = es.eigenvectors();
                D = es.eigenvalues().asDiagonal();
                for (int i = 0; i < k; i++)
                    D(i, i) = exp(D(i, i) * (t[j] - tmpt_t[0]));
                P = P * V * D * V.inverse();
                for (int i = 0; i < 4; i++)
                {
                    tmpt(i, i) = E(i, S[j] - 1);
                }
                T_j = T_j * P * tmpt;
                // Rprintf("T_j[1,1] = %f", T_j(0,0));
                // Rprintf("T_j[1,2] = %f", T_j(0,1));
            }

            tmpt_t[0] = t[m2 - 2];
            P = MatrixXd::Identity(k, k);
            while ((t[m2 - 1] - tmpt_t[0]) > (grid_pw + min_del))
            {
                Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
                Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(0, 0) = -(Q(0, 1) + Q(0, 3));

                Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
                Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
                Q(1, 1) = -(Q(1, 2) + Q(1, 3));

                EigenSolver<MatrixXd> es(Q);
                V = es.eigenvectors();
                D = es.eigenvalues().asDiagonal();
                for (int i = 0; i < k; i++)
                    D(i, i) = exp(D(i, i) * grid_pw);
                P = P * V * D * V.inverse();
                tmpt_t[0] += grid_pw;
            }

            Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
            Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(0, 0) = -(Q(0, 1) + Q(0, 3));

            Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
            Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(1, 1) = -(Q(1, 2) + Q(1, 3));

            EigenSolver<MatrixXd> es(Q);
            V = es.eigenvectors();
            D = es.eigenvalues().asDiagonal();
            for (int i = 0; i < k; i++)
                D(i, i) = exp(D(i, i) * (t[m2 - 1] - tmpt_t[0]));
            P = P * V * D * V.inverse();

            switch (S[m2 - 1])
            {
            case -2:
                tmpt(0, 0) = P(0, 0) + P(0, 1);
                tmpt(1, 1) = P(1, 0) + P(1, 1);
                tmpt(2, 2) = 0;
                tmpt(3, 3) = 0;

                T_j = T_j * tmpt;
                break;
            case 3:
                tmpt(0, 0) = 0;
                // tmpt(1, 1) = exp(a[2] + b[2] * t[m2 - 1]);
                tmpt(1, 1) = Q(1, 2);
                tmpt(2, 2) = 0;
                tmpt(3, 3) = 0;
                T_j = T_j * P * tmpt;
                break;
            default: // case 4
                // tmpt(0, 0) = exp(a[1] + b[1] * t[m2 - 1]);
                // tmpt(1, 1) = exp(a[1] + b[1] * t[m2 - 1]);
                tmpt(0, 0) = Q(0, 3);
                tmpt(1, 1) = Q(1, 3);
                tmpt(2, 2) = 0;
                tmpt(3, 3) = 0;
                T_j = T_j * P * tmpt;
                break;
            }
        }

        tmpt = MatrixXd::Zero(k, k);
        P = MatrixXd::Identity(k, k);
        tmpt_t[0] = t0;
        E(1, 1) = 1 / (1 + exp(c[0] + c[1] * t[m1]));
        E(1, 0) = 1 - E(1, 1);
        while ((t[m1] - tmpt_t[0]) > (grid_pw + min_del))
        {
            Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
            Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(0, 0) = -(Q(0, 1) + Q(0, 3));

            Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
            Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
            Q(1, 1) = -(Q(1, 2) + Q(1, 3));

            EigenSolver<MatrixXd> es(Q);
            V = es.eigenvectors();
            D = es.eigenvalues().asDiagonal();
            for (int i = 0; i < k; i++)
                D(i, i) = exp(D(i, i) * grid_pw);
            P = P * V * D * V.inverse();
            tmpt_t[0] += grid_pw;
        }

        Q(0, 1) = exp(a[0] + b[0] * tmpt_t[0]);
        Q(0, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(0, 0) = -(Q(0, 1) + Q(0, 3));

        Q(1, 2) = exp(a[2] + b[2] * tmpt_t[0]);
        Q(1, 3) = exp(a[1] + b[1] * tmpt_t[0]);
        Q(1, 1) = -(Q(1, 2) + Q(1, 3));
        EigenSolver<MatrixXd> n_es(Q);
        V = n_es.eigenvectors();
        D = n_es.eigenvalues().asDiagonal();
        for (int i = 0; i < k; i++)
            D(i, i) = exp(D(i, i) * (t[m1] - tmpt_t[0]));
        // P = P * V * D * V.inverse();
        P = P * V * D * V.inverse();
        switch (S[m1])
        {
        case -1:
            tmpt(0, 0) = P(0, 0) / (P(0, 0) + P(0, 1));
            tmpt(1, 1) = P(0, 1) / (P(0, 0) + P(0, 1));
            T_j = tmpt * T_j;
            break;
        case 1:
            tmpt(0, 0) = P(0, 0) * E(0, 0) / (P(0, 0) + P(0, 1));
            tmpt(1, 1) = P(0, 1) * E(1, 0) / (P(0, 0) + P(0, 1));
            T_j = tmpt * T_j;
            break;
        case 2:
            tmpt(0, 0) = 0;
            tmpt(1, 1) = P(0, 1) * E(1, 1) / (P(0, 0) + P(0, 1));
            T_j = tmpt * T_j;
            break;
        }

        L = real(T_j.array());
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < k; j++)
            {
                T_matrix[k * i + j] = L(i, j);
            }
        }
        delete[] tmpt_t;
    }
    SEXP lik_grid_4States(SEXP lambda, SEXP beta, SEXP nu, SEXP state, SEXP centered_age, SEXP no_subj, SEXP no_obs, SEXP screen_test, SEXP grid_piecewise, SEXP left_age, SEXP min_del)
    {
        // input
        double *a = REAL(lambda);
        double *b = REAL(beta);
        double *c = REAL(nu);
        int *S = INTEGER(state);
        double *t = REAL(centered_age);
        int *N = INTEGER(no_subj);
        int *m_cum = INTEGER(no_obs);       // cumsum no of obs for each subjects
        int *screen = INTEGER(screen_test); // taken screen yes or no
        double *grid_pw = REAL(grid_piecewise);
        double *t0 = REAL(left_age);
        double *p_min_del = REAL(min_del);
        // output
        SEXP ans;
        PROTECT(ans = allocVector(REALSXP, 1));
        double *rans = REAL(ans);
        rans[0] = 0;

        double *T_matrix = new double[16];

        for (int i = 0; i < N[0]; i++)
        {
            // screen = FALSE
            if (screen[m_cum[i]] == 0)
            {
                test_false_grid_4States(T_matrix, a, b, m_cum[i], m_cum[i + 1], t, S, 4, grid_pw[0], t0[0], p_min_del[0]);
            }
            else
            {
                test_true_grid_4States(T_matrix, a, b, c, m_cum[i], m_cum[i + 1], t, S, 4, grid_pw[0], t0[0], p_min_del[0]);
            }
            rans[0] += -2 * log(T_matrix[0] + T_matrix[1] + T_matrix[2] + T_matrix[3] + T_matrix[4] + T_matrix[5] + T_matrix[6] + T_matrix[7]); // add total of row 1 and 2.
        }

        delete[] T_matrix;
        UNPROTECT(1);
        return (ans);
    }
    void test_fixed_3States(double *T_matrix, double *a, double *b, int m1, int m2, double *t, int *S)
    {
        // ArrayXXd L(k,k);
        MatrixXd Q = MatrixXd::Zero(3, 3);
        // MatrixXcd tmpt = MatrixXd::Zero(k,k);

        switch (S[m2 - 1])
        {
        case 1:
            if ((b[0] != 0) & (b[1] != 0))
            {
                T_matrix[0] = exp(exp(a[0] + b[0] * t[m2 - 2]) / b[0] - exp(a[0] + b[0] * t[m2 - 1]) / b[0] +
                                  exp(a[1] + b[1] * t[m2 - 2]) / b[1] - exp(a[1] + b[1] * t[m2 - 1]) / b[1]);
            }
            else
            {
                if ((b[0] == 0) & (b[1] != 0))
                {
                    T_matrix[0] = exp(exp(a[0]) * (t[m2 - 2] - t[m2 - 1]) +
                                      exp(a[1] + b[1] * t[m2 - 2]) / b[1] - exp(a[1] + b[1] * t[m2 - 1]) / b[1]);
                }
                else
                {
                    if ((b[0] != 0) & (b[1] == 0))
                    {
                        T_matrix[0] = exp(exp(a[1]) * (t[m2 - 2] - t[m2 - 1]) +
                                          exp(a[0] + b[0] * t[m2 - 2]) / b[0] - exp(a[0] + b[0] * t[m2 - 1]) / b[0]);
                    }
                    else
                    {
                        if ((b[0] == 0) & (b[1] == 0))
                        {
                            T_matrix[0] = exp(exp(a[0]) * (t[m2 - 2] - t[m2 - 1]) +
                                              exp(a[1]) * (t[m2 - 2] - t[m2 - 1]));
                        }
                    }
                }
            }
            break;

        default:
            Q(0, 1) = exp(a[0] + b[0] * t[m2 - 1]);
            Q(0, 2) = exp(a[1] + b[1] * t[m2 - 1]);
            if ((b[0] != 0) & (b[1] != 0))
            {
                T_matrix[0] = exp(exp(a[0] + b[0] * t[m2 - 2]) / b[0] - exp(a[0] + b[0] * t[m2 - 1]) / b[0] +
                                  exp(a[1] + b[1] * t[m2 - 2]) / b[1] - exp(a[1] + b[1] * t[m2 - 1]) / b[1]) *
                              Q(0, S[m2 - 1] - 1);
            }
            else
            {
                if ((b[0] == 0) & (b[1] != 0))
                {
                    T_matrix[0] = exp(exp(a[0]) * (t[m2 - 2] - t[m2 - 1]) +
                                      exp(a[1] + b[1] * t[m2 - 2]) / b[1] - exp(a[1] + b[1] * t[m2 - 1]) / b[1]) *
                                  Q(0, S[m2 - 1] - 1);
                }
                else
                {
                    if ((b[0] != 0) & (b[1] == 0))
                    {
                        T_matrix[0] = exp(exp(a[1]) * (t[m2 - 2] - t[m2 - 1]) +
                                          exp(a[0] + b[0] * t[m2 - 2]) / b[0] - exp(a[0] + b[0] * t[m2 - 1]) / b[0]) *
                                      Q(0, S[m2 - 1] - 1);
                    }
                    else
                    {
                        if ((b[0] == 0) & (b[1] == 0))
                        {
                            T_matrix[0] = exp(exp(a[0]) * (t[m2 - 2] - t[m2 - 1]) +
                                              exp(a[1]) * (t[m2 - 2] - t[m2 - 1])) *
                                          Q(0, S[m2 - 1] - 1);
                        }
                    }
                }
            }
            break;
        }
    }
    SEXP lik_fixed_3States(SEXP lambda, SEXP beta, SEXP state, SEXP centered_age, SEXP no_subj, SEXP no_obs)
    {
        // input
        double *a = REAL(lambda);
        double *b = REAL(beta);
        int *S = INTEGER(state);
        double *t = REAL(centered_age);
        int *N = INTEGER(no_subj);
        int *m_cum = INTEGER(no_obs); // cumsum no of obs for each subjects
        // output
        SEXP ans;
        PROTECT(ans = allocVector(REALSXP, 1));
        double *rans = REAL(ans);
        rans[0] = 0;

        double *T_matrix = new double[1];

        for (int i = 0; i < N[0]; i++)
        {
            test_fixed_3States(T_matrix, a, b, m_cum[i], m_cum[i + 1], t, S);
            rans[0] += log(T_matrix[0]); // add total of row 1 and 2.
        }

        delete[] T_matrix;
        rans[0] = -2 * rans[0];
        UNPROTECT(1);
        return (ans);
    }
}
// open command promt= ctrl+shift+C
// compile : R CMD SHLIB Biometrical/C++/Lik_grid_Linux.cc -shared -o Biometrical/C++/Lik_grid_Linux.so
// in R
// .Call("myFunc",body(f),as.double(x),new.env())}