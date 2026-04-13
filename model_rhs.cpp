#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

inline bool mda_active_cpp(double t, const arma::vec &starts, double duration)
{
    for (arma::uword i = 0; i < starts.n_elem; ++i)
    {
        if (t >= starts[i] && t < starts[i] + duration)
            return true;
    }
    return false;
}

inline arma::vec ageing_term(const arma::vec &x, double age_rate)
{
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::zeros);
    if (n == 0)
        return out;

    out[0] = -age_rate * x[0];

    for (arma::uword i = 1; i < n - 1; ++i)
    {
        out[i] = age_rate * x[i - 1] - age_rate * x[i];
    }

    if (n > 1)
    {
        out[n - 1] = age_rate * x[n - 2];
    }

    return out;
}

// [[Rcpp::export]]
arma::vec rhs_mda_cpp(
    double t,
    const arma::vec &y,
    double betaS,
    double c,
    double uS,
    double uR,
    double uC,
    double k,
    double a,
    double aC,
    double tau,
    double theta,
    double amrd_rate,
    const arma::mat &m_contact,
    const arma::vec &mort,
    const arma::vec &azt,
    const arma::vec &mda_start_times,
    double mda_duration,
    bool use_mda)
{
    arma::uword n_age = mort.n_elem;
    double betaR = betaS * (1.0 - c);
    double age_rate = 1.0 / 365.25;

    // State layout: X, S, R, Sr, Rs, D, CumIncR, AMRD
    arma::vec X = y.subvec(0, n_age - 1);
    arma::vec S = y.subvec(n_age, 2 * n_age - 1);
    arma::vec R = y.subvec(2 * n_age, 3 * n_age - 1);
    arma::vec Sr = y.subvec(3 * n_age, 4 * n_age - 1);
    arma::vec Rs = y.subvec(4 * n_age, 5 * n_age - 1);

    arma::vec N = X + S + R + Sr + Rs;
    arma::vec denom = clamp(N, 1e-12, datum::inf);

    arma::vec S_tot = S + Sr;
    arma::vec R_tot = R + Rs;

    arma::vec lambdaS = betaS * (m_contact * (S_tot / denom));
    arma::vec lambdaR = betaR * (m_contact * (R_tot / denom));

    bool is_mda = use_mda && mda_active_cpp(t, mda_start_times, mda_duration);

    arma::vec a_t(n_age, arma::fill::value(tau));
    arma::vec aC_t(n_age, arma::fill::value(tau));

    if (is_mda)
    {
        a_t += a * azt;
        aC_t += aC * azt;
    }

    arma::vec mort_eff = mort;
    if (is_mda && n_age >= 5)
    {
        mort_eff.subvec(0, 4) *= (1.0 - theta);
    }

    // Stationary demography: births balance all live-population deaths
    arma::vec death_flow = mort_eff % N + amrd_rate * (R + Rs);
    arma::vec births(n_age, arma::fill::zeros);
    births[0] = sum(death_flow);

    arma::vec dX =
        births +
        (uS + a_t) % S +
        uR * R +
        uC * (Sr + Rs) -
        (lambdaS + lambdaR) % X +
        ageing_term(X, age_rate) -
        mort_eff % X;

    arma::vec dS =
        lambdaS % X -
        uS * S -
        k * lambdaR % S -
        a_t % S +
        ageing_term(S, age_rate) -
        mort_eff % S;

    arma::vec dR =
        lambdaR % X -
        uR * R -
        k * lambdaS % R +
        aC_t % (Sr + Rs) +
        ageing_term(R, age_rate) -
        mort_eff % R -
        amrd_rate * R;

    arma::vec dSr =
        k * lambdaR % S -
        uC * Sr -
        aC_t % Sr +
        ageing_term(Sr, age_rate) -
        mort_eff % Sr;

    arma::vec dRs =
        k * lambdaS % R -
        uC * Rs -
        aC_t % Rs +
        ageing_term(Rs, age_rate) -
        mort_eff % Rs -
        amrd_rate * Rs;

    arma::vec dD = mort_eff % (X + S + R + Sr + Rs) + amrd_rate * (R + Rs);
    arma::vec dCumIncR = lambdaR % X + k * lambdaR % S;
    arma::vec dAMRD = amrd_rate * (R + Rs);

    arma::vec dy(8 * n_age, arma::fill::zeros);
    dy.subvec(0, n_age - 1) = dX;
    dy.subvec(n_age, 2 * n_age - 1) = dS;
    dy.subvec(2 * n_age, 3 * n_age - 1) = dR;
    dy.subvec(3 * n_age, 4 * n_age - 1) = dSr;
    dy.subvec(4 * n_age, 5 * n_age - 1) = dRs;
    dy.subvec(5 * n_age, 6 * n_age - 1) = dD;
    dy.subvec(6 * n_age, 7 * n_age - 1) = dCumIncR;
    dy.subvec(7 * n_age, 8 * n_age - 1) = dAMRD;

    return dy;
}