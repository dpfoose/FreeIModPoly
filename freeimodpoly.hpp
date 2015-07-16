#include <armadillo>
using namespace std;
using namespace arma;

double IModPoly(const vec &spectrum, const vec &abscissa,
                const vec &baseline, const vec &corrected,
                int poly_order = 5, int max_it = 100, const double threshold = 0);


double CalcDEV(const vec &spectrum, const vec &fit);
uvec PeakRegions(const vec &spectrum, const double dev);
vec CalcPoly(const vec &coefs, const vec &x);
vec OrinaryLeastSquares(const mat &X, const vec &y);
mat DesignMatrix(const vec &x, const int poly_order);

double IModPoly(const vec &spectrum, const vec &abscissa,
                const vec &baseline, const vec &corrected,
                const int poly_order, const int max_it, const double threshold)
{
    double dev;
    double previous_dev;
    vec coefs(poly_order + 1); //coefs (0 is b, 1 thru poly_order correspond to each degree)
    uword i = 1; //we handle the i==0 iteration before the while loop to prevent
    //the evalutation of if(i==0) every iteration.
    double err;
    bool no_max_it = (max_it == 0);
    while (err > threshold && (i != max_it || no_max_it)){

        ++i;
    }
}


double CalcDev(const vec &spectrum, const vec &fit)
{
    vec R = spectrum - fit;
    double R_avg = mean(R);
    vec diff = R - R_avg*ones(R.n_rows);
    diff = pow(diff, 2.0);
    return std::sqrt(sum(diff)/diff.n_rows);
}

uvec PeakRegions(const vec &spectrum, const double dev)
{
    vec SUM = spectrum + dev * ones(spectrum.n_rows);
    return find(spectrum > SUM);
}

vec OrdinaryLeastSquares(const mat &X, const vec &y)
{
    return (vec) inv(X.t() * X) * (X.t() * y);
}

mat DesignMatrix(const vec &x, const int poly_order)
{
    mat X(x.n_rows, poly_order + 1);
    X.col(0) = ones(x.n_rows); //faster than pow(X, 0)
    for (uword i = 1; i < X.n_cols; ++i)
        X.col(i) = pow(X, i);
    return X;
}
