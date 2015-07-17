#include <armadillo>
using namespace std;
using namespace arma;

namespace FreeIModPoly{
    void IModPoly(const vec &spectrum, const vec &abscissa,
                    vec &baseline, vec &corrected,
                    int poly_order = 5, int max_it = 100,
                    const double threshold = 0);


    double CalcDEV(const vec &spectrum, const vec &fit);
    uvec NonPeakInd(const vec &spectrum, const double dev);
    vec CalcPoly(const vec &coefs, const vec &x);
    vec OrinaryLeastSquares(const mat &X, const vec &y);
    mat DesignMatrix(const vec &x, const int poly_order);
    double CalcErr(const double &dev, const double &prev_dev);

}

inline void FreeIModPoly::IModPoly(const vec &spectrum, const vec &abscissa,
                                     vec &baseline, vec &corrected,
                                     const int poly_order, const int max_it
                                     const double threshold)
{
    if (poly_order == 0)
        throw invalid_argument("Polynomial order must be 1 (linear) or greater.");
    if (threshold >= 1 || threshold <= 0)
        throw invalid_arument("Threshold value must be between 0 and 1.");
    if (spectrum.n_rows != abscissa.n_rows)
        throw invalid_argument("Spectrum and abscissa must be the same size.")

    using namespace FreeIModPoly;
    uword i = 1;
    bool no_max_it = (max_it == 0);

    //perform first regresion (on spectrum without removed major peaks)
    mat X = DesignMatrix(abscissa, poly_order);

    vec coefs = OrdinaryLeastSquares(X, spectrum);
    vec fit = CalcPoly(coefs, abscissa);
    double dev = CalcDev(specturm, fit);
    double prev_dev = dev; //used in while loop critera

    //find major peak areas to remove and remove them
    uvec non_peak_ind = NonPeakInd(spectrum, dev);
    vec new_abscissa = abscissa(non_peak_ind);

    vec prev_fit = spectrum(non_peak_ind);; //not a fit here, but in the loop.

    X = DesignMatrix(new_abscissa, poly_order);

    double err;
    uword rows = new_abscissa.n_rows;
    do{ //always perform at least one regression on the "shrunken" spectrum
        //Polynomial Fitting
        coefs = OrdinaryLeastSquares(X, prev_fit);
        fit = CalcPoly(coefs, new_abscissa);
        //Residual and dev calc (residual calculted inside CalcDev)
        dev = CalcDev(prev_fit, fit);
        //error criteria
        err = CalcErr(dev, prev_dev)
        //Reconstruction of model input
        fit += dev * ones(rows);
        //if a value in the previous fit is lower than this fit, take the previous
        uvec ind = find (prev_fit < fit);
        fit(ind) = prev_fit(ind);
        prev_fit = fit;
        prev_dev = dev;
        ++i;
    }while (err > threshold && (no_max_it || i != max_it));
    //calculate fit for all values in original abscissa
    baseline = CalcPoly(coefs, abscissa);
    corrected = spectrum - baseline;
}


inline double FreeIModPoly::CalcDev(const vec &spectrum, const vec &fit)
{
    vec R = spectrum - fit;
    double R_avg = mean(R);
    vec diff = R - R_avg*ones(R.n_rows);
    diff = pow(diff, 2.0);
    return std::sqrt(sum(diff)/diff.n_rows);
}

inline uvec FreeIModPoly::NonPeakInd(const vec &spectrum, const double dev)
{
    vec SUM = spectrum + dev * ones(spectrum.n_rows);
    return find(spectrum <= SUM);
}

inline vec CalcPoly(const vec &coefs, const vec &x)
{
    vec y = x*coefs(0) + x*coefs(1); //0th and 1st power of x
    //this loop only used for powers where pow(x, power) needs to be calculated
    if (coefs.n_rows > 1){
        for (uword i = 2; i < coefs.n_rows; ++i)
            y += coefs(i) * x;
    }
    return y;

}

inline vec FreeIModPoly::OrdinaryLeastSquares(const mat &X, const vec &y)
{
    return (vec) inv(X.t() * X) * (X.t() * y);
}

inline mat FreeIModPoly::DesignMatrix(const vec &x, const int poly_order)
{
    mat X(x.n_rows, poly_order + 1);
    X.col(0) = ones(x.n_rows); //faster than pow(X, 0)
    for (uword i = 1; i < X.n_cols; ++i)
        X.col(i) = pow(x, i);
    return X;
}

inline double FreeIModPoly::CalcErr(const double &dev, const double &prev_dev)
{
    return std::abs( (dev - prev_dev) / dev );
}
