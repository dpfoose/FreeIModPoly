// Copyright (C) 2015 Wright State University
// Author: Daniel P. Foose
// This file is part of FreeIModPoly.
//
// FreeIModPoly is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// FreeIModPoly is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file LICENSE.  If not, see
// <http://www.gnu.org/licenses/>.
//
// FreeIModPoly: A free software implementation of the Vancouver Raman Algorithm
// Please cite DOI: 10.1366/000370207782597003 and this project (see CITATION)
// The author of this implementation is not associated with the authors of the
// algorithm.

#include "freeimodpoly.h"
using namespace std;
using namespace arma;



///
/// \brief FreeIModPoly::IModPoly Perform the Vancouver Raman Algorithm to correct baseline
/// \param spectrum A vector containing the signal to be corrected
/// \param abscissa The x-values of spectrum
/// \param baseline Will contain the fitted baseline
/// \param corrected Will contain the corrected spectrum
/// \param poly_order Polynomial order for the baseline fits
/// \param max_it Maximum number of iterations to perform. if set to 0 there is no maximum
/// \param threshold A value specifying the upper limit of the error criterion
/// \param err The final error criterion
/// \return The number of iterations
///
uword FreeIModPoly::IModPoly(const vec &spectrum,
                             const vec &abscissa,
                             vec &baseline,
                             vec &corrected,
                             double &err,
                             const uword poly_order,
                             const uword max_it,
                             const double threshold)
{   
    if (poly_order == 0)
        throw invalid_argument("Polynomial order must be 1 (linear) or greater.");
    if (threshold >= 1 || threshold <= 0)
        throw invalid_argument("Threshold value must be between 0 and 1.");
    if (spectrum.n_rows != abscissa.n_rows)
        throw invalid_argument("Spectrum and abscissa must be the same size.");

    uword i = 1;
    bool no_max_it = (max_it == 0);
    mat X;
    vec coefs, fit;
    double dev;

    //perform first regresion (on spectrum without removed major peaks)
    X = Vandermonde(abscissa, poly_order);
    coefs = FreeIModPoly::OrdinaryLeastSquares(X, spectrum);
    fit = FreeIModPoly::CalcPoly(coefs, abscissa);
    dev = FreeIModPoly::CalcDev(spectrum, fit);
    cout << "CalcDev" << endl;

    double prev_dev = dev; //used in while loop critera

    //find major peak areas to remove and remove them
    uvec non_peak_ind = NonPeakInd(spectrum, dev);
    vec new_abscissa = abscissa(non_peak_ind);

    vec prev_fit = spectrum(non_peak_ind);; //not a fit here, but in the loop.

    X = Vandermonde(new_abscissa, poly_order);

    uword rows = new_abscissa.n_rows;
    do{ //always perform at least one regression on the "shrunken" spectrum
        //Polynomial Fitting
        coefs = FreeIModPoly::OrdinaryLeastSquares(X, prev_fit);
        fit = FreeIModPoly::CalcPoly(coefs, new_abscissa);
        //Residual and dev calc (residual calculted inside CalcDev)
        dev = FreeIModPoly::CalcDev(prev_fit, fit);
        //error criteria
        err = FreeIModPoly::CalcErr(dev, prev_dev);
        //Reconstruction of model input
        fit += dev * ones(rows);
        //if a value in the previous fit is lower than this fit, take the previous
        uvec ind = find (prev_fit < fit);
        fit(ind) = prev_fit(ind);
        prev_fit = fit;
        prev_dev = dev;
        ++i;
    }while (err > threshold && (no_max_it || i <= max_it));
    //calculate fit for all values in original abscissa
    baseline = FreeIModPoly::CalcPoly(coefs, abscissa);
    corrected = spectrum - baseline;
    return i;
}

///
/// \brief FreeIModPoly::CalcDev
/// \param spectrum
/// \param fit
/// \return
///
double FreeIModPoly::CalcDev(const vec &spectrum, const vec &fit)
{
    using namespace arma;
    vec R = spectrum - fit;
    double R_avg = mean(R);
    vec centered = R - R_avg*ones(R.n_rows);
    centered = pow(centered, 2.0);
    return std::sqrt(sum(centered)/centered.n_rows);
}

///
/// \brief FreeIModPoly::NonPeakInd
/// \param spectrum
/// \param dev
/// \return
///
uvec FreeIModPoly::NonPeakInd(const vec &spectrum, const double dev)
{
    using namespace arma;
    vec SUM = spectrum + dev * ones(spectrum.n_rows);
    return find(spectrum <= SUM);
}

///
/// \brief FreeIModPoly::CalcPoly Calculate the values of a polynomial
/// \param coefs The polynomial coefficients ordered from 0th order to nth order
/// \param x The values for which the polynomial is to be calculated
/// \return The calculated values
///
vec FreeIModPoly::CalcPoly(const vec &coefs, const vec &x)
{
    vec y = coefs(0) + x*coefs(1); //0th and 1st power of x
    //this loop only used for powers where pow(x, power) needs to be calculated
    if (coefs.n_rows > 1){
        for (uword i = 2; i < coefs.n_rows; ++i)
            y += coefs(i) * pow(x, i);
    }
    return y;

}

///
/// \brief FreeIModPoly::OrdinaryLeastSquares Perform Ordinary Least Squares
/// \param X The design matrix of the regression
/// \param y The response vector
/// \return
///
vec FreeIModPoly::OrdinaryLeastSquares(const mat &X, const vec &y)
{
    mat Q, R;
    qr(Q, R, X);
    return solve(R, Q.t()) * y;
}

///
/// \brief FreeIModPoly::Vandermonde Build a Vandermonde matrix for OLS
/// \param x A vector containing a signal
/// \param poly_order The polynomial order
/// \return A Vandermonde matrix of x for the polynomial of order poly_order
///
mat FreeIModPoly::Vandermonde(const vec &x, const int poly_order)
{
    mat X(x.n_rows, poly_order + 1);
    X.col(0) = ones(x.n_rows); //faster than pow(X, 0)
    X.col(1) = x;
    for (uword i = 2; i < X.n_cols; ++i)
        X.col(i) = pow(x, i);
    return X;
}

///
/// \brief FreeIModPoly::CalcErr Calculate the error criterion
/// \param dev
/// \param prev_dev
/// \return
///
double FreeIModPoly::CalcErr(const double &dev, const double &prev_dev)
{
    return std::abs( (dev - prev_dev) / dev );
}
