// Copyright (C) 2015 Wright State University
// Author: Daniel P. Foose
// This file is part of FreeIModPoly.
// FreeIModPoly is distributed under two licenses, the GNU GPL v3 and the MIT License
// Which license applies is up to your discretion.

// GPL Statement:
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

// MIT Statement:
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// FreeIModPoly: A free software implementation of the Vancouver Raman Algorithm
// Please cite DOI: 10.1366/000370207782597003 and this project (see CITATION)
// The author of this implementation is not associated with the authors of the
// algorithm.
#ifndef FREEIMODPOLY_H
#define FREEIMODPOLY_H

#include <armadillo>
using namespace std;
using namespace arma;

///
///
///
namespace FreeIModPoly{
    uword IModPoly(const vec &spectrum, const vec &abscissa,
                   vec &baseline, vec &corrected,
                   double &err,
                   const uword poly_order,
                   const uword max_it,
                   const double threshold);

    double CalcDev(const vec &spectrum, const vec &fit);
    uvec NonPeakInd(const vec &spectrum, const double dev);
    vec CalcPoly(const vec &coefs, const vec &x);
    vec OrdinaryLeastSquares(const mat &X, const vec &y);
    mat Vandermonde(const vec &x, const int poly_order);
    double CalcErr(const double &dev, const double &prev_dev);
}
#endif
