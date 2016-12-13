## Copyright (C) 2015 Wright State University
## Author: Daniel P. Foose
## This file is part of FreeIModPoly.

## FreeIModPoly is distributed under two licenses, the GNU General Public License
## v3 and the MIT License. Which license you use is left to your discretion

## GPL Statement:
## FreeIModPoly is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## FreeIModPoly is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## MIT License Statement:
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

## FreeIModPoly: A free software implementation of the Vancouver Raman Algorithm
## Please cite DOI: 10.1366/000370207782597003 and this project (see CITATION)
## The author of this implementation is not associated with the authors of the
## algorithm.

# abscissa and spectrum should be lists

"FreeIModPoly" <-
function (spectrum, abscissa, polyOrder=5, maxIt=100, threshold=0.05) {

  if (polyOrder < 1){
    warning("polyOrder must be 1 (linear) or greater")
    quit()
  }
  if (threshold >= 1 || threshold <= 0){
    warning("threshold must be between 0 and 1")
    quit()
  }
  if (length(spectrum) != length(abscissa)){
    warning("spectrum and abscissa must be of same size")
    quit()
  }
  i <- 2
  noMaxIt <- (maxIt == 0)
  model <- lm(unlist(spectrum) ~ poly(unlist(abscissa), polyOrder, raw=TRUE))
  coefs <- model$coefficients
  fit <- model$fitted.values
  dev <- CalcDev(model$residuals)
  prevDev <- dev;
  nonPeakInd <- NonPeakInd(spectrum, dev)
  newAbscissa <- abscissa[nonPeakInd]
  prevFit <- spectrum[nonPeakInd]
  err <- threshold

  repeat{
    #polynomial fitting
    model <- lm (unlist(prevFit) ~ poly(unlist(newAbscissa), polyOrder, raw=TRUE))
    fit <- model$fitted.values
    dev <- CalcDev(model$residuals)
    err <- CalcErr(dev, prevDev)
    #reconstruction of model input
    fit <- fit + dev
    ind <- which(prevFit < fit)
    fit[ind] <- prevFit[ind]
    prevFit <- fit
    prevDev <- dev
    i <- i + 1
    if((err < threshold) || ((noMaxIt != 0) && (i >= maxIt))){
      break
    }
  }
  coefs <- model$coefficients
  baseline <- CalcPoly(coefs, abscissa)
  corrected <- spectrum - baseline

  return(list("baseline" = baseline,
              "corrected" = corrected,
              "coefs" = coefs,
              "iterations" = i,
              "error" = err))
}


"CalcDev" <-
function(residual){
  averageResidual = mean(residual)
  centered = residual - averageResidual
  centered = centered ^ 2
  return(sqrt(sum(centered) / length(centered)))
}

"NonPeakInd" <-
function(spectrum, dev){
  SUM = spectrum + dev
  return(which(spectrum <= SUM))
}

"CalcPoly" <-
function(coefs, x){
  val = coefs[1] + x*coefs[2]
  if (length(coefs) > 1){
    for (i in 3:length(coefs)){
      val = val + coefs[i] * x^(i-1)
    }
  }
  return(val)
}

"CalcErr" <-
function(dev, prevDev){
  return (abs( (dev - prevDev) / dev))
}
