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

#include <iostream>
#include <string>
#include "freeimodpoly.h"
#include <algorithm>

using namespace std;
///
/// \brief main
/// \param argc
/// \param argv
/// \return
/// How to call: freeimodpoly order max_it threshold infile outfile outtype transpose
/// Valid outtypes:
/// "txt" - whitespace-delimited text
/// "csv" - comma-delimited text\
///
/// Default outfile is infile with _corrected appended to name (including old extension)
/// The file is structured in columns: abscissa baseline corrected for single spectra
/// Baseline and corrected are stored in separate files "_baseline, _corrected"
/// Default format is csv
///
/// Program crashes on invalid input, but I'm not going to do anything about it
/// If you need to transpose, you've got to specify all parameters, sorry.
int main(int argc, char* argv[])
{
    string infilename, outfilename, blfilename;
    unsigned int poly_order, max_it;
    double threshold;
    file_type outtype;
    string outtype_text;
    bool transpose = false;
    //parse commandline parameters:
    switch (argc){
      case 1: case 2: case 3: case 4:
        cerr << "No spectrum file specified" << endl;
        return 1;
      case 5:
        infilename = string(argv[4]);
        outfilename = infilename + "_corrected.csv";
        blfilename = infilename + "_baseline.csv";
        poly_order = atof(argv[1]);
        max_it = atof(argv[2]);
        threshold = atof(argv[3]);
        outtype = csv_ascii;
        outtype_text = "csv";
        break;
      case 6:
        infilename = string(argv[4]);
        outfilename = string(argv[5]);
        blfilename = outfilename + "_baseline.csv";
        poly_order = atof(argv[1]);
        max_it = atof(argv[2]);
        threshold = atof(argv[3]);
        outtype = csv_ascii;
        outtype_text = "csv";
        break;
      case 7:
        infilename = string(argv[4]);
        outfilename = string(argv[5]);
        blfilename = outfilename + "_baseline.csv";
        poly_order = atof(argv[1]);
        max_it = atof(argv[2]);
        threshold = atof(argv[3]);
        outtype_text = string(argv[6]);
        transform(outtype_text.begin(), outtype_text.end(), outtype_text.begin(), ::tolower);
        outtype = (outtype_text == "txt" ? raw_ascii : csv_ascii);
        break;
      case 8: default:
        infilename = string(argv[4]);
        outfilename = string(argv[5]);
        blfilename = outfilename + "_baseline.csv";
        poly_order = atof(argv[1]);
        max_it = atof(argv[2]);
        threshold = atof(argv[3]);
        outtype_text = string(argv[6]);
        transform(outtype_text.begin(), outtype_text.end(), outtype_text.begin(), ::tolower);
        outtype = (outtype_text == "txt" ? raw_ascii : csv_ascii);
        transpose = true;
        break;
    }
    using namespace arma;
    mat spec_file;
    spec_file.load(infilename);
    vec abscissa;
    mat spectra;
    try{
        abscissa = spec_file.col(0);
        spectra = spec_file.cols(1, spec_file.n_cols - 1);
    }catch(std::exception e){
        cerr << "Spectrum could not be loaded" << endl;
        return 1;
    }
    if(transpose){
        inplace_trans(spectra);
        inplace_trans(abscissa);
    }

    if (spectra.n_cols == 1){
        vec baseline, corrected;
        mat out_mat(spectra.n_rows, 3);
        double err;
        int iter;
        try{
            iter = FreeIModPoly::IModPoly(spectra.col(0),
                                          abscissa,
                                          baseline,
                                          corrected,
                                          err, poly_order, max_it, threshold);
        }catch(exception e){cerr << "exception thrown!" << endl; return 2;}
        out_mat.col(0) = abscissa;
        out_mat.col(1) = baseline;
        out_mat.col(2) = corrected;
        cout << "Final error: " << err << endl;
        if (iter >= max_it)
            cout << "Did not converge in " << iter << " iterations." << endl;
        else{
            cout << "Converged in " << iter << " iterations." << endl;
        }
        bool ok = out_mat.save(outfilename + "." + outtype_text, outtype);
        if (ok){return 0;}
        else{cout << "Save failed!" << endl; return 3;}
    }

    else{
        mat error_report (spectra.n_cols, 2);
        mat corrected_spectra(spectra.n_rows, spectra.n_cols + 1);
        mat baselines(spectra.n_rows, spectra.n_cols + 1);
        baselines.col(0) = abscissa;
        corrected_spectra.col(0) = abscissa;
        for (uword i = 0; i < spectra.n_cols; ++i){
            double err;
            vec spectrum = spectra.col(i);
            vec baseline, corrected;
            int iter = FreeIModPoly::IModPoly(spectrum, abscissa,
                                              baseline, corrected, err,
                                              poly_order, max_it, threshold);
            baselines.col(i+1) = baseline;
            corrected_spectra.col(i+1) = corrected;
            error_report(i, 0) = err;
            error_report(i, 1) = (double) (iter < max_it);
        }
        bool blok = baselines.save(outfilename + "_bl." + outtype_text, outtype);
        bool csok = corrected_spectra.save(outfilename + "_cs." + outtype_text, outtype);
        bool erok = error_report.save(outfilename + "_er." + outtype_text, outtype);
        cout << (blok ? "Saved baseline\n" : "Saving baslines failed\n")
             << (csok ? "Saved corrected\n" : "Saving spectra failed\n")
             << (erok ? "Saved error report\n" : "Saving error report failed\n");

        if(!(blok && csok && erok)){return 3;}
        else{return 0;}

    }
}

