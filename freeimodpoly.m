function baseline=FreeIModPoly(spectrum, abscissa, poly_order, max_it, threshold)
    if (poly_order < 1)
        exit("poly_order must be an integer greater than 0")
    endif
    if (threshold >= 1 || threshold <= 0)
        exit("threshold must be a value between 0 and 1")
    endif
    if(rows(spectrum) != rows(abscissa))
        exit("spectrum and abscissa must have same size")
    endif

    i = 2;
    no_max_it = (max_it == 0);
    X = DesignMatrix(abscissa, poly_order);
    coefs = OrdinaryLeastSquares(X, spectrum);
    fit = CalcPoly(coefs, abscissa);
    dev = CalcDev(spectrum, fit);
    prev_dev = dev;

    non_peak_ind = NonPeakInd(spectrum, dev);
    new_abscissa = abscissa(non_peak_ind);

    prev_fit = spectrum(non_peak_ind);

    X = DesignMatrix(new_abscissa, poly_order);

    err = threshold;
    row_ct = rows(new_abscissa);

    do
        %Polynomial fitting%
        coefs = OrdinaryLeastSquares(X, prev_fit);
        fit = CalcPoly(coefs, new_abscissa);
        dev = CalcDev(prev_fit, fit);
        err = CalcErr(dev, prev_dev);
        fit = fit + dev;
        ind = find(prev_fit < fit);
        fit(ind) = prev_fit(ind);
        prev_fit = fit;
        prev_dev = dev;
        i = i + 1;
    until (err < threshold || ((no_max_it == 0) && (i >= max_it)))
    baseline = CalcPoly(coefs, abscissa);
endfunction



function dev=CalcDev(specturm, fit)
    R = spectrum - fit;
    R_avg = mean(R);
    centered = R - R_avg;
    centered = centered .^ 2;
    return sqrt(sum(diff)/rows(diff));
endfunction




function ind=NonPeakInd(spectrum, dev)
    SUM = spectrum + dev;
    ind = find(spectrum <= SUM);
endfunction

function poly=CalcPoly(coefs, x)
    poly = x*coefs(1) + x*coefs(2);
    if (rows(coefs) > 1)
        for i = 3:rows(coefs)
            poly = poly + coefs(i) * x;
        endfor
    endif
endfunction

function fit=OrdinaryLeastSquares(X, y)
    fit = inv(X' * X) * (X' * y);
endfunction

function X=DesignMatrix(x, poly_order)
    X = zeros(rows(x), poly_order + 1);
    X(:, 1) = ones(rows(x));

    for 2:cols(X)
        X(:, i) = x .^ i
    endfor
endfunction

function err=CalcErr(dev, prev_dev)
    return abs( (dev - prev_dev) / dev);
endfunction
