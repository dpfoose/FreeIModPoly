function [baseline, corrected, coefs]=FreeIModPoly(spectrum, abscissa, poly_order=5, max_it=100, threshold=0.95)
    if (poly_order < 1)
        exit("poly_order must be an integer greater than 0");
    endif
    if (threshold >= 1 || threshold <= 0)
        exit("threshold must be a value between 0 and 1");
    endif
    if(rows(spectrum) != rows(abscissa))
        exit("spectrum and abscissa must have same size");
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
        %Calcualte residuals and dev%
        dev = CalcDev(prev_fit, fit);
        %error criteria%
        err = CalcErr(dev, prev_dev);
        %Reconstruction of model input
        fit = fit + dev;
        %if a value in the previous fit is lower than this fit, take previous
        ind = find(prev_fit < fit);
        fit(ind) = prev_fit(ind);
        prev_fit = fit;
        prev_dev = dev;
        i = i + 1;
    until (err < threshold || ((no_max_it == 0) && (i >= max_it)))
    baseline = CalcPoly(coefs, abscissa);
    corrected = spectrum - baseline;

    function dev=CalcDev(specturm, fit)
        R = spectrum - fit;
        R_avg = mean(R);
        centered = R - R_avg;
        centered = centered .^ 2;
        dev = sqrt(sum(centered)/rows(centered));
    endfunction




    function ind=NonPeakInd(spectrum, dev)
        SUM = spectrum + dev;
        ind = find(spectrum <= SUM);
    endfunction

    function poly=CalcPoly(coefs, x)
        poly = coefs(1) + x*coefs(2);
        if (rows(coefs) > 1)
           for i = 3:rows(coefs)
              poly = poly + coefs(i) * x .^ (i-1);
           endfor
        endif
    endfunction

    function fit=OrdinaryLeastSquares(X, y)
        fit = inv(X' * X) * (X' * y);
    endfunction

    function X=DesignMatrix(x, poly_order)
        X = zeros(rows(x), poly_order + 1);
        X(:, 1) = ones(rows(X), 1);
        X(:, 2) = x;
        for i = 2:columns(X)
            X(:, i) = x .^ (i-1);
        endfor
    endfunction

    function err=CalcErr(dev, prev_dev)
        err = abs( (dev - prev_dev) / dev);
    endfunction

endfunction
