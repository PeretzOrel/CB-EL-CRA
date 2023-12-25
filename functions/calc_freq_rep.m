function H = calc_freq_rep(coeff, params)
    % calc frequency response of FIR coeff
    H = [];
    for i = 1:size(coeff, 1)
        H = [H; freqz(coeff(i, :), 1, params.f_grid, params.Fs).'];
    end
end
