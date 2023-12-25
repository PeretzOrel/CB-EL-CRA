function wng = WNG(h, params)
    if params.K > 1
        wng = zeros(params.K, 1);
        for i = 1:params.K
            hi = h(:, i);
            wng(i) = (abs(hi'*params.d(params.f_grid(i), params.theta_d))^2)/abs(hi'*hi);
        end
    else
        wng = (abs(h'*params.d(params.f_grid, params.theta_d))^2)/abs(h'*h);
    end
end
