function df = DF(h, params)
    if params.K > 1
        df = zeros(params.K, 1);
        for i = 1:params.K
            hi = h(:, i);
            lambda_array = params.c/params.f_grid(i);
            Gamma = sinc(2*squareform(pdist(params.r))/lambda_array);
            df(i) = (abs(hi'*params.d(params.f_grid(i), params.theta_d))^2)/abs(hi'*Gamma*hi);
        end
    else
        lambda_array = params.c/params.f_grid;
        Gamma = sinc(2*squareform(pdist(params.r))/lambda_array);
        df = (abs(h'*params.d(params.f_grid, params.theta_d))^2)/abs(h'*Gamma*h);
    end
end
