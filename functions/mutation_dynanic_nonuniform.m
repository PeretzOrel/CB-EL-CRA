function z = mutation_dynanic_nonuniform(x, ga_params)
    z = x;
    n_genes = size(x, 2);

    % genens to be mutated
    p = rand(1, n_genes) <= ga_params.p_mu_genes;
    idx = find(p);
    idx = idx(randperm(length(idx)));

    for i=idx
        % lower and upper bounds of genes
        if i == 1
            LB = 0;
            UB = x(i+1);
        elseif i == n_genes
            LB = x(i-1);
            UB = 1;
        else
            LB = x(i-1);
            UB = x(i+1);
        end
        
        g_K = rand*(1-ga_params.iter/ga_params.iter_max)^ga_params.gamma;
        
        rb = mod(randi(2),2);
        if rb == 0
            z(i) = x(i) + (UB - x(i))*g_K;
        else % rb == 1
            z(i) = x(i) - (x(i) - LB)*g_K;
        end
    end
end
