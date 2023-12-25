function z = mutation_nonuniform(x, ga_params)
    z = x;
    n_genes = size(x, 2);
    
    % lower and upper bounds of genes
    LB = 0;
    UB = 1;
    
    % genens to be mutated
    p = rand(1, n_genes) <= ga_params.p_mu_genes;
    idx = find(p);
    idx = idx(randperm(length(idx)));

    for i=idx
        g_K = rand*(1-ga_params.iter/ga_params.iter_max)^ga_params.gamma;
        rb = randi(2)-1;
        if rb == 0
            z(i) = x(i) + (UB - x(i))*g_K;
        else % rb == 1
            z(i) = x(i) - (x(i) - LB)*g_K;
        end
    end
end