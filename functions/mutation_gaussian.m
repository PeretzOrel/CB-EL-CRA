function z = mutation_gaussian(x, ga_params)
    z = x;
    n_genes = size(x, 2);
    % genens to be mutated
    p = rand(1, n_genes) <= ga_params.p_mu_genes;
    idx = find(p);
    idx = idx(randperm(length(idx)));
    
    sigma = ga_params.sigma0*(1 - ga_params.iter/ga_params.iter_max)^ga_params.gamma;
    z(idx) = x(idx) + normrnd(0, sigma, [1, length(idx)]);
end