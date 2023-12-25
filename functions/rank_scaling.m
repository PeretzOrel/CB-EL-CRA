function fit_o = rank_scaling(fit)
    n_pop = length(fit);
    [~, fitness_sorted_idx] = sort(fit);
    rank = fitness_sorted_idx;
    min_fitness = min(fit);
    max_fitness = max(fit);
    s = min_fitness + (max_fitness - min_fitness)*(rank - 1)/(n_pop - 1);
    p = s/sum(s);
    fit_o = p*n_pop;
end