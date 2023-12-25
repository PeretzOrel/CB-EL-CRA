function NewChrIx = stochastic_universal_selection(fitness, n_sel)
    n_pop = length(fitness);
    cumfit = cumsum(fitness);
    trials = cumfit(n_pop) / n_sel * (rand + (0:n_sel-1)');
    Mf = cumfit(:, ones(1, n_sel));
    Mt = trials(:, ones(1, n_pop))';
    [NewChrIx, ~] = find(Mt < Mf & [zeros(1, n_sel); Mf(1:n_pop-1, :)] <= Mt);
    [~, shuf] = sort(rand(n_sel, 1));
    NewChrIx = NewChrIx(shuf);
end