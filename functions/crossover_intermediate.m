function [z] = crossover_intermediate(x1, x2)
    interval = [-0.25, 1.25];
    alpha = min(interval) + abs(diff(interval))*rand(1, length(x1));
    z = x1 + alpha.*(x2 - x1);
end