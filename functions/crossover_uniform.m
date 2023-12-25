function [z1, z2] = crossover_uniform(x1, x2)
    crossover_points = rand(1, length(x1)) > 0.5;
    z1 = x1.*crossover_points + x2.*(1 - crossover_points);
    z2 = x1.*(1 - crossover_points) + x2.*crossover_points;
end