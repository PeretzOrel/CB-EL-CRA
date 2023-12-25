function [z1, z2] = crossover_single_point(x1, x2)
    crossover_point = randi([1, length(x1)], 1);
    z1 = x1;
    z1(crossover_point+1:end) = x2(crossover_point+1:end);
    
    z2 = x2;
    z2(crossover_point+1:end) = x1(crossover_point+1:end);
end