function [child1, child2] = crossoverOnePoint(parent1, parent2)
    % crossoverOnePoint - One-point crossover for integer chromosomes
    n = numel(parent1);
    if n <= 1
        child1 = parent1;
        child2 = parent2;
        return
    end
    cp = randi([1, n-1], 1, 1);
    child1 = [parent1(1:cp), parent2(cp+1:end)];
    child2 = [parent2(1:cp), parent1(cp+1:end)];
end
