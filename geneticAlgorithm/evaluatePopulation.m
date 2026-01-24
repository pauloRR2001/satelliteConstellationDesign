function [fitness, order] = evaluatePopulation(population, objectiveFcn)
% evaluatePopulation - Evaluates fitness for all individuals and sorts
N = size(population,1);
fitness = zeros(N,1);
for i = 1:N
    fitness(i) = objectiveFcn(population(i,:));
end
[~, order] = sort(fitness, 'ascend');
end
