function [bestChrom, bestFitness, history] = runIntegerGA(objectiveFcn, initialPopulation, maxGenerations, options)
    % runIntegerGA - Simple integer-only genetic algorithm
    % Inputs:
    %   objectiveFcn      - function handle f(x) -> scalar fitness (minimize)
    %   initialPopulation - [Npop x Nvars] integer matrix (each row is a chromosome)
    %   maxGenerations    - scalar number of generations
    %   options           - struct with fields (optional):
    %       mutationRate   (default 0.1)
    %       crossoverRate  (default 0.8)
    %       eliteCount     (default 2)
    %       tournamentSize (default 3)
    %       mutationStep   (default 1)
    %       lowerBounds    (default [])
    %       upperBounds    (default [])
    %
    % Outputs:
    %   bestChrom  - best chromosome found
    %   bestFitness- fitness of bestChrom
    %   history    - struct with fields: fitnessTrace, bestPerGen
    
    if nargin < 4; options = struct(); end
    if ~isfield(options,'mutationRate');   options.mutationRate = 0.1; end
    if ~isfield(options,'crossoverRate');  options.crossoverRate = 0.8; end
    if ~isfield(options,'eliteCount');     options.eliteCount = 2; end
    if ~isfield(options,'tournamentSize'); options.tournamentSize = 3; end
    if ~isfield(options,'mutationStep');   options.mutationStep = 1; end
    if ~isfield(options,'lowerBounds');    options.lowerBounds = []; end
    if ~isfield(options,'upperBounds');    options.upperBounds = []; end
    
    pop = initialPopulation;
    [Npop, ~] = size(pop);
    
    [fitness, order] = evaluatePopulation(pop, objectiveFcn);
    pop = pop(order,:);
    fitness = fitness(order);
    
    bestChrom   = pop(1,:);
    bestFitness = fitness(1);
    
    history.fitnessTrace = zeros(maxGenerations+1,1);
    history.bestPerGen   = zeros(maxGenerations+1,1);
    history.fitnessTrace(1) = bestFitness;
    history.bestPerGen(1)   = bestFitness;
    
    for gen = 1:maxGenerations
        % Elitism
        elites = pop(1:options.eliteCount, :);
    
        % Build next population
        newPop = elites;
        while size(newPop,1) < Npop
            % Select parents
            p1 = selectTournament(pop, fitness, options.tournamentSize);
            p2 = selectTournament(pop, fitness, options.tournamentSize);
    
            child1 = pop(p1,:);
            child2 = pop(p2,:);
    
            % Crossover
            if rand < options.crossoverRate
                [child1, child2] = crossoverOnePoint(child1, child2);
            end
    
            % Mutation
            child1 = mutateRandomInt(child1, options.mutationRate, options.mutationStep, options.lowerBounds, options.upperBounds);
            child2 = mutateRandomInt(child2, options.mutationRate, options.mutationStep, options.lowerBounds, options.upperBounds);
    
            newPop = [newPop; child1; child2]; %#ok<AGROW>
        end
    
        % Trim to population size
        if size(newPop,1) > Npop
            newPop = newPop(1:Npop,:);
        end
    
        % Evaluate
        [fitness, order] = evaluatePopulation(newPop, objectiveFcn);
        pop = newPop(order,:);
        fitness = fitness(order);
    
        % Track best
        if fitness(1) < bestFitness
            bestFitness = fitness(1);
            bestChrom   = pop(1,:);
        end
    
        history.fitnessTrace(gen+1) = fitness(1);
        history.bestPerGen(gen+1)   = bestFitness;
    end

end
