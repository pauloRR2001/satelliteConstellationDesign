function idx = selectTournament(population, fitness, tournamentSize)
    % selectTournament - Picks one parent index via tournament selection
    % Minimization fitness assumed
    N = size(population,1);
    competitors = randi(N, tournamentSize, 1);
    [~, localOrder] = sort(fitness(competitors), 'ascend');
    idx = competitors(localOrder(1));
end
