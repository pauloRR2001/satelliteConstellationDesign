function chrom = mutateRandomInt(chrom, mutationRate, mutationStep, lowerBounds, upperBounds)
% mutateRandomInt - Random integer mutation per gene
% Adds +/- mutationStep with probability mutationRate per gene
n = numel(chrom);
for i = 1:n
    if rand < mutationRate
        if rand < 0.5
            chrom(i) = chrom(i) + mutationStep;
        else
            chrom(i) = chrom(i) - mutationStep;
        end
    end
end

% Optional bounds clamping
if ~isempty(lowerBounds)
    chrom = max(chrom, lowerBounds(:)');
end
if ~isempty(upperBounds)
    chrom = min(chrom, upperBounds(:)');
end
end
