clc
close all
clear

Nso = 1000;
NsoVec = 1:Nso;            % Total number of satellites
Ncvec = zeros(1, Nso);     % Number of unique Flower Constellations

for ns0 = 1:Nso
    for No = 1:ns0
        if mod(ns0, No) == 0
            Ncvec(ns0) = Ncvec(ns0) + No;
        end
    end
end

% Plotting
figure
scatter(NsoVec, Ncvec, 5, 'r', 'filled')
xlabel('Total Number of Satellites', 'Interpreter', 'latex')
ylabel('Number of Uniform Constellation Configurations', 'Interpreter', 'latex')
title('Uniform Constellation Configurations vs Satellite Count', 'Interpreter', 'latex')
grid minor
set(gca, 'TickLabelInterpreter', 'latex')
