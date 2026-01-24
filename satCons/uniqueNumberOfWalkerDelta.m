clc; clear; close all

N_max = 1000;
unique_constellation_count = zeros(N_max, 1);

for N = 1:N_max
    unique_keys = containers.Map();

    % Find all (P, S) pairs such that N = P * S
    for P = 1:N
        if mod(N, P) ~= 0
            continue
        end
        S = N / P;

        for f = 0:S-1
            % Generate Walker Delta (RAAN, Mean Anomaly) pairings
            raan = zeros(N, 1);
            M = zeros(N, 1);
            idx = 1;
            for i = 0:P-1
                for j = 0:S-1
                    raan(idx) = mod(2 * pi * i / P, 2*pi);
                    M(idx) = mod(2 * pi * j / S + 2 * pi * f * i / N, 2*pi);
                    idx = idx + 1;
                end
            end

            % Sort and flatten configuration for comparison
            key_matrix = sortrows([mod(rad2deg(raan), 360), mod(rad2deg(M), 360)]);
            key = mat2str(round(key_matrix, 4));  % rounding to 4 decimals for uniqueness

            if ~isKey(unique_keys, key)
                unique_keys(key) = true;
            end
        end
    end

    % Store number of unique configurations
    unique_constellation_count(N) = unique_keys.Count;
end

% Plotting
figure
scatter(1:N_max, unique_constellation_count, 10, 'b', 'filled')
xlabel('Total Number of Satellites (N)', 'Interpreter', 'latex')
ylabel('Number of Unique Uniform Constellations', 'Interpreter', 'latex')
title('Uniform Satellite Constellations (RAAN - Mean Anomaly Phasing)', 'Interpreter', 'latex')
grid on
set(gca, 'TickLabelInterpreter', 'latex')
