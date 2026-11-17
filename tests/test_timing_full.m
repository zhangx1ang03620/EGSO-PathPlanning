clc; clear;
[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

Algorithm_Names = {'EGSO', 'GWO', 'PSO', 'GA', 'DE', 'SSA', ...
                   'WOA', 'MFO', 'ALO', 'HHO', 'SCA', ...
                   'ABC', 'ACO', 'BBO', 'CS', 'FA', ...
                   'TLBO', 'BA', 'FPA', 'SMA', 'INFO', 'SBOA', 'ARO'};

total_time = 0;
for map_idx = 1:5
    G = Maps{map_idx};
    Xmax = size(G,1); X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);
    map_time = 0;
    for k = 1:numel(Algorithm_Names)
        fn = str2func(Algorithm_Names{k});
        tic;
        fn(50, 100, X_min, Xmax, dim, fobj, G);
        map_time = map_time + toc;
    end
    fprintf('Map %d (%dx%d): total time for 1 run of all %d algos = %.1f s\n', ...
        map_idx, Xmax, Xmax, numel(Algorithm_Names), map_time);
    total_time = total_time + map_time;
end
fprintf('\nTOTAL (1 run x 23 algos x 5 maps) = %.1f s\n', total_time);
fprintf('Estimated for 20 runs = %.1f s = %.1f min = %.2f hours\n', total_time*20, total_time*20/60, total_time*20/3600);
