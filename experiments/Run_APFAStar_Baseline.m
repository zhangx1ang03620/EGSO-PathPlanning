% Run_APFAStar_Baseline.m
% Guided-field hybrid baseline (APF-weighted A*) for the main comparison. Deterministic;
% one run per benchmark map (w_pot=1.0), with the same smoothing as A*/FMM-descent.
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

if ~exist('results/main/Classical_Baselines.mat', 'file')
    error('Missing results/main/Classical_Baselines.mat');
end
CB = load('results/main/Classical_Baselines.mat');

APF_AStar = struct();
fprintf('========================================\nAPF-A* Hybrid Baseline\n========================================\n\n');
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G,1);
    S = [1,1]; E = [M,M];

    tic;
    [path, ~] = APF_AStarPath(G, S, E, 1.0);
    t = toc;
    if ~isempty(path)
        path_s = GenerateSmoothPath(path, G);
        path_s = GenerateSmoothPath(path_s, G);
        len = 0;
        for i = 1:size(path_s,1)-1
            len = len + norm(path_s(i+1,:) - path_s(i,:));
        end
    else
        path_s = []; len = Inf;
    end
    APF_AStar(map_idx).Path = path_s;
    APF_AStar(map_idx).Length = len;
    APF_AStar(map_idx).Time = t;

    a_star_len = CB.Classical(map_idx).AStar.Length;
    gap_pct = 100*(len - a_star_len)/a_star_len;
    fprintf('Map %d: APF-A* length=%.2f (time=%.3fs), A*=%.2f, gap-to-A*=%+.2f%%\n', ...
        map_idx, len, t, a_star_len, gap_pct);
end

save('results/main/APFAStar_Baseline.mat', 'APF_AStar', 'Map_Info');
fprintf('\nSaved: results/main/APFAStar_Baseline.mat\nDone\n');
