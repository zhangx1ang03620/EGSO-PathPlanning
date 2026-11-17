% =========================================================================
% Run_EVFPrecompute_MainMaps.m
% EVF precomputation time on the five benchmark maps. The scalability study measures
% it only on the warehouse-style maps from M = 90 upwards. RRT*/Informed-RRT* have no
% precomputation step, so a fair end-to-end comparison must add the one-off EVF cost
% to EGSO's solve time.
%
%
% This script times T_nav + T_clear on each of the five benchmark maps (median of 10
% repetitions) and provides the EVF column of the end-to-end runtime table
% (EGSO = EVF precompute + solve time vs. RRT*/Informed-RRT*).
% =========================================================================
clc; clear;

if ~exist('results/main', 'dir')
    mkdir('results/main');
end

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

num_repeat = 10;
EVF_Precompute_Time = zeros(5, 1);
EVF_Precompute_Std   = zeros(5, 1);

fprintf('========================================\n');
fprintf('EVF (T_nav + T_clear) precomputation time on the five benchmark maps\n');
fprintf('========================================\n\n');

for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    goalMask = false(M, M); goalMask(M, M) = true;
    obstacleMask = (G == 1);

    times = zeros(num_repeat, 1);
    for r = 1:num_repeat
        tic;
        T_nav   = FastMarchingField(G, goalMask, 'navigate'); %#ok<NASGU>
        T_clear = FastMarchingField(G, obstacleMask, 'clearance'); %#ok<NASGU>
        times(r) = toc;
    end
    EVF_Precompute_Time(map_idx) = median(times);
    EVF_Precompute_Std(map_idx)  = std(times);
    fprintf('Map %d (%dx%d): EVF precompute = %.4f s (median of %d, std=%.4f)\n', ...
        map_idx, M, M, EVF_Precompute_Time(map_idx), num_repeat, EVF_Precompute_Std(map_idx));
end

save('results/main/EVFPrecompute_MainMaps.mat', 'EVF_Precompute_Time', 'EVF_Precompute_Std', 'Map_Info');

fid = fopen('results/main/EVFPrecompute_MainMaps.txt', 'w');
fprintf(fid, 'EVF (T_nav + T_clear) precomputation time on the 5 main benchmark maps (median of %d repeats)\n', num_repeat);
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, '%-8s %-10s %-16s %-12s\n', 'Map', 'Size', 'Precompute(s)', 'Std(s)');
for map_idx = 1:5
    fprintf(fid, '%-8d %-10s %-16.4f %-12.4f\n', map_idx, Map_Info(map_idx).size, ...
        EVF_Precompute_Time(map_idx), EVF_Precompute_Std(map_idx));
end
fclose(fid);
fprintf('\nSaved: results/main/EVFPrecompute_MainMaps.mat / .txt\n');

%% If the sampling-baseline results exist, build the end-to-end runtime table directly
if exist('results/sampling/Sampling_Results.mat', 'file') || exist('results/main/Sampling_Baselines.mat', 'file')
    cand = {'results/sampling/Sampling_Results.mat', 'results/main/Sampling_Baselines.mat'};
    for c = 1:numel(cand)
        if exist(cand{c}, 'file')
            SR = load(cand{c});
            fprintf('\n[note] found sampling-baseline results %s; the end-to-end table can be assembled from them\n', cand{c});
            break;
        end
    end
else
    fprintf('\n[note] no saved RRT*/Informed-RRT* results found; build the end-to-end table after running the sampling baselines\n');
end
