% Build_EndToEnd_Runtime_Table.m
% Adds the one-off EVF precomputation to EGSO's end-to-end time, for a fair
% comparison with RRT*/Informed-RRT* (which have no precomputation step).
clc; clear;

EV = load('results/main/EVFPrecompute_MainMaps.mat');
SB = load('results/sampling/Sampling_Baselines.mat');
CR = load('results/main/EGSO_Comparison_Results.mat');

disp(fieldnames(SB));
disp(fieldnames(CR));

Map_Info = EV.Map_Info;
fid = fopen('results/main/EndToEnd_Runtime_Table.txt', 'w');
fprintf(fid, 'End-to-End Runtime: EGSO (EVF precompute + per-run solve) vs RRT*/Informed-RRT* (no precompute)\n');
fprintf(fid, '=========================================================================\n\n');
header = sprintf('%-8s %-14s %-14s %-16s %-14s %-18s', 'Map', 'EVF-pre(s)', 'EGSO-solve(s)', 'EGSO-total(s)', 'RRT*(s)', 'Informed-RRT*(s)');
fprintf('%s\n', header); fprintf(fid, '%s\n', header);

for map_idx = 1:5
    evf_t = EV.EVF_Precompute_Time(map_idx);
    egso_solve_t = mean(CR.Results(map_idx).EGSO.Time);
    egso_total_t = evf_t + egso_solve_t;

    rrt_t = NaN; irrt_t = NaN;
    if isfield(SB, 'Sampling')
        S = SB.Sampling;
        if isfield(S, 'RRTStar')
            rrt_t = mean(S(map_idx).RRTStar.Time);
        end
        if isfield(S, 'InformedRRTStar')
            irrt_t = mean(S(map_idx).InformedRRTStar.Time);
        end
    end

    line = sprintf('%-8d %-14.4f %-14.4f %-16.4f %-14.4f %-18.4f', map_idx, evf_t, egso_solve_t, egso_total_t, rrt_t, irrt_t);
    fprintf('%s\n', line); fprintf(fid, '%s\n', line);
end
fclose(fid);
fprintf('\nSaved: results/main/EndToEnd_Runtime_Table.txt\n');
