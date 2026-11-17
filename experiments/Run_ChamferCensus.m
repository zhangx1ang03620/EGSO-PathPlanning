% =========================================================================
% Run_ChamferCensus.m
% A real chamfer distance transform (not just the
% Dijkstra/BFS proxy from the factorial ablation) run through the exact
% same exhaustive local-minima census protocol as Run_LocalMinimaCensus.m
% (Run_LocalMinimaCensus.m / GreedyDescentCensus.m), on the same 5 main
% maps and 5 stress-test maps, same seeds. Also reports BFS/Dijkstra in
% the same table for a complete three-way comparison: EVF(FMM) vs
% BFS/Dijkstra vs Chamfer vs the existing APF column.
% =========================================================================
clc; clear; close all;

if ~exist('results/localmin', 'dir')
    mkdir('results/localmin');
end

fprintf('========================================\n');
fprintf('Chamfer distance-field census (same protocol as the local-minima census table)\n');
fprintf('========================================\n\n');

katt = 1; apf_headline_rho0 = 5; apf_headline_krep = 50;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

MainResults = struct();
fprintf('---- five benchmark maps ----\n');
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    E = [M, M];
    goalMask = false(M, M); goalMask(E(1), E(2)) = true;

    T_nav = FastMarchingField(G, goalMask, 'navigate');
    T_bfs = BFSDistanceField(G, goalMask);
    T_cham = ChamferDistanceField(G, goalMask);

    [tr_evf, ns_evf]  = GreedyDescentCensus(T_nav, G, E);
    [tr_bfs, ns_bfs]  = GreedyDescentCensus(T_bfs, G, E);
    [tr_cham, ns_cham] = GreedyDescentCensus(T_cham, G, E);
    U = ArtificialPotentialField(G, E, katt, apf_headline_krep, apf_headline_rho0);
    [tr_apf, ns_apf]  = GreedyDescentCensus(U, G, E);

    MainResults(map_idx).EVF.trap_rate = tr_evf;   MainResults(map_idx).EVF.num_sinks = ns_evf;
    MainResults(map_idx).BFS.trap_rate = tr_bfs;   MainResults(map_idx).BFS.num_sinks = ns_bfs;
    MainResults(map_idx).Chamfer.trap_rate = tr_cham; MainResults(map_idx).Chamfer.num_sinks = ns_cham;
    MainResults(map_idx).APF.trap_rate = tr_apf;   MainResults(map_idx).APF.num_sinks = ns_apf;

    fprintf('Map%d: EVF trap=%6.2f%% sinks=%3d | BFS trap=%6.2f%% sinks=%3d | Chamfer trap=%6.2f%% sinks=%3d | APF trap=%6.2f%% sinks=%3d\n', ...
        map_idx, tr_evf, ns_evf, tr_bfs, ns_bfs, tr_cham, ns_cham, tr_apf, ns_apf);
end

fprintf('\n---- stress-test maps (decoy stages 1..5) ----\n');
StressResults = struct();
Mstress = 60;
for nstages = 1:5
    seed = 90000 + nstages;
    G = GenerateDecoyStressMap(Mstress, nstages, seed);
    E = [Mstress, Mstress];
    goalMask = false(Mstress, Mstress); goalMask(E(1), E(2)) = true;
    T = FastMarchingField(G, goalMask, 'navigate');
    attempt = 0;
    while ~isfinite(T(1,1)) && attempt < 5
        attempt = attempt + 1;
        seed = seed + 1000;
        G = GenerateDecoyStressMap(Mstress, nstages, seed);
        T = FastMarchingField(G, goalMask, 'navigate');
    end
    T_nav = T;
    T_bfs = BFSDistanceField(G, goalMask);
    T_cham = ChamferDistanceField(G, goalMask);

    [tr_evf, ns_evf]   = GreedyDescentCensus(T_nav, G, E);
    [tr_bfs, ns_bfs]   = GreedyDescentCensus(T_bfs, G, E);
    [tr_cham, ns_cham] = GreedyDescentCensus(T_cham, G, E);
    U = ArtificialPotentialField(G, E, katt, apf_headline_krep, apf_headline_rho0);
    [tr_apf, ns_apf]   = GreedyDescentCensus(U, G, E);

    StressResults(nstages).EVF.trap_rate = tr_evf;       StressResults(nstages).EVF.num_sinks = ns_evf;
    StressResults(nstages).BFS.trap_rate = tr_bfs;       StressResults(nstages).BFS.num_sinks = ns_bfs;
    StressResults(nstages).Chamfer.trap_rate = tr_cham;  StressResults(nstages).Chamfer.num_sinks = ns_cham;
    StressResults(nstages).APF.trap_rate = tr_apf;       StressResults(nstages).APF.num_sinks = ns_apf;

    fprintf('stages=%d: EVF trap=%6.2f%% sinks=%3d | BFS trap=%6.2f%% sinks=%3d | Chamfer trap=%6.2f%% sinks=%3d | APF trap=%6.2f%% sinks=%3d\n', ...
        nstages, tr_evf, ns_evf, tr_bfs, ns_bfs, tr_cham, ns_cham, tr_apf, ns_apf);
end

save('results/localmin/ChamferCensus_Data.mat', 'MainResults', 'StressResults', 'Map_Info');

fid = fopen('results/localmin/ChamferCensus_Summary.txt', 'w');
fprintf(fid, 'Local-Minima Census (same protocol as Run_LocalMinimaCensus.m): EVF(FMM) vs BFS/Dijkstra vs Chamfer vs APF\n');
fprintf(fid, 'Chamfer distance transform census\n');
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, '-- 5 main maps --\n');
fprintf(fid, '%-6s %14s %14s %14s %14s\n', 'Map', 'EVF', 'BFS', 'Chamfer', 'APF');
for map_idx = 1:5
    fprintf(fid, 'Map%-3d %6.2f%%/%-4d %6.2f%%/%-4d %6.2f%%/%-4d %6.2f%%/%-4d\n', map_idx, ...
        MainResults(map_idx).EVF.trap_rate, MainResults(map_idx).EVF.num_sinks, ...
        MainResults(map_idx).BFS.trap_rate, MainResults(map_idx).BFS.num_sinks, ...
        MainResults(map_idx).Chamfer.trap_rate, MainResults(map_idx).Chamfer.num_sinks, ...
        MainResults(map_idx).APF.trap_rate, MainResults(map_idx).APF.num_sinks);
end
fprintf(fid, '\n-- 5 stress-test maps (decoy stages 1-5) --\n');
fprintf(fid, '%-6s %14s %14s %14s %14s\n', 'Stage', 'EVF', 'BFS', 'Chamfer', 'APF');
for nstages = 1:5
    fprintf(fid, '%-6d %6.2f%%/%-4d %6.2f%%/%-4d %6.2f%%/%-4d %6.2f%%/%-4d\n', nstages, ...
        StressResults(nstages).EVF.trap_rate, StressResults(nstages).EVF.num_sinks, ...
        StressResults(nstages).BFS.trap_rate, StressResults(nstages).BFS.num_sinks, ...
        StressResults(nstages).Chamfer.trap_rate, StressResults(nstages).Chamfer.num_sinks, ...
        StressResults(nstages).APF.trap_rate, StressResults(nstages).APF.num_sinks);
end
fclose(fid);
fprintf('\nSaved: results/localmin/ChamferCensus_Data.mat / Summary.txt\n');
fprintf('\nDone\n');
