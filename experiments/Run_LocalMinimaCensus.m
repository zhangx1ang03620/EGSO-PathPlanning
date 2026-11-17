% =========================================================================
% Run_LocalMinimaCensus.m
%
% Turns the claim "the Eikonal viscosity field has no spurious local minimum (Sethian)"
% from a cited theorem into reproducible numbers:
%
%   on the five benchmark maps, exhaustively compares the EVF field (T_nav) with a
%   standard artificial potential field (APF, Khatib 1986 form, built only from the two
%   local quantities "Euclidean distance to the goal" and "distance to the nearest obstacle"):
%     1. num_sinks   -- number of local-minimum ("trap") cells other than the goal
%     2. trap_rate   -- fraction of free cells from which greedy 4-neighbour descent
%                        does not reach the goal
%
%   The APF is swept over a 3x3 grid of (rho0, krep) so the result does not hinge on a
%   poorly tuned APF; the EVF needs no tuning (in theory num_sinks=0, trap_rate=0%).
%
%   The census is repeated on the stress-test maps (decoy stages 1..5) to show the
%   APF trap count growing with topological complexity while the EVF stays at 0.
% =========================================================================
clc; clear; close all;

if ~exist('results/localmin', 'dir')
    mkdir('results/localmin');
end

fprintf('========================================\n');
fprintf('Local-minima census (APF vs EVF)\n');
fprintf('========================================\n\n');

katt = 1;
rho0_list = [3, 5, 8];
krep_list = [20, 50, 100];
apf_headline_rho0 = 5;
apf_headline_krep = 50;

%% -- Five benchmark maps -----------------------------------------------------
[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

MainResults = struct();
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    E = [M, M];
    goalMask = false(M, M); goalMask(E(1), E(2)) = true;
    T_nav = FastMarchingField(G, goalMask, 'navigate');

    [tr_evf, ns_evf, reach_evf] = GreedyDescentCensus(T_nav, G, E);
    fprintf('Map%d EVF          : trap_rate=%6.2f%%  num_sinks=%3d\n', map_idx, tr_evf, ns_evf);

    MainResults(map_idx).EVF.trap_rate = tr_evf;
    MainResults(map_idx).EVF.num_sinks = ns_evf;
    MainResults(map_idx).EVF.reach_mask = reach_evf;

    apf_tr = zeros(numel(rho0_list), numel(krep_list));
    apf_ns = zeros(numel(rho0_list), numel(krep_list));
    for ri = 1:numel(rho0_list)
        for ki = 1:numel(krep_list)
            rho0 = rho0_list(ri); krep = krep_list(ki);
            U = ArtificialPotentialField(G, E, katt, krep, rho0);
            [tr, ns] = GreedyDescentCensus(U, G, E);
            apf_tr(ri, ki) = tr;
            apf_ns(ri, ki) = ns;
            if rho0 == apf_headline_rho0 && krep == apf_headline_krep
                MainResults(map_idx).APF_headline.trap_rate = tr;
                MainResults(map_idx).APF_headline.num_sinks = ns;
                U_headline = ArtificialPotentialField(G, E, katt, krep, rho0);
                [~, ~, reach_apf] = GreedyDescentCensus(U_headline, G, E);
                MainResults(map_idx).APF_headline.reach_mask = reach_apf;
            end
            fprintf('Map%d APF rho0=%d krep=%-4d: trap_rate=%6.2f%%  num_sinks=%3d\n', ...
                    map_idx, rho0, krep, tr, ns);
        end
    end
    MainResults(map_idx).APF_sweep.rho0_list = rho0_list;
    MainResults(map_idx).APF_sweep.krep_list = krep_list;
    MainResults(map_idx).APF_sweep.trap_rate = apf_tr;
    MainResults(map_idx).APF_sweep.num_sinks = apf_ns;
    fprintf('\n');
end

%% -- Stress-test maps (decoy stages 1..5) ------------------------------------
fprintf('---- stress-test maps (decoy stages 1..5), APF with the headline parameters ----\n');
StressResults = struct();
Mstress = 60;
for nstages = 1:5
    seed = 90000 + nstages;
    G = GenerateDecoyStressMap(Mstress, nstages, seed);
    E = [Mstress, Mstress];
    goalMask = false(Mstress, Mstress); goalMask(E(1), E(2)) = true;
    T_nav = FastMarchingField(G, goalMask, 'navigate');
    T = FastMarchingField(G, goalMask, 'navigate');
    attempt = 0;
    while ~isfinite(T(1,1)) && attempt < 5
        attempt = attempt + 1;
        seed = seed + 1000;
        G = GenerateDecoyStressMap(Mstress, nstages, seed);
        T = FastMarchingField(G, goalMask, 'navigate');
    end
    T_nav = T;

    [tr_evf, ns_evf] = GreedyDescentCensus(T_nav, G, E);
    U = ArtificialPotentialField(G, E, katt, apf_headline_krep, apf_headline_rho0);
    [tr_apf, ns_apf] = GreedyDescentCensus(U, G, E);

    StressResults(nstages).EVF.trap_rate = tr_evf;
    StressResults(nstages).EVF.num_sinks = ns_evf;
    StressResults(nstages).APF.trap_rate = tr_apf;
    StressResults(nstages).APF.num_sinks = ns_apf;

    fprintf('stages=%d: EVF trap_rate=%6.2f%% sinks=%3d  |  APF trap_rate=%6.2f%% sinks=%3d\n', ...
            nstages, tr_evf, ns_evf, tr_apf, ns_apf);
end

save('results/localmin/LocalMinimaCensus_Data.mat', 'MainResults', 'StressResults', ...
     'rho0_list', 'krep_list', 'katt', 'apf_headline_rho0', 'apf_headline_krep');

%% -- Text summary tables -----------------------------------------------------
fid = fopen('results/localmin/LocalMinimaCensus_Summary.txt', 'w');
fprintf(fid, 'Local-Minima Census: EVF (T_nav) vs Artificial Potential Field (APF)\n');
fprintf(fid, '========================================================================\n\n');
fprintf(fid, 'Table 1: Five benchmark maps, EVF vs APF(headline: rho0=%d, krep=%d, katt=%d)\n', ...
        apf_headline_rho0, apf_headline_krep, katt);
fprintf(fid, '%-8s %14s %14s %14s %14s\n', 'Map', 'EVF trap(%)', 'EVF sinks', 'APF trap(%)', 'APF sinks');
for map_idx = 1:5
    fprintf(fid, '%-8s %14.2f %14d %14.2f %14d\n', sprintf('Map%d', map_idx), ...
        MainResults(map_idx).EVF.trap_rate, MainResults(map_idx).EVF.num_sinks, ...
        MainResults(map_idx).APF_headline.trap_rate, MainResults(map_idx).APF_headline.num_sinks);
end

fprintf(fid, '\nTable 2: APF (rho0 x krep) parameter sweep -- trap_rate(%%) / num_sinks, all 5 maps\n');
for map_idx = 1:5
    fprintf(fid, '\n-- Map%d --\n', map_idx);
    fprintf(fid, '%10s', 'rho0\\krep');
    for ki = 1:numel(krep_list)
        fprintf(fid, ' %14d', krep_list(ki));
    end
    fprintf(fid, '\n');
    for ri = 1:numel(rho0_list)
        fprintf(fid, '%10d', rho0_list(ri));
        for ki = 1:numel(krep_list)
            fprintf(fid, ' %8.2f%%/%-4d', MainResults(map_idx).APF_sweep.trap_rate(ri,ki), ...
                    MainResults(map_idx).APF_sweep.num_sinks(ri,ki));
        end
        fprintf(fid, '\n');
    end
end

fprintf(fid, '\nTable 3: Stress-test maps (decoy-trap stages 1..5), EVF vs APF(headline params)\n');
fprintf(fid, '%-8s %14s %14s %14s %14s\n', 'Stages', 'EVF trap(%)', 'EVF sinks', 'APF trap(%)', 'APF sinks');
for nstages = 1:5
    fprintf(fid, '%-8d %14.2f %14d %14.2f %14d\n', nstages, ...
        StressResults(nstages).EVF.trap_rate, StressResults(nstages).EVF.num_sinks, ...
        StressResults(nstages).APF.trap_rate, StressResults(nstages).APF.num_sinks);
end
fclose(fid);
fprintf('\nSaved: results/localmin/LocalMinimaCensus_Summary.txt\n');

%% -- LaTeX table ------------------------------------------------------------
fid = fopen('results/localmin/LocalMinimaCensus_LaTeX.txt', 'w');
fprintf(fid, '%% Table: Local-minima census, EVF vs APF (headline rho0=%d, krep=%d)\n', ...
        apf_headline_rho0, apf_headline_krep);
fprintf(fid, '\\begin{table}[h]\n\\caption{Exhaustive local-minima census: EVF ($T_{nav}$) vs. a standard artificial potential field (APF).}\\label{tab:localmin}\n');
fprintf(fid, '\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lcccc@{}}\n\\toprule\n');
fprintf(fid, 'Map & EVF Trap Rate (\\%%) & EVF \\#Sinks & APF Trap Rate (\\%%) & APF \\#Sinks\\\\\n\\midrule\n');
for map_idx = 1:5
    fprintf(fid, 'Map %d & %.2f & %d & %.2f & %d\\\\\n', map_idx, ...
        MainResults(map_idx).EVF.trap_rate, MainResults(map_idx).EVF.num_sinks, ...
        MainResults(map_idx).APF_headline.trap_rate, MainResults(map_idx).APF_headline.num_sinks);
end
fprintf(fid, '\\botrule\n\\end{tabular*}\n\\end{table}\n');
fclose(fid);
fprintf('Saved: results/localmin/LocalMinimaCensus_LaTeX.txt\n');

%% -- Figure 1: trap basins on Map 1 (EVF vs headline APF) --------------------
G1 = Maps{1};
M1 = size(G1,1);
fig = figure('Position', [80,80,1300,650], 'Visible', 'off', 'Color', 'w');

subplot(1,2,1);
img = ones(M1, M1, 3);
reach = MainResults(1).EVF.reach_mask;
for r = 1:M1
    for c = 1:M1
        if G1(r,c) == 1
            img(r,c,:) = [0.25 0.25 0.25];
        elseif reach(r,c)
            img(r,c,:) = [0.75 0.92 0.75];
        else
            img(r,c,:) = [0.95 0.35 0.35];
        end
    end
end
image(img); axis equal tight; set(gca,'YDir','normal');
hold on; plot(1,1,'go','MarkerFaceColor','g','MarkerSize',8); plot(M1,M1,'bs','MarkerFaceColor','b','MarkerSize',8);
title(sprintf('EVF (T_{nav}): trap rate = %.2f%%, sinks = %d', ...
      MainResults(1).EVF.trap_rate, MainResults(1).EVF.num_sinks), 'FontSize', 11, 'FontWeight','bold');

subplot(1,2,2);
img2 = ones(M1, M1, 3);
reach2 = MainResults(1).APF_headline.reach_mask;
for r = 1:M1
    for c = 1:M1
        if G1(r,c) == 1
            img2(r,c,:) = [0.25 0.25 0.25];
        elseif reach2(r,c)
            img2(r,c,:) = [0.75 0.92 0.75];
        else
            img2(r,c,:) = [0.95 0.35 0.35];
        end
    end
end
image(img2); axis equal tight; set(gca,'YDir','normal');
hold on; plot(1,1,'go','MarkerFaceColor','g','MarkerSize',8); plot(M1,M1,'bs','MarkerFaceColor','b','MarkerSize',8);
title(sprintf('APF (\\rho_0=%d, k_{rep}=%d): trap rate = %.2f%%, sinks = %d', ...
      apf_headline_rho0, apf_headline_krep, MainResults(1).APF_headline.trap_rate, MainResults(1).APF_headline.num_sinks), ...
      'FontSize', 11, 'FontWeight','bold');

sgtitle('Map1: Cells whose greedy descent reaches the goal (green) vs. gets trapped (red)', ...
        'FontSize', 13, 'FontWeight', 'bold');
print(fig, 'results/localmin/Fig_TrapBasins_Map1', '-dpng', '-r300');
close(fig);
fprintf('Saved: results/localmin/Fig_TrapBasins_Map1.png\n');

%% -- Figure 2: trap rate on the five maps (bar chart) ------------------------
fig2 = figure('Position', [100,100,800,500], 'Visible', 'off', 'Color', 'w');
tr_evf_all = arrayfun(@(x) x.EVF.trap_rate, MainResults);
tr_apf_all = arrayfun(@(x) x.APF_headline.trap_rate, MainResults);
b = bar([tr_evf_all; tr_apf_all]');
b(1).FaceColor = [0.85 0.10 0.10]; b(2).FaceColor = [0.20 0.45 0.75];
set(gca, 'XTickLabel', {'Map1','Map2','Map3','Map4','Map5'}, 'FontSize', 11);
ylabel('Trap Rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'EVF (T_{nav})', 'APF (headline)'}, 'Location', 'northwest', 'FontSize', 11);
title('Fraction of free cells whose greedy descent fails to reach the goal', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
print(fig2, 'results/localmin/Fig_TrapRate_Bar', '-dpng', '-r300');
close(fig2);
fprintf('Saved: results/localmin/Fig_TrapRate_Bar.png\n');

%% -- Figure 3: number of sinks versus stress-test complexity -----------------
fig3 = figure('Position', [100,100,800,500], 'Visible', 'off', 'Color', 'w');
ns_evf_stress = arrayfun(@(x) x.EVF.num_sinks, StressResults);
ns_apf_stress = arrayfun(@(x) x.APF.num_sinks, StressResults);
plot(1:5, ns_evf_stress, '-o', 'Color', [0.85 0.10 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.10 0.10]);
hold on;
plot(1:5, ns_apf_stress, '-s', 'Color', [0.20 0.45 0.75], 'LineWidth', 2, 'MarkerFaceColor', [0.20 0.45 0.75]);
xlabel('Decoy-trap stages (topological complexity level)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of non-goal local minima (sinks)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'EVF (T_{nav})', 'APF (headline)'}, 'Location', 'northwest', 'FontSize', 11);
title('Local-minima count vs. topological complexity', 'FontSize', 13, 'FontWeight', 'bold');
xticks(1:5); grid on;
print(fig3, 'results/localmin/Fig_SinksVsComplexity', '-dpng', '-r300');
close(fig3);
fprintf('Saved: results/localmin/Fig_SinksVsComplexity.png\n');

fprintf('\n========================================\n');
fprintf('Local-minima census finished\n');
fprintf('========================================\n');
