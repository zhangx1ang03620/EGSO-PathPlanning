function Run_Comparison_SingleMap(map_idx)
% Run_Comparison_SingleMap  Runs the full 23-algorithm comparison for ONE
% map only (paper-standard N=50, T=100, 20 runs) and saves a partial
% result file. Designed to be launched as 5 independent parallel MATLAB
% batch processes (one per map) to cut wall-clock time roughly 5x versus
% running all 5 maps sequentially in one process, since each map's
% experiment is fully independent.
%
% Usage (from shell): matlab -batch "Run_Comparison_SingleMap(1)"
% Aggregate all 5 partial files afterwards with Aggregate_Comparison_Results.m

fprintf('========================================\n');
fprintf('Running full comparison for Map %d only (paper-standard params)\n', map_idx);
fprintf('========================================\n\n');

if ~exist('results/main', 'dir')
    mkdir('results/main');
end

[Maps, Map_Info_full] = GenerateEGSOMaps();
close(gcf);

max_gen        = 100;
num_population = 50;
num_runs       = 20;
collision_threshold = 1e4;

Algorithm_Names = {'EGSO', 'GWO', 'PSO', 'GA', 'DE', 'SSA', ...
                   'WOA', 'MFO', 'ALO', 'HHO', 'SCA', ...
                   'ABC', 'ACO', 'BBO', 'CS', 'FA', ...
                   'TLBO', 'BA', 'FPA', 'SMA', 'INFO', 'SBOA', 'ARO'};
num_algorithms = length(Algorithm_Names);

G = Maps{map_idx};
Map_Info_this = Map_Info_full(map_idx);
[Xmax, ~] = size(G);
X_min = 1;
dimensions = Xmax - 2;
S = [1, 1];
E = [Xmax, Xmax];
fobj = @(x) fitness(x, G);

MapResult = struct();
for algo_idx = 1:num_algorithms
    algo_name = Algorithm_Names{algo_idx};
    MapResult.(algo_name).Best   = zeros(num_runs, 1);
    MapResult.(algo_name).Curves = zeros(num_runs, max_gen);
    MapResult.(algo_name).Paths  = cell(num_runs, 1);
    MapResult.(algo_name).Time   = zeros(num_runs, 1);
end

for run = 1:num_runs
    fprintf('[Map %d] Run %d/%d ...\n', map_idx, run, num_runs);

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        try
            tic;
            switch algo_name
                case 'EGSO', [Best_score, Best_pos, Curve] = EGSO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'GWO',  [Best_score, Best_pos, Curve] = GWO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'PSO',  [Best_score, Best_pos, Curve] = PSO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'GA',   [Best_score, Best_pos, Curve] = GA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'DE',   [Best_score, Best_pos, Curve] = DE(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'SSA',  [Best_score, Best_pos, Curve] = SSA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'WOA',  [Best_score, Best_pos, Curve] = WOA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'MFO',  [Best_score, Best_pos, Curve] = MFO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'ALO',  [Best_score, Best_pos, Curve] = ALO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'HHO',  [Best_score, Best_pos, Curve] = HHO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'SCA',  [Best_score, Best_pos, Curve] = SCA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'ABC',  [Best_score, Best_pos, Curve] = ABC(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'ACO',  [Best_score, Best_pos, Curve] = ACO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'BBO',  [Best_score, Best_pos, Curve] = BBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'CS',   [Best_score, Best_pos, Curve] = CS(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'FA',   [Best_score, Best_pos, Curve] = FA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'TLBO', [Best_score, Best_pos, Curve] = TLBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'BA',   [Best_score, Best_pos, Curve] = BA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'FPA',  [Best_score, Best_pos, Curve] = FPA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'SMA',  [Best_score, Best_pos, Curve] = SMA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'INFO', [Best_score, Best_pos, Curve] = INFO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'SBOA', [Best_score, Best_pos, Curve] = SBOA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 'ARO',  [Best_score, Best_pos, Curve] = ARO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
            end
            elapsed_time = toc;

            Best_pos = round(Best_pos);
            Best_pos = max(Best_pos, 1);
            Best_pos = min(Best_pos, Xmax);
            route = [S(1) Best_pos E(1)];
            path = generateContinuousRoute(route, G);
            if ~isempty(path)
                path = GenerateSmoothPath(path, G);
                path = GenerateSmoothPath(path, G);
            end

            MapResult.(algo_name).Best(run)   = Best_score;
            MapResult.(algo_name).Curves(run,:) = Curve;
            MapResult.(algo_name).Paths{run}  = path;
            MapResult.(algo_name).Time(run)   = elapsed_time;

            status = 'ok';
            if Best_score >= collision_threshold
                status = 'COLLISION';
            end
            fprintf('  [Map %d][%s] %s: %.4f (%.2fs)\n', map_idx, status, algo_name, Best_score, elapsed_time);
        catch ME
            fprintf('  [Map %d][ERROR] %s: %s\n', map_idx, algo_name, ME.message);
            MapResult.(algo_name).Best(run) = inf;
        end
    end
end

save(sprintf('results/main/partial_map%d.mat', map_idx), 'MapResult', 'Map_Info_this', ...
     'Algorithm_Names', 'num_runs', 'max_gen', 'num_population', 'map_idx');
fprintf('\n[Map %d] DONE. Saved results/main/partial_map%d.mat\n', map_idx, map_idx);

end
