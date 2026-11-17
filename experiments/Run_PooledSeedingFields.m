function Run_PooledSeedingFields()
% Run_PooledSeedingFields  Comparison of seeding fields (Sec. 4.5 of the paper).
% The three elite-field cells of each initialization in the crossed factorial
% ablation (results/factorial2) are runs of the same algorithm, because the
% elite field-directed step of the earlier implementation never changed a
% candidate (removing it leaves every result bit-identical). Pool them
% into 30 runs per initialization: exact pooled mean and sample SD from
% the per-cell means and sample SDs (n = 10 each, all runs valid).
% Run from the repository root.
F = load(fullfile('results', 'factorial2', 'FactorialAblation2_Data.mat'));
out_dir = fullfile('results', 'factorial_pooled');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
assert(all(F.SR_data(:) == 100), 'all cells must have 100% SR for exact pooling');
n = F.num_runs;
labels = cellfun(@(c) c{4}, F.Cells, 'UniformOutput', false);
groups = {'GSI(eikonal)', 'GSI(apf)', 'GSI(dijkstra)', 'Random'};
pm = nan(4, 3); ps = nan(4, 3); members = cell(4, 1);
for g = 1:4
    idx = find(startsWith(labels, groups{g}));
    assert(numel(idx) == 3, 'group %s', groups{g});
    members{g} = labels(idx);
    for mi = 1:3
        m = F.Mean_data(idx, mi); s = F.Std_data(idx, mi);
        N = n * numel(idx);
        mu = sum(n * m) / N;
        ss = sum((n - 1) * s.^2 + n * (m - mu).^2);
        pm(g, mi) = mu; ps(g, mi) = sqrt(ss / (N - 1));
    end
end
fid = fopen(fullfile(out_dir, 'PooledSeedingFields_Summary.txt'), 'w');
for out = [1 fid]
    fprintf(out, 'Pooled seeding-field comparison (Maps %s, 30 runs per row, V.Mean +/- SD)\n', mat2str(F.map_indices));
    for g = 1:4
        fprintf(out, '%-14s', groups{g});
        for mi = 1:3, fprintf(out, '  %8.2f +/- %5.2f', pm(g, mi), ps(g, mi)); end
        fprintf(out, '   [%s]\n', strjoin(members{g}, '; '));
    end
    fprintf(out, 'Average over maps: %s\n', mat2str(round(mean(pm, 2)', 2)));
    fprintf(out, 'Map 1 APF vs Eikonal: %+.1f%%; Random vs GSI(eikonal) average: %+.1f%%\n', ...
        100 * (pm(2, 1) / pm(1, 1) - 1), 100 * (mean(pm(4, :) ./ pm(1, :)) - 1));
end
fclose(fid);
save(fullfile(out_dir, 'PooledSeedingFields_Data.mat'), 'pm', 'ps', 'groups', 'members');
end
