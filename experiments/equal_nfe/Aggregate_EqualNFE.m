function Aggregate_EqualNFE(resdir)
% Aggregate_EqualNFE  Summarise the equal-NFE experiment (Run_EqualNFE_Job).
%
% Produces, for each MaxFEs budget:
%   * per-map SR / V.Best / V.Mean / V.Worst / V.Std / NFE actually used,
%   * per-map ranks by V.Mean (0%-SR algorithms tie for last place, average
%     rank) and the average rank, side by side with the T = 100 protocol,
%   * Wilcoxon rank-sum W/T/L of EGSO vs. every competitor (valid-only and
%     all-runs incl. collision penalties, same protocol as the manuscript),
%   * the Friedman test on the 5 x 23 V.Mean matrix,
% plus the performance-versus-NFE figure (15,000 budget) and LaTeX snippets.
% Outputs go to resdir (default: results/equal_nfe).

    here = fileparts(mfilename('fullpath'));
    project_root = fullfile(here, '..', '..');
    addpath(project_root); setup_egso();
    cd(here);

    if nargin < 1, resdir = fullfile(project_root, 'results', 'equal_nfe'); end
    budgets = [5000 15000];
    collision_threshold = 1e4;
    alpha = 0.05;
    astar = [77.11 132.19 86.45 117.44 139.68];   % A* path lengths (results/main/Classical_vs_EGSO_Summary.txt)

    % ---------- load ----------
    D = struct();
    for bi = 1:numel(budgets)
        fe = budgets(bi);
        for m = 1:5
            files = dir(fullfile(resdir, sprintf('EqNFE_map%d_fe%d_r*.mat', m, fe)));
            assert(~isempty(files), 'no results for map %d, MaxFEs %d', m, fe);
            B = []; NF = []; TR = []; TM = []; ER = {};
            for f = 1:numel(files)
                S = load(fullfile(files(f).folder, files(f).name));
                B  = [B;  S.Best];        %#ok<AGROW>
                NF = [NF; S.NFEused];     %#ok<AGROW>
                TR = cat(1, TR, S.Trace);
                TM = [TM; S.Time];        %#ok<AGROW>
                ER = [ER; S.ErrMsg];      %#ok<AGROW>
            end
            Algos = S.Algos;
            D(bi, m).Best = B; D(bi, m).NFE = NF; D(bi, m).Trace = TR;
            D(bi, m).Time = TM; D(bi, m).ckpt = S.ckpt;
            D(bi, m).nErr = sum(~cellfun(@isempty, ER(:)));
        end
    end
    nA = numel(Algos);
    egso = find(strcmp(Algos, 'EGSO'));
    others = setdiff(1:nA, egso);

    % T = 100 protocol reference (main manuscript results)
    R100 = load(fullfile(project_root, 'results', 'main', 'EGSO_Comparison_Results.mat'));
    VM100 = nan(5, nA);
    for m = 1:5
        for a = 1:nA
            v = R100.Results(m).(Algos{a}).Mean;
            if isfinite(v), VM100(m, a) = v; end
        end
    end
    Rank100 = rank_rows(VM100);

    fid = fopen(fullfile(resdir, 'EqualNFE_Summary.txt'), 'w');
    out = @(varargin) fprintf_both(fid, varargin{:});
    out('Equal-NFE comparison, N = 50, 23 algorithms, 5 maps\n');
    out('Result of a run = best fitness evaluated within MaxFEs (anytime best); valid if < %g\n', collision_threshold);
    out('Ranks: by V.Mean per map; 0%%-SR algorithms tie for last (average rank)\n\n');

    Summary = struct();
    for bi = 1:numel(budgets)
        fe = budgets(bi);
        VM = nan(5, nA); SR = zeros(5, nA); VB = nan(5, nA); VW = nan(5, nA); VS = nan(5, nA);
        NFEm = zeros(5, nA); NFEmin = zeros(5, nA); TimeM = zeros(5, nA);
        WTLv = zeros(nA, 3); WTLa = zeros(nA, 3); detail = cell(5, nA);
        for m = 1:5
            B = D(bi, m).Best;
            for a = 1:nA
                v = B(:, a); val = v(v < collision_threshold);
                SR(m, a) = 100 * numel(val) / numel(v);
                if ~isempty(val)
                    VM(m, a) = mean(val); VB(m, a) = min(val); VW(m, a) = max(val); VS(m, a) = std(val);
                end
                NFEm(m, a) = mean(D(bi, m).NFE(:, a)); NFEmin(m, a) = min(D(bi, m).NFE(:, a));
                TimeM(m, a) = mean(D(bi, m).Time(:, a));
            end
            % Wilcoxon, same protocol as Run_SignificanceTests_v2.m
            e_all = B(:, egso); e_val = e_all(e_all < collision_threshold);
            for a = others
                o_all = B(:, a); o_val = o_all(o_all < collision_threshold);
                if isempty(o_val)
                    WTLv(a, 1) = WTLv(a, 1) + 1; dv = '+(0%SR)';
                elseif isempty(e_val)
                    dv = 'N/A';
                else
                    p = ManualRankSum(e_val, o_val);
                    dv = classify(p, mean(e_val), mean(o_val), alpha);
                    WTLv(a, :) = WTLv(a, :) + tally(dv);
                end
                pa = ManualRankSum(e_all, o_all);
                da = classify(pa, mean(e_all), mean(o_all), alpha);
                WTLa(a, :) = WTLa(a, :) + tally(da);
                detail{m, a} = sprintf('%s / %s', dv, da);
            end
        end
        Rk = rank_rows(VM);
        avgRk = mean(Rk, 1);
        VMf = VM;
        for m = 1:5                                   % 0%-SR -> worse than any valid mean
            row = VMf(m, :); row(isnan(row)) = max(row(isfinite(row))) * 1.5; VMf(m, :) = row;
        end
        [chi2, df, pF] = ManualFriedman(VMf);

        Summary(bi).fe = fe; Summary(bi).VM = VM; Summary(bi).SR = SR; Summary(bi).VB = VB;
        Summary(bi).VW = VW; Summary(bi).VS = VS; Summary(bi).Rank = Rk; Summary(bi).avgRank = avgRk;
        Summary(bi).NFEmean = NFEm; Summary(bi).NFEmin = NFEmin; Summary(bi).Time = TimeM;
        Summary(bi).WTLvalid = WTLv; Summary(bi).WTLall = WTLa; Summary(bi).detail = detail;
        Summary(bi).friedman = [chi2 df pF];

        out('==================== MaxFEs = %d ====================\n', fe);
        for m = 1:5
            out('\n--- Map %d (A* = %.2f) ---  errors: %d\n', m, astar(m), D(bi, m).nErr);
            out('%-6s %6s %8s %8s %8s %7s %5s %9s %9s %8s\n', 'Algo', 'SR%', 'V.Best', 'V.Mean', 'V.Worst', 'V.Std', 'Rank', 'NFEmean', 'NFEmin', 'Time(s)');
            [~, ord] = sort(Rk(m, :));
            for a = ord
                out('%-6s %6.1f %8.2f %8.2f %8.2f %7.2f %5.1f %9.0f %9.0f %8.2f\n', Algos{a}, SR(m, a), VB(m, a), VM(m, a), VW(m, a), VS(m, a), Rk(m, a), NFEm(m, a), NFEmin(m, a), TimeM(m, a));
            end
        end
        out('\nAverage rank (MaxFEs %d) vs. T = 100 protocol:\n', fe);
        [~, ord] = sort(avgRk);
        out('%-6s %s   %s\n', 'Algo', 'ranks M1..M5 (avg)  [MaxFEs]', 'ranks M1..M5 (avg)  [T=100]');
        for a = ord
            out('%-6s %4.1f %4.1f %4.1f %4.1f %4.1f (%5.2f)   %4.1f %4.1f %4.1f %4.1f %4.1f (%5.2f)\n', Algos{a}, Rk(:, a), avgRk(a), Rank100(:, a), mean(Rank100(:, a)));
        end
        n_comp = numel(others) * 5;
        out('\nWilcoxon EGSO vs others (%d comparisons):\n', n_comp);
        out('  valid-only: W=%d T=%d L=%d\n', sum(WTLv(:, 1)), sum(WTLv(:, 2)), sum(WTLv(:, 3)));
        out('  all-runs  : W=%d T=%d L=%d\n', sum(WTLa(:, 1)), sum(WTLa(:, 2)), sum(WTLa(:, 3)));
        for a = others
            if any(WTLv(a, 2:3)) || any(WTLa(a, 2:3))
                out('  %-6s valid %d/%d/%d  all %d/%d/%d   per map: %s\n', Algos{a}, WTLv(a, :), WTLa(a, :), strjoin(detail(:, a)', ' | '));
            end
        end
        out('Friedman: chi2 = %.2f, df = %d, p = %.3g\n\n', chi2, df, pF);
        out('EGSO vs best competitor by V.Mean:\n');
        for m = 1:5
            vm = VM(m, :); vm(egso) = Inf; [bv, bk] = min(vm);
            out('  Map %d: EGSO %.2f (SR %.0f%%) | best other %s %.2f | EGSO gap to A* %+.2f%%\n', m, VM(m, egso), SR(m, egso), Algos{bk}, bv, 100 * (VM(m, egso) / astar(m) - 1));
        end
        out('\n');
    end
    fclose(fid);
    save(fullfile(resdir, 'EqualNFE_Summary.mat'), 'Summary', 'Algos', 'budgets', 'Rank100', 'VM100', 'astar');

    plot_perf_vs_nfe(resdir, D, Summary, Algos, budgets, astar, egso);
    write_latex(resdir, Summary, Algos, Rank100, egso, others);
    fprintf('Aggregation done.\n');
end

% ---------------------------------------------------------------------
function d = classify(p, me, mo, alpha)
    if p < alpha
        if me < mo, d = '+'; else, d = '-'; end
    else
        d = '~';
    end
end

function t = tally(d)
    t = [strcmp(d, '+') strcmp(d, '~') strcmp(d, '-')];
end

function R = rank_rows(X)
% Rank each row ascending with average ranks for ties; NaN (0% SR) entries
% tie for last place.
    [n, k] = size(X);
    R = zeros(n, k);
    for i = 1:n
        row = X(i, :);
        nanm = isnan(row);
        row(nanm) = Inf;
        [s, idx] = sort(row);
        r = zeros(1, k);
        j = 1;
        while j <= k
            jj = j;
            while jj < k && s(jj + 1) == s(j), jj = jj + 1; end
            r(idx(j:jj)) = (j + jj) / 2;
            j = jj + 1;
        end
        R(i, :) = r;
    end
end

function fprintf_both(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end

function plot_perf_vs_nfe(resdir, D, Summary, Algos, budgets, astar, egso)
    bi = find(budgets == 15000);
    avgRk = Summary(bi).avgRank;
    [~, ord] = sort(avgRk);
    ord = ord(ord ~= egso);
    show = [egso ord(1:8)];
    cols = [0 0 0; lines(8)];
    styles = {'-', '-', '--', '-.', ':', '-', '--', '-.', ':'};
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 860]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for m = 1:5
        nexttile;
        TR = D(bi, m).Trace;                 % runs x algos x checkpoints
        ck = D(bi, m).ckpt;
        x = (1:size(TR, 3)) * ck;
        hold on;
        hs = gobjects(1, numel(show));
        for s = 1:numel(show)
            a = show(s);
            y = median(squeeze(TR(:, a, :)), 1, 'omitnan');
            y(y >= 1e4) = NaN;               % not yet feasible in the median run
            lw = 1.3; if a == egso, lw = 2.6; end
            hs(s) = plot(x, y, styles{s}, 'Color', cols(s, :), 'LineWidth', lw);
        end
        hy = yline(astar(m), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4);
        hx = xline(5000, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
        yl = ylim; ylim([astar(m) * 0.97, min(yl(2), astar(m) * 1.8)]);
        xlim([0 15000]);
        set(gca, 'FontSize', 11, 'Box', 'on');
        title(sprintf('(%c) Map %d', 'a' + m - 1, m), 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Number of function evaluations (NFE)');
        ylabel('Median best-so-far path length');
        grid on;
    end
    lgd = legend([hs hy hx], [Algos(show), {'A* reference length', 'NFE = 5,000'}], 'FontSize', 12, 'Box', 'off');
    lgd.Layout.Tile = 6;
    exportgraphics(fig, fullfile(resdir, 'Fig_PerfVsNFE_15000.png'), 'Resolution', 300);
    close(fig);
end

function write_latex(resdir, Summary, Algos, Rank100, egso, others)
    fid = fopen(fullfile(resdir, 'EqualNFE_LaTeX.txt'), 'w');
    b5 = find([Summary.fe] == 5000); b15 = find([Summary.fe] == 15000);
    % ranking table: avg rank under T=100, MaxFEs=5000, MaxFEs=15000 + per-map ranks at 15000
    [~, ord] = sort(Summary(b15).avgRank);
    fprintf(fid, '%% Ranking under the equal-NFE protocol\n');
    fprintf(fid, 'Algorithm & Map~1 & Map~2 & Map~3 & Map~4 & Map~5 & Avg. (15k) & Avg. (5k) & Avg. ($T=100$)\\\\\n\\midrule\n');
    for a = ord
        r = Summary(b15).Rank(:, a);
        fprintf(fid, '%s & %s & %s & %s & %s & %s & %.2f & %.2f & %.2f\\\\\n', Algos{a}, rk(r(1)), rk(r(2)), rk(r(3)), rk(r(4)), rk(r(5)), ...
            Summary(b15).avgRank(a), Summary(b5).avgRank(a), mean(Rank100(:, a)));
    end
    fprintf(fid, '\n%% EGSO vs best competitor (V.Mean, SR) per map and budget\n');
    for bi = [b5 b15]
        S = Summary(bi);
        for m = 1:5
            vm = S.VM(m, :); vm(egso) = Inf; [bv, bk] = min(vm);
            fprintf(fid, 'MaxFEs %d & Map~%d & %.2f & %.0f & %s & %.2f & %.0f\\\\\n', S.fe, m, S.VM(m, egso), S.SR(m, egso), Algos{bk}, bv, S.SR(m, bk));
        end
    end
    fprintf(fid, '\n%% Wilcoxon / Friedman\n');
    for bi = [b5 b15]
        S = Summary(bi);
        fprintf(fid, 'MaxFEs %d: valid-only W/T/L = %d/%d/%d; all-runs W/T/L = %d/%d/%d; Friedman chi2 = %.2f, df = %d, p = %.3g\n', S.fe, ...
            sum(S.WTLvalid(others, :), 1), sum(S.WTLall(others, :), 1), S.friedman);
    end
    fclose(fid);
end

function s = rk(r)
    if r == round(r), s = sprintf('%d', r); else, s = sprintf('%.1f', r); end
end
