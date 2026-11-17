function run_all_tests()
% RUN_ALL_TESTS  Run the smoke and unit tests in tests/ (about 2-3 min).
%
%   >> setup_egso; run_all_tests
%
% Each test is a script; a test passes if it runs to completion without an
% error. The timing scripts (test_timing_*) are excluded: they only print
% wall-clock measurements.

root = setup_egso();
tests = {'test_fmm', 'test_generate_maps', 'test_egso_smoke', ...
         'test_egso_variant_smoke', 'test_ablation_smoke', ...
         'test_classical_baselines', 'test_new_baselines'};
old = pwd; cleanup = onCleanup(@() cd(old));
cd(root);
failed = {};
for i = 1:numel(tests)
    fprintf('\n=== %s ===\n', tests{i});
    try
        run_isolated(fullfile(root, 'tests', [tests{i} '.m']));
        close all;
        fprintf('PASS %s\n', tests{i});
    catch ME
        fprintf(2, 'FAIL %s: %s\n', tests{i}, ME.message);
        failed{end + 1} = tests{i}; %#ok<AGROW>
    end
end
fprintf('\n%d/%d tests passed\n', numel(tests) - numel(failed), numel(tests));
if ~isempty(failed)
    error('EGSO:tests', 'failed: %s', strjoin(failed, ', '));
end
end

function run_isolated(file)
% Runs a test script in its own workspace (scripts may define a variable
% named 'run' or call clear).
run(file);
end
